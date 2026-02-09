#!/usr/bin/env python3
"""
Generalized taxonomy reconciliation utilities that build on the heuristics used
for Hantaviridae classification. The script reads pre-computed lookup tables
per virus family and applies standardized taxonomy labels to an input dataset.

Example:
    python TaxonomyReconciler.py \
        --input Data/TreeGrafting/.../AllViruses.tsv \
        --output Data/TreeGrafting/.../AllViruses_taxonomy.tsv \
        --lookup Hantaviridae=hantaviridae_ncbi_taxonomy_lookup.tsv \
        --lookup Arenaviridae=arenaviridae_ncbi_taxonomy_lookup.tsv
"""
from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

import pandas as pd
import re
import unicodedata


RANK_PRIORITY = {
    "species": 0,
    "subspecies": 1,
    "varietas": 2,
    "forma": 3,
    "genus": 4,
    "no rank": 5,
}


def norm_name(value: Optional[str]) -> str:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return ""
    text = unicodedata.normalize("NFKC", str(value))
    text = text.lower().strip()
    text = re.sub(r"\s+", " ", text)
    return text


def norm_token(value: Optional[str]) -> str:
    if value is None:
        return ""
    return unicodedata.normalize("NFKC", str(value)).strip().lower()


def parse_lineage(field: str) -> List[str]:
    if not isinstance(field, str):
        return []
    return [chunk.strip() for chunk in field.split(";") if chunk.strip()]


def clean_species_value(value: Optional[str]) -> str:
    if not isinstance(value, str):
        return ""
    text = value.strip().strip('"').strip()
    if not text:
        return ""
    if "," in text:
        parts = [part.strip().strip('"') for part in text.split(",") if part.strip()]
        ortho = next((p for p in parts if p.lower().startswith("orthohantavirus")), None)
        text = ortho or parts[0]
    text = text.rstrip(",; ").strip()
    return text


class TaxonomyLookup:
    """Lookup helper that loads an ETE-based TSV table."""

    def __init__(self, lookup_path: Path):
        self.lookup_path = lookup_path
        self._name_map, self._genus_whitelist = self._build_lookup()

    @property
    def genus_whitelist(self) -> Iterable[str]:
        return self._genus_whitelist

    def _build_lookup(self) -> Tuple[Dict[str, List[Dict[str, object]]], List[str]]:
        df = pd.read_csv(self.lookup_path, sep="\t")
        records = df.to_dict("records")
        name_map: Dict[str, List[Dict[str, object]]] = {}
        genus_whitelist: List[str] = []
        seen_genera: Set[str] = set()
        for rec in records:
            sci_name = rec.get("sci_name")
            rank = (rec.get("rank") or "").lower()
            if isinstance(sci_name, str) and rank == "species":
                genus = sci_name.split()[0].lower()
                if genus not in seen_genera:
                    genus_whitelist.append(genus)
                    seen_genera.add(genus)
            key = norm_name(rec.get("name_txt") or sci_name)
            if not key:
                continue
            name_map.setdefault(key, []).append(rec)
            alt_key = norm_name(rec.get("name_norm"))
            if alt_key and alt_key != key:
                name_map.setdefault(alt_key, []).append(rec)
        return name_map, genus_whitelist

    def find(self, raw_name: Optional[str]) -> Optional[List[Dict[str, object]]]:
        key = norm_name(raw_name)
        if not key:
            return None
        return self._name_map.get(key)


class BaseFamilyAssigner:
    """Shared assignment logic across virus families."""

    def __init__(
        self,
        family: str,
        lookup: TaxonomyLookup,
        placeholder_tokens: Sequence[str],
        recognized_genera: Sequence[str],
        default_placeholder: str,
        primary_name_fields: Optional[Sequence[str]] = None,
        variant_markers: Optional[Sequence[str]] = None,
    ):
        self.family = family
        self.lookup = lookup
        self.placeholder_tokens = {token.lower() for token in placeholder_tokens}
        self.recognized_genera = {norm_token(genus): genus for genus in recognized_genera}
        self.default_placeholder = default_placeholder
        self.lookup_genera = {g.lower() for g in lookup.genus_whitelist}
        self.primary_name_fields = self._build_primary_fields(primary_name_fields)
        self.variant_markers = self._build_variant_markers(variant_markers)

    def assign_row(self, row: pd.Series) -> Dict[str, object]:
        organism = self._extract_primary_name(row)
        assignment: Optional[Dict[str, object]] = None
        candidates = self.lookup.find(organism)
        if candidates:
            assignment = self._pick_assignment(candidates)
        meta_species, meta_source = self.resolve_metadata_species(row)
        species_valid = bool(meta_species and meta_source)
        if assignment and assignment.get("assignment_source") == "fallback" and species_valid:
            assignment = self._format_assignment(
                taxid=None,
                rank="species",
                name=meta_species,
                source=meta_source,
            )
        if assignment is None and species_valid:
            assignment = self._format_assignment(
                taxid=None,
                rank="species",
                name=meta_species,
                source=meta_source,
            )
        if assignment is None:
            heuristic = self.heuristic_assignment(str(organism or ""))
            if heuristic:
                assignment = heuristic
        if assignment is None:
            assignment = self._format_assignment(
                taxid=None,
                rank="unassigned",
                name=organism,
                source="unmatched",
            )
        assignment = self.post_process_assignment(row, assignment)
        assignment["assigned_name"] = clean_species_value(assignment.get("assigned_name"))
        assignment["assigned_family"] = self.family
        return assignment
    
    def _build_primary_fields(self, provided: Optional[Sequence[str]]) -> List[str]:
        fields: List[str] = []
        default_order = ["Organism_Name", "Species", "Organism", "Virus_name"]
        for value in provided or []:
            if not value:
                continue
            if value not in fields:
                fields.append(value)
        for value in default_order:
            if value not in fields:
                fields.append(value)
        return fields

    def _build_variant_markers(self, markers: Optional[Sequence[str]]) -> Set[str]:
        result: Set[str] = set()
        for marker in markers or []:
            token = norm_token(marker)
            if token:
                result.add(token)
        return result
    
    def _extract_primary_name(self, row: pd.Series) -> str:
        for column in self.primary_name_fields:
            value = row.get(column)
            if isinstance(value, str) and value.strip():
                return clean_species_value(value)
        return ""

    def _should_apply_variant_placeholder(self, stripped: str) -> bool:
        if not self.variant_markers:
            return False
        if not stripped:
            return False
        if not (
            any(ch.isdigit() for ch in stripped)
            or "/" in stripped
            or "_" in stripped
        ):
            return False
        lowered = stripped.lower()
        return any(marker in lowered for marker in self.variant_markers)

    def post_process_assignment(self, row: pd.Series, assignment: Dict[str, object]) -> Dict[str, object]:
        return assignment

    def resolve_metadata_species(self, row: pd.Series) -> Tuple[Optional[str], str]:
        raw = row.get("Species")
        if not isinstance(raw, str) or not raw.strip():
            return None, ""
        genus_raw = row.get("Genus")
        genus_norm = norm_token(genus_raw)
        species = clean_species_value(raw)
        if self._is_placeholder_species(species):
            placeholder = self.placeholder_for_genus(genus_raw)
            return placeholder, "metadata_placeholder"
        if self._is_generic_species(species):
            return None, ""
        if genus_norm in self.recognized_genera and self.recognized_genera[genus_norm].lower() not in species.lower():
            placeholder = self.placeholder_for_genus(genus_raw)
            return placeholder, "metadata_unclassified"
        return species, "metadata_species"

    def heuristic_assignment(self, raw_name: str) -> Optional[Dict[str, object]]:
        if not isinstance(raw_name, str):
            return None
        stripped = raw_name.strip()
        if not stripped:
            return None
        if self._is_placeholder_species(stripped):
            return self._placeholder_assignment(None, source="heuristic_placeholder")
        parts = stripped.split()
        first = parts[0].lower()
        if first in self.placeholder_tokens:
            return self._placeholder_assignment(first, source="heuristic_unclassified")
        if first in self.lookup_genera and len(parts) >= 2:
            canonical = f"{parts[0]} {parts[1]}"
            return self._format_assignment(
                taxid=None,
                rank="species",
                name=canonical,
                source="heuristic_genus",
            )
        if first.endswith("virus") and len(parts) >= 2:
            canonical = f"{parts[0]} {parts[1]}"
            return self._format_assignment(
                taxid=None,
                rank="species",
                name=canonical,
                source="heuristic_virus",
            )
        if self._should_apply_variant_placeholder(stripped):
            return self._placeholder_assignment(None, source="heuristic_variant")
        return None

    def placeholder_for_genus(self, genus_value: Optional[str]) -> str:
        canon = self.recognized_genera.get(norm_token(genus_value))
        if canon:
            return f"Unclassified {canon}"
        if isinstance(genus_value, str) and genus_value.strip():
            return f"Unclassified {genus_value.strip().title()}"
        return self.default_placeholder

    def _placeholder_assignment(self, genus_token: Optional[str], source: str) -> Dict[str, object]:
        canon = None
        if genus_token:
            canon = self.recognized_genera.get(genus_token.lower())
        label = self.placeholder_for_genus(canon)
        return self._format_assignment(
            taxid=None,
            rank="no rank",
            name=label,
            source=source,
        )

    def _pick_assignment(self, candidates: List[Dict[str, object]]) -> Optional[Dict[str, object]]:
        sorted_candidates = sorted(
            candidates,
            key=lambda r: RANK_PRIORITY.get((r.get("rank") or "no rank").lower(), 99),
        )
        for rec in sorted_candidates:
            sci_name = rec.get("sci_name") or rec.get("name_txt")
            rank = (rec.get("rank") or "no rank").lower()
            taxid = rec.get("taxid")
            if isinstance(taxid, float) and math.isnan(taxid):
                taxid = None
            if rank == "species" and sci_name:
                if self._is_generic_species(sci_name):
                    continue
                return self._format_assignment(
                    taxid=taxid,
                    rank="species",
                    name=sci_name,
                    source="species_exact",
                )
            lineage_pick = self._select_lineage_ancestor(rec)
            if lineage_pick:
                if lineage_pick["assigned_taxid"] is None and taxid is not None:
                    lineage_pick["assigned_taxid"] = int(taxid)
                return lineage_pick
            if sci_name:
                return self._format_assignment(
                    taxid=taxid,
                    rank=rec.get("rank") or "no rank",
                    name=sci_name,
                    source="fallback",
                )
        return None

    def _select_lineage_ancestor(self, record: Dict[str, object]) -> Optional[Dict[str, object]]:
        names = parse_lineage(record.get("lineage_names", ""))
        taxids = parse_lineage(record.get("lineage_taxids", ""))
        ranks = parse_lineage(record.get("lineage_ranks", ""))
        lineage = list(zip(names, taxids, ranks))
        if not lineage:
            return None
        for name, taxid, rank in reversed(lineage[:-1]):
            lname = name.lower()
            if "strain" in lname:
                continue
            if lname == "root":
                continue
            if lname.startswith("unclassified"):
                parts = lname.split(maxsplit=1)
                suffix = parts[1] if len(parts) > 1 else ""
                clean = f"Unclassified {suffix.title()}".strip()
            else:
                clean = name
            try:
                tid = int(taxid)
            except (TypeError, ValueError):
                tid = None
            return self._format_assignment(
                taxid=tid,
                rank=rank if rank else "no rank",
                name=clean,
                source="lineage_rollup",
            )
        return None

    def _format_assignment(
        self,
        taxid: Optional[object],
        rank: str,
        name: Optional[str],
        source: str,
    ) -> Dict[str, object]:
        if isinstance(taxid, float) and math.isnan(taxid):
            taxid = None
        if isinstance(taxid, int):
            fmt_taxid: Optional[int] = taxid
        elif isinstance(taxid, str) and taxid.isdigit():
            fmt_taxid = int(taxid)
        elif taxid is None:
            fmt_taxid = None
        else:
            fmt_taxid = None
        return {
            "assigned_taxid": fmt_taxid,
            "assigned_rank": rank,
            "assigned_name": name,
            "assignment_source": source,
        }

    def _is_generic_species(self, name: str) -> bool:
        lname = name.lower()
        if lname.startswith("unclassified"):
            return True
        for prefix in self.placeholder_tokens:
            prefix_with_space = f"{prefix} "
            if lname.startswith(prefix_with_space):
                return True
        return False

    @staticmethod
    def _is_placeholder_species(text: str) -> bool:
        if not isinstance(text, str):
            return False
        lowered = text.strip().lower()
        if not lowered:
            return False
        if lowered.startswith("species:") or lowered.startswith("species "):
            return True
        if lowered.startswith("species"):
            head = lowered.split(maxsplit=1)[0]
            if any(sep in lowered[:10] for sep in (":", "-", ";")):
                return True
        return False


class HantaviridaeAssigner(BaseFamilyAssigner):
    EDGE_CASE_OVERRIDES = {
        "gou virus": "Orthohantavirus seoulense",
        "araucaria virus": "Orthohantavirus andesense",
        "jabora virus": "Orthohantavirus andesense",
        "juquitiba virus": "Orthohantavirus andesense",
        "mayotte virus": "Unclassified Orthohantavirus",
        "jemez virus": "Orthohantavirus sinnombreense",
        "isla vista hantavirus": "Orthohantavirus sinnombreense",
        "monongahela hantavirus": "Orthohantavirus sinnombreense",
        "kielder hantavirus": "Orthohantavirus tulaense",
        "amur virus": "Orthohantavirus hantanense",
        "azagny virus": "Unclassified Orthohantavirus",
        "ape aime-itapua": "Orthohantavirus andesense",
        "boginia virus": "Orthohantavirus tulaense",
        "calabazo virus": "Orthohantavirus andesense",
        "lohja virus": "Orthohantavirus puumalaense",
        "dahonggou creek virus": "Orthohantavirus hantanense",
        "alto paraguay": "Orthohantavirus andesense",
        "alto paraguay hantavirus": "Orthohantavirus andesense",
        "muleshoe hantavirus": "Unclassified Orthohantavirus",
        "rusne orthohantavirus": "Unclassified Orthohantavirus",
        "jemez springs virus": "Unclassified Orthohantavirus",
        "orthohantavirus negraense": "Orthohantavirus mamorense",

    }

    def __init__(self, lookup: TaxonomyLookup, primary_name_fields: Optional[Sequence[str]] = None):
        super().__init__(
            family="Hantaviridae",
            lookup=lookup,
            placeholder_tokens=["hantavirus"],
            recognized_genera=["Orthohantavirus"],
            default_placeholder="Unclassified Orthohantavirus",
            primary_name_fields=primary_name_fields,
            variant_markers=["hv", "hantavirus", "orthohantavirus"],
        )
        self.edge_case_overrides = {norm_name(k): v for k, v in self.EDGE_CASE_OVERRIDES.items()}

    def post_process_assignment(self, row: pd.Series, assignment: Dict[str, object]) -> Dict[str, object]:
        override_target = self._lookup_edge_override(row, assignment)
        if override_target:
            assignment["assigned_name"] = override_target
            assignment["assignment_source"] = "edge_case_override"
            return assignment
        current_norm = norm_name(assignment.get("assigned_name"))
        search_name = assignment.get("assigned_name") or self._extract_primary_name(row)
        parent = self._find_orthospecies_parent(search_name)
        if parent and norm_name(parent) != current_norm:
            assignment["assigned_name"] = parent
            assignment["assignment_source"] = "lineage_orthospecies"
            return assignment
        if assignment.get("assignment_source") == "unmatched":
            return self._placeholder_assignment(None, source="heuristic_unclassified")
        return assignment

    def _lookup_edge_override(self, row: pd.Series, assignment: Dict[str, object]) -> Optional[str]:
        raw_species = self._extract_primary_name(row)
        candidates = [
            norm_name(raw_species),
            norm_name(assignment.get("assigned_name")),
            norm_name(row.get("pathogen_species_cleaned")),
        ]
        for key in candidates:
            if key and key in self.edge_case_overrides:
                return self.edge_case_overrides[key]
        return None

    def _find_orthospecies_parent(self, assigned_name: Optional[str]) -> Optional[str]:
        if not assigned_name:
            return None
        records = self.lookup.find(assigned_name)
        if not records:
            return None
        for rec in records:
            parent = self._extract_ortho_lineage_parent(rec)
            if parent:
                return parent
        return None

    def _extract_ortho_lineage_parent(self, record: Dict[str, object]) -> Optional[str]:
        names = parse_lineage(record.get("lineage_names", ""))
        ranks = parse_lineage(record.get("lineage_ranks", ""))
        for name, rank in zip(reversed(names), reversed(ranks)):
            lname = name.lower()
            if not lname.startswith("orthohantavirus "):
                continue
            if lname.startswith("orthohantavirus sp"):
                continue
            return name
        return None


class ArenaviridaeAssigner(BaseFamilyAssigner):
    EDGE_CASE_OVERRIDES = {
        "morogoro mammarenavirus": "Mammarenavirus mopeiaense",
        "catarina virus": "Mammarenavirus whitewaterense",
        "patawa virus": "Unclassified Arenaviridae",
        "pinhal virus": "Unclassified Arenaviridae",
        "middle pease river virus": "Unclassified Arenaviridae",
        "kodoko virus": "Unclassified Arenaviridae",
    }

    RECOGNIZED_GENERA = [
        "Mammarenavirus",
        "Antennavirus",
        "Hartmanivirus",
        "Innmovirus",
        "Reptanavirus",
    ]

    def __init__(self, lookup: TaxonomyLookup, primary_name_fields: Optional[Sequence[str]] = None):
        super().__init__(
            family="Arenaviridae",
            lookup=lookup,
            placeholder_tokens=["arenavirus"],
            recognized_genera=self.RECOGNIZED_GENERA,
            default_placeholder="Unclassified Arenaviridae",
            primary_name_fields=primary_name_fields,
            variant_markers=(
                ["arenavirus", "arev", "mamv"]
                + [g.lower() for g in self.RECOGNIZED_GENERA]
            ),
        )
        self.edge_case_overrides = {norm_name(k): v for k, v in self.EDGE_CASE_OVERRIDES.items()}

    def post_process_assignment(self, row: pd.Series, assignment: Dict[str, object]) -> Dict[str, object]:
        override_target = self._lookup_edge_override(row, assignment)
        if override_target:
            assignment["assigned_name"] = override_target
            assignment["assignment_source"] = "edge_case_override"
            return assignment
        if assignment.get("assignment_source") == "unmatched":
            return self._placeholder_assignment(None, source="heuristic_unclassified")
        return assignment

    def _lookup_edge_override(self, row: pd.Series, assignment: Dict[str, object]) -> Optional[str]:
        candidates = [
            norm_name(self._extract_primary_name(row)),
            norm_name(assignment.get("assigned_name")),
            norm_name(row.get("pathogen_species_cleaned")),
        ]
        for key in candidates:
            if key and key in self.edge_case_overrides:
                return self.edge_case_overrides[key]
        return None


class TaxonomyReconciler:
    def __init__(self, assigners: Dict[str, BaseFamilyAssigner]):
        self.assigners = assigners

    def annotate_file(self, input_path: Path, output_path: Path, sep: str = "\t") -> None:
        df = pd.read_csv(input_path, sep=sep, encoding="windows-1252", low_memory=False)
        annotated = self.annotate_frame(df)
        annotated.to_csv(output_path, sep=sep, index=False)
        total = len(annotated)
        resolved = (annotated["assignment_source"] != "unmatched").sum()
        supported = total - (annotated["assignment_source"] == "family_not_supported").sum()
        print(
            f"Annotated {resolved}/{total} rows "
            f"(families supported for {supported} rows). "
            f"Output written to {output_path}"
        )

    def annotate_frame(self, df: pd.DataFrame) -> pd.DataFrame:
        assignments: List[Dict[str, object]] = []
        for _, row in df.iterrows():
            family = row.get("Family") or row.get("pathogen_family")
            family_key = norm_token(family)
            assigner = None
            for name, handler in self.assigners.items():
                if norm_token(name) == family_key:
                    assigner = handler
                    break
            if assigner is None:
                assignments.append(
                    {
                        "assigned_taxid": None,
                        "assigned_rank": "unassigned",
                        "assigned_name": row.get("Species") or row.get("Organism_Name"),
                        "assignment_source": "family_not_supported",
                        "assigned_family": row.get("Family"),
                    }
                )
                continue
            assignments.append(assigner.assign_row(row))
        assign_df = pd.DataFrame(assignments)
        return pd.concat([df, assign_df], axis=1)


def parse_lookup_args(values: Sequence[str]) -> Dict[str, Path]:
    paths: Dict[str, Path] = {}
    for raw in values:
        if "=" not in raw:
            raise ValueError(f"Lookup argument '{raw}' must be of the form Family=PATH")
        family, path_str = raw.split("=", 1)
        fam = family.strip()
        if not fam:
            raise ValueError(f"Empty family label in lookup argument '{raw}'")
        path = Path(path_str).expanduser()
        paths[fam] = path
    return paths


def build_assigners(lookup_map: Dict[str, Path], species_columns: Sequence[str]) -> Dict[str, BaseFamilyAssigner]:
    assigners: Dict[str, BaseFamilyAssigner] = {}
    if "Hantaviridae" in lookup_map:
        assigners["Hantaviridae"] = HantaviridaeAssigner(TaxonomyLookup(lookup_map["Hantaviridae"]), primary_name_fields=species_columns)
    if "Arenaviridae" in lookup_map:
        assigners["Arenaviridae"] = ArenaviridaeAssigner(TaxonomyLookup(lookup_map["Arenaviridae"]), primary_name_fields=species_columns)
    return assigners


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Apply taxonomy lookups for Hantaviridae and Arenaviridae.")
    parser.add_argument("--input", type=Path, required=True, help="Path to the TSV containing genome metadata.")
    parser.add_argument("--output", type=Path, required=True, help="Destination file for the annotated TSV.")
    parser.add_argument("--species-column", 
                        action="append",
                        dest="species_columns",
                        default=["species_cleaned"], 
                        help=("Column name to prioritize when extracting organism/species labels. Repeat the argument to specify multiple columns in priority order."),
    )
    parser.add_argument(
        "--lookup",
        action="append",
        default=[],
        help="Mapping of Family=PATH for lookup tables (one argument per family). e.g. --lookup Hantaviridae=path/to/hantaviridae_lookup.tsv --lookup Arenaviridae=path/to/arenaviridae_lookup.tsv",
    )
    parser.add_argument(
        "--sep",
        default="\t",
        help="Delimiter used for both input and output files (default: tab).",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if not args.lookup:
        raise SystemExit("At least one --lookup Family=PATH argument is required.")
    lookup_map = parse_lookup_args(args.lookup)
    assigners = build_assigners(lookup_map, args.species_columns)
    if not assigners:
        raise SystemExit("No supported families found in --lookup arguments.")
    reconciler = TaxonomyReconciler(assigners)
    reconciler.annotate_file(args.input, args.output, sep=args.sep)


if __name__ == "__main__":
    main()
