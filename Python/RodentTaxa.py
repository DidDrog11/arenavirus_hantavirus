# arha_tools.py
from __future__ import annotations
import subprocess
import shutil
import re, time, unicodedata
from pathlib import Path
from typing import Optional, Dict, Any, List, Union

import pandas as pd
from pygbif import species as gbif
import requests_cache

# Cache GBIF calls for speed + politeness
requests_cache.install_cache("gbif_cache", expire_after=7*24*3600)  # 7 days

class ArHaPipeline:
    """
    Utilities for (1) extracting tables from an .rds bundle via Rscript and
    (2) first-pass cleaning of host species taxonomy using GBIF.
    """

    FLAG_TOKENS = {"cf.", "aff.", "sp.", "spp.", "nr.", "gr.", "sensu", "sl.", "ss."}

    def __init__(self,
                 project_root: Union[str, Path],
                 r_script: Union[str, Path],
                 rscript_bin: Optional[str] = None):
        """
        Parameters
        ----------
        project_root : Path to repo root (e.g., .../arenavirus_hantavirus).
        r_script     : Path to the R exporter script (e.g., R/RSD2CSV.R).
        rscript_bin  : Optional path to Rscript binary; if None, uses PATH.
        """
        self.root = Path(project_root).resolve()
        self.r_script = Path(r_script).resolve()
        self.rscript = rscript_bin or shutil.which("Rscript")
        if self.rscript is None:
            raise RuntimeError("Rscript not found on PATH. Provide rscript_bin or fix PATH.")

    # --------------------- R EXPORTER ---------------------

    def list_tables(self, input_rds: Union[str, Path]) -> List[str]:
        """Return available table names inside the .rds bundle (via --list)."""
        cmd = [self.rscript, str(self.r_script), str(Path(input_rds).resolve()), "--list"]
        res = subprocess.run(cmd, text=True, capture_output=True)
        if res.returncode != 0:
            raise subprocess.CalledProcessError(res.returncode, cmd, res.stdout, res.stderr)
        # Parse list lines like: " - table_name"
        lines = [ln.strip().lstrip("-").strip() for ln in res.stdout.splitlines()
                 if ln.strip() and not ln.strip().lower().startswith("available")]
        return lines

    def extract_table(self,
                      input_rds: Union[str, Path],
                      db_name: str,
                      output_csv: Union[str, Path]) -> Path:
        """Extract a named table from .rds to UTF-8 CSV."""
        out = Path(output_csv).resolve()
        out.parent.mkdir(parents=True, exist_ok=True)
        cmd = [self.rscript, str(self.r_script),
               str(Path(input_rds).resolve()), db_name, str(out)]
        res = subprocess.run(cmd, text=True, capture_output=True)
        if res.returncode != 0:
            # Surface R error message
            msg = f"Rscript failed (code {res.returncode}).\nSTDOUT:\n{res.stdout}\nSTDERR:\n{res.stderr}"
            raise RuntimeError(msg)
        return out

    # --------------------- TAXONOMY CLEANING ---------------------

    @staticmethod
    def _normalize_name(name: str) -> Optional[str]:
        """Unicode-normalize, strip authorships/flags, keep Genus species."""
        if not isinstance(name, str):
            return None
        x = unicodedata.normalize("NFKC", name).strip()
        x = re.sub(r"\s+", " ", x)
        if x == "" or x.lower() in {"na", "n/a", "unknown", "undetermined"}:
            return None
        x = re.sub(r"\([^)]*\)", "", x).strip()  # remove (authors/years)
        x = re.sub(r"\b[A-Z][a-zA-Z-]+(?:\s*&\s*[A-Z][a-zA-Z-]+)?(?:\s*,\s*\d{4})?$", "", x).strip()
        tokens = x.split(" ")
        cleaned = []
        for t in tokens:
            if t.lower() in ArHaPipeline.FLAG_TOKENS:
                break
            cleaned.append(t)
        if len(cleaned) >= 2:
            genus = cleaned[0].capitalize()
            species = cleaned[1].lower()
            if not re.match(r"^[A-Z][a-zA-Z-]+$", genus):
                return None
            if not re.match(r"^[a-z-]+$", species):
                return None
            return f"{genus} {species}"
        return None

    @staticmethod
    def _gbif_resolve(name: str,
                      min_confidence: int = 85,
                      retries: int = 3,
                      sleep: float = 0.3) -> Dict[str, Any]:
        last_err = None
        for i in range(retries):
            try:
                res = gbif.name_backbone(name=name, strict=False)
                if not res or "matchType" not in res:
                    return {"best_match": None, "confidence": None, "status": "no_match", "note": "empty_response"}
                conf = res.get("confidence", 0)
                can = res.get("canonicalName") or res.get("scientificName")
                status = res.get("status")  # ACCEPTED, SYNONYM, etc.
                accepted_usage_key = res.get("acceptedUsageKey") or res.get("usageKey")
                accepted_name = None
                if status == "SYNONYM" and res.get("acceptedUsageKey"):
                    acc = gbif.name_usage(key=res["acceptedUsageKey"])
                    accepted_name = acc.get("canonicalName") or acc.get("scientificName")
                out = {
                    "query": name,
                    "best_match": can,
                    "accepted_name": accepted_name,
                    "usageKey": res.get("usageKey"),
                    "acceptedUsageKey": accepted_usage_key,
                    "rank": res.get("rank"),
                    "status": status,
                    "matchType": res.get("matchType"),
                    "confidence": conf,
                    "kingdom": res.get("kingdom"),
                    "phylum": res.get("phylum"),
                    "class": res.get("class"),
                    "order": res.get("order"),
                    "family": res.get("family"),
                    "genus": res.get("genus"),
                    "note": None,
                    "source": "GBIF",
                }
                if conf is None or conf < min_confidence:
                    out["note"] = f"low_confidence({conf})"
                    out["status"] = "low_confidence"
                return out
            except Exception as e:
                last_err = str(e)
                time.sleep(sleep * (i + 1))
        return {"query": name, "best_match": None, "confidence": None, "status": "error", "note": last_err, "source": "GBIF"}

    def firstpass_cleaning(self,
                           species_list: List[str],
                           min_confidence: int = 85) -> pd.DataFrame:
        """
        Returns mapping table:
        original, normalized, cleaned_name, match_status, usageKey, acceptedUsageKey, confidence, note, ...
        """
        records = []
        for raw in species_list:
            normalized = self._normalize_name(raw)
            if not normalized:
                records.append({
                    "original": raw, "normalized": None, "cleaned_name": None,
                    "match_status": "unusable_input", "source": None,
                    "usageKey": None, "acceptedUsageKey": None, "confidence": None, "note": "no_binomial"
                })
                continue
            m = self._gbif_resolve(normalized, min_confidence=min_confidence)
            cleaned = None
            status = m.get("status")
            if status in {"ACCEPTED", "low_confidence"} and m.get("best_match"):
                cleaned = m["best_match"]
            elif status == "SYNONYM" and (m.get("accepted_name") or m.get("best_match")):
                cleaned = m.get("accepted_name") or m.get("best_match")
            elif status == "no_match":
                cleaned = normalized
            records.append({
                "original": raw,
                "normalized": normalized,
                "cleaned_name": cleaned,
                "match_status": status,
                "source": m.get("source"),
                "usageKey": m.get("usageKey"),
                "acceptedUsageKey": m.get("acceptedUsageKey"),
                "confidence": m.get("confidence"),
                "note": m.get("note"),
                "order": m.get("order"),
                "family": m.get("family"),
                "genus": m.get("genus"),
            })
        return pd.DataFrame.from_records(records)

    # --------------------- MERGE/IO & METRICS ---------------------

    @staticmethod
    def merge_cleaned(host_df: pd.DataFrame, mapping: pd.DataFrame) -> pd.DataFrame:
        """Merge mapping into host_df on host_species -> original."""
        out = (host_df.merge(mapping[["original", "cleaned_name", "match_status", "confidence"]],
                             left_on="host_species", right_on="original", how="left")
                      .drop(columns=["original"])
                      .rename(columns={"cleaned_name": "species_clean"}))
        return out

    @staticmethod
    def quantify_corrections(host_df_with_clean: pd.DataFrame) -> Dict[str, Any]:
        """Counts corrected/unchanged/missing."""
        total = len(host_df_with_clean)
        valid = host_df_with_clean.dropna(subset=["species_clean"])
        corrected = (valid["host_species"] != valid["species_clean"]).sum()
        unchanged = (valid["host_species"] == valid["species_clean"]).sum()
        missing = host_df_with_clean["species_clean"].isna().sum()
        return {
            "total_rows": total,
            "with_cleaned_match": len(valid),
            "corrected": corrected,
            "unchanged": unchanged,
            "missing_or_unmatched": missing,
            "correction_rate": corrected / total if total > 0 else 0.0
        }

    @staticmethod
    def list_unmatched_originals(host_df_with_clean: pd.DataFrame) -> pd.Series:
        """Return value counts of original host_species where species_clean is NaN."""
        return host_df_with_clean.loc[host_df_with_clean["species_clean"].isna(), "host_species"].value_counts()
    

