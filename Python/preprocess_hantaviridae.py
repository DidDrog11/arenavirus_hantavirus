from Bio import SeqIO
import subprocess
import os
import sys
import re

TARGET_SPECIES = {
    "Mobatvirus novaense",
    "Orthohantavirus andesense",
    "Orthohantavirus asamaense",
    "Orthohantavirus asikkalaense",
    "Orthohantavirus bayoui",
    "Orthohantavirus boweense",
    "Orthohantavirus caobangense",
    "Orthohantavirus chocloense",
    "Orthohantavirus dabieshanense",
    "Orthohantavirus dobravaense",
    "Orthohantavirus hantanense",
    "Orthohantavirus jejuense",
    "Orthohantavirus khabarovskense",
    "Orthohantavirus mamorense",
    "Orthohantavirus maporalense",
    "Orthohantavirus montanoense",
    "Orthohantavirus moroense", #now carrizalense
    "Orthohantavirus necocliense", #sp. 'necocliense'
    "Orthohantavirus nigrorivense",
    "Orthohantavirus oxbowense", #sp. 'oxbowense'
    "Orthohantavirus prospectense", 
    "Orthohantavirus puumalaense",
    "Orthohantavirus rockportense",
    "Orthohantavirus sangassouense", #Orthohantavirus sangassouense
    "Orthohantavirus seewisense", #Orthohantavirus sp. 'seewisense'
    "Orthohantavirus seoulense",
    "Orthohantavirus sinnombreense",
    "Orthohantavirus thailandense",
    "Orthohantavirus tigrayense",
    "Orthohantavirus tulaense",
    "Thottimvirus imjinense",
    "Thottimvirus thottapalayamense",
}

# Alternative names encountered in headers. Keys are canonical names to keep
# (as in pathogen_taxonomy_pcr.csv), values are substrings to match.
ALIAS_MAP = {
    "Orthohantavirus moroense": ["carrizalense"],
    "Orthohantavirus necocliense": ["necocliense"],
    "Orthohantavirus oxbowense": ["oxbowense"],
    "Orthohantavirus mamorense": ["marmorense", "mamorense", "negraense", "nc_038505.1", "d1r27_ssgp1"],
    "Orthohantavirus sangassouense": ["sangassouense"],
    "Orthohantavirus seewisense": ["seewisense"],
    "Orthohantavirus nigrorivense": ["nigrorivense", "black creek canal", "nc_043075.1"],
    "Orthohantavirus prospectense": ["prospect hill", "nc_038938.1"],
    "Orthohantavirus tigrayense": ["tigrayense", "ku934010.1"],
}


def capture_nucleoprotein(fasta_file, output_file):
    sequences = SeqIO.parse(fasta_file, "fasta")
    kept = 0
    with open(output_file, "w") as out_f:
        for record in sequences:
            desc = record.description.lower()
            has_np = re.search(r"\bnp\b", desc) is not None
            has_old_tag = "htnvssgp1" in desc
            has_n_protein = "n protein" in desc
            special_include = any(
                key in desc
                for key in (
                    "black creek canal",
                    "prospect hill",
                    "nc_043075.1",
                    "nc_038938.1",
                )
            )
            if "nucleo" in desc or has_np or has_old_tag or has_n_protein or special_include:
                SeqIO.write(record, out_f, "fasta")
                kept += 1
    return kept


def filter_trimmed_alignment(input_fasta, output_fasta, target_species):
    targets = {s.lower() for s in target_species}
    kept = 0
    with open(output_fasta, "w") as out_f:
        for record in SeqIO.parse(input_fasta, "fasta"):
            desc = record.description
            desc_lower = desc.lower()
            canonical = resolve_canonical_species(desc_lower, targets)
            if canonical:
                cleaned = normalize_header(desc)
                if cleaned and "|" in cleaned:
                    cleaned = apply_canonical_species(cleaned, canonical)
                if cleaned:
                    first_token = cleaned.split()[0]
                    record.id = first_token
                    record.name = first_token
                    record.description = cleaned
                SeqIO.write(record, out_f, "fasta")
                kept += 1
    return kept


def normalize_header(description: str) -> str:
    cleaned = description
    cleaned = re.sub(r"complement\(([^)]*)\)", r"\1", cleaned, flags=re.IGNORECASE)
    cleaned = cleaned.replace("complement", "")
    cleaned = " ".join(cleaned.split())
    if "|" not in cleaned:
        return re.sub(r"\s+", "_", cleaned)
    parts = [p.strip() for p in cleaned.split("|")]
    date_raw = parts[-1]
    date_norm = normalize_date(date_raw)
    if date_norm:
        parts[-1] = date_norm
    parts = [re.sub(r"\s+", "_", p) for p in parts]
    return "|".join(parts)


def resolve_canonical_species(desc_lower: str, targets_lower: set[str]) -> str | None:
    # Prefer exact canonical matches in the header
    for canon in TARGET_SPECIES:
        if canon.lower() in desc_lower:
            return canon
    # Otherwise try aliases mapped to canonical names
    for canon, aliases in ALIAS_MAP.items():
        for alias in aliases:
            if alias in desc_lower:
                return canon
    return None


def apply_canonical_species(cleaned: str, canonical: str) -> str:
    parts = [p.strip() for p in cleaned.split("|")]
    if len(parts) >= 3:
        parts[2] = canonical.replace(" ", "_")
    return "|".join(parts)


def normalize_date(date_raw: str) -> str:
    text = (date_raw or "").strip()
    if not text:
        return ""
    # yyyy-mm-dd
    if len(text) == 10 and text[4] == "-" and text[7] == "-":
        return text
    # yyyy-mm
    if len(text) == 7 and text[4] == "-":
        return f"{text}-01"
    # yyyy
    if len(text) == 4 and text.isdigit():
        return f"{text}-06-01"
    return text


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Preprocess Hantaviridae sequences.")
    parser.add_argument("--input-fasta", required=True, type=str, help="Input FASTA file with viral sequences.")
    parser.add_argument("--output-fasta", required=True, type=str, help="Output FASTA file for nucleoprotein sequences.")
    args = parser.parse_args()

    kept = capture_nucleoprotein(args.input_fasta, args.output_fasta)
    if kept == 0:
        sys.exit("No sequences matched 'nucleo' in header; nothing to align.")

    output_aligned = f"{os.path.splitext(args.output_fasta)[0]}_aligned.fasta"
    align_cmd = ["augur", "align", "--sequences", args.output_fasta, "--output", output_aligned]
    trimmed_path = f"{os.path.splitext(output_aligned)[0]}_trimmed.fasta"
    filtered_path = f"{os.path.splitext(output_aligned)[0]}_trimmed_filtered.fasta"
    trimal_cmd = ["trimal", "-in", output_aligned, "-out", trimmed_path, "-automated1", "-keepheader"]

    subprocess.run(align_cmd, check=True)
    subprocess.run(trimal_cmd, check=True)

    kept_filtered = filter_trimmed_alignment(trimmed_path, filtered_path, TARGET_SPECIES)
    if kept_filtered == 0:
        sys.exit("No sequences matched target Hantaviridae species in trimmed alignment.")
