from Bio import SeqIO
import subprocess
import os
import sys
import re

TARGET_SPECIES = {
    "Mammarenavirus aporeense",
    "Mammarenavirus bearense",
    "Mammarenavirus chapareense",
    "Mammarenavirus choriomeningitidis",
    "Mammarenavirus gairoense",
    "Mammarenavirus guanaritoense",
    "Mammarenavirus juninense",
    "Mammarenavirus lassaense",
    "Mammarenavirus latinum",
    "Mammarenavirus lunaense",
    "Mammarenavirus lunkense",
    "Mammarenavirus marientalense",
    "Mammarenavirus merinoense",
    "Mammarenavirus mopeiaense",
    "Mammarenavirus oliverosense",
    "Mammarenavirus piritalense",
    "Mammarenavirus tacaribeense",
    "Mammarenavirus wenzhouense",
    "Mammarenavirus whitewaterense",
    "Mammarenavirus xapuriense",
}

def capture_nucleoprotein(fasta_file, output_file):
    sequences = SeqIO.parse(fasta_file, "fasta")
    kept = 0
    with open(output_file, "w") as out_f:
        for record in sequences:
            desc = record.description.lower()
            has_np = re.search(r"\bnp\b", desc) is not None
            if "nucleo" in desc or has_np:
                SeqIO.write(record, out_f, "fasta")
                kept += 1
    return kept


def filter_trimmed_alignment(input_fasta, output_fasta, target_species):
    targets = {s.lower() for s in target_species}
    kept = 0
    with open(output_fasta, "w") as out_f:
        for record in SeqIO.parse(input_fasta, "fasta"):
            desc = record.description
            if any(t in desc.lower() for t in targets):
                cleaned = normalize_header(desc)
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
    cleaned = re.sub(r"complement\\(([^)]*)\\)", r"\\1", cleaned, flags=re.IGNORECASE)
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
    parser = argparse.ArgumentParser(description="Preprocess Arenaviridae sequences.")
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
        sys.exit("No sequences matched target Arenaviridae species in trimmed alignment.")
