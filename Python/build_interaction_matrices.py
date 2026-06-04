import argparse
import pandas as pd
from pathlib import Path

EXCLUDE_NAMES = {
    "unclassified orthohantavirus",
    "unclassified hantaviridae",
    "unclassified arenaviridae",
    "unclassified mammarenavirus",
    "hantaviridae",
    "arenaviridae",
}


def normalize_name(value: str) -> str:
    if pd.isna(value):
        return ""
    return str(value).strip().lower()


def load_inputs(pathogen_path: Path, host_path: Path) -> pd.DataFrame:
    pathogens = pd.read_csv(pathogen_path, low_memory=False)
    hosts = pd.read_csv(host_path, low_memory=False)

    hosts = hosts[["host_record_id", "host_species"]].dropna(subset=["host_record_id"])
    merged = pathogens.merge(hosts, on="host_record_id", how="left")
    return merged


def filter_records(df: pd.DataFrame, family: str) -> pd.DataFrame:
    out = df.copy()
    out = out[out["assigned_family"].astype(str).str.lower() == family.lower()]
    out = out.dropna(subset=["assigned_name", "host_species"])
    out = out[~out["assigned_name"].map(normalize_name).isin(EXCLUDE_NAMES)]
    out = out[~out["host_species"].map(normalize_name).isin(EXCLUDE_NAMES)]
    out = out[["assigned_name", "host_species"]].drop_duplicates()
    return out


def build_matrix(pairs: pd.DataFrame) -> pd.DataFrame:
    pairs = pairs.copy()
    pairs["value"] = 1
    matrix = (
        pairs.pivot_table(
            index="host_species",
            columns="assigned_name",
            values="value",
            aggfunc="max",
            fill_value=0,
        )
        .sort_index()
    )
    matrix = matrix.reindex(sorted(matrix.columns), axis=1)
    return matrix


def main() -> None:
    parser = argparse.ArgumentParser(description="Build host-pathogen interaction matrices.")
    parser.add_argument("--pathogen", required=True, type=Path)
    parser.add_argument("--host", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    merged = load_inputs(args.pathogen, args.host)

    for family, tag in [("Hantaviridae", "hanta"), ("Arenaviridae", "arena")]:
        pairs = filter_records(merged, family)
        pairs_path = args.out_dir / f"interaction_pairs_{tag}.csv"
        pairs.to_csv(pairs_path, index=False)

        matrix = build_matrix(pairs)
        matrix_path = args.out_dir / f"interaction_matrix_{tag}.csv"
        matrix.to_csv(matrix_path)


if __name__ == "__main__":
    main()
