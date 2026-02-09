import argparse
from pathlib import Path
import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser(description="Prepare cophylogeny inputs for arenavirus.")
    parser.add_argument("--matrix", required=True, type=Path, help="Path to interaction_matrix_arena.csv")
    parser.add_argument("--out-dir", required=True, type=Path)
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.matrix, index_col=0)
    df.index = df.index.map(lambda x: str(x).strip())
    df.columns = [str(c).strip() for c in df.columns]
    df = df.loc[df.sum(axis=1) > 0, df.sum(axis=0) > 0]

    matrix_path = args.out_dir / "arena_assoc_matrix.csv"
    df.to_csv(matrix_path)

    pairs = df.stack().reset_index()
    pairs.columns = ["host", "virus", "link"]
    pairs = pairs[pairs["link"] > 0]
    pairs_path = args.out_dir / "arena_assoc_pairs.csv"
    pairs.to_csv(pairs_path, index=False)

    (args.out_dir / "arena_hosts.txt").write_text("\n".join(df.index) + "\n")
    (args.out_dir / "arena_viruses.txt").write_text("\n".join(df.columns) + "\n")


if __name__ == "__main__":
    main()
