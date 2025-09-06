# main.py
from pathlib import Path
import pandas as pd
from RodentTaxa import ArHaPipeline

pipeline = ArHaPipeline(
    project_root="/Users/ricardorivero/Documents/GitHub/arenavirus_hantavirus",
    r_script="/Users/ricardorivero/Documents/GitHub/arenavirus_hantavirus/R/RSD2CSV.R",
)

# 1) Inspect tables in the .rds
rds = Path(pipeline.root / "data/database/Project_ArHa_database_2025-08-20.rds")
tables = pipeline.list_tables(rds)
print("Tables:", tables)

# 2) Extract one table to CSV (UTF-8)
out_csv = pipeline.extract_table(rds, "pathogen",
                                 pipeline.root / "data/test_data/paper_data.csv")
print("CSV written:", out_csv)

# 3) Load host data and clean species
host_df = pd.read_csv(pipeline.root / "data/test_data/host_data.csv")
species = host_df["host_species"].astype(str).unique().tolist()
mapping = pipeline.firstpass_cleaning(species, min_confidence=85)

# 4) Merge and write outputs
host_clean = pipeline.merge_cleaned(host_df, mapping)
host_clean.to_csv(pipeline.root / "data/test_data/host_data.cleaned.csv",
                  index=False, encoding="utf-8")
mapping.to_csv(pipeline.root / "data/test_data/host_species_cleaning_map.csv",
               index=False, encoding="utf-8")

# 5) Metrics
summary = pipeline.quantify_corrections(host_clean)
print("Summary:", summary)

# 6) Unmatched originals (for manual curation)
print(pipeline.list_unmatched_originals(host_clean))