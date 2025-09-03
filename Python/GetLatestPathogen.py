from pathlib import Path
import pandas as pd
from RodentTaxa import ArHaPipeline

pipeline = ArHaPipeline(
    project_root="/Users/ricardorivero/Documents/GitHub/arenavirus_hantavirus",
    r_script="/Users/ricardorivero/Documents/GitHub/arenavirus_hantavirus/R/RSD2CSV.R",
)

rds = Path(pipeline.root / "data/database/Project_ArHa_database_2025-08-20.rds")
tables = pipeline.list_tables(rds)
print("Tables:", tables)


out_csv = pipeline.extract_table(rds, "pathogen",
                                 pipeline.root / "data/test_data/pathogen_data.csv")
print("CSV written:", out_csv)