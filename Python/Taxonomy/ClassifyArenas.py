#!/usr/bin/env python3

from ete3 import NCBITaxa
import pandas as pd
import re
import unicodedata

ROOT_TAXID = 11617

def norm_name(s: str) -> str:
    """Normalize names for lookup"""
    if s is None:
        return ""
    s = unicodedata.normalize("NFKC", s)
    s = s.strip().lower()
    s = re.sub(r"\s+", " ", s)
    s = re.sub(r"[‘’´`]", "'", s)
    return s


ncbi = NCBITaxa()
# Build local taxonomy db on first run

# get subtree taxids (descendants + root)
desc_taxids = ncbi.get_descendant_taxa(ROOT_TAXID, collapse_subspecies=False) + [ROOT_TAXID]

# For canonical metadata, rank + scientific name + lineage
taxid2rank = ncbi.get_rank(desc_taxids)
taxid2name = ncbi.get_taxid_translator(desc_taxids)

rows = []
for taxid in desc_taxids:
    lineage = ncbi.get_lineage(taxid)
    lin_names = ncbi.get_taxid_translator(lineage)
    rows.append({
        "taxid": taxid,
        "sci_name": taxid2name.get(taxid),
        'rank': taxid2rank.get(taxid),
        'lineage_taxids': ';'.join(map(str,lineage)),
        'lineage_names': ';'.join([lin_names[t] for t in lineage if t in lin_names]),
        'lineage_ranks': ';'.join([taxid2rank.get(t, '') or '' for t in lineage]),
    })

meta_df = pd.DataFrame(rows)

#add names for each taxid

import sqlite3
db = sqlite3.connect(ncbi.dbfile)

taxid_set = set(desc_taxids)
names_rows = []

schema = pd.read_sql("SELECT name FROM sqlite_master WHERE type='table'", db)
tables = set(schema['name'].tolist())

if "names" in tables:
    q = """SELECT taxid, name, name_class FROM names"""
    for taxid, name, name_class in db.execute(q):
        if taxid in taxid_set:
            names_rows.append({
                'taxid': taxid,
                'name_txt': name,
                'name_class': name_class,
                'name_norm': norm_name(name),
            })
elif {"species", "synonym"}.issubset(tables):
    # Newer ETE taxonomy DB schema (>=3.1) exposes species/common names in `species`
    # and alternate spellings in `synonym`.
    for taxid, spname, common in db.execute("SELECT taxid, spname, common FROM species"):
        if taxid not in taxid_set:
            continue
        if spname:
            names_rows.append({
                'taxid': taxid,
                'name_txt': spname,
                'name_class': 'scientific name',
                'name_norm': norm_name(spname),
            })
        if common:
            names_rows.append({
                'taxid': taxid,
                'name_txt': common,
                'name_class': 'common name',
                'name_norm': norm_name(common),
            })
    for taxid, syn in db.execute("SELECT taxid, spname FROM synonym"):
        if taxid in taxid_set and syn:
            names_rows.append({
                'taxid': taxid,
                'name_txt': syn,
                'name_class': 'synonym',
                'name_norm': norm_name(syn),
            })
else:
    raise RuntimeError(
        f"Could not locate expected taxonomy tables ('names' or 'species'/'synonym'). "
        f"Found tables: {sorted(tables)}"
    )

names_df = pd.DataFrame(names_rows)

# Join names <-> metadata
lut_df = names_df.merge(meta_df, on="taxid", how="left")

#Build lookup dicts
amb = (lut_df.groupby("name_norm")["taxid"].nunique().sort_values(ascending=False))
ambiguous_keys = set(amb[amb > 1].index)

lookup = {}
lookup_ambiguous = {}

for _, r in lut_df.iterrows():
    k = r['name_norm']
    if k in ambiguous_keys:
        lookup_ambiguous.setdefault(k, set()).add(int(r['taxid']))
    else:
        lookup[k] = int(r['taxid'])

# Save table
lut_df.to_csv("arenaviridae_ncbi_taxonomy_lookup.tsv", sep="\t", index=False)
              
print("rows:", len(lut_df))
print("unique lookup keys:", len(lookup))
print("ambiguous keys:", len(lookup_ambiguous))
