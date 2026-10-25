# ArHa: arenavirus and hantavirus sampling in wild small mammals

ArHa is a database of published sampling of wild small mammals for arenaviruses and hantaviruses. It was compiled through a systematic review of the literature and links each study to its host records, assay results and GenBank sequences, with harmonised taxonomy and a record of how precisely each location and date was reported. The data are aligned to Darwin Core.

## Getting the data

- **Zenodo** (latest version): https://doi.org/10.5281/zenodo.21978982. The release contains the four relational tables (`sampling_events`, `host_occurrences`, `pathogen_mof`, `resource_relationships`) as Parquet and CSV files, a Darwin Core Archive, and the data dictionary.
- **ArHa Database Explorer**, to filter, map and download subsets: https://diddrog11.shinyapps.io/arenavirus_hantavirus_app/
- **This repository**: the same files are written to `data/data_release/` when the database is built.

## Documentation

- **Data paper**: Simons et al., *ArHa: a globally standardised, relational database of arenavirus and hantavirus sampling in wild small mammal populations* (in preparation). It describes the database, how it was compiled and checked, and how to use it.
- **Protocol**: Simons D, Rivero R, Martinez-Checa Guiote A, Gordon H, Milne GC, Rickard G, Redding DW, Seifert SN (2025). Protocol to produce a systematic Arenavirus and Hantavirus host-pathogen database: Project ArHa. *Wellcome Open Research* 10:227. https://doi.org/10.12688/wellcomeopenres.24037.2
- **Field definitions**: `data/data_release/data_dictionary.csv`.

## Rebuilding the database

The scripts in `R/` run in numbered order, from `00_load_data.R` to `08_build_dwca.R`; `09_figures.R` and `99_prisma_fig.R` draw the figures for the data paper. `00_load_data.R` reads the extraction data from a CSV snapshot in `data/raw_data/` (`raw_data_source <- "csv"`), which is included in the Zenodo archive. Corrections to the extracted records are listed, each with its source in the paper, in `data/corrections/`, and are applied by script; the extraction sheets themselves are not edited.

## Outputs using ArHa

1. [Viral reservoir status in small mammals emerges as a predictable life-history trait after correcting for surveillance bias, bioRxiv](https://doi.org/10.64898/2026.06.03.729761)
2. [Molecular and ecological determinants of effective reassortment in orthohantaviruses, bioRxiv](https://doi.org/10.64898/2026.06.10.731004)
3. [Open ecological, genomic and digital behavioural data resolve the 2026 MV Hondius Andes virus outbreak to the southern-Andean reservoir zone, medRxiv](https://doi.org/10.64898/2026.09.28.26364179)
4. [The Socio-economic Shield Limits Lassa Virus Spillover in Urban West Africa, Epidemiology and Infection](https://doi.org/10.1017/s0950268826101824)
5. [A protocol for biodiversity-informed wildlife disease surveillance, PLOS Ecosystems](https://doi.org/10.1371/journal.pesy.0000023)

## Citation and licence

Please cite the data through the Zenodo DOI and the data paper; `CITATION.cff` gives the details. The data and code are released under CC0 1.0 (`LICENSE`).
