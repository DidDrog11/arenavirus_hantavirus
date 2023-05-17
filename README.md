# Synthesising Arena- and Hantavirus data from rodents to understand current known host distributions and viral pathogens

David Simons and Steph Seifert

Rodents are global hosts of zoonotic pathogens and potential hosts of novel pathogens of epidemic potential. Existing efforts to catalogue host-pathogen associations in these species are limited by global datasets which lack temporal and geographic specificity. Current research is hindered by spatial-, host taxa- and temporal biases within these datasets that are challenging to quantify. Recent work within the Verena consortium has synthesised publicly available data on Alpha- and Beta-coronaviruses in bats which provides a template for this work to be expanded to rodents. In this proposal, I outline a project to expand on a database produced during my PhD studies on rodent-pathogen associations in West Africa to a global scale and incorporating additional virus-host features which will support future hypothesis-testing. The database produced will provide a novel and accessible resource to explore a range of questions at a global scale about rodents and their pathogens. The produced data will be deposited within the Global Biodiversity Information Facility to support wider re-use as they continue to expand their GBIF health programme.

The project is proposed to be limited to two important, globally distributed, rodent-associated virus families in the order Bunyavirales, Hantaviridae and Arenaviridae, which contain several known zoonoses. These viruses have similarities in genomic architecture, segmented and negative sense, which allow us to ask questions about viral evolution, geography and niche-overlap. A literature search of peer-reviewed and pre-printed articles, alongside ecological reports and “grey” literature will be conducted with data extraction. Extracted data will contain information on host presence, the number of individuals identified, trapping or sampling effort, temporo-spatial information, and results from pathogen assays. Where possible, database records will be associated with pathogen sequences stored on repositories including GenBank and pathogenesis studies from experimental challenges in laboratory models. An initial search and data extraction from a sample of studies have been conducted to construct the database aligned to Darwin Core Terms and to pilot data abstraction steps (see this repository). We will use these data to test the hypothesis that niche overlap facilitates viral reassortment between species and that reassortment is more likely to be detected between closely related viruses. 

Recombination is dampened in negative sense RNA viruses relative to positive sense RNA viruses which undergo homologous recombination. The segmented genomic architecture of Bunyavirales allows for reassortment of viral genomic segments that facilitates rapid acquisition of mutations that may contribute to host breadth, replication dynamics, transmission and stability, or pathogenesis. Viral reassortment is driven by both ecological factors (e.g. host overlap) and molecular compatibility (e.g. coinfection of a single cell with two virus particles, functional compatibility of reassorted segments). Using a Bayesian stochastic search variable selection and GLM framework implemented in BEAST v2, we will determine the relative contributions of ecological, geographic, and molecular factors contributing to cross-species transmission and reassortment in rodent-associated arenaviruses and hantaviruses. We will perform co-phylogenetic analysis on available pathogen and host sequences to investigate adaptation of these viruses to rodent hosts. 
   
This project will provide an additional source of data for researchers across the Verena consortium and beyond. Future research could include exploration of changes in host or pathogen presence and the potential for spatio-temporal analysis of pathogen sharing. Based on prior work, I anticipate that data synthesis can be completed within a 3-month timeline, with initial searches and screening of the literature in month 1, and data extraction in months 2-3. Phylogenomic and generalized linear models can be implemented in months 4-6 including computational time. 

Participating in this fellowship programme will allow me to contribute to and access a network of researchers interested in endemic and emerging zoonoses in a changing environment will be of immense support in the final stages of my PhD as I continue to develop my research career and think about next steps. To complete the project, I will further develop my skills in producing and maintaining a dynamic database and gain new skills in mixed data modelling including genomic data. Finally, I share many of the core values of the network and hope to continue to promote open science and equitable access to data and research outputs in my future career. 

## Useful resources

https://www.bv-brc.org/ - BACTERIAL AND VIRAL BIOINFORMATICS RESOURCE CENTER  
https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/ - NCBI VIRUS   
https://www.neonscience.org/ - NEON 

# Data abstraction

Data from studies meeting inclusion criteria will be abstracted into a Google sheet. The current version (v2) is available at the following [link](https://docs.google.com/spreadsheets/d/14eW_YwSP6EWWuDrnsvDX-vwTi-7KnVyRdYL8FDqYivk/edit?usp=sharing), the sheet can be commented on but not modified. Five sheets have been created to capture the required data. Column titles where appropriate have been harmonised to the [Darwin Core Terms](https://dwc.tdwg.org/terms/) for subsequent submission to the [Global Biodiversity Information Facility](https://www.gbif.org/publishing-data) database. The column names, accepted values and DWC reference (where applicable), are shown below.

<details>
   
   <summary>Descriptive sheet :books:</summary>

This sheet will contain the general information about included studies. Data contributing to the dataset will have a single row, if there are studies which are complementary and the data has been pooled from the same sampling activity additional articles will be contained in the `linked_manuscripts` field. 

| Variable name | Description | Values | DWC term link |
| --- | -- | -- | -- |
| study_id | An identifier for the included study | numeric | NA |
| pubmed_id | The Pubmed ID of the included study | numeric | NA |
| bibliographicCitation | Full link to the DOI or webpage hosting the study | string | [link](https://dwc.tdwg.org/list/#dcterms_bibliographicCitation) |
| identifiedBy | The first author of the included study | string | [link](https://dwc.tdwg.org/list/#dwc_identifiedBy) |
| datasetName | The title of the included study | string | [link](https://dwc.tdwg.org/list/#dwc_datasetName) |
| journal | The journal of publication | string | NA |
| study_design | The method of sampling of small rodents and pathogens | protected string currently `prospective rodent sampling`, `purposeful rodent sampling`, `sampling around cases and control locations`, `not described` | NA |
| sampling_effort | A measure of the intensity of sampling effort | string. Ideally capturing trap-nights | NA |
| data_access | A descriptor of the availability of included data | protected string currently `summarised data` | NA |
| linked_manuscripts | An identifier for manuscripts where data is separated across publications to prevent duplication of records | NA |
| rightsHolder | A description of whether the copyright remains with the author or publisher | string | [link](https://dwc.tdwg.org/list/#dcterms_rightsHolder) |
| license | A description of the permissions associated with the included study | string | [link](https://dwc.tdwg.org/list/#dcterms_license) |

</details>

<details>
   
   <summary>Rodent sheet :mouse:</summary>

Each occurrence of a species will have a single record. If the study breaks down occurrences by date and location this will be replicated here. If summary data only is available across the locations of sampling or the entire study the data will be grouped as reported.

| Variable name | Description | Values | DWC term link |
| --- | -- | -- | -- |
| rodent_record_id | A unique identifier of the rodent species, trapped at a specific location (decimalLatitude and decimalLongitude) and time point (eventDate). Will be converted to a unique occurrenceID  | numeric | [occurrenceID](https://dwc.tdwg.org/list/#dwc_occurrenceID) |
| study_id | A linking value to the `study_id` in the `Descriptive sheet` | numeric | NA |
| eventDate | The date-time or interval where the occurrence or sampling occurred. Dates conform to ISO 8601-1:2019 | date | [link](https://dwc.tdwg.org/list/#dwc_eventDate) |
| basisOfRecord | The nature of the data record | protected string, most likely `Human observation` for live trapping | [link](https://dwc.tdwg.org/list/#dwc:basisOfRecord) |
| taxonRank | The rank of the identification of the occurrence | protected string, typically `species`, `genus` or `family` | [link](https://dwc.tdwg.org/list/#dwc:taxonRank) |
| genus | The name of the genus in which the occurrence is from | string | [link](https://dwc.tdwg.org/list/#dwc_genus) |
| scientificName | The full scientific name at the lowest level taxonomic rank that is identified | string | [link](https://dwc.tdwg.org/list/#dwc_scientificName) |
| locality | The specific description of the sampling location in the included study, this will be the highest resolution at which occurrences are identified | string | [link](https://dwc.tdwg.org/list/#dwc_locality) |
| country | The name of the country in which the locality occurs | string | [link](https://dwc.tdwg.org/list/#dwc_country) |
| verbatimLocality | This will contain a description of the locality that may be helpful for additional analysis but does not provide location data to higher resolution (i.e., the habitat type sampled). | string | [link](https://dwc.tdwg.org/list/#dwc_verbatimLocality) |
| coordinate_resolution | A description of the resolution of the subsequent coordinates, this will be used to estimate a `coordinateUncertaintyInMeters` | string | [coordinateUncertaintyInMeters](https://dwc.tdwg.org/list/#dwc_coordinateUncertaintyInMeters) |
| decimalLatitude | Latitude of sampling site to the highest resolution provided in a study. Non decimal coordinates will be converted at the point of extraction to EPSG:4326 | numeric | [link](https://dwc.tdwg.org/list/#dwc_decimalLatitude) |
| decimalLongitude | Longitude of sampling site to the highest resolution provided in a study. Non decimal coordinates will be converted at the point of extraction to EPSG:4326 | numeric | [link](https://dwc.tdwg.org/list/#dwc_decimalLongitude) |
| occurrenceStatus | Whether the occurrence was detected as present or absent at the sampling location, for the purposes here this is more correctly termed detection and non-detection | protected string, `Present` or `Absent` | [link](https://dwc.tdwg.org/list/#dwc_occurrenceStatus) |
| individualCount | The number of individuals present at the location and time of the occurrence | numeric | [link](https://dwc.tdwg.org/list/#dwc_individualCount) |

</details>

<details>
   
   <summary>Pathogen sheet :biohazard:</summary>

Pathogen data will be extracted as presented in included studies. If multiple pathogens or multiple assays for the same pathogen have been reported each of these will be given an occurrence record. These data will be linked back to the sampled rodents where possible.

| Variable name | Description | Values | DWC term link |
| --- | -- | -- | -- |
| pathogen_record_id | A unique identifier of the pathogen species, using a specific detection method, from rodents of the same species, trapped at the same location at the same time, will be converted to a unique occurrenceID in conjuction with the `rodent_record_id` | numeric |  [occurrenceID](https://dwc.tdwg.org/list/#dwc_occurrenceID) |
| associated_rodent_record_id | A link to the `rodent_record_id` from which the pathogens were sampled in the `Rodent sheet` | numeric | NA |
| grouping | A descriptor of how these pathogens are grouped. (e.g.,  at the `study`, `habitat` or `host` level) | string | NA |
| study_id | A linking value to the `study_id` in the `Descriptive sheet` | numeric | NA |
| eventDate | The date-time or interval where the occurrence or sampling occurred. Dates conform to ISO 8601-1:2019 | date | [link](https://dwc.tdwg.org/list/#dwc_eventDate) |
| pathway | A term used to describe how an organism came to be in a given place at a given time, will generally be `parasites on animal` | protected string | [link](https://dwc.tdwg.org/list/#dwc_pathway) |
| basisOfRecord | The nature of the data record | protected string, most likely `Material sample` for pathogen assays | [link](https://dwc.tdwg.org/list/#dwc:basisOfRecord) |
| host_genus | The genus of the host sampled for the pathogen | string | NA |
| associatedTaxa | A term to link the pathogen to its host organism | string | [link](https://dwc.tdwg.org/list/#dwc_associatedTaxa) |
| locality | The specific description of the sampling location in the included study, this will be the highest resolution at which occurrences are identified | string | [link](https://dwc.tdwg.org/list/#dwc_locality) |
| country | The name of the country in which the locality occurs | string | [link](https://dwc.tdwg.org/list/#dwc_country) |
| verbatimLocality | This will contain a description of the locality that may be helpful for additional analysis but does not provide location data to higher resolution (i.e., the habitat type sampled). | string | [link](https://dwc.tdwg.org/list/#dwc_verbatimLocality) |
| coordinate_resolution | A description of the resolution of the subsequent coordinates, this will be used to estimate a `coordinateUncertaintyInMeters` | string | [coordinateUncertaintyInMeters](https://dwc.tdwg.org/list/#dwc_coordinateUncertaintyInMeters) |
| decimalLatitude | Latitude of sampling site to the highest resolution provided in a study. Non decimal coordinates will be converted at the point of extraction to EPSG:4326 | numeric | [link](https://dwc.tdwg.org/list/#dwc_decimalLatitude) |
| decimalLongitude | Longitude of sampling site to the highest resolution provided in a study. Non decimal coordinates will be converted at the point of extraction to EPSG:4326 | numeric | [link](https://dwc.tdwg.org/list/#dwc_decimalLongitude) |
| taxonRank | The rank of the identification of the pathogen occurrence | protected string, typically `species`, `genus` or `family` | [link](https://dwc.tdwg.org/list/#dwc:taxonRank) |
| family | The name of the family in which the pathogen occurrence is from | string | [link](https://dwc.tdwg.org/list/#dwc_family) |
| scientificName | The full scientific name at the lowest level taxonomic rank that is identified | string | [link](https://dwc.tdwg.org/list/#dwc_scientificName) |
| identificationRemarks | A description of the method of identification of the pathogen (i.e., `antibody`, `pcr`) | string | [link](https://dwc.tdwg.org/list/#dwc_identificationRemarks) |
| occurrenceStatus | Whether the pathogen occurrence was detected as present or absent at the sampling location, for the purposes here this is more correctly termed detection and non-detection | protected string, `Present` or `Absent` | [link](https://dwc.tdwg.org/list/#dwc_occurrenceStatus) |
| occurrenceRemarks | Further information on the pathogen occurrence, here it will be limited to the number of individuals assayed using the `identificationRemarks` method for the pathogen | numeric | [link](https://dwc.tdwg.org/list/#dwc_occurrenceRemarks) |
| number_negative | To allow for inconclusive results the number determined to be negative is extracted | numeric | NA |
| organismQuantity | The number of individuals assayed that were found to be positive for the pathogen using the assay | numeric | [link](https://dwc.tdwg.org/list/#dwc_organismQuantity) |
| number_inconclusive | If the number of individuals returning an inconclusive assay is returned this will be extracted | numeric | NA |

</details>

<details>

   <summary>Pathogen sequences sheet :dna:</summary>  
  
For subsequent analysis we will obtain the sequences from NCBI where available. This sheet is to collate the accession numbers of these sequences and relate them back to the rodents from which they were sampled. It may be that more accurate geographic and temporal data is available on NCBI than reported within the manuscripts. The highest resolution data that is available will be used for subsequent analysis.

| Variable name | Description | Values | DWC term link |
| --- | -- | -- | -- |
| sequence_record_id | A unique identifier for the sequence of a pathogen | numeric | NA |
| associated_rodent_record_id | A link to the `rodent_record_id` from which the pathogens were sampled in the `Rodent sheet` | numeric | NA |
| study_id | A linking value to the `study_id` in the `Descriptive sheet` | numeric | NA |
| host_genus | The genus of the host sampled for the pathogen | string | NA |
| associatedTaxa | A term to link the pathogen sequence to its host organism | string | [link](https://dwc.tdwg.org/list/#dwc_associatedTaxa) |
| scientificName | The full scientific name | string | [link](https://dwc.tdwg.org/list/#dwc_scientificName) |
| accession_number | The record number in NCBI of the pathogen sequence | string | NA |

</details>

<details>

   <summary>Known zoonoses sheet :stethoscope:</summary>  

All included pathogens will be designated as known or not-known zoonoses. Currently this will be performed by a manual search of all the recorded micro-organsisms at species level. There may be additional resources that can contribute to this.

| Variable name | Description | Values |
| --- | -- | -- |
| pathogen_species_id | A unique identifier of the pathogen species | numeric |
| pathogen_family | The family of the pathogen species | string |
| pathogen_species | The name of the pathogen species | string |
| known_zoonosis | Whether there are records of this pathogen causing disease in humans from exposure to an animal host | logical |
| disease_name | The name(s) of diseases associated with infection from this zoonosis | string |
| icd_10 | The ICD-10 code associated with this/these disease(s) | string |
| disease_reference | A link, prefereably a DOI to an article demonstrating whether this pathogen is classed as a zoonosis | string |

</details>

# RShiny Application :computer:

A companion application is being developed to display the contained information. This interactive application will allow the dataset to be explored, subset and downloaded for subsequent reuse. This is being developed in a separate Github project called [arenavirus_hantavirus_app](https://github.com/DidDrog11/arenavirus_hantavirus_app)

