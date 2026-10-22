# Project ArHa: 08_build_dwca.R
#
# Purpose: package the release tables written by 07 as a Darwin Core Archive
# (Sampling Event dataset): an Event core (one event per study) with
# Occurrence, Extended MeasurementOrFact and ResourceRelationship extensions.
# Only Darwin Core terms go into the archive; the parquet/csv release keeps the
# ArHa-specific columns. The archive is not deposited with GBIF.
#
# Output: data/data_release/dwca/ (data files, meta.xml, eml.xml) and
# data/data_release/ArHa_dwca.zip

library(here)
library(dplyr)
library(arrow)
library(countrycode)

release_dir <- here("data", "data_release")
dwca_dir <- file.path(release_dir, "dwca")
dir.create(dwca_dir, showWarnings = FALSE)

events    <- read_parquet(file.path(release_dir, "parquet", "sampling_events.parquet"))
hosts     <- read_parquet(file.path(release_dir, "parquet", "host_occurrences.parquet"))
pathogens <- read_parquet(file.path(release_dir, "parquet", "pathogen_mof.parquet"))
sequences <- read_parquet(file.path(release_dir, "parquet", "resource_relationships.parquet"))

# --- 1. Event core -------------------------------------------------------------
# The release's publication_year is not the sampling year, so
# it is left out. The event's eventDate is the span of its host records.
host_dates <- hosts |>
  filter(!is.na(eventDate)) |>
  mutate(start = sub("/.*", "", eventDate), end = sub(".*/", "", eventDate)) |>
  group_by(eventID) |>
  summarise(first = min(start), last = max(end), .groups = "drop") |>
  mutate(eventDate = if_else(first == last, first, paste0(first, "/", last)))

# Duplicate-group flags and data access level go into eventRemarks
duplicate_note <- if_else(!is.na(events$duplicate_group),
                          paste0("Same animals as another study: duplicate group ", events$duplicate_group,
                                 " (", gsub("_", " ", events$duplicate_role), ")"),
                          NA_character_)
access_note <- if_else(!is.na(events$dataAccessLevel), paste("Data access level:", events$dataAccessLevel), NA_character_)
remarks <- cbind(events$eventRemarks, duplicate_note, access_note)
event_remarks <- apply(remarks, 1, \(x) if (all(is.na(x))) NA_character_ else paste(x[!is.na(x)], collapse = ". "))

event_core <- events |>
  mutate(eventRemarks = event_remarks,
         sampleSizeUnit = if_else(!is.na(sampleSizeValue), "trap nights", NA_character_),
         samplingEffort = if_else(!is.na(samplingEffort), paste(samplingEffort, "sampling sessions"), NA_character_)) |>
  left_join(host_dates |> select(eventID, eventDate), by = "eventID") |>
  select(eventID, eventDate, sampleSizeValue, sampleSizeUnit, samplingEffort,
         bibliographicCitation = associatedReferences, eventRemarks)

# --- 2. Occurrence extension ---------------------------------------------------
# Records come from published tables, so basisOfRecord is MaterialCitation. A
# host record with 0 animals is a recorded absence in that sample.

# Coordinates outside the GADM 4.1 outline of their stated country were checked against
# the source papers (data/corrections/coordinate_reviews.csv). The reason each was kept goes
# into georeferenceRemarks, so that tools re-running a country check can see it.
coordinate_key <- function(iso3c, lat, lon) paste(iso3c, sprintf("%.5f", lat), sprintf("%.5f", lon))
coordinate_reviews <- read.csv(here("data", "corrections", "coordinate_reviews.csv")) |>
  filter(outcome == "accepted", confirmed) |>
  transmute(review_key = coordinate_key(iso3c, decimalLatitude, decimalLongitude), review_reason = reason) |>
  distinct(review_key, .keep_all = TRUE)
hosts <- hosts |>
  mutate(review_key = coordinate_key(countryCode, decimalLatitude, decimalLongitude)) |>
  left_join(coordinate_reviews, by = "review_key")
n_unreviewed <- sum(!is.na(hosts$dist_from_expected) & is.na(hosts$review_reason))
if (n_unreviewed > 0) warning(n_unreviewed, " records outside their stated country have no coordinate review", call. = FALSE)

occurrence_ext <- hosts |>
  mutate(basisOfRecord = "MaterialCitation",
         occurrenceStatus = case_when(individualCount > 0 ~ "present", individualCount == 0 ~ "absent"),
         kingdom = "Animalia",
         taxonID = if_else(!is.na(taxonID), paste0("https://www.gbif.org/species/", format(taxonID, scientific = FALSE, trim = TRUE)),
                           NA_character_),
         countryCode = countrycode(countryCode, origin = "iso3c", destination = "iso2c"),
         geodeticDatum = if_else(!is.na(decimalLatitude), "EPSG:4326", NA_character_),
         resolution_note = if_else(!is.na(coordinate_resolution_processed),
                                   paste0("Coordinate resolution: ", coordinate_resolution_processed, "."), NA_character_),
         review_note = if_else(!is.na(dist_from_expected),
                               paste0("Outside the GADM 4.1 outline of ", country, " (", dist_from_expected,
                                      "); checked against the source paper: ", coalesce(review_reason, "not yet reviewed.")),
                               NA_character_),
         georeferenceRemarks = case_when(!is.na(review_note) ~ paste(coalesce(resolution_note, ""), review_note) |> trimws(),
                                         TRUE ~ resolution_note)) |>
  select(eventID, occurrenceID, basisOfRecord, occurrenceStatus, scientificName, kingdom, class, order, family, genus,
         taxonRank, taxonID, verbatimIdentification, individualCount, organismQuantity, organismQuantityType,
         eventDate, verbatimEventDate, country, countryCode, stateProvince, county, municipality, higherGeographyID,
         verbatimLocality, decimalLatitude, decimalLongitude, geodeticDatum, georeferenceRemarks)

# --- 3. Extended MeasurementOrFact extension -----------------------------------
# Two measurements per pathogen record: number tested and number positive.
# Records without a host record keep their eventID only.
pathogen_name <- sub(" detection$", "", pathogens$measurementType)
mof_tested <- pathogens |>
  transmute(eventID, occurrenceID, measurementID = paste0(measurementID, "_tested"),
            measurementType = paste(pathogen_name, "- number tested"), measurementValue = number_tested,
            measurementUnit, measurementMethod, measurementRemarks)
mof_positive <- pathogens |>
  transmute(eventID, occurrenceID, measurementID = paste0(measurementID, "_positive"),
            measurementType = paste(pathogen_name, "- number positive"), measurementValue,
            measurementUnit, measurementMethod, measurementRemarks)
mof_ext <- bind_rows(mof_tested, mof_positive) |>
  filter(!is.na(measurementValue)) |>
  arrange(measurementID)

# --- 4. ResourceRelationship extension -----------------------------------------
# resourceID is a host occurrenceID or a pathogen measurementID; the core id is
# the event that record belongs to.
record_event <- bind_rows(hosts |> select(resourceID = occurrenceID, eventID),
                          pathogens |> select(resourceID = measurementID, eventID))
relationship_ext <- sequences |>
  left_join(record_event, by = "resourceID") |>
  mutate(relationshipRemarks = paste0("GenBank accession ", accession_primary,
                                      if_else(!is.na(gene_name), paste0("; gene: ", gene_name), ""),
                                      "; ", tolower(sequence_type), " sequence")) |>
  select(eventID, resourceRelationshipID, resourceID, relatedResourceID, relationshipOfResource, relationshipRemarks)

# Identifiers unique and every row linked to an event
stopifnot(!anyDuplicated(event_core$eventID), !anyDuplicated(occurrence_ext$occurrenceID),
          !anyDuplicated(mof_ext$measurementID), !anyDuplicated(relationship_ext$resourceRelationshipID),
          all(occurrence_ext$eventID %in% event_core$eventID), all(mof_ext$eventID %in% event_core$eventID),
          all(relationship_ext$eventID %in% event_core$eventID))

# --- 5. Write the data files ---------------------------------------------------
# Tab-delimited UTF-8 without quoting, so tabs and line breaks inside text
# fields are replaced with spaces.
write_dwc <- function(df, file) {
  df <- df |> mutate(across(everything(), \(x) gsub("[\t\r\n]+", " ", as.character(x))))
  write.table(df, file.path(dwca_dir, file), sep = "\t", quote = FALSE, row.names = FALSE, na = "", fileEncoding = "UTF-8")
}
write_dwc(event_core, "event.txt")
write_dwc(occurrence_ext, "occurrence.txt")
write_dwc(mof_ext, "extendedmeasurementorfact.txt")
write_dwc(relationship_ext, "resourcerelationship.txt")

# --- 6. meta.xml ---------------------------------------------------------------
# One field per column with its term URI; column 0 of each file is the id
# (core) or coreid (extensions).
dwc <- "http://rs.tdwg.org/dwc/terms/"
term_uri <- function(col) if_else(col == "bibliographicCitation", "http://purl.org/dc/terms/bibliographicCitation", paste0(dwc, col))

file_block <- function(tag, row_type, file, df) {
  id_tag <- if (tag == "core") "id" else "coreid"
  fields <- paste0('    <field index="', seq_along(df) - 1, '" term="', term_uri(names(df)), '"/>', collapse = "\n")
  paste0('  <', tag, ' encoding="UTF-8" fieldsTerminatedBy="\t" linesTerminatedBy="\n" fieldsEnclosedBy="" ',
         'ignoreHeaderLines="1" rowType="', row_type, '">\n',
         '    <files><location>', file, '</location></files>\n',
         '    <', id_tag, ' index="0"/>\n', fields, '\n  </', tag, '>')
}

meta_xml <- c('<?xml version="1.0" encoding="UTF-8"?>',
              '<archive xmlns="http://rs.tdwg.org/dwc/text/" metadata="eml.xml">',
              file_block("core", paste0(dwc, "Event"), "event.txt", event_core),
              file_block("extension", paste0(dwc, "Occurrence"), "occurrence.txt", occurrence_ext),
              file_block("extension", "http://rs.iobis.org/obis/terms/ExtendedMeasurementOrFact",
                         "extendedmeasurementorfact.txt", mof_ext),
              file_block("extension", paste0(dwc, "ResourceRelationship"), "resourcerelationship.txt", relationship_ext),
              '</archive>')
writeLines(meta_xml, file.path(dwca_dir, "meta.xml"))

# --- 7. eml.xml: people and coverage -------------------------------------------
esc <- function(x) gsub(">", "&gt;", gsub("<", "&lt;", gsub("&", "&amp;", x, fixed = TRUE), fixed = TRUE), fixed = TRUE)

authors <- tribble(
  ~given,         ~surname,         ~orcid,                ~affiliation,
  "David",        "Simons",         "0000-0001-9655-1656", "Department of Anthropology, Pennsylvania State University",
  "Ricardo",      "Rivero",         NA,                    "Paul G. Allen School of Global Health, Washington State University, Pullman, Washington, USA",
  "Grant",        "Rickard",        "0009-0005-4620-5911", "Paul G. Allen School of Global Health, Washington State University, Pullman, Washington, USA",
  "Ana",          "Martinez-Checa", "0009-0003-4263-444X", "Science Department, Natural History Museum, London, UK",
  "Harry",        "Gordon",         NA,                    "School of Biological Sciences, University of East Anglia, Norwich, UK",
  "Gregory C.",   "Milne",          "0000-0002-9594-3321", "Science Department, Natural History Museum, London, UK",
  "David W.",     "Redding",        "0000-0001-8615-1798", "Science Department, Natural History Museum, London, UK",
  "Stephanie N.", "Seifert",        "0000-0002-4397-6156", "Paul G. Allen School of Global Health, Washington State University, Pullman, Washington, USA")

person_xml <- function(tag, given, surname, orcid, affiliation) {
  org <- if (is.na(affiliation)) "      <!-- TODO: affiliation -->" else paste0("      <organizationName>", esc(affiliation), "</organizationName>")
  # userId comes last in an EML party; left out for authors without an ORCID
  user_id <- if (is.na(orcid)) NULL else paste0("      <userId directory=\"https://orcid.org/\">", orcid, "</userId>\n")
  paste0("    <", tag, ">\n",
         "      <individualName><givenName>", esc(given), "</givenName><surName>", esc(surname), "</surName></individualName>\n",
         org, "\n",
         user_id,
         "    </", tag, ">")
}
creators_xml <- paste(mapply(person_xml, "creator", authors$given, authors$surname, authors$orcid, authors$affiliation), collapse = "\n")
contact_xml <- person_xml("contact", authors$given[1], authors$surname[1], authors$orcid[1], authors$affiliation[1])

# Geographic, temporal and taxonomic coverage from the data
lat <- round(range(hosts$decimalLatitude, na.rm = TRUE), 4)
lon <- round(range(hosts$decimalLongitude, na.rm = TRUE), 4)
all_dates <- unlist(strsplit(na.omit(hosts$eventDate), "/"))
first_date <- min(all_dates); last_date <- max(all_dates)
taxon_xml <- function(rank, names) {
  paste0("        <taxonomicClassification><taxonRankName>", rank, "</taxonRankName><taxonRankValue>",
         esc(names), "</taxonRankValue></taxonomicClassification>", collapse = "\n")
}
host_orders <- sort(unique(na.omit(hosts$order)))
pathogen_families <- sort(unique(na.omit(pathogens$family_pathogen)))

# --- 8. eml.xml: assemble ------------------------------------------------------
dataset_title <- "ArHa: A globally standardised, relational database of arenavirus and hantavirus sampling in wild small mammal populations"
concept_doi <- "10.5281/zenodo.21978982"
n_fmt <- function(x) format(x, big.mark = ",")

abstract <- paste0(
  "ArHa is a database of arenavirus and hantavirus testing in wild rodents and other small mammals, extracted from the published literature ",
  "through a systematic review. This archive holds ", n_fmt(nrow(event_core)), " sampling events (one per study), ",
  n_fmt(nrow(occurrence_ext)), " host occurrence records, ", n_fmt(nrow(pathogens)), " pathogen assay records ",
  "(as number tested and number positive) and ", n_fmt(nrow(relationship_ext)), " links to GenBank sequences. ",
  "Host taxonomy is matched to GBIF and pathogen taxonomy to NCBI; each record carries the resolution of its coordinates.")

methods <- paste0(
  "Searches of PubMed, Web of Science and Google Scholar on 24 August 2023 gave 2,448 citations after de-duplication. ",
  "Included studies were extracted into structured sheets, keeping the verbatim text. An R pipeline cleans dates, host and pathogen ",
  "names, coordinates and trapping effort, matches records to GADM 4.1 administrative units, enriches GenBank accessions from NCBI, ",
  "applies documented corrections and resolves duplicate publication of the same field data. Code: https://github.com/DidDrog11/arenavirus_hantavirus")

eml_xml <- paste0(
  '<?xml version="1.0" encoding="UTF-8"?>\n',
  '<eml:eml xmlns:eml="https://eml.ecoinformatics.org/eml-2.2.0" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"\n',
  '  xsi:schemaLocation="https://eml.ecoinformatics.org/eml-2.2.0 https://rs.gbif.org/schema/eml-gbif-profile/1.3/eml.xsd"\n',
  '  packageId="https://doi.org/', concept_doi, '" system="https://doi.org" xml:lang="en">\n',
  '  <dataset>\n',
  '    <alternateIdentifier>https://doi.org/', concept_doi, '</alternateIdentifier>\n',
  '    <title xml:lang="en">', esc(dataset_title), '</title>\n',
  creators_xml, '\n',
  '    <pubDate>', format(Sys.Date()), '</pubDate>\n',
  '    <language>en</language>\n',
  '    <abstract><para>', esc(abstract), '</para></abstract>\n',
  '    <keywordSet>\n',
  paste0('      <keyword>', c("Hantaviridae", "Arenaviridae", "rodent-borne zoonoses", "wildlife disease surveillance",
                            "systematic review"), '</keyword>', collapse = "\n"), '\n',
  '      <keywordThesaurus>N/A</keywordThesaurus>\n',
  '    </keywordSet>\n',
  '    <keywordSet><keyword>Samplingevent</keyword>',
  '<keywordThesaurus>GBIF Dataset Type Vocabulary: http://rs.gbif.org/vocabulary/gbif/dataset_type.xml</keywordThesaurus></keywordSet>\n',
  '    <intellectualRights><para>This work is licensed under a <ulink url="http://creativecommons.org/publicdomain/zero/1.0/legalcode">',
  '<citetitle>Creative Commons Zero (CC0) 1.0 License</citetitle></ulink>.</para></intellectualRights>\n',
  '    <distribution><online><url function="information">https://doi.org/', concept_doi, '</url></online></distribution>\n',
  '    <coverage>\n',
  '      <geographicCoverage>\n',
  '        <geographicDescription>Global; records from ', length(unique(na.omit(hosts$countryCode))), ' countries</geographicDescription>\n',
  '        <boundingCoordinates><westBoundingCoordinate>', lon[1], '</westBoundingCoordinate><eastBoundingCoordinate>', lon[2],
  '</eastBoundingCoordinate><northBoundingCoordinate>', lat[2], '</northBoundingCoordinate><southBoundingCoordinate>', lat[1],
  '</southBoundingCoordinate></boundingCoordinates>\n',
  '      </geographicCoverage>\n',
  '      <temporalCoverage><rangeOfDates><beginDate><calendarDate>', first_date, '</calendarDate></beginDate>',
  '<endDate><calendarDate>', last_date, '</calendarDate></endDate></rangeOfDates></temporalCoverage>\n',
  '      <taxonomicCoverage>\n',
  '        <generalTaxonomicCoverage>Small mammal hosts (mainly Rodentia) and the arenaviruses and hantaviruses tested for in them</generalTaxonomicCoverage>\n',
  taxon_xml("order", host_orders), '\n',
  taxon_xml("family", pathogen_families), '\n',
  '      </taxonomicCoverage>\n',
  '    </coverage>\n',
  '    <maintenance><description><para>Versioned releases on Zenodo; development on GitHub.</para></description>',
  # "unkown" is the EML vocabulary's own spelling
  '<maintenanceUpdateFrequency>unkown</maintenanceUpdateFrequency></maintenance>\n',
  contact_xml, '\n',
  '    <methods>\n',
  '      <methodStep><description><para>', esc(methods), '</para></description></methodStep>\n',
  '    </methods>\n',
  '  </dataset>\n',
  '</eml:eml>\n')
writeLines(eml_xml, file.path(dwca_dir, "eml.xml"), useBytes = TRUE)

# --- 9. Check and zip ------------------------------------------------------------
# Both XML files must parse, and each data file must have as many columns as
# meta.xml declares.
meta_doc <- xml2::read_xml(file.path(dwca_dir, "meta.xml"))
invisible(xml2::read_xml(file.path(dwca_dir, "eml.xml")))
for (f in c("event.txt", "occurrence.txt", "extendedmeasurementorfact.txt", "resourcerelationship.txt")) {
  header <- strsplit(readLines(file.path(dwca_dir, f), n = 1, encoding = "UTF-8"), "\t")[[1]]
  declared <- xml2::xml_find_all(meta_doc, paste0("//*[d1:files/d1:location='", f, "']/d1:field"), xml2::xml_ns(meta_doc))
  stopifnot(length(header) == length(declared))
}

archive_files <- c("meta.xml", "eml.xml", "event.txt", "occurrence.txt", "extendedmeasurementorfact.txt", "resourcerelationship.txt")
zip_path <- file.path(release_dir, "ArHa_dwca.zip")
if (file.exists(zip_path)) file.remove(zip_path)
zip::zip(zip_path, files = archive_files, root = dwca_dir)
message("Darwin Core Archive written: ", zip_path, " (", n_fmt(nrow(event_core)), " events, ", n_fmt(nrow(occurrence_ext)),
        " occurrences, ", n_fmt(nrow(mof_ext)), " measurements, ", n_fmt(nrow(relationship_ext)), " relationships)")
