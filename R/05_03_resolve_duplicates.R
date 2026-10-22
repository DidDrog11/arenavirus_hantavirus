# Project ArHa: Resolve cross-study duplicates
# 05_03_resolve_duplicates.R
# Purpose: Some animals appear in ArHa under more than one study, because the same
# field data were published in several papers, a compilation reproduces primary
# studies, or the wrong PDF was filed under a citation and extracted again.
# Pairs were found with R/99_02_cross_study_duplicates.R and checked against the
# papers (notes in data/duplicate_check/, not tracked: collaborator material).
#
# Decisions are listed in the tables below, one row per action. Only rows with
# confirmed = TRUE are applied; unconfirmed rows are reported as a dry run so
# they can be checked before confirming. Every action works on host records and
# carries the change through to their pathogen and sequence records, so no
# foreign keys are orphaned.
#
# Rules agreed for confirmed duplicates (2026-10-20):
# - Keep the most complete source data, reshaping in script where two papers
#   hold complementary breakdowns of the same animals.
# - Same animals at a breakdown that cannot be merged: keep both, flag them with
#   a shared duplicate_group, primary = finer spatial breakdown.
# - Compilations: keep the compilation's extra data, exclude overlapping
#   primary-paper records where separable, otherwise keep the compilation.
# - A review or a wrongly filed PDF is not a source: its records go.

combined_data <- read_rds(here("data", "data_cleaning", "05_02_output.rds"))

# --- 1. Record-level actions ---------------------------------------------------
# Host records are matched on study_id and then either on an explicit list of
# host_ids (rodent_record_id, separated by ";") or on species (NA = all species;
# species_exclude removes one species) and verbatimLocality_raw (NA = any).
# action:
#   drop     - remove the host records with their pathogen and sequence records
#   move     - reassign the host records (and linked records) to to_study_id
# Relocations (new coordinates, locality and coordinate resolution) are field
# corrections in data/corrections/field_corrections.csv, applied in 00_01, so that
# the border check (03_04) and the admin units (03_05) see the final coordinates.
# The lat, lon, locality and resolution columns below are no longer used.
record_actions <- tibble::tribble(
  ~case,               ~study_id,   ~host_ids, ~species,            ~species_exclude,    ~verbatim_locality, ~action,    ~to_study_id, ~lat,     ~lon,     ~locality,         ~resolution, ~confirmed,
  # Klein 2012 (david_147) reports A. agrarius serology by firing point and
  # year (Table 1); field corrections place those rows on the O'Guinn 2008
  # (david_7) firing points (cases klein_oguinn_fp10 and _fp60).
  # Klein gives the other species only pooled across both sites; O'Guinn gives
  # them by site with Klein's serology attached, so the pooled rows go.
  "klein_oguinn_pool", "david_147", NA,        NA,                  "Apodemus agrarius", NA,                 "drop",     NA,           NA,       NA,       NA,                NA,          TRUE,
  # O'Guinn has no virus data of its own (it defers serology to later papers);
  # its A. agrarius site totals duplicate Klein's site x year rows.
  "klein_oguinn_aa",   "david_7",   NA,        "Apodemus agrarius", NA,                  NA,                 "drop",     NA,           NA,       NA,       NA,                NA,          TRUE,
  # The remaining O'Guinn rows (other species by site) join the Klein study.
  "klein_oguinn_move", "david_7",   NA,        NA,                  NA,                  NA,                 "move",     "david_147",  NA,       NA,       NA,                NA,          TRUE,
  # Abramov 2011 (grant_259) Myodes records at the 12 sites shared with Yashina
  # 2015 (david_175) are identical counts; Yashina has year per site, RT-PCR and
  # sequences. Abramov's other species and three Myodes records stay.
  "yashina_abramov",   "grant_259", "grant_2459;grant_2460;grant_2461;grant_2462;grant_2463;grant_2464;grant_2466;grant_2467;grant_2468;grant_2470;grant_2471;grant_2472;grant_2473;grant_2474;grant_2475;grant_2476;grant_2477;grant_2478;grant_2479;grant_2480;grant_2481;grant_2482;grant_2483;grant_2484;grant_2486;grant_2487;grant_2488;grant_2489;grant_2490;grant_2491",
                                               NA,                  NA,                  NA,                 "drop",     NA,           NA,       NA,       NA,                NA,          TRUE,
  # Luis 2018 (david_79) compiles the long-term Sin Nombre webs. Douglass 2001
  # (ricardo_2) reports seven of the same Montana grids for 1994-99 (matched on
  # Table 1 average MNA): C11, CMR18, GC8, GC9, P5, P6 and one Cutbank grid. Only
  # the six certain grids are dropped; which Cutbank grid (ricardo_39 or
  # ricardo_40) matches Luis Cut15 is unresolved, so both stay. Luis is kept as
  # the continuous series; Douglass's other grids and species stay.
  "luis_douglass",     "ricardo_2", "ricardo_42;ricardo_31;ricardo_34;ricardo_44;ricardo_46;ricardo_32",
                                               NA,                  NA,                  NA,                 "drop",     NA,           NA,       NA,       NA,                NA,          TRUE,
  # Calisher 1999 (david_263): its Fort Lewis and Molina webs are Luis's Hesperus
  # and Molina sites; the deer mouse rows duplicate Luis, other species stay.
  "luis_calisher",     "david_263", NA,        "Peromyscus maniculatus", NA,             NA,                 "drop",     NA,           NA,       NA,       NA,                NA,          TRUE,
  # Olsson 2003 (grant_174) is kept over Olsson 2002; field corrections move its
  # point to the centre of the study area (case olsson_2003_coords).
  # Holsomback 2009 (ana_316) vs McIntyre 2005 (david_160): only the rows that
  # are certainly the same animals (March 2002 trapline pilot, O. palustris
  # May 2004). The rest of ana_316 is trapline animals that cannot be separated.
  "holsomback_pilot",  "ana_316",   "ana_472;ana_479;ana_480;ana_488;ana_496;ana_504",
                                               NA,                  NA,                  NA,                 "drop",     NA,           NA,       NA,       NA,                NA,          TRUE,
  # Jeske 2022 (david_21) Microtus rows marked in Table S2 as taken from
  # Drewes 2021 (ana_110), which has both assays and the 2015 season
  "drewes_jeske",      "david_21",  "david_4181;david_4190;david_4191;david_4199;david_4226;david_4227;david_4270;david_4280;david_4281;david_4360;david_4362;david_4369;david_4370;david_4371;david_4378;david_4380;david_4388;david_4389;david_4398;david_4406;david_4414;david_4415;david_4416;david_4424;david_4442;david_4443",
                                               NA,                  NA,                  NA,                 "drop",     NA,           NA,       NA,       NA,                NA,          TRUE,
  # Kalpokai 2018: Jeske's own TULV positive (david_4262) on what are probably
  # the voles Drewes reports negative (ana_44); Jeske's record is kept
  "drewes_kalpokai",   "ana_110",   "ana_44",  NA,                  NA,                  NA,                 "drop",     NA,           NA,       NA,       NA,                NA,          TRUE,
  # Ahlm 1997 (david_258) reports autumn 1995 at the Olsson sites; Olsson's
  # pooled bank vole total includes these animals. Other species are not in Olsson.
  "ahlm_1997",         "david_258", NA,        "Myodes glareolus",  NA,                  NA,                 "drop",     NA,           NA,       NA,       NA,                NA,          TRUE,
  # Escutenaire 2002 (david_163) splits the 2000 paper's CMR grids (david_210)
  # into dense and low ground cover; the two strata sum exactly to the 2000 grid
  # rows from autumn 1996 to autumn 1998 (after the spring-date correction).
  # ArHa does not carry ground cover, so the 2002 strata go for those seasons.
  "escutenaire_96_98", "david_163", "david_36972;david_36973;david_36976;david_36977;david_36978;david_36979;david_36982;david_36983;david_36988;david_36989;david_36990;david_36991;david_36992;david_36993;david_36996;david_36997;david_37000;david_37001;david_37002;david_37003;david_37004;david_37005;david_37006;david_37007;david_37010;david_37011;david_37014;david_37015;david_37016;david_37017;david_37018;david_37019;david_37020;david_37021;david_37024;david_37025",
                                               NA,                  NA,                  NA,                 "drop",     NA,           NA,       NA,       NA,                NA,          TRUE,
  # Spring 1999: the 2002 paper has 3-4 more animals per grid (tissue was still
  # being processed for the 2000 paper), so its rows are kept and the 2000 grid
  # rows go. Autumn 1999 is only in the 2002 paper.
  "escutenaire_s99",   "david_210", "david_37787;david_37799;david_37812;david_37828",
                                               NA,                  NA,                  NA,                 "drop",     NA,           NA,       NA,       NA,                NA,          TRUE,
  # Tagliapietra 2009 (david_246) is pooled to one row per species; only its
  # M. arvalis row (1/5, Val Cembra 2002) is separable and equals Kallio-Kokko
  # 2006 (ana_30). Site-level re-extraction of david_246 is for ArHa v2.
  "kallio_tagliapietra","david_246", "david_38732", NA,              NA,                  NA,                 "drop",     NA,           NA,       NA,       NA,                NA,          TRUE
)

# --- 1a. Notes added to kept records -------------------------------------------
# Appended to the pathogen note of every assay record of the host
record_notes <- tibble::tribble(
  ~case,             ~host_id,      ~note,                                                                                                         ~confirmed,
  "drewes_kalpokai", "david_4262", "Same site and year as Drewes et al. 2021, which reports these voles TULV negative (0/7); Jeske et al. 2022 detected one.", TRUE
)

# --- 1c. Assay merges ------------------------------------------------------------
# Where two papers report the same animals tested with different assays or for
# different viruses, ArHa keeps one set of host records with all the assay
# records attached. The host-by-host mapping is in data/corrections/assay_merges.csv
# (case, from_host_id, to_host_id, confirmed, ...): the pathogen and sequence
# records of each from_host move to its to_host, and the from_host is removed.
assay_merges_file <- here("data", "corrections", "assay_merges.csv")
assay_merges <- if (file.exists(assay_merges_file)) {
  read_csv(assay_merges_file, col_types = cols(.default = col_character())) |>
    mutate(confirmed = as.logical(confirmed))
} else {
  tibble(case = character(), from_host_id = character(), to_host_id = character(), confirmed = logical())
}

# --- 1b. Sequences moved to the kept study's records --------------------------
# Where a dropped study carries GenBank accessions for animals that the kept
# study also holds, the sequences are re-pointed before the drop.
sequence_moves <- tibble::tribble(
  ~case,      ~from_host_id, ~to_host_id,   ~to_pathogen_id, ~confirmed,
  # Oliveira 2011 (david_353) expeditions are the first seven of Oliveira 2014
  # (david_334); same Jabora session, March 2004, positive PCR records
  "oliveira", "david_48466", "david_47829", "david_60058",   TRUE,  # Akodon montensis, 5 accessions
  "oliveira", "david_48471", "david_47835", "david_60065",   TRUE   # Oligoryzomys nigripes, 4 accessions
)

# --- 2. Study-level actions ----------------------------------------------------
# drop_study removes every remaining record of the study and its sampling event.
# duplicate_group / duplicate_role link events that keep the same animals at
# breakdowns that cannot be merged (role: primary or alternative_stratification);
# event_remarks records what was done.
study_actions <- tibble::tribble(
  ~case,                ~study_id,   ~drop_study, ~duplicate_group, ~duplicate_role,             ~event_remarks, ~confirmed,
  # Klein 2012 / O'Guinn 2008: merged above into one study
  "klein_oguinn",       "david_147", FALSE,       NA,               NA,                          "Combines Klein et al. 2012 (Hantaan virus serology by firing point and year) with O'Guinn et al. 2008 (captures by firing point); both papers report the same trapping at Firing Points 10 and 60, 2001-2005.", TRUE,
  "klein_oguinn",       "david_7",   TRUE,        NA,               NA,                          NA, TRUE,
  # Wrongly filed PDFs: the file under the citation is another paper, which is
  # extracted in full elsewhere (md5 or title checked)
  # One paper entered twice
  "seijo_twice",        "grant_31",  TRUE,        NA,               NA,                          NA, TRUE,  # same Medicina 2003 article as ana_687 (identical file)
  "kozuch_twice",       "grant_270", TRUE,        NA,               NA,                          NA, TRUE,  # ft_773 is a second scan of ft_674 (david_356)
  "arai_twice",         "grant_126", TRUE,        NA,               NA,                          NA, TRUE,  # ft_327 and ft_836 are the same file; counts identical to david_377
  "bennett_twice",      "grant_36",  TRUE,        NA,               NA,                          NA, TRUE,  # ft_654 = ft_246; abstract totals mislabelled as PCR; david_2 matches Table 1
  "jakab_twice",        "grant_129", TRUE,        NA,               NA,                          NA, TRUE,  # ft_338 = ft_736; RT-PCR 3/11 (grant_56, with sequences); 8 was 11 - 3 negatives
  # Same animals published again, nothing new
  "olsson_2002",        "grant_32",  TRUE,        NA,               NA,                          NA, TRUE,  # same 1,568 voles as Olsson 2003 (grant_174), which is kept
  "salazar_review",     "grant_189", TRUE,        NA,               NA,                          NA, TRUE,  # review chapter; 703/3,166 is Borucki 2000 Table 1 (grant_272)
  "adesina_in_oyeyiola","grant_276", TRUE,        NA,               NA,                          NA, TRUE,  # Oyeyiola 2025 (david_401) states it reuses Adesina's PCR data
  "torres_2004",        "david_197", TRUE,        NA,               NA,                          NA, TRUE,  # all 20 sites re-reported in Torres-Perez 2019 (david_10)
  "blasdell_2016",      "harry_144", TRUE,        NA,               NA,                          NA, TRUE,  # 2016 paper states it uses the 2011 dataset (david_17)
  "fulhorst_is_tesh",   "ana_375",   TRUE,        NA,               NA,                          NA, TRUE,  # ArHa holds only Tesh 1993's data (grant_2) under Fulhorst 1999
  "hjelle_1996",        "david_374", TRUE,        NA,               NA,                          NA, TRUE,  # coarser totals of Bharadwaj 1997 (david_265)
  # Raboni 2009 (david_151) is the November 2006 expedition of Raboni 2012
  # (david_186); the corrections step re-extracts david_186 by expedition and
  # moves david_151's accessions onto it, leaving this event empty
  "raboni_2009",        "david_151", TRUE,        NA,               NA,                          NA, TRUE,
  # Same animals at breakdowns that cannot be merged: keep both, flag
  "zou_zhang",          "david_257", FALSE,       "dup_zou_zhang",  "primary",                   NA, TRUE,
  "zou_zhang",          "david_277", FALSE,       "dup_zou_zhang",  "alternative_stratification", NA, TRUE,
  # Oliveira 2014 (david_334) contains all of Oliveira 2011 and revises its
  # counts; 2011's sequences are moved first (1b)
  "oliveira",           "david_353", TRUE,        NA,               NA,                          NA, TRUE,
  # Madhav 2007 pools May sessions across 10 grids, 7 of which are Luis grids;
  # the 3 others cannot be separated, so both stay, Luis primary
  "luis_madhav",        "david_79",  FALSE,       "dup_luis_madhav","primary",                   NA, TRUE,
  "luis_madhav",        "david_157", FALSE,       "dup_luis_madhav","alternative_stratification", NA, TRUE,
  # Mills 1991 (grant_252) shares the six Junín grid populations, March 1988 -
  # August 1989, with Mills 1992 (grant_111, grid rows re-added in the
  # corrections step); its other localities are unique, so both stay
  "mills_1991_1992",    "grant_111", FALSE,       "dup_mills",      "primary",                   NA, TRUE,
  "mills_1991_1992",    "grant_252", FALSE,       "dup_mills",      "alternative_stratification", NA, TRUE
)

# --- 3. Apply ------------------------------------------------------------------
# study_id is a factor in the sequence table; character throughout so that
# moved records can take a new study_id and filters do not drop NA levels
host     <- combined_data$host |> mutate(study_id = as.character(study_id))
pathogen <- combined_data$pathogen |> mutate(study_id = as.character(study_id))
sequence <- combined_data$sequence |> mutate(study_id = as.character(study_id))

# One row per applied action, written to data/duplicate_resolution.csv for the paper
action_log <- tibble(case = character(), action = character(), study_id = character(),
                     host_records = integer(), pathogen_records = integer(), sequence_records = integer())
log_action <- function(log, case, action, study_id, host = 0L, pathogen = 0L, sequence = 0L) {
  bind_rows(log, tibble(case = case, action = action, study_id = study_id,
                        host_records = as.integer(host), pathogen_records = as.integer(pathogen),
                        sequence_records = as.integer(sequence)))
}

# Host records matched by one action row
match_hosts <- function(a) {
  hit <- host$study_id == a$study_id
  if (!is.na(a$host_ids)) {
    hit <- hit & host$rodent_record_id %in% strsplit(a$host_ids, ";")[[1]]
  } else {
    if (!is.na(a$species))           hit <- hit & host$resolved_species %in% a$species
    if (!is.na(a$species_exclude))   hit <- hit & !host$resolved_species %in% a$species_exclude
    if (!is.na(a$verbatim_locality)) hit <- hit & host$verbatimLocality_raw %in% a$verbatim_locality
  }
  host$rodent_record_id[hit]
}

for (i in seq_len(nrow(record_actions))) {
  a   <- record_actions[i, ]
  ids <- match_hosts(a)
  # Sequences can be linked through the pathogen record without a host id
  path_ids <- pathogen$pathogen_record_id[pathogen$associated_rodent_record_id %in% ids]
  seq_hit  <- sequence$associated_rodent_record_id %in% ids | sequence$associated_pathogen_record_id %in% path_ids
  message(if (a$confirmed) "APPLY " else "DRY RUN ", a$case, ": ", a$action, " ", a$study_id, " - ",
          length(ids), " host, ", length(path_ids), " pathogen, ", sum(seq_hit), " sequence records")
  if (!is.na(a$host_ids) && length(ids) != length(strsplit(a$host_ids, ";")[[1]])) {
    warning("Not all listed host_ids found for ", a$case, call. = FALSE)
  }
  if (length(ids) == 0) warning("No host records matched for ", a$case, call. = FALSE)
  if (!a$confirmed) next
  action_log <- log_action(action_log, a$case, paste(a$action, "records"), a$study_id, length(ids), length(path_ids), sum(seq_hit))

  if (a$action == "drop") {
    host     <- host |> filter(!rodent_record_id %in% ids)
    pathogen <- pathogen |> filter(!pathogen_record_id %in% path_ids)
    sequence <- sequence[!seq_hit, ]
  }
  if (a$action == "move") {
    host     <- host |> mutate(study_id = if_else(rodent_record_id %in% ids, a$to_study_id, study_id))
    pathogen <- pathogen |> mutate(study_id = if_else(pathogen_record_id %in% path_ids, a$to_study_id, study_id))
    sequence$study_id[seq_hit] <- a$to_study_id
  }
}

for (i in seq_len(nrow(record_notes))) {
  n <- record_notes[i, ]
  rows <- pathogen$associated_rodent_record_id %in% n$host_id
  message(if (n$confirmed) "APPLY " else "DRY RUN ", n$case, ": note on ", sum(rows), " pathogen records of ", n$host_id)
  if (sum(rows) == 0) warning("No pathogen records for ", n$host_id, call. = FALSE)
  if (!n$confirmed) next
  pathogen$note[rows] <- if_else(is.na(pathogen$note[rows]) | pathogen$note[rows] == "",
                                 n$note, paste(pathogen$note[rows], n$note, sep = " | "))
}

for (cs in unique(assay_merges$case)) {
  am <- assay_merges |> filter(case == cs)
  ok <- am$from_host_id %in% host$rodent_record_id & am$to_host_id %in% host$rodent_record_id
  message(if (all(am$confirmed)) "APPLY " else "DRY RUN ", cs, ": merge ", sum(ok), " of ", nrow(am), " host records")
  if (!all(ok)) warning(cs, ": host records not found: ", paste(am$from_host_id[!ok], collapse = ", "), call. = FALSE)
  if (!all(am$confirmed)) next
  am <- am[ok, ]
  to_study <- host$study_id[match(am$to_host_id, host$rodent_record_id)]
  action_log <- log_action(action_log, cs, "merge host records", paste(unique(to_study), collapse = ";"), nrow(am),
                           sum(pathogen$associated_rodent_record_id %in% am$from_host_id),
                           sum(sequence$associated_rodent_record_id %in% am$from_host_id))
  # Keep the from-record's verbatim locality where the to-record has none
  # (Nitatpattana's village names exist only on the hantavirus records)
  from_vl <- host$verbatimLocality_raw[match(am$from_host_id, host$rodent_record_id)]
  to_rows <- match(am$to_host_id, host$rodent_record_id)
  fill    <- is.na(host$verbatimLocality_raw[to_rows]) | host$verbatimLocality_raw[to_rows] == ""
  host$verbatimLocality_raw[to_rows[fill]] <- from_vl[fill]
  p_rows <- match(pathogen$associated_rodent_record_id, am$from_host_id)
  s_rows <- match(sequence$associated_rodent_record_id, am$from_host_id)
  pathogen$study_id[!is.na(p_rows)] <- to_study[p_rows[!is.na(p_rows)]]
  pathogen$associated_rodent_record_id[!is.na(p_rows)] <- am$to_host_id[p_rows[!is.na(p_rows)]]
  sequence$study_id[!is.na(s_rows)] <- to_study[s_rows[!is.na(s_rows)]]
  sequence$associated_rodent_record_id[!is.na(s_rows)] <- am$to_host_id[s_rows[!is.na(s_rows)]]
  host <- host |> filter(!rodent_record_id %in% am$from_host_id)
}

for (i in seq_len(nrow(sequence_moves))) {
  m <- sequence_moves[i, ]
  rows <- sequence$associated_rodent_record_id %in% m$from_host_id
  target <- host |> filter(rodent_record_id == m$to_host_id)
  message(if (m$confirmed) "APPLY " else "DRY RUN ", m$case, ": move ", sum(rows), " sequence records from ",
          m$from_host_id, " to ", m$to_host_id, " / ", m$to_pathogen_id)
  if (nrow(target) != 1 || !m$to_pathogen_id %in% pathogen$pathogen_record_id) {
    stop("Target host or pathogen record not found for ", m$case)
  }
  if (!m$confirmed) next
  action_log <- log_action(action_log, m$case, "move sequences", target$study_id, sequence = sum(rows))
  sequence$associated_rodent_record_id[rows]   <- m$to_host_id
  sequence$associated_pathogen_record_id[rows] <- m$to_pathogen_id
  sequence$study_id[rows]                      <- target$study_id
}

# Study-level actions: new columns are always present so 07 can select them
descriptives <- combined_data$descriptives_cleaned |>
  mutate(duplicate_group = NA_character_, duplicate_role = NA_character_, event_remarks = NA_character_)

for (i in seq_len(nrow(study_actions))) {
  s <- study_actions[i, ]
  if (s$drop_study) {
    message(if (s$confirmed) "APPLY " else "DRY RUN ", s$case, ": drop study ", s$study_id, " - ",
            sum(host$study_id == s$study_id), " host, ", sum(pathogen$study_id == s$study_id), " pathogen, ",
            sum(sequence$study_id == s$study_id), " sequence records")
  } else {
    message(if (s$confirmed) "APPLY " else "DRY RUN ", s$case, ": annotate study ", s$study_id)
  }
  if (!s$study_id %in% descriptives$study_id) warning("No sampling event for ", s$study_id, call. = FALSE)
  if (!s$confirmed) next

  if (s$drop_study) {
    action_log <- log_action(action_log, s$case, "drop study", s$study_id, sum(host$study_id == s$study_id),
                             sum(pathogen$study_id == s$study_id), sum(sequence$study_id == s$study_id))
    host         <- host |> filter(study_id != s$study_id)
    pathogen     <- pathogen |> filter(study_id != s$study_id)
    sequence     <- sequence |> filter(study_id != s$study_id)
    descriptives <- descriptives |> filter(study_id != s$study_id)
  } else {
    action_log <- log_action(action_log, s$case, if (!is.na(s$duplicate_group)) "link as duplicate group" else "annotate study", s$study_id)
    rows <- descriptives$study_id == s$study_id
    descriptives$duplicate_group[rows] <- s$duplicate_group
    descriptives$duplicate_role[rows]  <- s$duplicate_role
    descriptives$event_remarks[rows]   <- s$event_remarks
  }
}

# --- 4. Checks and save --------------------------------------------------------
# Pathogen and sequence records pointing to a host or pathogen record that does
# not exist. Compared with the input so that only links broken here stop the run.
count_orphans <- function(h, p, s) {
  c(pathogen = sum(!is.na(p$associated_rodent_record_id) & !p$associated_rodent_record_id %in% h$rodent_record_id),
    sequence = sum((!is.na(s$associated_rodent_record_id) & !s$associated_rodent_record_id %in% h$rodent_record_id) |
                   (!is.na(s$associated_pathogen_record_id) & !s$associated_pathogen_record_id %in% p$pathogen_record_id)))
}
orphans_before <- count_orphans(combined_data$host, combined_data$pathogen, combined_data$sequence)
orphans_after  <- count_orphans(host, pathogen, sequence)
message("Orphaned records before / after: pathogen ", orphans_before["pathogen"], " / ", orphans_after["pathogen"],
        ", sequence ", orphans_before["sequence"], " / ", orphans_after["sequence"])
stopifnot(all(orphans_after <= orphans_before))

# Each duplicate group needs exactly one primary
groups <- descriptives |> filter(!is.na(duplicate_group)) |> count(duplicate_group, primary = duplicate_role == "primary") |>
  filter(primary)
stopifnot(all(groups$n == 1))

message("Host: ", nrow(combined_data$host), " -> ", nrow(host),
        " | pathogen: ", nrow(combined_data$pathogen), " -> ", nrow(pathogen),
        " | sequence: ", nrow(combined_data$sequence), " -> ", nrow(sequence),
        " | events: ", nrow(combined_data$descriptives_cleaned), " -> ", nrow(descriptives))

combined_data$host     <- host
combined_data$pathogen <- pathogen
combined_data$sequence <- sequence
combined_data$descriptives_cleaned <- descriptives

write_rds(combined_data, here("data", "data_cleaning", "05_03_output.rds"))
write_csv(action_log, here("data", "duplicate_resolution.csv"))
