source(here::here("R", "00_load_data.R"))

# Cleaned data will be consolidated as combined_data
combined_data <- list()

# Clean citation sheet ----------------------------------------------------
clean_citations <- combined_data_v3$inclusion_full_text %>%
  full_join(bind_rows(
    combined_data_v2$studies %>%
      dplyr::select(study_id, full_text_id),
    combined_data_v3$studies %>%
      dplyr::select(study_id, full_text_id)), by = "full_text_id") %>%
  select(-any_of(c("study_id.x", "study_id.y"))) %>%
  mutate(full_text_id = fct_inorder(full_text_id)) %>%
  relocate(study_id, .before = 1)

extractions <- clean_citations %>%
  select(full_text_id, study_id, processed, decision, reason) %>%
  left_join(bind_rows(combined_data_v2$studies, combined_data_v3$studies) %>%
              mutate(full_text_id = fct(full_text_id)),
            by = c("full_text_id")) %>%
  drop_na(full_text_id) %>%
  mutate(study_id = coalesce(study_id.x, study_id.y)) %>%
  distinct(full_text_id, study_id.x, study_id.y, .keep_all = TRUE) %>%
  group_by(full_text_id) %>%
  mutate(n_records = n())
  
# These have been extracted multiple times, they can be used to assess accuracy of extraction
multiple_extractions <- extractions %>%
  filter(n_records > 1)

single_extractions <- extractions %>%
  filter(n_records == 1) %>%
  filter(!is.na(study_id)) %>%
  select(-any_of(c("study_id.x", "study_id.y"))) %>%
  relocate(study_id, .before = 1)

excluded <- extractions %>%
  filter(str_detect(decision, "Exclude")) %>%
  distinct(full_text_id, .keep_all = TRUE) %>%
  select(full_text_id, decision, reason) %>%
  left_join(clean_citations)

no_extractions <- extractions %>%
  filter(!study_id %in% c(multiple_extractions$study_id, single_extractions$study_id)) %>%
  filter(!full_text_id %in% excluded$full_text_id)


# Clean Descriptive -------------------------------------------------------

descriptives_v2 <- combined_data_v2$studies

descriptives_v3 <- combined_data_v3$studies

# Clean Rodent sheet ------------------------------------------------------

# Check rodent IDs are unique
rodent_uid <- tibble(source = c(rep("v2", times = nrow(combined_data_v2$host)), rep("v3", times = nrow(combined_data_v3$host))),
                      row_number = c(seq(from = 1, to = nrow(combined_data_v2$host)), seq(from = 1, to = nrow(combined_data_v3$host))),
                      rodent_id = c(combined_data_v2$host$rodent_record_id, combined_data_v3$host$rodent_record_id)) %>%
  group_by(rodent_id) %>%
  mutate(n = n())

table(id_count = rodent_uid$n)

#v2
v2_raw_length = nrow(combined_data_v2$host)

#v3
v3_raw_length = nrow(combined_data_v3$host)


## Clean rodent names ------------------------------------------------------

rodent_names <- tibble(rodent_names = str_to_sentence(str_squish(sort(c(unique(combined_data_v2$host$scientificName), unique(combined_data_v3$host$scientificName))))),
                       # Correct entries where Genus has been abbreviated as these are not handled well by taxize
                       # Some spelling mistakes are too far for fuzzy matching so these are also being corrected
                       # Remove cf. from species matching, will retain the uncertainty in the extracted name
                       cleaned_rodent_names = case_when(
                         str_detect(rodent_names, "\"Oriental vole\"") ~ "Eothenomys eleusis",
                         str_detect(rodent_names, "\"Asian house shrew\"|Asian house shrew") ~ "Suncus murinus",
                         str_detect(rodent_names, "Abrothrix longipilus|A longipilis|Abrotrix longipilis") ~ "Abrothrix longipilis",
                         str_detect(rodent_names, "Abrothrix olivacea|A olivaceus") ~ "Abrothrix olivaceus",
                         str_detect(rodent_names, "Aterelix albiventrix") ~ "Atelerix albiventrix",
                         str_detect(rodent_names, "A\\.? agrarius|Aoidemus agrarius|Apodeums agrarius|Apodumus agrarius") ~ "Apodemus agrarius",
                         str_detect(rodent_names, "A\\.? niloticus|Arvicanthis dembeensis") ~ "Arvicanthis niloticus",
                         str_detect(rodent_names, "A\\.? amphibius") ~ "Arvicola amphibius",
                         str_detect(rodent_names, "A\\.? flavicollis|A flavocilis|Yellow-necked mouse|Apodemus flavicolis|Adodemus flavicollis|Apodeums flavicollis") ~ "Apodemus flavicollis",
                         str_detect(rodent_names, "Aoidemus peninsulae|Aoidmus peninsulae|Apodeums peninsulae|A peninsulae|Adodemus peninsulae") ~ "Apodemus peninsulae",
                         str_detect(rodent_names, "A\\.? sylvaticus|Apodemus syhaticus|Wood mouse|Apodemus sylvaticus s.l.|Apodeums sylvaticus|Apodemus sylvatcus") ~ "Apodemus sylvaticus",
                         str_detect(rodent_names, "Apodeums witherby") ~ "Apodemus witherby",
                         str_detect(rodent_names, "A\\.? montensis|Akodon arviculoides montensis") ~ "Akodon montensis",
                         str_detect(rodent_names, "A\\.? serrensis") ~ "Akodon serrensis",
                         str_detect(rodent_names, "A\\.? azarae|Akodon azare") ~ "Akodon azarae",
                         str_detect(rodent_names, "A\\.? caenosus") ~ "Akodon caenosus",
                         str_detect(rodent_names, "A\\.? cursor") ~ "Akodon cursor",
                         str_detect(rodent_names, "A\\.? nigrita") ~ "Akodon nigrita",
                         str_detect(rodent_names, "A\\.? toba") ~ "Akodon toba",
                         str_detect(rodent_names, "A\\.? varius") ~ "Akodon varius",
                         str_detect(rodent_names, "A\\.? scherman") ~ "Arvicola scherman",
                         str_detect(rodent_names, "B\\.? taylori") ~ "Baimoys taylori",
                         str_detect(rodent_names, "Bolomys lasiururs") ~ "Bolomys lasiurus",
                         str_detect(rodent_names, "B\\.? obscurus|B obscurus") ~ "Necromys obscurus",
                         str_detect(rodent_names, "\"Yellow-breasted rat\"") ~ "Rattus flavipectus",
                         str_detect(rodent_names, "B\\.? indica|Bandicota indicia") ~ "Bandicota indica",
                         str_detect(rodent_names, "Bandicota savelei|B savilei") ~ "Bandicota savilei",
                         str_detect(rodent_names, "Berykmys bowersi") ~ "Berylmys bowersi",
                         str_detect(rodent_names, "Bolomys lasirus|Bolomys lasiururs|B lasiurus") ~ "Bolomys lasiurus",
                         str_detect(rodent_names, "Brush mouse|Brush mice") ~ "Peromyscus boylii",
                         str_detect(rodent_names, "Crocidura cf. macmillani") ~ "Crocidura macmillani",
                         str_detect(rodent_names, "Crocidura cf. yaldeni") ~ "Crocidura yaldeni",
                         str_detect(rodent_names, "C\\.? buettikoferi") ~ "Crocidura buettikoferi",
                         str_detect(rodent_names, "C\\.? glareolus|Bank vole|Clethryonomys glareolus|C glareoulus|C glareous|Cl glareolus") ~ "Clethrionomys glareolus",
                         str_detect(rodent_names, "C\\.? laucha") ~ "Calomys laucha",
                         str_detect(rodent_names, "C\\.? musculinus") ~ "Calomys musculinus",
                         str_detect(rodent_names, "C\\.? callosus") ~ "Calomys callosus",
                         str_detect(rodent_names, "M\\.? nelsoni") ~ "Chaetodipus nelsoni",
                         str_detect(rodent_names, "Crocidura barabensis") ~ "Cricetulus barabensis",
                         str_detect(rodent_names, "Crocidura oliveri") ~ "Crocidura olivieri",
                         str_detect(rodent_names, "Crocidura rusla") ~ "Crocidura russula",
                         str_detect(rodent_names, "Crocidura shauntungensis") ~ "Crocidura shantungensis",
                         str_detect(rodent_names, "Crocidura subaveolens|C suaveolens") ~ "Crocidura suaveolens",
                         str_detect(rodent_names, "Corcidura somalica") ~ "Crocidura somalica",
                         str_detect(rodent_names, "Dendromus angolensis") ~ "Dendromus sp.",
                         str_detect(rodent_names, "Dideiphis marsupialis") ~ "Didelphis marsupialis",
                         str_detect(rodent_names, "Dipodomys ordi") ~ "Dipodomys ordii",
                         str_detect(rodent_names, "Diplodidae sagitta") ~ "Dipus sagitta",
                         str_detect(rodent_names, "D gliroides") ~ "Dromiciops gliroides",
                         str_detect(rodent_names, "Eliomys tunetae") ~ "Eliomys quercinus",
                         str_detect(rodent_names, "Echimyidae proechimys") ~ "Proechimys sp.",
                         str_detect(rodent_names, "Fossorial water vole") ~ "Arvicola scherman",
                         str_detect(rodent_names, "Citellus erthyrogenys|Citellus erthrogenys") ~ "Spermophilus erythrogenys",
                         str_detect(rodent_names, "Tristriatus trisrtriatus") ~ "Funambulus tristriatus",
                         str_detect(rodent_names, "Gerbilliscus angolae|Gerbilliscus taborae") ~ "Gerbilliscus sp.",
                         str_detect(rodent_names, "Gerbilliscus humpatensis") ~ "Gerbilliscus brantsii",
                         str_detect(rodent_names, "G glis") ~ "Glis glis",
                         str_detect(rodent_names, "Grammomys surdaster") ~ "Grammomys dolichurus",
                         str_detect(rodent_names, "Heteromiydae heteromys") ~ "Heteromys sp.",
                         str_detect(rodent_names, "H\\.? chacarius") ~ "Holochilus chacarius",
                         str_detect(rodent_names, "H\\.? brasiliensis") ~ "Holochilus brasiliensis",
                         str_detect(rodent_names, "L\\.? linulus") ~ "Lemniscomys linulus",
                         str_detect(rodent_names, "L\\.? sibiricus") ~ "Lemmus sibiricus",
                         str_detect(rodent_names, "L micropus") ~ "Loxodontomys micropus",
                         str_detect(rodent_names, "M\\.? agrestis|M\\.? agretis") ~ "Microtus agrestis",
                         str_detect(rodent_names, "M\\.? californicus") ~ "Microtus californicus",
                         str_detect(rodent_names, "M\\.? glareolus|M\\.? glareoulus|Cletherionyms glareolus|Myodes glareolus|Myodes glareoulus|Myodes glareous") ~ "Clethrionomys glareolus",
                         str_detect(rodent_names, "M\\.? minutoides") ~ "Mus minutoides",
                         str_detect(rodent_names, "M\\.? musculus") ~ "Mus musculus",
                         str_detect(rodent_names, "M\\.? natalensis") ~ "Mastomys natalensis",
                         str_detect(rodent_names, "Mastomys erythrileucus|M\\.? erythroleucus") ~ "Mastomys erythroleucus",
                         str_detect(rodent_names, "Microtus arvalis obscurus|M\\.? arvalis") ~ "Microtus arvalis",
                         str_detect(rodent_names, "Microtus fortis buchner") ~ "Microtus fortis",
                         str_detect(rodent_names, "Microtus guentheri lydius") ~ "Microtus guentheri",
                         str_detect(rodent_names, "Microtus liechtensterini") ~ "Microtus liechtensteini",
                         str_detect(rodent_names, "Montane vole|M\\.? montanus") ~ "Microtus montanus",
                         str_detect(rodent_names, "Microtus obscuruc") ~ "Microtus obscurus",
                         str_detect(rodent_names, "M oeconomus") ~ "Microtus oeconomus",
                         str_detect(rodent_names, "Microtus rossiaemeridonalis") ~ "Microtus rossiaemeridionalis",
                         str_detect(rodent_names, "Mictrotus subterraneus") ~ "Microtus subterraneus",
                         str_detect(rodent_names, "Rattus meltada") ~ "Millardia meltada",
                         str_detect(rodent_names, "Montemys delectorum") ~ "Praomys delectorum",
                         str_detect(rodent_names, "Mus nannomys") ~ "Mus sp.",
                         str_detect(rodent_names, "Mus (nannomys) mahomet") ~ "Mus mahomet",
                         str_detect(rodent_names, "Mus musculus castaneus|Mus musculus domesticus|Mus musculus musculus|Mus domesticus|Mus musculis|M domesticus") ~ "Mus musculus",
                         str_detect(rodent_names, "Muscardinus aranus") ~ "Muscardinus sp.",
                         str_detect(rodent_names, "Myodes rufocanus bedfordiae|C rufocanus|M\\.? rufocanus") ~ "Myodes rufocanus",
                         str_detect(rodent_names, "Myodes rutilus mikado|C rutilus") ~ "Myodes rutilus",
                         str_detect(rodent_names, "Nannomys minutoides") ~ "Mus minutoides",
                         str_detect(rodent_names, "Nannomys setulosus|Nannomys seulosus") ~ "Mus setulosus",
                         str_detect(rodent_names, "Muscardinus avellinarius") ~ "Muscardinus avellanarius",
                         str_detect(rodent_names, "Mustela nivales|Mustela nivalus") ~ "Mustela nivalis",
                         str_detect(rodent_names, "Neoromicia africanus") ~ "Pipistrellus nanus",
                         str_detect(rodent_names, "Neotoma montanus") ~ "Neotoma sp.",
                         str_detect(rodent_names, "N\\.? stephensi") ~ "Neotoma stephensi",
                         str_detect(rodent_names, "N fodiens") ~ "Neomys fodiens",
                         str_detect(rodent_names, "Niviventer cf. confucianus|N\\.? confucianus") ~ "Niviventer confucianus",
                         str_detect(rodent_names, "Necromys lasirurus|N lasiurus") ~ "Necromys lasiurus",
                         str_detect(rodent_names, "N\\.? indica") ~ "Nesokia indica",
                         str_detect(rodent_names, "T\\.? quadrivittatus") ~ "Neotamias quadrivittatus",
                         str_detect(rodent_names, "N\\.? albigula") ~ "Neotoma albigula",
                         str_detect(rodent_names, "Neotoma cinera|Neotoma cineria") ~ "Neotoma cinerea",
                         str_detect(rodent_names, "Dusky-footed|Neotoma fiscipes") ~ "Neotoma fuscipes",
                         str_detect(rodent_names, "N\\.? lepida") ~ "Neotoma lepida",
                         str_detect(rodent_names, "N\\.? leucodon") ~ "Neotoma leucodon",
                         str_detect(rodent_names, "N\\.? micropus|N\\.? microtus") ~ "Neotoma micropus",
                         str_detect(rodent_names, "N\\.? mexican") ~ "Neotoma mexicana",
                         str_detect(rodent_names, "Rattus fulvescens|R\\.? fulvescens") ~ "Niviventer fulvescens",
                         str_detect(rodent_names, "Rattus niviventer") ~ "Niviventer niviventer",
                         str_detect(rodent_names, "Oecomys speciousus") ~ "Oecomys speciosus",
                         str_detect(rodent_names, "Oligoryzomys albigularis") ~ "Nephelomys albigularis",
                         str_detect(rodent_names, "O\\.? delticola") ~ "Oligoryzomys delticola",
                         str_detect(rodent_names, "O\\.? fornesi") ~ "Oligoryzomys fornesi",
                         str_detect(rodent_names, "O\\.? leucogaster") ~ "Onychomys leucogaster",
                         str_detect(rodent_names, "O\\.? judex|O judex|Oxymycterus judex") ~ "Oxymycterus quaestor",
                         str_detect(rodent_names, "Oximycterus rutilans") ~ "Oxymycterus rutilans",
                         str_detect(rodent_names, "O\\.? nigripes|Ol. nigripes") ~ "Oligoryzomys nigripes",
                         str_detect(rodent_names, "Oligoyzomys microtis|Oligoryzomys mictrotis") ~ "Oligoryzomys microtis",
                         str_detect(rodent_names, "Oryzomys grupo nitidus") ~ "Oryzomys capito",
                         str_detect(rodent_names, "O\\.? palustris|Rice rat|Oryzomys paulustris") ~ "Oryzomys palustris",
                         str_detect(rodent_names, "Oligoryzomys capito") ~ "Oryzomys nitidus",
                         str_detect(rodent_names, "Oligoryzomys yunganus") ~ "Oryzomys yunganus",
                         str_detect(rodent_names, "O\\.? flavescens|Oi. flavescens|Oligoryzomysfulvescens|Oligoryzomys favescens|Oligoryomys falvescens") ~ "Oligoryzomys flavescens",
                         str_detect(rodent_names, "O\\.? chacoensis") ~ "Oligoryzomys chacoensis",
                         str_detect(rodent_names, "O longicaudatus") ~ "Oligoryzomys longicaudatus",
                         str_detect(rodent_names, "Ochromyscus niveiventris") ~ "Mus niveiventris",
                         str_detect(rodent_names, "Oryzomys cous|O\\.? couesi") ~ "Oryzomys couesi",
                         str_detect(rodent_names, "Onychomys torridue") ~ "Onychomys torridus",
                         str_detect(rodent_names, "P\\.? daltoni") ~ "Praomys daltoni",
                         str_detect(rodent_names, "Deer mouse|Deer mice") ~ "Peromyscus sp.",
                         str_detect(rodent_names, "P\\.? leicopus|P\\.? leucopus") ~ "Peromyscus leucopus",
                         str_detect(rodent_names, "P\\.? boylii") ~ "Peromyscus boylii",
                         str_detect(rodent_names, "P\\.? megalops") ~ "Peromyscus megalops",
                         str_detect(rodent_names, "P\\.? melanophrys") ~ "Peromyscus melanophrys",
                         str_detect(rodent_names, "P\\.? melanotis") ~ "Peromyscus melanotis",
                         str_detect(rodent_names, "P\\.? maniculatus|P\\.? manicultas|Peromyscus maniculatus bairdii|Permyscus maniculatus|Peromuscus maniculatus|	
Peromycus manicalatus|P maiculatus|P maniculatis") ~ "Peromyscus maniculatus",
                         str_detect(rodent_names, "Peromyscus banderanus|P\\.? mexicanus") ~ "Peromyscus mexicanus",
                         str_detect(rodent_names, "Peromyscus evides") ~ "Peromyscus aztecus",
                         str_detect(rodent_names, "Peromyscus eremicus") ~ "Peromycus eremicus",
                         str_detect(rodent_names, "P\\.? truei|Pinyon mouse") ~ "Peromyscus truei",
                         str_detect(rodent_names, "P\\.? californicus|Peromycus californicus|Peromyscus calofornicus|Peromysucs californicus") ~ "Peromyscus californicus",
                         str_detect(rodent_names, "Peromyscus leucopus noveboracensis") ~ "Peromyscus leucopus",
                         str_detect(rodent_names, "Peromyscus flavus") ~ "Perognathus flavus",
                         str_detect(rodent_names, "Praomys missonei") ~ "Praomys misonnei",
                         str_detect(rodent_names, "Phyllotic darwini") ~ "Phyllotis darwini",
                         str_detect(rodent_names, "R\\.? gracilis") ~ "Reithrodontomys gracilis",
                         str_detect(rodent_names, "R\\.? megalotis") ~ "Reithrodontomys megalotis",
                         str_detect(rodent_names, "R\\.? mexicanus") ~ "Reithrodontomys mexicanus",
                         str_detect(rodent_names, "R\\.? mantanus") ~ "Reithrodontomys montanus",
                         str_detect(rodent_names, "R\\.? sumichrasti") ~ "Reithrodontomys sumichrasti",
                         str_detect(rodent_names, "Rattus edwardsi") ~ "Leopoldamys edwardsi",
                         str_detect(rodent_names, "O\\.? rufus") ~ "Oxymycterus rufus",
                         str_detect(rodent_names, "R\\.? argentiventer|R argentiventer") ~ "Rattus argentiventer",
                         str_detect(rodent_names, "R\\.? exulans") ~ "Rattus exulans",
                         str_detect(rodent_names, "R\\.? losea|R losea") ~ "Rattus losea",
                         str_detect(rodent_names, "Rattus nitidius") ~ "Rattus nitidus",
                         str_detect(rodent_names, "R\\.? norvegicus|R vorvegicus") ~ "Rattus norvegicus",
                         str_detect(rodent_names, "R\\.? rattus|Rattus rattus alexandrinus|Rattus rattus frugivorus|Rattus rattus species complex|Rattus arboreus|Rattus rufescens|Rattus wroughtoni") ~ "Rattus rattus",
                         str_detect(rodent_names, "R\\.? tanezumi|Rattus tanezumi species complex|R\\.? flavipectus|R flavipectus") ~ "Rattus tanezumi",
                         str_detect(rodent_names, "R tiomanicus") ~ "Rattus tiomanicus",
                         str_detect(rodent_names, "Reihtrodon spp.") ~ "Reithrodon spp.",
                         str_detect(rodent_names, "Rattus yunnanensis") ~ "Hadromys yunnanensis",
                         str_detect(rodent_names, "Scapteromys acquaticus") ~ "Scapteromys aquaticus",
                         str_detect(rodent_names, "Scapteromys tumidis") ~ "Scapteromys tumidus",
                         str_detect(rodent_names, "S\\.? hispidus|S\\.? hispidis") ~ "Sigmodon hispidus",
                         str_detect(rodent_names, "S\\.? toltecus") ~ "Sigmodon toltecus",
                         str_detect(rodent_names, "U.u. Soricidae") ~ "Soricidae",
                         str_detect(rodent_names, "Sorex ananeus|S araneus") ~ "Sorex araneus",
                         str_detect(rodent_names, "S\\.? alpinus") ~ "Sorex alpinus",
                         str_detect(rodent_names, "Sorex gacillimus") ~ "Sorex gracillimus",
                         str_detect(rodent_names, "Serengetimys pernanus") ~ "Mastomys pernanus",
                         str_detect(rodent_names, "Spermophilus beecheyii") ~ "Spermophilus beecheyi",
                         str_detect(rodent_names, "Stenocephalemys sp. A") ~ "Stenocephalemys sp.",
                         str_detect(rodent_names, "Sylvilagus auduboni") ~ "Sylvilagus audubonii",
                         str_detect(rodent_names, "T\\.? dorsalis") ~ "Tamias dorsalis",
                         str_detect(rodent_names, "T\\.? minimus") ~ "Tamias minimus",
                         str_detect(rodent_names, "T\\.? indica") ~ "Tatera indica",
                         str_detect(rodent_names, "T\\.? nigrita") ~ "Thaptomys nigrita",
                         str_detect(rodent_names, "T\\.? elegans") ~ "Thylamys elegans",
                         str_detect(rodent_names, "Tamiascurus douglasii") ~ "Tamiasciurus douglasii",
                         str_detect(rodent_names, "Western harvest mouse") ~ "Reithrodontomys megalotis",
                         str_detect(rodent_names, "Citellus undulatus") ~ "Urocitellus undulatus",
                         str_detect(rodent_names, "Urotrichis talpoides") ~ "Urotrichus talpoides",
                         str_detect(rodent_names, "Oligorysomys|Oligorysomys spp.") ~ "Oligoryzomys",
                         str_detect(rodent_names, "Zygogontomy s") ~ "Zygodontomys",
                         str_detect(rodent_names, "Zygodontomys cherriei|Zygodontomys brevicauda cherriei") ~ "Zygodontomys brevicauda",
                         str_detect(rodent_names, "Other") ~ "Mammalia",
                         TRUE ~  rodent_names),
                       species_level = case_when(str_detect(cleaned_rodent_names, "sp\\.$|sp\\.?$|spp\\.?$|species") ~ FALSE, # if sp./spp. in name then not identified to species
                                                 lengths(str_split(cleaned_rodent_names, "\\s+")) != 2 ~ FALSE, # if no binomial also not identified to species
                                                 TRUE ~ TRUE)) %>%
  distinct()

# These ones didn't seem to work within the case_when
rodent_names$cleaned_rodent_names[rodent_names$rodent_names == "Mus (nannomys) mahomet"] <- "Mus mahomet"
rodent_names$cleaned_rodent_names[rodent_names$rodent_names == "Dendromus angolensis"] <- "Dendromus sp."
rodent_names$cleaned_rodent_names[rodent_names$rodent_names == "Gerbilliscus angolae"] <- "Gerbilliscus sp."
rodent_names$cleaned_rodent_names[rodent_names$rodent_names == "Mus Nannomys"] <- "Mus sp."
rodent_names$cleaned_rodent_names[rodent_names$rodent_names == "Gerbilliscus taborae"] <- "Gerbilliscus sp."
rodent_names$species_level[rodent_names$rodent_names == "Mus (nannomys) mahomet"] <- TRUE
rodent_names$species_level[rodent_names$rodent_names == "Dendromus angolensis"] <- FALSE
rodent_names$species_level[rodent_names$rodent_names == "Gerbilliscus angolae"] <- FALSE
rodent_names$species_level[rodent_names$rodent_names == "Mus Nannomys"] <- FALSE
rodent_names$species_level[rodent_names$rodent_names == "Gerbilliscus taborae"] <- FALSE


## Clean rodent genera -----------------------------------------------------

# Clean rodent genus names

rodent_genus <- tibble(rodent_genus = c(str_split(unique(rodent_names$cleaned_rodent_names), " ", simplify = TRUE)[, 1], str_to_sentence(combined_data_v2$host$genus)),
                       species_level = FALSE) %>%
  distinct() %>%
  mutate(cleaned_rodent_genus = case_when(str_detect(rodent_genus, "Aeothomys") ~ "Aethomys",
                                  str_detect(rodent_genus, "Ammonospermophilus") ~ "Ammospermophilus",
                                  str_detect(rodent_genus, "Apododemus|Apodemdus|Apodmus") ~ "Apodemus",
                                  str_detect(rodent_genus, "Arvincanthis") ~ "Arvicanthis",
                                  str_detect(rodent_genus, "Aterelix") ~ "Atelerix",
                                  str_detect(rodent_genus, "Baoimys|Biomys") ~ "Baiomys",
                                  str_detect(rodent_genus, "Bucepattersonius") ~ "Brucepattersonius",
                                  str_detect(rodent_genus, "Caomys") ~ "Calomys",
                                  str_detect(rodent_genus, "Cerradiomys") ~ "Cerradomys",
                                  str_detect(rodent_genus, "Chaetodipus") ~ "Chaeotdipus",
                                  str_detect(rodent_genus, "Cletherionomys") ~ "Clethrionomys",
                                  str_detect(rodent_genus, "Corcidura") ~ "Crocidura",
                                  str_detect(rodent_genus, "Lopuromys") ~ "Lophuromys",
                                  str_detect(rodent_genus, "Serengetimys") ~ "Mastomys",
                                  str_detect(rodent_genus, "Maxomus") ~ "Maxomys",
                                  str_detect(rodent_genus, "Monodelphys") ~ "Monodelphis",
                                  str_detect(rodent_genus, "Mictorus") ~ "Microtus",
                                  str_detect(rodent_genus, "Mydoes") ~ "Myodes",
                                  str_detect(rodent_genus, "Ochromyscus") ~ "Myomyscus",
                                  str_detect(rodent_genus, "Neotomas") ~ "Neotoma",
                                  str_detect(rodent_genus, "Oeonomys") ~ "Oenomys",
                                  str_detect(rodent_genus, "Oligorysomys|Oligoyzomys|Oligozomys|Oiigoryzomys|Oligorymys") ~ "Oligoryzomys",
                                  str_detect(rodent_genus, "Onchomys") ~ "Onychomys",
                                  str_detect(rodent_genus, "Peromycus") ~ "Peromyscus",
                                  str_detect(rodent_genus, "Pitimys") ~ "Pitymys",
                                  str_detect(rodent_genus, "Montemys") ~ "Praomys",
                                  str_detect(rodent_genus, "Scuirdae|Sciuridae") ~ "Sciuridae",
                                  str_detect(rodent_genus, "Sigmondon") ~ "Sigmodon",
                                  str_detect(rodent_genus, "Sorext") ~ "Sorex",
                                  str_detect(rodent_genus, "Suneus") ~ "Suncus",
                                  str_detect(rodent_genus, "Zygogontomys") ~ "Zygodontomys",
                                  TRUE ~ rodent_genus)) %>%
  distinct()


## Match to GBIF species taxa -------------------------------------------------

# Species level resolution
species_names <- read_rds(here("data", "raw_data", "gbif_species_names.rds"))
species_hierarchy <- read_rds(here("data", "raw_data", "gbif_species_hierarchy.rds"))
# Compare the species names previously resolved with names in the data
# If these do not match it will access GBIF and ask for those with multiple matches to match to the closest species
if(!length(unique(species_names$cleaned_rodent_names)) == length(unique(rodent_names$cleaned_rodent_names))) {
  
  gbif_species <- get_gbifid(unique(rodent_names$cleaned_rodent_names[rodent_names$species_level == TRUE]), rank = "species")
  resolve_species <- classification(gbif_species, db = "gbif", return_id = TRUE)
  species_hierarchy <- rbind(resolve_species) %>% 
    group_by(query) %>% 
    distinct() %>% 
    ungroup() %>% 
    pivot_wider(id_cols = "query", names_from = rank, values_from = name) %>% 
    left_join(rbind(resolve_species) %>%
                filter(rank == "species") %>%
                select(query, species_id = id),
              by = "query") %>%
    arrange(kingdom, phylum, class, order, family, genus, species, species_id) %>%
    drop_na(species)
  species_names <- tibble(cleaned_rodent_names = unique(rodent_names$cleaned_rodent_names[rodent_names$species_level == TRUE]),
                          query = as.character(gbif_species)) %>%
    # match to raw data names
    right_join(rodent_names,
               by = "cleaned_rodent_names") %>%
    left_join(species_hierarchy %>%
                distinct(species, query, genus, family, order, class, species_id), by = c("query")) %>%
    rename("resolved_name" = species,
           "gbif_id" = species_id)
  write_rds(species_names, here("data", "raw_data", "gbif_species_names.rds"))
  write_rds(species_hierarchy, here("data", "raw_data", "gbif_species_hierarchy.rds"))
  
}

# Some species names remain unmatched
# Need to think about the processes to incorporate or rename these
unmatched_species_names <- species_names %>%
  filter(species_level == TRUE) %>%
  filter(is.na(gbif_id)) %>%
  pull(rodent_names)
#gnr_species <- gnr_resolve(unmatched_species_names)


## Join GBIF species taxa to raw data --------------------------------------

# Join these classifications to the raw data where species level data is available
# v2

v2_host_species <- combined_data_v2$host %>%
  filter(!is.na(scientificName)) %>%
  filter(taxonRank == "species") %>%
  mutate(rodent_names = str_to_sentence(str_squish(scientificName))) %>%
  select(-genus) %>%
  left_join(species_names %>%
              select(-cleaned_rodent_names),
            by = c("rodent_names")) %>%
  drop_na(gbif_id)

# v3
 
v3_host_species <- combined_data_v3$host %>%
  mutate(rodent_names = str_to_sentence(str_squish(scientificName))) %>%
  left_join(rodent_names, by = "rodent_names") %>%
  left_join(species_names %>%
              select(-cleaned_rodent_names),
            by = c("rodent_names")) %>%
  drop_na(gbif_id)

## Match to GBIF genus taxa ------------------------------------------------

# Genus level resolution
genus_names <- read_rds(here("data", "raw_data", "gbif_genus_names.rds"))
genus_hierarchy <- read_rds(here("data", "raw_data", "gbif_genus_hierarchy.rds"))
# Compare the genus names previously resolved with names in the data
# If these do not match it will access GBIF and ask for those with multiple matches to match to the closest genus
# Remove common names before matching
if(!length(unique(genus_names$cleaned_rodent_genus)) == length(unique(rodent_genus$cleaned_rodent_genus))) {
  
  genera_to_match <- rodent_genus %>%
    filter(!str_detect(rodent_genus, "\"Social")) %>%
    distinct(cleaned_rodent_genus) %>%
    arrange() %>%
    pull()
  
  gbif_genus <- get_gbifid(genera_to_match, rank = "genus")
  resolve_genus <- classification(gbif_genus, db = "gbif")
  genus_hierarchy <- rbind(resolve_genus) %>%
    group_by(query) %>%
    distinct() %>% 
    ungroup() %>% 
    pivot_wider(id_cols = "query", names_from = rank, values_from = name) %>% 
    left_join(rbind(resolve_genus) %>% 
                group_by(query) %>% 
                distinct() %>%
                summarise(gbif_id = id[which.max(rank == "genus")]), by = "query") %>% 
    arrange(kingdom, phylum, class, order, family, genus) %>%
    drop_na(genus)
  
  genus_names <- tibble(cleaned_rodent_genus = genera_to_match,
                        query = as.character(gbif_genus)) %>%
    right_join(rodent_genus,
               by = "cleaned_rodent_genus") %>%
    left_join(genus_hierarchy %>%
                distinct(genus, query, family, order, class, gbif_id), by = c("query")) %>%
    distinct(rodent_genus, genus, .keep_all = TRUE)
  
  write_rds(genus_names, here("data", "raw_data", "gbif_genus_names.rds"))
  write_rds(genus_hierarchy, here("data", "raw_data", "gbif_genus_hierarchy.rds"))
  
}

# v2

v2_host_genus <- combined_data_v2$host %>%
  filter(!rodent_record_id %in% v2_host_species$rodent_record_id) %>%
  mutate(rodent_genus = str_to_sentence(str_squish(genus))) %>%
  select(-genus) %>%
  left_join(genus_names %>%
              select(-cleaned_rodent_genus) %>%
              distinct(),
            by = "rodent_genus") %>%
  drop_na(gbif_id)

# v3

v3_host_genus <- combined_data_v3$host %>%
  filter(!rodent_record_id %in% v3_host_species$rodent_record_id) %>%
  mutate(scientificName = str_to_sentence(str_squish(scientificName))) %>%
  left_join(rodent_names, by = c("scientificName" = "rodent_names")) %>%
  mutate(cleaned_rodent_names = str_trim(str_remove_all(cleaned_rodent_names, "sp\\.$|sp\\.?$|spp\\.?$|species"))) %>%
  left_join(genus_names %>%
              select(-rodent_genus) %>%
              distinct(),
            by = c("cleaned_rodent_names" = "cleaned_rodent_genus")) %>%
  drop_na(gbif_id)


# Higher taxonomy ---------------------------------------------------------
# Some records are only associated to rodentia, sciuridae and mammalia

unmatched_genera <- combined_data_v3$host %>%
  filter(!rodent_record_id %in% v3_host_species$rodent_record_id) %>%
  filter(!rodent_record_id %in% v3_host_genus$rodent_record_id) %>%
  select(-any_of(c("cleaned_rodent_names", "species_level.x", "query", "species_level.y", "genus", "family", "order", "class", "gbif_id")))

higher_taxa <- unmatched_genera %>%
  distinct(scientificName) %>%
  mutate(clean_higher_taxa = case_when(str_detect(scientificName, "Rodentia|Rodentia sp.|\"Social Rat\"|Unknown") ~ "Rodentia",
                                       str_detect(scientificName, "Other") ~ "Mammalia",
                                       str_detect(scientificName, "Neotominae") ~ "Cricetidae",
                                       str_detect(scientificName, "Sciuridae|Scuirdae") ~ "Sciuridae",
                                       str_detect(scientificName, "Soricidae") ~ "Soricidae",
                                       TRUE ~ scientificName))

gbif_higher <- get_gbifid(unique(higher_taxa$clean_higher_taxa))
resolve_higher <- classification(gbif_higher, db = "gbif")
higher_hierarchy <- rbind(resolve_higher) %>% 
  group_by(query) %>% 
  distinct() %>% 
  ungroup() %>% 
  pivot_wider(id_cols = "query", names_from = rank, values_from = name) %>%
  mutate(gbif_id = as.integer(query))
higher_taxa_table <- tibble(clean_higher_taxa = unique(higher_taxa$clean_higher_taxa),
                            query = as.character(gbif_higher)) %>%
  left_join(higher_hierarchy, by = "query") %>%
  left_join(higher_taxa, by = "clean_higher_taxa")

v3_unmatched_higher <- unmatched_genera %>%
  left_join(higher_taxa_table) %>%
  rename("cleaned_rodent_names" = "clean_higher_taxa")

# Combine species and genus level matching --------------------------------

# v2

if(nrow(v2_host_species) + nrow(v2_host_genus) == v2_raw_length) {
  
  v2_species_matched <- v2_host_species
  
  v2_genus_matched <- v2_host_genus
  
  v2_rodent_cleaned <- bind_rows(v2_species_matched,
                                 v2_genus_matched) %>%
    arrange(rodent_record_id) %>%
    rename("species" = "resolved_name",
           "extracted_name" = scientificName) %>%
    mutate(taxa_resolution = fct(case_when(!is.na(species) ~ "species",
                                           !is.na(genus) ~ "genus",
                                           !is.na(family) ~ "family",
                                           !is.na(order) ~ "order",
                                           !is.na(class) ~ "class",
                                           TRUE ~ "unresolved"), # Default for rows with no taxonomic information
                                 levels = c("species", "genus", "family", "order", "class")),
           individualCount = if_else(is.na(individualCount), 0, individualCount),
           scientificName = coalesce(species, genus, family, order, class)) %>%
    select(-any_of(c("row_n", "rodent_genus", "taxonRank", "rodent_names"))) %>%
    select(rodent_record_id, study_id, eventDate, extracted_name, scientificName, locality, country, verbatimLocality, coordinate_resolution, decimalLatitude, decimalLongitude, individualCount,
           trapEffort, trapEffortResolution,
           query, taxa_resolution, species, genus, family, order, class, gbif_id)
  
  message("Species and genus matching appears complete (v2).")
  
} else {
  
  message("Species and genus matching leads to a change in the amount of raw data (v2). This will need manual checking")
  
}

# v3

if(nrow(v3_host_species) + nrow(v3_host_genus) + nrow(v3_unmatched_higher) == v3_raw_length) {
  
  v3_species_matched <- v3_host_species %>%
    mutate(row_n = as.numeric(str_extract(rodent_record_id, "\\d+")),
           study_id_n = as.numeric(str_extract(study_id, "\\d+")))
  
  v3_genus_matched <- v3_host_genus %>%
    mutate(row_n = as.numeric(str_extract(rodent_record_id, "\\d+")),
           study_id_n = as.numeric(str_extract(study_id, "\\d+")))
  
  v3_higher_taxa_matched <- v3_unmatched_higher %>%
    mutate(row_n = as.numeric(str_extract(rodent_record_id, "\\d+")),
           study_id_n = as.numeric(str_extract(study_id, "\\d+")))
  
  v3_rodent_cleaned <- bind_rows(v3_species_matched,
                                 v3_genus_matched,
                                 v3_higher_taxa_matched) %>%
    arrange(rodent_record_id) %>%
    rename("species" = "resolved_name",
           "extracted_name" = scientificName) %>%
    mutate(taxa_resolution = fct(case_when(!is.na(species) ~ "species",
                                           !is.na(genus) ~ "genus",
                                           !is.na(family) ~ "family",
                                           !is.na(order) ~ "order",
                                           !is.na(class) ~ "class",
                                           TRUE ~ "unresolved"), # Default for rows with no taxonomic information
                                 levels = c("species", "genus", "family", "order", "class", "unresolved")),
           individualCount = if_else(is.na(individualCount), 0, individualCount),
           scientificName = coalesce(species, genus, family, order, class)) %>%
    select(rodent_record_id, study_id, eventDate, extracted_name, scientificName, locality, country, verbatimLocality, coordinate_resolution, decimalLatitude, decimalLongitude, individualCount,
           trapEffort, trapEffortResolution,
           query, taxa_resolution, species, genus, family, order, class, gbif_id)
  
  message("Species and genus matching appears complete (v3).")
  
} else {
  
  message("Species and genus matching leads to a change in the amount of raw data (v3). This will need manual checking")
  
}


# Clean coordinates -------------------------------------------------------
v2_missing_coords <- v2_rodent_cleaned %>%
  filter(is.na(decimalLatitude) | is.na(decimalLongitude))

v3_missing_coords <- v3_rodent_cleaned %>%
  filter(is.na(decimalLatitude) | is.na(decimalLongitude))

# Check coordinates fall within the -180, 180
v2_illogical_coords <- v2_rodent_cleaned %>%
  filter(!rodent_record_id %in% v2_missing_coords$rodent_record_id) %>%
  mutate(logical_coords = decimalLatitude >= -180 & decimalLatitude <= 180 &
           decimalLongitude >= -180 & decimalLongitude <= 180) %>%
  filter(logical_coords != TRUE)

nrow(v2_illogical_coords)

v3_illogical_coords <- v3_rodent_cleaned %>%
  filter(!rodent_record_id %in% v2_missing_coords$rodent_record_id) %>%
  mutate(logical_coords = decimalLatitude >= -180 & decimalLatitude <= 180 &
           decimalLongitude >= -180 & decimalLongitude <= 180) %>%
  filter(logical_coords != TRUE)

nrow(v3_illogical_coords)

# Check and clean countrynames
# v2

v2_country_cleaned <- v2_rodent_cleaned %>%
  filter(!rodent_record_id %in% v2_missing_coords$rodent_record_id) %>%
  distinct(rodent_record_id, study_id, country, decimalLatitude, decimalLongitude) %>%
  mutate(country = case_when(str_detect(country, "England|Scotland|Wales") ~ "United Kingdom",
                             TRUE ~ country),
         iso3 = if_else(str_detect(country, ", "), 
                        NA_character_, 
                        countrycode(country, "country.name", "iso3c")))

v2_single_countries <- v2_country_cleaned %>%
  distinct(study_id, country, decimalLatitude, decimalLongitude, iso3, .keep_all = TRUE) %>%
  drop_na() %>%
  vect(geom = c("decimalLongitude", "decimalLatitude"), crs = "EPSG:4326")

v2_single_countries_mapping <- terra::extract(world_shapefile, v2_single_countries)

v2_single_countries_status <- v2_single_countries %>%
  cbind(v2_single_countries_mapping) %>%
  mutate(coords_in_country = if_else(iso3 == GID_0, TRUE, FALSE)) %>%
  as_tibble(geom = "XY") %>%
  rename("decimalLongitude" = x,
         "decimalLatitude" = y)

v2_single_countries_flagged <- v2_single_countries_status %>%
  filter(coords_in_country == FALSE)

v2_multiple_countries <- v2_country_cleaned %>%
  filter(is.na(iso3) & str_detect(country, ", ")) %>%
  mutate(iso3_composite = str_split(country, ", ") %>% 
           map(~ countrycode(.x, "country.name", "iso3c")) %>% 
           map_chr(~ paste0(.x, collapse = ", ")))

check_multiple_countries <- function(lat, lon, iso3_list, world_shapefile) {
  # If iso3_list is NA or empty, use coordinates to determine the country
  if (is.null(iso3_list) || all(is.na(iso3_list)) || all(iso3_list == "NA") || length(iso3_list) == 0) {
    message("iso3_list is NA or empty, checking by coordinates")
    
    point <- vect(data.frame(lon = lon, lat = lat), geom = c("lon", "lat"), crs = "EPSG:4326")
    extracted <- terra::extract(world_shapefile, point)
    message("Extracted country: ", extracted$GID_0)
    # Return both a logical for inside-country check and the extracted iso3 code
    return(list(coords_in_country = !is.na(extracted$GID_0), iso3_code = extracted$GID_0))
  }
  
  # Subset world_shapefile for the composite countries
  composite_boundary <- world_shapefile[world_shapefile$GID_0 %in% iso3_list, ]
  if (nrow(composite_boundary) == 0) return(list(coords_in_country = NA, iso3_code = NA)) # If no boundary, return NA for both
  
  # Create a single composite polygon
  composite_union <- aggregate(composite_boundary, dissolve = TRUE)
  
  # Create a point object
  point <- vect(data.frame(lon = lon, lat = lat), geom = c("lon", "lat"), crs = "EPSG:4326")
  # Check if the point falls within the composite boundary
  in_country <- relate(point, composite_union, "intersects")
  
  return(list(coords_in_country = as.logical(in_country), iso3_code = NA))
}

v2_multiple_countries <- v2_multiple_countries %>%
  rowwise() %>%
  mutate(country_check_results = list(check_multiple_countries(
    decimalLatitude, 
    decimalLongitude, 
    str_split(iso3_composite, ", ") %>% unlist(), 
    world_shapefile)
    )) %>%
  ungroup()  %>%
  mutate(
    coords_in_country = map_lgl(country_check_results, ~ .x$coords_in_country),  # Extract logical values
    iso3_code = map_chr(country_check_results, ~ .x$iso3_code)                     # Extract country codes
  ) %>%
  select(-country_check_results)

v2_multiple_countries_flagged <- v2_multiple_countries %>%
  filter(coords_in_country == FALSE)

v2_coords_requiring_check <- bind_rows(v2_single_countries_flagged,
                                       v2_multiple_countries_flagged)

v2_coords_checked <- bind_rows(v2_single_countries_status,
                               v2_multiple_countries) %>%
  select(rodent_record_id, study_id, country, iso3, iso3_composite, NAME_0, GID_0, country_flag = coords_in_country) %>%
  right_join(v2_rodent_cleaned, by = c("rodent_record_id", "study_id", "country")) %>%
  arrange(rodent_record_id) %>%
  group_by(study_id, country, decimalLatitude, decimalLongitude) %>%
  fill(iso3, .direction = "down") %>%
  fill(NAME_0, .direction = "down") %>%
  fill(GID_0, .direction = "down") %>%
  fill(country_flag, .direction = "down")

# Based on the coordinate resolution we will match locations to administrative areas

v2_coords_checked %>%
  select(rodent_record_id, study_id, country, locality, coordinate_resolution, decimalLatitude, decimalLongitude) %>%
  distinct(country, locality, coordinate_resolution, decimalLatitude, decimalLongitude, .keep_all = TRUE)

# v3

v3_country_cleaned <- v3_rodent_cleaned %>%
  filter(!rodent_record_id %in% v3_missing_coords$rodent_record_id) %>%
  distinct(rodent_record_id, study_id, country, decimalLatitude, decimalLongitude) %>%
  mutate(country = case_when(str_detect(country, "England|Scotland|Wales") ~ "United Kingdom",
                             str_detect(country, "California U.S.") ~ "USA",
                             str_detect(country, "Finald") ~ "Finland",
                             str_detect(country, "Phillipines") ~ "Philippines",
                             str_detect(country, "Solvenia") ~ "Slovenia",
                             str_detect(country, "Uruaguay") ~ "Uruguay",
                             str_detect(country, "Venezuala") ~ "Venezuela",
                             TRUE ~ country),
         iso3 = if_else(str_detect(country, ", |Multiple|Czechoslovakia"), 
                        NA_character_, 
                        countrycode(country, "country.name", "iso3c")))

v3_single_countries <- v3_country_cleaned  %>%
  distinct(study_id, country, decimalLatitude, decimalLongitude, iso3, .keep_all = TRUE) %>%
  drop_na() %>%
  vect(geom = c("decimalLongitude", "decimalLatitude"), crs = "EPSG:4326")

v3_single_countries_mapping <- terra::extract(world_shapefile, v3_single_countries)

v3_single_countries_status <- v3_single_countries %>%
  cbind(v3_single_countries_mapping) %>%
  mutate(coords_in_country = if_else(iso3 == GID_0, TRUE, FALSE)) %>%
  as_tibble(geom = "XY") %>%
  rename("decimalLongitude" = x,
         "decimalLatitude" = y)

v3_single_countries_flagged <- v3_single_countries_status %>%
  filter(coords_in_country == FALSE)

v3_multiple_countries <- v3_country_cleaned %>%
  filter(is.na(iso3) | str_detect(country, ", ")) %>%
  mutate(iso3_composite = str_split(country, ", ") %>% 
           map(~ countrycode(.x, "country.name", "iso3c")) %>% 
           map_chr(~ paste0(.x, collapse = ", ")))

v3_multiple_countries_status <- v3_multiple_countries %>%
  rowwise() %>%
  mutate(country_check_results = list(check_multiple_countries(
    decimalLatitude, 
    decimalLongitude, 
    str_split(iso3_composite, ", ") %>% unlist(), 
    world_shapefile)
  )) %>%
  ungroup() %>%
  mutate(
    coords_in_country = map_lgl(country_check_results, ~ .x$coords_in_country),  # Extract logical values
    iso3_code = map_chr(country_check_results, ~ .x$iso3_code)                     # Extract country codes
  ) %>%
  select(-country_check_results)

v3_multiple_countries_flagged <- v3_multiple_countries_status %>%
  filter(coords_in_country == FALSE)

v3_coords_requiring_check <- bind_rows(v3_single_countries_flagged,
                                       v3_multiple_countries_flagged)


v3_coords_checked <- bind_rows(v3_single_countries_status,
                               v3_multiple_countries) %>%
  select(rodent_record_id, study_id, country, iso3, iso3_composite, NAME_0, GID_0, country_flag = coords_in_country) %>%
  right_join(v3_rodent_cleaned, by = c("rodent_record_id", "study_id", "country")) %>%
  arrange(rodent_record_id) %>%
  group_by(study_id, country, decimalLatitude, decimalLongitude) %>%
  fill(iso3, .direction = "down") %>%
  fill(NAME_0, .direction = "down") %>%
  fill(GID_0, .direction = "down") %>%
  fill(country_flag, .direction = "down")

# Add location hierarchy --------------------------------------------------

# Clean coordinate resolution level to site (locations of traps, trap grids, trap lines etc), village (), town, city, adm3 (lowest level of administrative unit the country), adm2 (level 2 administrative unit), adm1 (level 1 administrative unit)

# Some locations may span multiple regions at the same level need to be used with caution (maybe set up a flag for these)

## Join cleaned coordinate data -------------------------------------------

clean_host <- bind_rows(v2_coords_checked,
                        v3_coords_checked)

# Impute non-detections ---------------------------------------------------
# Leave to the end
# impute non-detected species for studies with summarised detections
# we will identify those studies that report multiple sites and are not explicitly labelled as individual level data
# this will only be performed for v3

studies_to_impute <- combined_data_v3$studies %>%
  filter(!str_detect(data_access, "individual")) %>%
  arrange(study_id) %>%
  pull(study_id)

# Finalise clean rodent data ----------------------------------------------

# Clean Pathogen sheet ----------------------------------------------------
virus_mapping <- read_csv(here("data", "matching", "virus_names_matching.csv"))

virus_matched_long <- virus_mapping %>%
  pivot_longer(cols = starts_with("virus_"),
               names_to = "virus_number",
               values_to = "matched_virus",
               values_drop_na = TRUE) %>%
  select(-virus_number) %>%
  distinct(virus, clean_name, matched_virus, taxonomic_level)   %>%
  group_by(virus) %>%
  summarise(
    clean_name = str_c(unique(clean_name), collapse = ", "),
    matched_virus = str_c(unique(matched_virus), collapse = ", "),
    taxonomic_level = str_c(unique(taxonomic_level), collapse = ", "),
    .groups = "drop"
  ) %>%
  rename("virus_clean" = clean_name)

# To-do generate taxonomy from looking up on NCBI or equivalent for pathogens

assay_clean <- tibble(assay = unique(c(unique(combined_data_v2$pathogen$identificationRemarks), unique(combined_data_v3$pathogen$assay)))) %>%
  mutate(assay_clean = case_when(str_detect(tolower(assay), "elisa|antibody|antigen|ifa|immuno|serology|frnt|western|eia|rapid field test|igg|complement|elsia|ib|hdp") ~ "Serology",
                                 str_detect(tolower(assay), "pcr|rt-pcr") ~ "PCR",
                                 str_detect(tolower(assay), "culture") ~ "Culture",
                                 str_detect(tolower(assay), "sequencing|illumina") ~ "Sequencing",))
# Check pathogen IDs are unique
pathogen_uid <- tibble(source = c(rep("v2", times = nrow(combined_data_v2$pathogen)), rep("v3", times = nrow(combined_data_v3$pathogen))),
                     row_number = c(seq(from = 1, to = nrow(combined_data_v2$pathogen)), seq(from = 1, to = nrow(combined_data_v3$pathogen))),
                     pathogen_id = c(combined_data_v2$pathogen$pathogen_record_id, combined_data_v3$pathogen$pathogen_record_id)) %>%
  group_by(pathogen_id) %>%
  mutate(n = n())

table(id_count = pathogen_uid$n)
# manually recode but has been updated in grant_v3
levels(combined_data_v3$pathogen$pathogen_record_id) <- c(levels(combined_data_v3$pathogen$pathogen_record_id), "grant_3038")
combined_data_v3$pathogen$pathogen_record_id[combined_data_v3$pathogen$pathogen_record_id == "grant_2409"] <- c("grant_2409", "grant_3038")

# v2 pathogen processing --------------------------------------------------

clean_v2_path <- combined_data_v2$pathogen %>%
  left_join(virus_matched_long, by = c("scientificName" = "virus")) %>%
  mutate(virus_clean = case_when(is.na(scientificName) & str_detect(tolower(family), "hanta") ~ "Orthohantavirus",
                                 is.na(scientificName) & str_detect(tolower(family), "arena") ~ "Mammarenavirus",
                                 TRUE ~ virus_clean),
         taxonomic_level = case_when(is.na(scientificName) & str_detect(tolower(family), "hanta|arena") ~ "genus",
                                     TRUE ~ taxonomic_level)) %>%
  rename("original_name" = scientificName) %>%
  mutate(original_name = coalesce(original_name, family)) %>%
  left_join(assay_clean, by = c("identificationRemarks" = "assay")) %>%
  rename("original_assay" = identificationRemarks) %>%
  mutate(n_assayed = as.numeric(occurrenceRemarks),
         n_negative = as.numeric(number_negative),
         n_positive = as.numeric(organismQuantity),
         n_inconclusive = as.numeric(number_inconclusive)) %>%
  mutate(n_negative = replace_na(n_negative, 0),
         n_positive = replace_na(n_positive, 0),
         n_inconclusive = replace_na(n_inconclusive, 0)) %>%
  dplyr::select(pathogen_record_id, associated_rodent_record_id, study_id, family, virus_clean, matched_virus, taxonomic_level, assay_clean,
         n_assayed, n_positive, n_negative, n_inconclusive, original_name, original_assay, note) %>%
  left_join(clean_host %>%
              dplyr::select(rodent_record_id, study_id, eventDate, host_species = species, host_genus = genus, host_family = family, host_order = order,
                            locality, country, verbatimLocality, coordinate_resolution, decimalLatitude, decimalLongitude),
            by = c("associated_rodent_record_id" = "rodent_record_id", "study_id"))

# Some records are missing an associated rodent record id
orphaned_pathogen_v2 <- clean_v2_path %>%
  filter(is.na(associated_rodent_record_id))

# v3 pathogen processing --------------------------------------------------

clean_v3_path <- combined_data_v3$pathogen %>%
  left_join(virus_matched_long, by = c("scientificName" = "virus")) %>%
  mutate(virus_clean = case_when(is.na(virus_clean) & str_detect(tolower(scientificName), "hanta") ~ "Orthohantavirus",
                                 is.na(virus_clean) & str_detect(tolower(scientificName), "arena") ~ "Mammarenavirus",
                                 is.na(virus_clean) & str_detect(tolower(family), "hanta") ~ "Orthohantavirus",
                                 is.na(virus_clean) & str_detect(tolower(family), "arena") ~ "Mammarenavirus",
                                 TRUE ~ virus_clean),
         taxonomic_level = case_when(is.na(taxonomic_level) & str_detect(tolower(family), "hanta|arena") ~ "genus",
                                     TRUE ~ taxonomic_level)) %>%
  rename("original_name" = scientificName) %>%
  mutate(original_name = coalesce(original_name, family)) %>%
  left_join(assay_clean, by = "assay") %>%
  mutate(n_assayed = as.numeric(tested),
         n_negative = as.numeric(negative),
         n_positive = as.numeric(positive),
         n_inconclusive = as.numeric(number_inconclusive)) %>%
  mutate(n_negative = replace_na(n_negative, 0),
         n_positive = replace_na(n_positive, 0),
         n_inconclusive = replace_na(n_inconclusive, 0)) %>%
  dplyr::select(pathogen_record_id, associated_rodent_record_id, study_id, family, virus_clean, matched_virus, taxonomic_level, assay_clean,
         n_assayed, n_positive, n_negative, n_inconclusive, original_name, note) %>%
  left_join(clean_host %>%
              dplyr::select(rodent_record_id, study_id, eventDate, host_species = species, host_genus = genus, host_family = family, host_order = order,
                     locality, country, verbatimLocality, coordinate_resolution, decimalLatitude, decimalLongitude),
            by = c("associated_rodent_record_id" = "rodent_record_id", "study_id"))

orphaned_pathogen_v3 <- clean_v3_path %>%
  filter(is.na(associated_rodent_record_id))

orphaned_pathogen <- bind_rows(orphaned_pathogen_v2, orphaned_pathogen_v3)

clean_pathogen <- clean_v2_path %>%
  mutate(version = "v2") %>%
  bind_rows(clean_v3_path %>%
              mutate(version = "v3"))

# Clean sequence_data sheet -----------------------------------------------

pathogen_sequences_v2 <- combined_data_v2$sequence_data %>%
  filter(str_detect(sequenceType, "Pathogen")) %>%
  left_join(virus_matched_long, by = c("scientificName" = "virus")) %>%
  dplyr::select(sequence_record_id, associated_pathogen_record_id, associated_rodent_record_id, study_id, sequenceType, virus_clean, matched_virus, original_name = scientificName, accession_number) %>%
  left_join(clean_host %>%
              dplyr::select(rodent_record_id, study_id, eventDate, species, host_genus = genus, gbif_id, coordinate_resolution, decimalLatitude, decimalLongitude),
            by = c("associated_rodent_record_id" = "rodent_record_id", "study_id")) %>%
  left_join(clean_pathogen %>%
              dplyr::select(pathogen_record_id, study_id, n_assayed, n_positive, n_negative),
            by = c("associated_pathogen_record_id" = "pathogen_record_id", "study_id"))

host_sequences_v2 <-  combined_data_v2$sequence_data %>%
  filter(str_detect(sequenceType, "Host")) %>%
  dplyr::select(sequence_record_id, associated_pathogen_record_id, associated_rodent_record_id, study_id, sequenceType, accession_number) %>%
  left_join(clean_host %>%
              dplyr::select(rodent_record_id, study_id, eventDate, species, host_genus = genus, gbif_id, coordinate_resolution, decimalLatitude, decimalLongitude),
            by = c("associated_rodent_record_id" = "rodent_record_id", "study_id"))

clean_sequences_v2 <- bind_rows(pathogen_sequences_v2,
                                host_sequences_v2) %>%
  mutate(version = "v2")

pathogen_sequences_v3 <- combined_data_v3$sequence_data %>%
  filter(str_detect(sequenceType, "Pathogen")) %>%
  left_join(virus_matched_long, by = c("scientificName" = "virus")) %>%
  dplyr::select(sequence_record_id, associated_pathogen_record_id, associated_rodent_record_id, study_id, sequenceType, virus_clean, matched_virus, original_name = scientificName, accession_number) %>%
  left_join(clean_host %>%
              dplyr::select(rodent_record_id, study_id, eventDate, species, host_genus = genus, gbif_id, coordinate_resolution, decimalLatitude, decimalLongitude),
            by = c("associated_rodent_record_id" = "rodent_record_id", "study_id")) %>%
  left_join(clean_pathogen %>%
              dplyr::select(pathogen_record_id, study_id, n_assayed, n_positive, n_negative),
            by = c("associated_pathogen_record_id" = "pathogen_record_id", "study_id"))

host_sequences_v3 <-  combined_data_v3$sequence_data %>%
  filter(str_detect(sequenceType, "Host")) %>%
  dplyr::select(sequence_record_id, associated_pathogen_record_id, associated_rodent_record_id, study_id, sequenceType, accession_number) %>%
  left_join(clean_host %>%
              dplyr::select(rodent_record_id, study_id, eventDate, species, host_genus = genus, gbif_id, coordinate_resolution, decimalLatitude, decimalLongitude),
            by = c("associated_rodent_record_id" = "rodent_record_id", "study_id"))

clean_sequences_v3 <- bind_rows(pathogen_sequences_v3,
                                host_sequences_v3) %>%
  mutate(version = "v3")

clean_sequences <- bind_rows(clean_sequences_v2,
                             clean_sequences_v3)

# Several Accessions are entered more than once
table(clean_sequences %>%
        group_by(accession_number) %>%
        summarise(n = n()) %>%
        pull(n))
# Add extraction descriptives ---------------------------------------------
group_descriptives <- c("study_id", "full_text_id", "identifiedBy", "datasetName", "publication_year", "data",
                        "sampling_effort", "data_access", "linked_manuscripts", "data_extractor", "data_checker",
                        "notes", "source", "data_resolution")

combined_descriptive <- bind_rows(descriptives_v2 %>%
                                    mutate(source = "v2"),
                                  descriptives_v3 %>%
                                    mutate(source = "v3")) %>%
  mutate(full_text_id = fct(full_text_id, levels = levels(clean_citations$full_text_id))) %>%
  left_join(clean_host %>%
              ungroup() %>%
              select(rodent_record_id, study_id, species, individualCount), by = "study_id") %>%
  group_by(across(any_of(group_descriptives))) %>%
  summarise(n_rodent_records = n(),
            n_rodent_species = length(unique(species)),
            n_rodent_individuals = sum(individualCount, na.rm = TRUE),
            .groups = "keep") %>%
  left_join(clean_pathogen %>%
              ungroup() %>%
              select(pathogen_record_id, study_id, virus_clean, n_assayed, n_positive), by = "study_id") %>%
  group_by(n_rodent_records, n_rodent_species, n_rodent_individuals, .add = TRUE) %>%
  summarise(n_pathogen_records = n(),
            n_pathogen_species = length(unique(virus_clean)),
            n_assays = sum(n_assayed, na.rm = TRUE),
            n_positive = sum(n_positive, na.rm = TRUE),
            .groups = "keep") %>%
  left_join(clean_sequences %>%
              select(sequence_record_id, study_id, accession_number, sequenceType)) %>%
  group_by(n_pathogen_records, n_pathogen_species, n_assays, .add = TRUE) %>%
  summarise(n_sequence_records = sum(!is.na(accession_number)),
            n_pathogen_sequences = sum(sequenceType == "Pathogen", na.rm = TRUE),
            n_host_sequences = sum(sequenceType == "Host", na.rm = TRUE),
            .groups = "keep")

# Clean data --------------------------------------------------------------
combined_data <- list(citations = clean_citations,
                      descriptive = combined_descriptive,
                      host = clean_host,
                      pathogen = clean_pathogen,
                      sequences = clean_sequences)

write_rds(combined_data, here("data", "clean_data", paste0(analysis_date, "_data.rds")))
