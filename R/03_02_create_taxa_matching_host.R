# Project ArHa: Create Manual Taxonomic Matching Tables
# 03_02_create_taxa_matching_host.R
# Purpose: This script loads the raw rodent data, identifies all unique
# scientific names and genera, and applies a set of manual cleaning rules
# to create standardized lookup tables.
# 

combined_data_v2 <- read_rds(here::here("data", "raw_data", paste0(analysis_date, "_v2_data.rds")))
combined_data_v3 <- read_rds(here::here("data", "raw_data", paste0(analysis_date, "_v3_data.rds")))


# Clean Rodent Species ----------------------------------------------------

rodent_names <- tibble(rodent_names = str_to_sentence(str_squish(sort(c(unique(combined_data_v2$rodent$scientificName), unique(combined_data_v3$rodent$scientificName))))),
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
                         str_detect(rodent_names, "Neotoma alstoni") ~ "Neotomadon alstoni",
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


# Clean Rodent Genera -----------------------------------------------------

rodent_genus <- tibble(rodent_genus = c(str_split(unique(rodent_names$cleaned_rodent_names), " ", simplify = TRUE)[, 1], str_to_sentence(unique(combined_data_v2$rodent$genus))),
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


# Write Dictionaries ------------------------------------------------------

write_rds(rodent_names, here("data", "matching", "rodent_names_manual.rds"))
write_rds(rodent_genus, here("data", "matching", "rodent_genus_manual.rds"))


# Match to GBIF Species Taxa ----------------------------------------------

# Check if the GBIF dictionary needs to be updated
species_names_file <- here("data", "matching", "gbif_species_names.rds")
if (file.exists(species_names_file)) {
  gbif_species_names_old <- read_rds(species_names_file)
  if (length(unique(gbif_species_names_old$cleaned_rodent_names)) == length(unique(rodent_names$cleaned_rodent_names))) {
    message("GBIF species dictionary is up-to-date. Skipping API call.")
    species_names <- gbif_species_names_old
  } else {
    message("New species found. Updating GBIF species dictionary.")
    # Run the GBIF API calls
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
    write_rds(species_names, species_names_file)
    write_rds(species_hierarchy, here("data", "matching", "gbif_species_hierarchy.rds"))
  } 
  } else {
  message("GBIF species dictionary not found. Creating new dictionary.")
  # Run the GBIF API calls
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
  
  # Save the new files
  write_rds(species_names, species_names_file)
  write_rds(species_hierarchy, here("data", "matching", "gbif_species_hierarchy.rds"))
  }


# Some species names remain unmatched
# Need to think about the processes to incorporate or rename these
unmatched_species_names <- species_names %>%
  filter(species_level == TRUE) %>%
  filter(is.na(gbif_id)) %>%
  pull(rodent_names)

# Match to GBIF Genus Taxa ------------------------------------------------

genus_names_file <- here("data", "matching", "gbif_genus_names.rds")
genus_hierarchy_file <- here("data", "matching", "gbif_genus_hierarchy.rds")

if (file.exists(genus_names_file) && file.exists(genus_hierarchy_file)) {
  gbif_genus_names_old <- read_rds(genus_names_file)
  gbif_genus_hierarchy_old <- read_rds(genus_hierarchy_file)
  
  # Check if the manual list of genera has expanded
  if (length(unique(gbif_genus_names_old$cleaned_rodent_genus)) == length(unique(rodent_genus$cleaned_rodent_genus))) {
    message("GBIF genus dictionary is up-to-date. Skipping API call.")
    genus_names <- gbif_genus_names_old
    genus_hierarchy <- gbif_genus_hierarchy_old
  } else {
    message("New genera found. Updating GBIF genus dictionary.")
    
    # Identify genera already resolved from the species hierarchy
    resolved_genera_from_species <- unique(gbif_species_hierarchy_old$genus)
    
    # Filter the manual genus list to only include those that are not yet resolved
    genera_to_match <- rodent_genus %>%
      filter(!str_detect(rodent_genus, "sp\\.$|sp\\.?$|spp\\.?$|species")) %>% # Exclude ambiguous names
      distinct(cleaned_rodent_genus) %>%
      filter(!cleaned_rodent_genus %in% resolved_genera_from_species) %>%
      pull(cleaned_rodent_genus)
    
    # Perform the API call only on the new genera
    if (length(genera_to_match) > 0) {
      gbif_genus <- get_gbifid(genera_to_match, rank = "genus")
      resolve_genus <- classification(gbif_genus, db = "gbif")
      
      new_genus_hierarchy <- bind_rows(resolve_genus) %>%
        group_by(query) %>%
        distinct() %>%
        ungroup() %>%
        pivot_wider(id_cols = "query", names_from = rank, values_from = name) %>%
        left_join(bind_rows(resolve_genus) %>%
                    group_by(query) %>%
                    distinct() %>%
                    summarise(gbif_id = id[which.max(rank == "genus")]), by = "query") %>%
        arrange(kingdom, phylum, class, order, family, genus) %>%
        drop_na(genus)
      
      new_genus_names <- tibble(
        cleaned_rodent_genus = genera_to_match,
        query = as.character(gbif_genus)
      ) %>%
        left_join(new_genus_hierarchy %>%
                    distinct(genus, query, family, order, class, gbif_id), by = "query")
      
      # Combine old and new genus data
      genus_hierarchy <- bind_rows(gbif_genus_hierarchy_old, new_genus_hierarchy)
      genus_names <- bind_rows(gbif_genus_names_old, new_genus_names)
      
    } else {
      # No new genera to match, so keep the old data
      genus_hierarchy <- gbif_genus_hierarchy_old
      genus_names <- gbif_genus_names_old
    }
    
    # Save the updated files
    write_rds(genus_names, genus_names_file)
    write_rds(genus_hierarchy, genus_hierarchy_file)
  }
} else {
  message("GBIF genus dictionary not found. Creating new dictionary.")
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
  
  # Save the new files
  write_rds(genus_names, genus_names_file)
  write_rds(genus_hierarchy, genus_hierarchy_file)
}


# Match to Higher GBIF Taxa -----------------------------------------------
# Step 1: Find all names from the manual list that were not matched at the species level.
unmatched_at_species_level <- rodent_names %>%
  anti_join(species_names %>% filter(!is.na(gbif_id)), by = "cleaned_rodent_names")

# Step 2: From that list, find all names that were also not matched at the genus level.
unmatched_at_genus_level <- unmatched_at_species_level %>%
  # Create a temporary column with just the genus name for joining
  mutate(cleaned_genus_for_join = word(cleaned_rodent_names, 1)) %>%
  anti_join(genus_names %>% filter(!is.na(gbif_id)), by = c("cleaned_genus_for_join" = "cleaned_rodent_genus")) %>%
  select(-cleaned_genus_for_join)

# The final list of higher taxa to match is the combination of these unresolved names
higher_taxa_to_match <- unique(unmatched_at_genus_level$cleaned_rodent_names)

higher_taxa_table_file <- here("data", "matching", "gbif_higher_taxa.rds")
if (file.exists(higher_taxa_table_file) && length(higher_taxa_to_match) == 0) {
  message("GBIF higher taxa dictionary is up-to-date. Skipping API call.")
  higher_taxa_table <- read_rds(higher_taxa_table_file)
} else {
  message("New higher taxa found. Updating GBIF higher taxa dictionary.")
  
  # Manual cleaning for higher taxa
  higher_taxa_manual <- tibble(scientificName = higher_taxa_to_match) %>%
    mutate(clean_higher_taxa = case_when(
      # Rodentia and its variants
      str_detect(tolower(scientificName), "rodentia|rodentia sp.|social rat|unknown") ~ "Rodentia",
      # Specific families and orders
      str_detect(scientificName, "Muridae|Murinae") ~ "Muridae",
      str_detect(scientificName, "Cricetidae|Neotominae") ~ "Cricetidae",
      str_detect(scientificName, "Soricidae|Soricidae sp.") ~ "Soricidae",
      str_detect(scientificName, "Sciuridae|Scuirdae") ~ "Sciuridae",
      str_detect(scientificName, "Didelphimorphia") ~ "Didelphimorphia",
      str_detect(scientificName, "Dipodidae") ~ "Dipodidae",
      str_detect(scientificName, "Ochotonidae") ~ "Ochotonidae",
      # Mammalia and other general categories
      str_detect(scientificName, "Other") ~ "Mammalia",
      # Fallback to the original name for any remaining specific taxa
      TRUE ~ scientificName
    ))
  
  # GBIF API calls
  gbif_higher <- get_gbifid(unique(higher_taxa_manual$clean_higher_taxa))
  resolve_higher <- classification(gbif_higher, db = "gbif")
  
  higher_hierarchy <- rbind(resolve_higher) %>%
    group_by(query) %>%
    distinct() %>%
    ungroup() %>%
    pivot_wider(id_cols = "query", names_from = rank, values_from = name) %>%
    mutate(gbif_id = as.integer(query))
  
  higher_taxa_table <- tibble(
    clean_higher_taxa = unique(higher_taxa_manual$clean_higher_taxa),
    query = as.character(gbif_higher)
  ) %>%
    left_join(higher_hierarchy, by = "query") %>%
    left_join(higher_taxa_manual, by = "clean_higher_taxa") %>%
    rename("extracted_name" = "scientificName")
  
  # Save the new file
  write_rds(higher_taxa_table, higher_taxa_table_file)
}

