# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SURVEYS ADDED IN PR #5 (ELSOC 2016-2023, GSS 2014): loaders, recodes and the
# per-survey configuration of the calibrated networks.
# source("scripts/01_survey_data.R") from the repo root (02, 02c, 02_ELSOC_homophily,
# 03_ELSOC_MUR_calculation, 03_GSS14_MUR_calculation, 06e, 07, 09, 10).
#
# Raw data live outside the repo, next to GSS 2004 and ATP 2014:
#   <NDFU_SURVEYS>/ELSOC 2016-2023/ELSOC_Long_2016_2023.RData   (20,761 rows, waves 1-7 = 2016-2023)
#   <NDFU_SURVEYS>/GSS 2014/General Social Survey 2014 Cross-Section and Panel Combined.DTA
# NDFU_SURVEYS defaults to ../A - Proyecto de Tesis/A - Trabajo 1/B - Surveys Data
# (relative to the repo root); set the environment variable to override.
#
# Surveys of the sweep (SURVEY_CONFIG):
#   GSS    GSS 2004 nodes, collective action (9 ISSP items); GSS 2004 race targets
#   ATP    ATP 2014 nodes, new products (6 METECH items); GSS 2004 race targets
#   GSS14  GSS 2014 nodes, collective action (8 of the 9 items: no INTERPOL); GSS 2004
#          race targets (no ego networks in 2014), POLVIEWS as attribute 'ideol'
#   ELSOC  ELSOC 2018 nodes, willingness to act (c09, 4 items); targets and homophily
#          coefficients from the ELSOC confidants (r13, waves 2017 and 2019),
#          calibrated by EDUCATION (5 groups); Blau space without race, with ideology
#
# ELSOC recodes (alters have coarser categories than egos; egos are mapped to them):
#   education  ego m01 (10 levels) and alter r13_educ (5) -> Primary (basica or less),
#              SecondaryInc, Secondary, Technical (complete or not), University (complete
#              or not, and graduate); educ_num = approximate years 6, 10, 12, 14, 16
#   religion   ego m38 (9) and alter r13_relig (5) -> Catholic, Evangelical (ego:
#              evangelical or historical protestant), None (ego: none or believer without
#              affiliation), AtheistAgnostic, Other
#   ideology   ego c15 (0 left ... 10 right; 11 independent, 12 none) and alter r13_ideol
#              (right, centre-right, centre, centre-left, left, none) -> Left (ego 0-4;
#              alter left or centre-left), Right (6-10; right or centre-right), CentreNone
#              (5, 11, 12; centre or none). The ego scale and the perceived ideology of
#              the alters do not line up at the centre (23.5% of egos at 5, 6% of alters
#              "centre"); with three categories the marginals match (egos 21/64/15%,
#              alters 19/67/14%). Decided with the user, 2026-11-11.
#   alter age  < 12 set to missing (132 alters)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
suppressWarnings(suppressMessages({ library(data.table); library(haven) }))

NDFU_SURVEYS <- Sys.getenv("NDFU_SURVEYS", file.path("..", "A - Proyecto de Tesis", "A - Trabajo 1", "B - Surveys Data"))
ELSOC_FILE <- file.path(NDFU_SURVEYS, "ELSOC 2016-2023", "ELSOC_Long_2016_2023.RData")
GSS14_FILE <- file.path(NDFU_SURVEYS, "GSS 2014", "General Social Survey 2014 Cross-Section and Panel Combined.DTA")
GSS04_FILE <- "data/02_GSS_network_ergm/GSS_2004_NORC.dta"

# ---- per-survey configuration of the networks ----------------------------------------------------
#   group      vertex attribute calibrated by group (activity + homophily, Newton)
#   groups     its levels; cal = the calibrated ones (the rest: relative degree 1 and the
#              common coefficient coef_group, as 'Other' in GSS)
#   cats       categorical attributes kept EXACT by the x2 rewiring (mixing matrix of the
#              group; same-attribute counts of the others)
#   nums       continuous attributes kept within 1% by the rewiring
#   blau       inputs of the Gower distance d_ij of the gate (05_unified_engine.R)
#   cell       attributes whose combinations are the circles of the theta arm
#   ideol      TRUE if the nodes carry 'ideol' (adds the theta arm with ideology circles)
RACES <- c("White", "Black", "Hispanic", "Asian", "Other")
EDU5  <- c("Primary", "SecondaryInc", "Secondary", "Technical", "University")
EDU_YEARS <- c(Primary = 6, SecondaryInc = 10, Secondary = 12, Technical = 14, University = 16)
REL5_CL <- c("Catholic", "Evangelical", "None", "AtheistAgnostic", "Other")
IDE3  <- c("Left", "CentreNone", "Right")
US_CONFIG <- list(group = "race", groups = RACES, cal = 1:4, cats = c("race", "sex", "relig"), nums = c("age", "educ_num"),
                  blau = c("age", "educ_num", "race", "relig", "sex"), cell = c("race", "relig"), ideol = FALSE,
                  targets = "GSS 2004")
SURVEY_CONFIG <- list(
  GSS   = US_CONFIG,
  ATP   = US_CONFIG,
  GSS14 = modifyList(US_CONFIG, list(ideol = TRUE)),
  ELSOC = list(group = "educ", groups = EDU5, cal = 1:5, cats = c("educ", "sex", "relig", "ideol"), nums = c("age", "educ_num"),
               blau = c("age", "educ_num", "sex", "relig", "ideol"), cell = c("educ", "relig"), ideol = TRUE,
               targets = "ELSOC 2017+2019"))
survey_config <- function(S) { cf <- SURVEY_CONFIG[[S]]; if (is.null(cf)) stop("unknown survey ", S, "; known: ",
                                                                           paste(names(SURVEY_CONFIG), collapse = ", ")); cf }
SURVEY_LABEL <- c(GSS = "GSS 2004", ATP = "ATP 2014", GSS14 = "GSS 2014", ELSOC = "ELSOC 2018")

# ---- ELSOC ------------------------------------------------------------------------------------------
.elsoc_cache <- new.env()
elsoc_load <- function() {
  if (is.null(.elsoc_cache$df)) {
    if (!file.exists(ELSOC_FILE)) stop("ELSOC not found at ", ELSOC_FILE, " (set NDFU_SURVEYS)")
    e <- new.env(); load(ELSOC_FILE, envir = e); .elsoc_cache$df <- as.data.table(e$elsoc_long_2016_2023)
  }
  .elsoc_cache$df
}
ELSOC_YEAR <- c(`1` = 2016, `2` = 2017, `3` = 2018, `4` = 2019, `5` = 2021, `6` = 2022, `7` = 2023)
num_pos <- function(x) { x <- as.numeric(x); x[x < 0] <- NA; x }      # ELSOC: -999 ... -666 are missing
elsoc_educ_ego   <- function(m01) EDU5[c(1, 1, 1, 2, 3, 4, 4, 5, 5, 5)[num_pos(m01)]]
elsoc_educ_alter <- function(x) EDU5[num_pos(x)]
elsoc_relig_ego  <- function(m38) REL5_CL[c(1, 2, 2, 5, 3, 5, 4, 4, 3)[num_pos(m38)]]
elsoc_relig_alter <- function(x) REL5_CL[num_pos(x)]
elsoc_ideol_ego  <- function(c15) { x <- num_pos(c15)
  ifelse(is.na(x), NA_character_, ifelse(x <= 4, "Left", ifelse(x >= 6 & x <= 10, "Right", "CentreNone"))) }
elsoc_ideol_alter <- function(x) c("Right", "Right", "CentreNone", "Left", "Left", "CentreNone")[num_pos(x)]
elsoc_sex <- function(x) c("Male", "Female")[num_pos(x)]
ELSOC_C09 <- sprintf("c09_%02d", 1:4)             # willingness: petition, march, strike, social media (1-5)
ELSOC_C08 <- sprintf("c08_%02d", 1:4)             # frequency of the same four acts (1-5); c08_05 cacerolazo from 2019

# egos of the network waves (r13: 2 = 2017, 4 = 2019, 6 = 2022): one row per respondent and wave
elsoc_egos <- function(waves = c(2, 4)) {
  df <- elsoc_load()
  d <- df[ola %in% waves & !is.na(num_pos(r13_sexo_01))]
  d[, .(idencuesta = as.numeric(idencuesta), ola = as.integer(ola), year = ELSOC_YEAR[as.character(ola)],
        sex = elsoc_sex(m0_sexo), age = num_pos(m0_edad), educ = elsoc_educ_ego(m01), relig = elsoc_relig_ego(m38),
        ideol = elsoc_ideol_ego(c15), nredes = num_pos(r13_nredes))][, educ_num := unname(EDU_YEARS[educ])][]
}
# confidants (up to 5 per ego); rel: 1 spouse/partner, 2 child, 3 other relative, 4 friend, 5 other acquaintance
elsoc_alters <- function(waves = c(2, 4)) {
  df <- elsoc_load()
  d <- df[ola %in% waves]
  a <- rbindlist(lapply(1:5, function(k) data.table(
    idencuesta = as.numeric(d$idencuesta), ola = as.integer(d$ola), k = k,
    rel = num_pos(d[[sprintf("r13_relacion_%02d", k)]]), sex = elsoc_sex(d[[sprintf("r13_sexo_%02d", k)]]),
    age = num_pos(d[[sprintf("r13_edad_%02d", k)]]), educ = elsoc_educ_alter(d[[sprintf("r13_educ_%02d", k)]]),
    relig = elsoc_relig_alter(d[[sprintf("r13_relig_%02d", k)]]), ideol = elsoc_ideol_alter(d[[sprintf("r13_ideol_%02d", k)]]))))
  a <- a[!is.na(sex)]
  a[age < 12, age := NA]
  a[, `:=`(kin = rel %in% 1:3, educ_num = unname(EDU_YEARS[educ]))][]
}
# node pool of the networks: one wave (default 2018), willingness c09 complete, attributes of that wave
#   propensity_score = mean of the four c09 items rescaled to [0, 1] (high = willing)
#   mur_score = 1 - propensity_score (q: aversion)
elsoc_nodes <- function(wave = 3) {
  df <- elsoc_load(); d <- df[ola == wave]
  X <- sapply(ELSOC_C09, function(v) num_pos(d[[v]]))
  out <- data.table(idencuesta = as.numeric(d$idencuesta), sex = elsoc_sex(d$m0_sexo), age = num_pos(d$m0_edad),
                    educ = elsoc_educ_ego(d$m01), relig = elsoc_relig_ego(d$m38), ideol = elsoc_ideol_ego(d$c15))
  out <- cbind(out, as.data.table(X))
  out <- out[complete.cases(out)]
  out[, educ_num := unname(EDU_YEARS[educ])]
  out[, propensity_score := (rowSums(.SD) - 4) / 16, .SDcols = ELSOC_C09]
  out[, mur_score := 1 - propensity_score][]
}

# ---- GSS 2004 and 2014 ---------------------------------------------------------------------------------
race_gss <- function(d) with(d, dplyr::case_when(hispanic %in% c(2:6, 10, 15, 20:24, 30, 40, 41) ~ "Hispanic",   # as in 02
                                                 racecen1 == 1 ~ "White", racecen1 == 2 ~ "Black",
                                                 racecen1 %in% c(3, 11, 15) ~ "Other", racecen1 %in% c(4:10) ~ "Asian",
                                                 racecen1 == 16 ~ "Hispanic", TRUE ~ NA_character_))
GSS_CA9 <- c("signdpet", "avoidbuy", "joindem", "attrally", "cntctgov", "polfunds", "usemedia", "interpol", "actlaw")
GSS_CA8 <- setdiff(GSS_CA9, "interpol")           # INTERPOL was not asked in 2014
GSS_SCI <- c("nextgen", "toofast", "advfront", "scibnfts")
.gss_cache <- new.env()
gss_load <- function(year = 2014) {
  key <- paste0("y", year)
  if (is.null(.gss_cache[[key]])) {
    f <- if (year == 2014) GSS14_FILE else GSS04_FILE
    if (!file.exists(f)) stop("GSS ", year, " not found at ", f)
    g <- read_dta(f); names(g) <- tolower(names(g)); .gss_cache[[key]] <- g
  }
  .gss_cache[[key]]
}
# respondent attributes coded as in 02 (sex, race, educ in years, age, relig in 5) + polviews -> ideol (3)
gss_attributes <- function(g) {
  nn <- function(x) as.numeric(zap_labels(x))
  data.table(gss_id = nn(g$id), sex = c("Male", "Female")[nn(g$sex)], race = race_gss(g),
             educ_num = ifelse(nn(g$educ) %in% c(98, 99), NA_real_, nn(g$educ)),
             age = ifelse(nn(g$age) %in% c(98, 99), NA_real_, nn(g$age)),
             relig = dplyr::case_when(nn(g$relig) == 1 ~ "Protestant", nn(g$relig) == 2 ~ "Catholic", nn(g$relig) == 3 ~ "Jewish",
                                      nn(g$relig) == 4 ~ "None", nn(g$relig) %in% 5:13 ~ "OtherRelig", TRUE ~ NA_character_),
             polviews = ifelse(nn(g$polviews) %in% 1:7, nn(g$polviews), NA_real_))[
    , ideol := ifelse(is.na(polviews), NA_character_, ifelse(polviews <= 3, "Liberal", ifelse(polviews == 4, "Moderate", "Conservative")))][]
}
# item matrix: ISSP action items 1-4 (1 done in past year ... 4 never) and actlaw 1-4; 0, 8, 9 missing
gss_items <- function(g, items) { X <- sapply(items, function(v) { x <- as.numeric(zap_labels(g[[v]])); x[!x %in% 1:4] <- NA; x })
  if (!is.matrix(X)) X <- matrix(X, ncol = length(items), dimnames = list(NULL, items)); X }
# collective action propensity = sum(4 - x) / (3 x items): 0 = would never do any, 1 = did all in the past year
ca_propensity <- function(X) rowSums(4 - X) / (3 * ncol(X))
# science attitudes: nextgen and advfront agree = pro (1-4), toofast agree = con (1-4), scibnfts 1 benefits
# greater, 2 about equal, 3 harmful greater; each rescaled to [0, 1] toward pro-science, then averaged
sci_items <- function(g) { nn <- function(v, ok) { x <- as.numeric(zap_labels(g[[v]])); x[!x %in% ok] <- NA; x }
  cbind(nextgen = (4 - nn("nextgen", 1:4)) / 3, toofast = (nn("toofast", 1:4) - 1) / 3,
        advfront = (4 - nn("advfront", 1:4)) / 3, scibnfts = (3 - nn("scibnfts", 1:3)) / 2) }
# GSS 2014 node pool: the 2014 cross-section (samptype 2014; the ISSP module was asked only there)
#   index "ca"   collective action, 8 items; "sci" science attitudes, 4 items
gss14_nodes <- function(index = c("ca", "sci")) {
  index <- match.arg(index); g <- gss_load(2014); cs <- as.numeric(g$samptype) == 2014
  A <- gss_attributes(g)
  if (index == "ca") { X <- gss_items(g, GSS_CA8); p <- ca_propensity(X) } else { X <- sci_items(g); p <- rowMeans(X) }
  out <- cbind(A, as.data.table(X))[, propensity_score := p][cs]
  out <- out[complete.cases(out[, c("age", "sex", "educ_num", "race", "relig", "ideol", "propensity_score"), with = FALSE])]
  out[, mur_score := 1 - propensity_score][]
}
