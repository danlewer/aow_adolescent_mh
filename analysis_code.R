setwd("H:/AoW_SP795_MH_inequalities")

# ====================================
# AOW mental health question locations
# ------------------------------------

# MODULE 1: https://borninbradford.nhs.uk/our-data/questionnaires/age-of-wonder-module-1-survey-231/
# - EDE-QS (eating problems) - can you score boys and girls the same?

# MODULE 2: https://borninbradford.nhs.uk/our-data/questionnaires/age-of-wonder-module-2-survey-232/ 
# - RCADS25 (anxiety & depression)
# - SDQ (emotional & behavioual problems)
# - SWEMWBS (wellbeing)
# - Self harm question
# - ULS-4 (loneliness) - unsure of provenance
# - PLIKS-8 (unusal experiences/ psychosis) - unsure of provenance

# RECRUITMENT FILE
# - age
# - sex
# - ethnicity
# - IMD
# - FSM

# =======================
# functions and libraries
# -----------------------

library(haven) # for read_dta (reading Stata files)
library(data.table) # for processing data
library(RColorBrewer) # for plot colours
library(sandwich) # for sandwich estimator
library(lmtest) # for sandwich estimator
library(Amelia) # for multiple imputation

read_dta2 <- function (...) data.table(read_dta(...)) # read Stata file and convert to data.table

ecdf2 <- function (v) { # simple ecdf function (retruns x/y values rather than callable function)
  v <- v[!is.na(v)]
  ux <- sort(unique(v))
  idx <- match(v, ux)
  counts <- tabulate(idx, nbins = length(ux))
  y <- cumsum(counts) / length(v)
  list(x = ux, y = y)
}

miqr <- function (var = 'age_mod2_m', digs = 1, level = '', DAT) {
  v <- DAT[, get(var)]
  v <- quantile(v, na.rm = T, probs = c(0.5, 0.25, 0.75))
  v <- format(round(v, digs), nsmall = digs, big.mark = ',')
  v <- paste0(v[1], ' [', v[2], '-', v[3], ']')
  data.frame(var = var, level = level, value = v)
}

cat <- function (var = 'eth7', digs = 1, ord = T, DAT) {
  tab <- DAT[, .N, get(var)]
  names(tab) <- c('level', 'value')
  if (ord) {
    tab <- tab[order(value, decreasing = T)]
  } else {
    tab <- tab[order(level, decreasing = F)]
  }
  pc <- tab$value / sum(tab$value) * 100
  pc <- format(round(pc, digs), digits = digs, nsmall = digs)
  tab$value <- formatC(tab$value, big.mark = ',')
  tab$value <- paste0(tab$value, ' (', pc, ')')
  tab$value <- gsub('\\(', ' (', gsub(' ', '', tab$value))
  cbind(var = var, tab)
}

propci <- function (X, N, form = T, digs = 1) { # vectorized confidence interval for proportions
  pt <- function (x, n) c(x/n, binom.test(x, n)$conf.int[1:2])
  r <- mapply(pt, x = X, n = N)
  r <- `colnames<-`(t(r), c('prop', 'lower', 'upper'))
  if (form) {
    r <- format(round(r * 100, digs), nsmall = digs)
    paste0(r[,1], ' (', r[,2], '-', r[,3], ')')
  } else {
    r
  }
}

# ==========================================
# read Stata files and convert to data.table
# ------------------------------------------

rec   <- read_dta2("AoWRecsummary.dta", col_select = c('aow_person_id', 'age_recruitment_m', 'recruitment_era', 'recruitment_date', 'year_group', 'gender', 'school_id', 'ethnicity_1', 'ethnicity_2', 'fsm', 'sen', 'IMD_2019_decile', 'LSOA11CD'))
mod1  <- read_dta2("AoW_MH_SP795_BiB_AgeOfWonder.survey_mod231_main_dr24.dta")
mod1$year_group <- NULL # use recruitment year group
mod1d <- read_dta2("AoW_MH_SP795_BiB_AgeOfWonder.survey_mod231_derived_dr24.dta")
mod2  <- read_dta2("AoW_MH_SP795_BiB_AgeOfWonder.survey_mod232_main_dr24.dta")
sdqs  <- read_dta2("AoWRawSDQv2.dta")[, -c('BiBPersonID', 'is_bib', 'age_recruitment_y', 'age_recruitment_m', 'school_id', 'year_group')]
pliks_addition <- read_dta2("Additional_BiB_AgeOfWonder.survey_mod232_main_dr24.dta")
table(pliks_addition$aow_person_id == mod2$aow_person_id) # should all be TRUE - OK to bind this new variable to mod2
mod2$awb2_11_psychosis_10_r4 <- pliks_addition$awb2_11_psychosis_10_r4

# ================
# limit to 2023/24
# ----------------

start_2324 <- as.Date('2023-09-01', origin = '1970-01-01')

mod1 <- mod1[survey_date >= start_2324]
mod2 <- mod2[survey_date >= start_2324]
rec <- rec[recruitment_era == '2023-24']
sdqs <- sdqs[recruitment_era == '2023-24']

# ===========
# deduplicate
# -----------

# mod1: 5 individuals completed more than once, eg. due to moving schools - take first response
# ---------------------------------------------------------------------------------------------

table(table(mod1$aow_person_id))
mod1 <- mod1[, .(minsurvy = min(survey_date)), aow_person_id][mod1, on = 'aow_person_id']
mod1 <- mod1[survey_date == minsurvy, -'minsurvy']

# mod1: 4 individuals completed more than once, eg. due to moving schools - take first response
# ---------------------------------------------------------------------------------------------

table(table(mod2$aow_person_id))
mod2 <- mod2[, .(minsurvy = min(survey_date)), aow_person_id][mod2, on = 'aow_person_id']
mod2 <- mod2[survey_date == minsurvy, -'minsurvy']

# SDQs: 4 individuals completed more than once - take first response
# ------------------------------------------------------------------

table(table(sdqs$aow_person_id))
sdqs <- sdqs[, .(minsurvy = min(survey_date)), aow_person_id][sdqs, on = 'aow_person_id']
sdqs <- sdqs[survey_date == minsurvy, -'minsurvy']

# recruitment file - 59 individuals recorded more than once - take the first
# --------------------------------------------------------------------------

table(table(rec$aow_person_id))
rec <- rec[, .(minsurvy = min(recruitment_date)), aow_person_id][rec, on = 'aow_person_id']
rec <- rec[recruitment_date == minsurvy, -'minsurvy']

# ====================================
# create dataset of unique individuals
# ------------------------------------

d <- data.table(aow_person_id = unique(c(mod1$aow_person_id, mod2$aow_person_id)))
d[, mod1 := aow_person_id %in% mod1$aow_person_id]
d[, mod2 := aow_person_id %in% mod2$aow_person_id]

setnames(mod1, c('age_survey_m', 'survey_date'), c('age_mod1_m', 'mod1_date'))
setnames(mod2, c('age_survey_m', 'survey_date'), c('age_mod2_m', 'mod2_date'))

d <- mod1[d, on = 'aow_person_id']
d <- mod2[d, on = 'aow_person_id']

# ================================
# summarise mental health outcomes
# --------------------------------

# RCADS25 (module 2)
# ------------------

# 1 = never, 2 = sometimes, 3 = often, 4 = always; scores 0-3 respectively
# https://psycnet.apa.org/doiLanding?doi=10.1037%2Fa0027283

rcads_vars <- paste0('awb2_1_illhealth_', 1:25)
rcads_depression_vars <- paste0('awb2_1_illhealth_', c(1, 4, 8, 10, 13, 15, 16, 19, 21, 24))
rcads_anxiety_vars <- paste0('awb2_1_illhealth_', c(2, 3, 5, 6, 7, 9, 11, 12, 14, 17, 18, 20, 22, 23, 25))

d[, (rcads_vars) := lapply(.SD, function (x) x-1), .SDcols = rcads_vars]
d[, rcads25 := rowSums(d[, rcads_vars, with = F])]
d[, rcads25_depression := rowSums(d[, rcads_depression_vars, with = F])]
d[, rcads25_anxiety := rowSums(d[, rcads_anxiety_vars, with = F])]

# RCADS T-scores
# --------------

tscores <- fread("https://raw.githubusercontent.com/danlewer/aow_adolescent_mh/refs/heads/main/rcads_tscores.csv", col.names = c('raw', 'gender', 'rcads25_depression_tscore', 'rcads25_anxiety_tscore', 'rcads25_tscore'))
tscores$gender <- ifelse(tscores$gender == 'Boys', 'Male', 'Female')
d$gender <- ifelse(d$gender == 1, "Female", "Male")
tscores_rcads25 <- tscores[, .(gender = gender, rcads25 = raw, rcads25_tscore = rcads25_tscore)]
d <- tscores_rcads25[d, on = c('gender', 'rcads25')]
tscores_rcads25_anxiety <- tscores[, .(gender = gender, rcads25_anxiety = raw, rcads25_anxiety_tscore = rcads25_anxiety_tscore)]
d <- tscores_rcads25_anxiety[d, on = c('gender', 'rcads25_anxiety')]
tscores_rcads25_depression <- tscores[, .(gender = gender, rcads25_depression = raw, rcads25_depression_tscore = rcads25_depression_tscore)]
d <- tscores_rcads25_depression[d, on = c('gender', 'rcads25_depression')]

# SDQ (module 2)
# --------------

# 0 = not true; 1  = somewhat true; 2 = certainly true
# https://www.sdqinfo.org/py/sdqinfo/c0.py

sdq_vars <- paste0('awb2_1_sdq_', 1:25, '_a10')

# reverse score positive questions
sdqs$awb2_1_sdq_7_a10  <- 2 - sdqs$awb2_1_sdq_7_a10
sdqs$awb2_1_sdq_21_a10 <- 2 - sdqs$awb2_1_sdq_21_a10
sdqs$awb2_1_sdq_25_a10 <- 2 - sdqs$awb2_1_sdq_25_a10
sdqs$awb2_1_sdq_11_a10 <- 2 - sdqs$awb2_1_sdq_11_a10
sdqs$awb2_1_sdq_14_a10 <- 2 - sdqs$awb2_1_sdq_14_a10

sdq_emotional_vars <- c('awb2_1_sdq_3_a10', 'awb2_1_sdq_8_a10', 'awb2_1_sdq_13_a10', 'awb2_1_sdq_16_a10', 'awb2_1_sdq_24_a10')
sdq_conduct_vars <- c('awb2_1_sdq_5_a10', 'awb2_1_sdq_7_a10', 'awb2_1_sdq_12_a10', 'awb2_1_sdq_18_a10', 'awb2_1_sdq_22_a10')
sdq_hyper_vars <- c('awb2_1_sdq_2_a10', 'awb2_1_sdq_10_a10', 'awb2_1_sdq_15_a10', 'awb2_1_sdq_21_a10', 'awb2_1_sdq_25_a10')
sdq_peer_vars <- c('awb2_1_sdq_6_a10', 'awb2_1_sdq_11_a10', 'awb2_1_sdq_14_a10', 'awb2_1_sdq_19_a10', 'awb2_1_sdq_23_a10')

sdqs[, sdq_emotional := rowSums(sdqs[, sdq_emotional_vars, with = F])]
sdqs[, sdq_conduct := rowSums(sdqs[, sdq_conduct_vars, with = F])]
sdqs[, sdq_hyper := rowSums(sdqs[, sdq_hyper_vars, with = F])]
sdqs[, sdq_peer := rowSums(sdqs[, sdq_peer_vars, with = F])]
sdqs$sdq_internal  <- sdqs$sdq_emotional + sdqs$sdq_peer
sdqs$sdq_external  <- sdqs$sdq_conduct + sdqs$sdq_hyper
sdqs$sdq_total     <- sdqs$sdq_internal + sdqs$sdq_external

d <- sdqs[, c('aow_person_id', sdq_vars, 'sdq_total', 'sdq_emotional', 'sdq_conduct', 'sdq_hyper', 'sdq_peer', 'sdq_internal', 'sdq_external'), with = F][d, on = 'aow_person_id']

# SWEMWBS (module 2)
# ------------------

# 1 = none of the time, 2 = rarely, 3 = some of the time, 4 = often , 5 = all of the time
# https://www.corc.uk.net/outcome-experience-measures/short-warwick-edinburgh-mental-wellbeing-scale-swemwbs/
# they want you to translate it to a 'metric score' (not sure why) - https://hqlo.biomedcentral.com/articles/10.1186/1477-7525-7-15

swem_vars <- c('awb2_2_optmstc_1_a4', 'awb2_2_useful_2_a4', 'awb2_2_relxed_3_a4', 'awb2_2_problems_4_a4', 'awb2_2_think_clr_5_a4', 'awb2_2_close_othrs_6_a4', 'awb2_2_own_mnd_7_a4')
d[, swemwbs := rowSums(d[, swem_vars, with = F])]
metric_lookup <- fread('https://raw.githubusercontent.com/danlewer/aow_adolescent_mh/refs/heads/main/swemwbs_metric_lookup.csv', col.names = c('swemwbs', 'swemwbs_metric'))
d <- metric_lookup[d, on = 'swemwbs']

# EDE-QS (module 1)
# -----------------

# https://www.corc.uk.net/outcome-experience-measures/eating-disorder-examination-questionnaire-ede-q/
# eat vars: 1 = 0 days, 2 = 1-2 days, 3 = 3-5 days, 4 = 6-7 days; scores 0-3
# weight vars: 1 = not at all, 2 = slightly, 3 = moderately, 4 = markedly, scores 0 -3 
# EDE-QS cutoff (15): https://bmcpsychiatry.biomedcentral.com/articles/10.1186/s12888-020-02565-5
# note piped question awb2_12_eat_hbt_10_a5 - scores 0 if not answered

eat_vars <- paste0('awb2_12_eat_hbt_', 1:10, '_a5')
weight_vars <- c('awb2_12_wght_1_a5', 'awb2_12_wght_2_a5')
edeqs_vars <- c(eat_vars, weight_vars)
d[, (edeqs_vars) := lapply(.SD, function (x) x-1), .SDcols = edeqs_vars]
d$awb2_12_eat_hbt_10_a5[is.na(d$awb2_12_eat_hbt_10_a5) & !is.na(is.na(d$awb2_12_eat_hbt_10_a5))] <- 0
d[, edeqs := rowSums(d[, edeqs_vars, with = F])]

# self-harm (module 2)
# --------------------

d[, self_harm := as.integer(awb2_9_seek_hurt_self_a5)]

# ULS-4 (module 2) - UCLA loneliness scale
# ----------------------------------------

# https://pubmed.ncbi.nlm.nih.gov/7431205/
# 1 = hardly ever, 2 = some of the time, 3 = often
# unsure exactly where these questions came from

uls_vars <- paste0('awb2_4_loneliness_', 1:4)
d[, uls4 := rowSums(d[, uls_vars, with = F])]

# PLIKS-8
# -------

# unsure of the provenance of these questions
# 1 = yes, definitely; 2 = yes, maybe; 3 = no, never
# think 'awb2_11_psychosis_10_r4' is missing

pliks_vars <- c(paste0('awb2_11_psychosis_', c(1:5, 9, 10), '_r4'), 'awb2_11_pwrs_read_a4')
d[, (pliks_vars) := lapply(.SD, function (x) factor(x, c(-4, -2, 1:3), c(NA_integer_, NA_integer_, 2:0))), .SDcols = pliks_vars]
d[, (pliks_vars) := lapply(.SD, function (x) as.integer(as.character(x))), .SDcols = pliks_vars]
d[, pliks := rowSums(d[, pliks_vars, with = F])]

# no strong association with schizphrenia risk factors and pliks - https://pubmed.ncbi.nlm.nih.gov/18562177/

# ===================================================
# characteristics from school data / recruitment file
# ---------------------------------------------------

# population density
# note AoW has LSOA11CD, including LSOAs E01010677 and E01010835 which are inactive in LSOA22CD
popdens <- fread('https://raw.githubusercontent.com/danlewer/aow_adolescent_mh/refs/heads/main/geographical_data/lsoa_pop_density_census_2021.csv', select = c('LSOA21CD', 'calc_density'), col.names = c('LSOA11CD', 'popdens'))
rec <- popdens[rec, on = 'LSOA11CD']

# index of multiple deprivation
imdrank <- fread("https://raw.githubusercontent.com/danlewer/aow_adolescent_mh/refs/heads/main/imd2019rank3.csv")
rec <- imdrank[rec, on = 'LSOA11CD']
rec[, imd5 := factor(IMD_2019_decile, 1:10, ceiling(1:10 / 2))]

# special educational needs
rec$SEN <- as.integer(rec$sen)
rec$SEN[is.na(rec$SEN)] <- 3
rec[, SEN := factor(SEN, 0:3, c('none', 'SEN', 'EHCP', 'missing'))]

# free school meals
rec$fsm <- as.integer(rec$fsm)

# ethnicity
eth_lookup <- fread("https://raw.githubusercontent.com/danlewer/aow_adolescent_mh/refs/heads/main/eth_lookup.csv", col.names = c('ethnicity_2', 'eth_desc', 'eth5'))
rec <- eth_lookup[rec, on = 'ethnicity_2']
rec$eth5[is.na(rec$eth5)] <- 'Unknown'
rec$eth6 <- rec$eth5
rec$eth6[rec$ethnicity_2 == 41] <- 'White British'
rec$eth7 <- rec$eth6
rec$eth7[rec$ethnicity_2 == 12] <- 'Pakistani'
rec[, eth7 := factor(eth7, c('Pakistani', 'White British', 'Asian', 'Black', 'Mixed', 'Other', 'Unknown', 'White'))]
rec$eth7_2 <- ifelse(rec$eth7 == 'Unknown', NA_character_, as.character(rec$eth7)) # to allow imputation of missing ethnicity
eth_levs <- table(rec$eth7_2)
eth_levs <- names(eth_levs)[order(eth_levs, decreasing = T)]
rec$eth7_2 <- factor(rec$eth7_2, eth_levs)
rec$eth7_2 <- relevel(rec$eth7_2, ref = 'White British')

# year group
rec$year_group <- as.integer(rec$year_group)

# sex
rec[, sex := factor(gender, 1:2, c('Female', 'Male'))]

# add to main dataset
d <- rec[, c('aow_person_id', 'school_id', 'imd5', 'popdens', 'sex', 'fsm', 'SEN', 'eth7_2', 'year_group')][d, on = 'aow_person_id']

# =================================
# process variables in main dataset
# ---------------------------------

d$age <- ifelse(is.na(d$age_mod2_m), d$age_mod1_m, d$age_mod2_m)
d[, age := age/12]
d$season <- ifelse(is.na(d$mod2_season), d$mod1_season, d$mod2_season)
d[, season := factor(season, 1:4, c('winter', 'spring', 'summer', 'autumn'))]
d$survey_date <- as.Date(ifelse(is.na(d$mod2_date), d$mod1_date, d$mod2_date), origin = '1970-01-01')
d$season <- month(d$survey_date)
d$season <- factor(d$season, 1:12, c('winter', 'winter', 'spring', 'spring', 'spring', 'summer', 'summer', 'summer', 'autumn', 'autumn', 'autumn', 'winter'))

# ========================
# missing data in measures
# ------------------------

edeqs_vars2 <- setdiff(edeqs_vars, "awb2_12_eat_hbt_10_a5")
measures <- c('rcads25', 'rcads25_anxiety', 'rcads25_depression', 'swemwbs', 'uls4', 'sdq_total', 'sdq_emotional', 'sdq_peer', 'sdq_conduct', 'sdq_hyper', 'edeqs', 'pliks')

missf <- function (v, years_included = 8:10) {
  n <- length(v)
  x <- d[year_group %in% years_included, v, with = F]
  x <- rowSums(!is.na(x))
  responses <- c(sum(x == 0), sum(x > 0 & x < n), sum(x == n), sum(x > 0))
  pc <- format(round(responses / length(x) * 100, 1), nsmall = 1, digits = 1)
  responses <- formatC(responses, big.mark = ',')
  responses <- paste0(responses, '(', pc, ')')
  responses <- gsub('\\(', ' (', gsub(' ', '', responses))
  `names<-`(c(formatC(length(x), big.mark = ','), n, responses), c('n', 'items', 'none', 'some', 'all', 'any'))
}

miss <- mapply(missf,
               v = list(rcads_vars, rcads_anxiety_vars, rcads_depression_vars, swem_vars, uls_vars, sdq_vars, sdq_emotional_vars, sdq_peer_vars, sdq_conduct_vars, sdq_hyper_vars, edeqs_vars2, pliks_vars),
               years_included = list(c(8, 10), c(8, 10), c(8, 10), 8:10, 8:10, 9, 9, 9, 9, 9, 8:10, 10))
colnames(miss) <- measures
miss <- t(miss)

write.csv(miss, 'missing.csv')

# FAW missing data

d$missingRCADS <- rowSums(!is.na(d[, rcads_vars, with = F])) == 0
tmp <- d[year_group != 9, table(missingRCADS, fsm)]
t(t(tmp) / colSums(tmp))
tmp <- d[year_group != 9, table(missingRCADS, eth7_2)]
t(t(tmp) / colSums(tmp))
tmp <- d[year_group != 9, table(missingRCADS, mod2)]
t(t(tmp) / colSums(tmp))
tmp <- d[year_group != 9, table(missingRCADS, school_id)]
t(t(tmp) / colSums(tmp))
tmp <- d[year_group != 9, table(missingRCADS, imd5)]
t(t(tmp) / colSums(tmp))
tmp <- d[year_group != 9, table(missingRCADS, school_id)]
t(t(tmp) / colSums(tmp))
  
# ===============================
# describe sample characteristics 
# -------------------------------

desc_tab_sample <- rbind(miqr('age', DAT = d),
                         cat('sex', DAT = d),
                         cat('year_group', DAT = d),
                         miqr('popdens', digs = 0, DAT = d),
                         cat('season', DAT = d),
                         cat('SEN', ord = F, DAT = d),
                         cat('fsm', DAT = d),
                         cat('eth7_2', DAT = d),
                         cat('imd5', ord = F, DAT = d))
desc_tab_rec <- rbind(cat('sex', DAT = rec),
                      cat('year_group', DAT = rec),
                      miqr('popdens', digs = 0, DAT = rec),
                      cat('SEN', ord = F, DAT = rec),
                      cat('fsm', DAT = rec),
                      cat('eth7_2', DAT = rec),
                      cat('imd5', ord = F, DAT = rec))

names(desc_tab_sample)[3] <- 'sample'
names(desc_tab_rec)[3] <- 'recruitmemt'
desc_tab <- merge(desc_tab_rec, desc_tab_sample, all.x = T, sort = F)

write.csv(desc_tab, 'desc.csv')

# ===============================
# describe mental health measures
# -------------------------------

yax <- function(x, tickabove = F, ntick = 5) { # y-axis tick points
  l <- c(c(1, 2, 4, 5, 25) %o% 10^(0:8))
  d <- l[which.min(abs(x/ntick - l))]
  d <- 0:(ntick+1) * d
  i <- findInterval(x, d)
  if (tickabove) {i <- i + 1}
  d[seq_len(i)]
}

# histograms with cutoffs

par(mfrow = c(4, 3), mar = c(1, 1, 1, 1))
lapply(measures, function (x) hist(d[, get(x)], main = x))

p <- function (var = 'rcads25', cols = brewer.pal(11, 'Spectral')[7], cutoff = 70, TITLE = 'Depression and Anxiety\nRCADS-25') {
  x <- d[, .N, .(v = get(var))]
  x <- x[!is.na(v)]
  ymax <- max(x$N) * 1.3
  xmin <- min(x$v)
  xmax <- max(x$v) + 1
  plot(1, type = 'n', xlim = c(xmin, xmax), ylim = c(0, ymax), axes = F, xlab = NA, ylab = NA)
  rect(xmin, 0, xmax, ymax, col = 'grey98')
  rect(x$v, 0, x$v + 1, x$N, col = cols[1])
  if (cutoff <= max(x$v)) {with(x[v >= cutoff], rect(v, 0, v+1, N, col = cols[3]))}
  axis(1, pos = 0, tck = -0.03)
  axis(1, min(x$v):max(x$v), pos = 0, labels = F, tck = -0.01)
  axis(2, yax(ymax), pos = min(x$v), las = 2)
  text(xmax/2 + xmin/2, max(x$N) * 1.15, TITLE)
}

titles <- c('Anxiety & depression\nRCADS-25',
            'Anxiety\nRCADS-25 subscale',
            'Depression\nRCADS-25 subscale',
            'Wellbeing\nSWEMWBS',
            'Loneliness\nULS-4',
            'Behavioural & emotional problems\nSDQ',
            'Emotional problems\nSDQ subscale',
            'Peer problems\nSDQ subscale',
            'Conduct problems\nSDQ subscale',
            'Hyperactivity problems\nSDQ subscale',
            'Eating disorders\nEDE-QS',
            'Psychosis-like experiences\nPLIKS')

# guess the "non T-score" cutoff by estimating the same quantile in the raw score
quantile(d$rcads25, 1 - mean(d$rcads25_tscore > 70, na.rm = T), na.rm = T) # 36
quantile(d$rcads25_anxiety, 1 - mean(d$rcads25_anxiety_tscore > 70, na.rm = T), na.rm = T) #22
quantile(d$rcads25_depression, 1 - mean(d$rcads25_depression_tscore > 70, na.rm = T), na.rm = T) # 17

cutoffs <- c(36, 22, 17, Inf, Inf, 17, 8, 11, 15, Inf)
cutoffs <- rep(Inf, 12)

png('hists.png', height = 9.3, width = 7, units = 'in', res = 300)
par(mfrow = c(4, 3), mar = c(1, 1, 3, 2), oma = c(4, 4, 0, 0))
mapply(p, var = measures, TITLE = titles, cutoff = cutoffs)
title(xlab = 'Score', outer = T, line = 2)
title(ylab = 'Number of participants', outer = T)
dev.off()

# eCDFs by ethnicity

p2 <- function (var = 'rcads25', exs = 'eth7_2', levs = c('White British', 'Pakistani'), cols = brewer.pal(5, 'Set1'), TITLE = 'RCADS25', cutoff = 20) {
  #x <- d[eth7 %in% c('White British', 'Pakistani'), .(eth7 = eth7, var = get(var))]
  x <- d[get(exs) %in% levs, .(exs = get(exs), var = get(var))]
  x <- droplevels(x)
  iqr <- aggregate(var ~ exs, data = x, FUN = quantile, probs = c(0.25, 0.5, 0.75))
  xmax <- max(x$var, na.rm = T)
  x <- split(x, f = x$exs)
  x <- lapply(x, function (y) ecdf2(y$var))
  plot(1, type = 'n', xlim = c(0, xmax), ylim = c(0, 1), axes = F, xlab = NA, ylab = NA)
  rect(cutoff, 0, xmax, 1, col = 'grey92', border = NA)
  rect(0, 0, xmax, 1)
  mapply(lines, x, col = cols[1:2], type = 's')
  mapply(points, x, col = cols[1:2], cex = 0.7, pch = 19)
  segments(0, c(0.25, 0.5, 0.75), x1 = xmax, lty = 3)
  ys <- -0.075 * 1:4 + 0.05
  segments(iqr[1,-1], ys[1], y1 = c(0.25, 0.5, 0.75), col = cols[1], lty = 3)
  segments(iqr[2,-1], ys[2], y1 = c(0.25, 0.5, 0.75), col = cols[2], lty = 3)
  text(iqr[1,-1], ys[3], iqr[1,-1], col = cols[1])
  text(iqr[2,-1], ys[4], iqr[2,-1], col = cols[2])
  axis(1, pos = 0)
  axis(2, 0:4/4, pos = 0, las = 2)
  text(xmax * 0.05, 0.90, TITLE, adj = 0)
}

png('ecdfs_eth.png', height = 9.3, width = 7, units = 'in', res = 300)
par(mfrow = c(4, 3), mar = c(1, 1, 3, 2), oma = c(4, 4, 0, 0), xpd = NA)
mapply(p2, var = measures, TITLE = titles, cutoff = cutoffs)
title(xlab = 'Score', outer = T)
title(ylab = 'Empirical cumulative distribution', outer = T)
dev.off()

png('ecdfs_fsm.png', height = 9.3, width = 7, units = 'in', res = 300)
par(mfrow = c(4, 3), mar = c(1, 1, 3, 2), oma = c(4, 4, 0, 0), xpd = NA)
mapply(p2, var = measures, TITLE = titles, cutoff = cutoffs, exs = 'fsm', levs = list(c(0, 1)))
title(xlab = 'Score', outer = T)
title(ylab = 'Empirical cumulative distribution', outer = T)
dev.off()

# ===================
# multiple imputation
# -------------------

nonMHvars <- c('aow_person_id', 'imd5', 'popdens', 'age', 'sex', 'fsm', 'SEN', 'eth7_2', 'year_group', 'season')
dm <- d[, c(nonMHvars, rcads_vars, swem_vars, uls_vars, sdq_vars, edeqs_vars, pliks_vars), with = F]

# set.seed(2906)
# a <- amelia(dm, m = 20,
#             idvars = 'aow_person_id',
#             noms = c('imd5', 'sex', 'fsm', 'SEN', 'eth7_2', 'year_group', 'season'),
#             ords = c(rcads_vars, swem_vars, uls_vars, sdq_vars, edeqs_vars, pliks_vars))
# save(a, file = 'imputations_23sept2025b.RDS')
load("imputations_23sept2025b.RDS")

# make scales

scales <- function (x) {
  x$rcads25 <- rowSums(x[, rcads_vars, with = F])
  x$rcads25_anxiety <- rowSums(x[, rcads_anxiety_vars, with = F])
  x$rcads25_depression <- rowSums(x[, rcads_depression_vars, with = F])
  x$swemwbs <- rowSums(x[, swem_vars, with = F])
  x$uls4 <- rowSums(x[, uls_vars, with = F])
  x$sdq_total <- rowSums(x[, sdq_vars, with = F])
  x$sdq_emotional <- rowSums(x[, sdq_emotional_vars, with = F])
  x$sdq_conduct <- rowSums(x[, sdq_conduct_vars, with = F])
  x$sdq_hyper <- rowSums(x[, sdq_hyper_vars, with = F])
  x$sdq_peer <- rowSums(x[, sdq_peer_vars, with = F])
  x$edeqs <- rowSums(x[, edeqs_vars, with = F])
  x$pliks <- rowSums(x[, pliks_vars, with = F])
  x
}

imps <- lapply(a$imputations, scales)

# ==========================
# examine heteroskedasticity
# --------------------------

par(mfrow = c(3, 4), mar = c(2, 2, 4, 1), oma = c(4, 4, 0, 0))
het <- lapply(seq_along(measures), function (x) {
  f <- paste0(measures[x], '~eth7_2 + age + sex + season')
  m <- lm(f, data = d)
  plot(m$fitted.values, m$residuals, xlab = NA, ylab = NA)
  pm <- lm(m$residuals ~ poly(m$fitted.values, 3))
  pm <- cbind(m$fitted.values, predict(pm))
  pm <- pm[order(pm[,1]),]
  lines(pm[,1], pm[,2], col = 'red')
  abline(h = 0, col = 'blue')
  title(main = titles[x])
})
title(xlab = 'Fitted values', outer = T, line = 2)
title(ylab = 'Residuals', outer = T, line = 2)

# ===========================
# effect of ethnicity and SES
# ---------------------------

# smds with robust sandwich estimator

smd_sandwich <- function (DATA, exposure = 'eth7_2', outcome = 'rcads25', adj = c('age', 'sex', 'season'), type = 'HC1', years = c(8, 10), limit = F) {
  DATA <- DATA[year_group %in% years] # drop year groups
  form <- paste0(outcome, '~', exposure, '+',  paste0(adj, collapse = '+')) # formula
  m <- lm(form, data = DATA) # model
  coef <- coeftest(m, vcov. = vcovHC(m, type = type)) # extract coefficents
  coef <- coef[, c('Estimate', 'Std. Error')]
  sd_y <- sd(m$model[, outcome])
  coef[grepl(exposure, row.names(coef)),] / sd_y # smd
}

years_included <- list(
  rcads25 = c(8, 10),
  rcads25_anxiety = c(8, 10),
  rcads25_depression = c(8, 10),
  swemwbs = 8:10,
  uls4 = 8:10,
  sdq_total = 9,
  sdq_internal = 9,
  sdq_external = 9,
  sdq_emotional = 9,
  sdq_conduct = 9,
  sdq_hyper = 9,
  sdq_peer = 9,
  edeqs = 8:10,
  pliks = 10
)
years_included <- years_included[measures]

# (a) SMD for ethnicity
# ---------------------

r_eth1 <- lapply(imps, function (x) {
  r <- mapply(smd_sandwich, DATA = list(x), outcome = measures, exposure = 'eth7_2', years = years_included, SIMPLIFY = F)
  names(r) <- measures
  r
})

r_eth2 <- lapply(imps, function (x) {
  r <- mapply(smd_sandwich, DATA = list(x), outcome = measures, exposure = 'eth7_2', years = years_included, adj = list(c('age', 'sex', 'season', 'fsm', 'imd5')), SIMPLIFY = F)
  names(r) <- measures
  r
})

# (b) FSM
# -------

r_fsm <- lapply(imps, function (x) {
  r <- mapply(smd_sandwich, DATA = list(x), outcome = measures, exposure = 'fsm', years = years_included, SIMPLIFY = F)
  names(r) <- measures
  lapply(lapply(r, as.matrix), t)
})

# combine imputations

extract_amelia <- function (amres, alpha = 0.05, trans = F) {
  q <- lapply(measures, function (y) t(sapply(amres, function (x) x[[y]][,1])))
  names(q) <- measures
  se <- lapply(measures, function (y) t(sapply(amres, function (x) x[[y]][,2])))
  names(se) <- measures
  if (trans) {
    q <- lapply(q, t)
    se <- lapply(se, t)
  }
  rubin <- lapply(measures, function (x) {
    mi.meld(q = q[[x]], se = se[[x]])
  })
  names(rubin) <- measures
  rubin <- lapply(rubin, function (x) `colnames<-`(as.data.frame.matrix(t(rbind(x$q.mi, x$se.mi))), c('q', 'se')))
  z <- qnorm(1 - alpha/2)
  lapply(rubin, function (x) cbind(x, lower = x$q - x$se * z, upper = x$q + x$se * z))
}

r_eth1 <- extract_amelia(r_eth1)
r_eth2 <- extract_amelia(r_eth2)
r_fsm <- extract_amelia(r_fsm, trans = T)

# =====
# plots
# -----

offs <- 0.2
yoff <- seq(offs, -offs, length.out = nrow(r_eth1[[1]]))
ymids <- rev(seq_along(measures))
ys <- lapply(ymids, `+`, y = yoff)
ys <- c(outer(yoff, ymids, `+`))
x <- do.call(rbind, r_eth1)
x2 <- do.call(rbind, r_fsm)
gap <- 0.2
main_measures <- c(0, 1, 6, 7, 8, 11) + 0.5

png('smd_plot_v3.png', height = 11, width = 11.5, units = 'in', res = 300)

par(mar = c(11, 13, 4, 11), xpd = NA)

cols <- c('black', brewer.pal(5, 'Paired'))
plot(1, type = 'n', xlim= c(-1, 1.5), ylim = c(0, 13), axes = F, xlab = NA, ylab = NA)
rect(-1, main_measures, 0.5, main_measures + 1, col = 'grey96', border = NA)
rect(-1, 0.5, 0.5, 12.5)
segments(-1, 0:12 + 0.5, x1 = 0.5, lty = 3, lwd = 0.5)
segments(0, 0.5, y1 = 12.5)
points(x$q, ys, pch = 19, col = cols, cex = c(1.15, 0.7, 0.7, 0.7, 0.7, 0.7))
arrows(x$lower, ys, x1 = x$upper, length = 0.04, code = 3, angle = 90, col = cols)
axis(1, at = seq(-1, 0.5, 0.5), pos = 0.5)
arrows(-0.1, -1, x1 = -0.5, length = 0.13)
arrows(0.1, -1, x1 = 0.5, length = 0.13)
text(-0.1, -1.8, 'Higher score in\nWhite British\nadolescents', adj = 1)
text(0.1, -1.8, 'Higher score in\nminority ethnic\nadolescents', adj = 0)

rect(-0.5 + 1+gap, main_measures, 0.5 + 1+gap, main_measures + 1, col = 'grey96', border = NA)
rect(-0.5 + 1+gap, 0.5, 0.5 + 1+gap, 12.5)
segments(-0.5 +1+gap, 0:12 + 0.5, x1 = 0.5 +1+gap, lty = 3, lwd = 0.5)
segments(0 + 1+gap, 0.5, y1 = 12.5)
points(x2$q + 1+gap, ymids, pch = 19 , cex = 1.15)
arrows(x2$lower +1+gap, ymids, x1 = x2$upper +1+gap, length = 0.04, code = 3, angle = 90)
axis(1, at = seq(-0.5 +1+gap, 0.5 +1+gap, 0.5), labels = c(-0.5, 0, 0.5), pos = 0.5)
arrows(-0.1 +1+gap, -1, x1 = -0.5 +1+gap, length = 0.15)
arrows(0.1 +1+gap, -1, x1 = 0.5 +1+gap, length = 0.15)

text(-0.1 +1+gap, -1.8, 'Higher score in\nineligible\nadolescents', adj = 1)
text(0.1 +1+gap, -1.8, 'Higher score in\neligible\nadolescents', adj = 0)

text(-1.1, ymids, titles, adj = 1)
text(0, 13.5, 'Associations with\nethnicity\n(ref = White British: see colour key)')
text(0 +1+gap, 13.5, 'Associations with\nFree School Meal\neligibility')
text(0.5 +gap/2, -3, 'Standardised mean difference')

ysl <- seq(12, 7.5, length.out = 6)
points(rep(0.7 +1+gap, 6), ysl, col = cols, pch = 19, cex = c(1.15, 0.7, 0.7, 0.7, 0.7, 0.7))
arrows(0.6 +1+gap, ysl, x1 = 0.8 +1+gap, col = cols, angle = 90, code = 3, length = 0.05)
text(0.82 +1+gap, ysl, c('British\nPakistani', 'Other\nAsian', 'Mixed\nethnicities', 'Other White', 'Black', 'Other\nethnicities'), adj = 0)
text(0.6 +1+gap, 13.5, 'Key for\nethnicity\nresults', adj = 0)

dev.off()

# =====
# table
# -----

ethtab <- do.call(rbind, r_eth1)
ethtab <- format(round(ethtab, digits = 2), nsmall = 2, digits = 2)
ethtab <- cbind(measure = sub(".eth7_2.*", "", rownames(ethtab)), ethtab)
ethtab <- cbind(ethnicity = rep(levels(d$eth7_2)[-1], length(measures)), ethtab)
ethtab$smd <- paste0(ethtab$q, '(', ethtab$lower, ',', ethtab$upper, ')')
ethtab$smd <- gsub(' ', '', ethtab$smd)
ethtab$smd <- gsub('\\(', ' (', ethtab$smd)
ethtab$smd <- gsub(',', ', ', ethtab$smd)
ethtab$ethnicity <- factor(ethtab$ethnicity, c('Pakistani', 'Asian', 'Mixed', 'White', 'Black', 'Other'))
setDT(ethtab)
ethtab <- dcast(ethtab, measure ~ ethnicity, value.var = 'smd')

fsmtab <- do.call(rbind, r_fsm)
fsmtab <- format(round(fsmtab, digits = 2), nsmall = 2, digits = 2)
fsmtab <- cbind(measure = rownames(fsmtab), fsmtab)
fsmtab$smd <- paste0(fsmtab$q, '(', fsmtab$lower, ',', fsmtab$upper, ')')
fsmtab$smd <- gsub(' ', '', fsmtab$smd)
fsmtab$smd <- gsub('\\(', ' (', fsmtab$smd)
fsmtab$smd <- gsub(',', ', ', fsmtab$smd)

tab <- cbind(ethtab, fsm = fsmtab$smd)
fwrite(tab, 'smd_table.csv')

# compare complete case and multiple imputation results
# -----------------------------------------------------

cc <- function (x, var = 'eth7_2') {
  f <- paste0(x, '~', var, '+age+sex+season')
  sd_y <- sd(d[, get(x)], na.rm = T)
  m <- lm(f, data = d)
  r <- cbind(coef(m), confint(m))
  r[grepl(var, row.names(r)),] / sd_y
}

eth_cc <- sapply(measures, cc, simplify = F)
fsm_cc <- sapply(measures, cc, var = 'fsm', simplify = F)

x$y <- rep(12:1, each = 5) + rep(seq(0.2, -0.2, length.out = 5), 12)
x2 <- do.call(rbind, eth_cc)
x2 <- as.data.frame.matrix(x2)
colnames(x2) <- c('q', 'lower', 'upper')
off <- 0.05

png('mi_vs_cc.png', height = 15, width = 9, units = 'in', res = 300)

par(xpd = NA, mar = c(3, 10, 0, 0))
plot(1, type = 'n', xlim = c(-1, 0.5), ylim = c(0, 13), axes = F, xlab = NA, ylab = NA)
axis(1, pos = 0.5)

rect(-1, 0.5, 0.5, 12.5)
segments(0, 0.5, y1 = 12.5)
segments(-1, 0:12 + 0.5, x1 = 0.5, lty = 3)

points(x$q, x$y, pch = 19, col = cols)
arrows(x$lower, x$y, x1 = x$upper, angle = 90, code = 3, length = 0.05, col = cols)

points(x2$q, x$y - off, pch = 15, col = cols)
arrows(x2$lower, x$y - off, x1 = x2$upper, angle = 90, code = 3, length = 0.05, col = cols, lty = 3)

text(-1.05, 12:1, titles, adj = 1)

dev.off()
