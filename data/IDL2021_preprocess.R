## This script takes the IDL data and creates
## two data files: one with the French records
## (who died aged 105 and above between 1987 
## and 2017) and one with the rest of the IDL
## supercentenarians
## 

library(lubridate) # data manipulation
library(tidyverse)
# TODO set working directory to data location
# setwd("")
# The following CSV file can be downloaded
# from supercentenarians.org upon registering
# - see Supplementary material for instructions
idl_raw <- read.csv(file = "idl_complete.csv",
                    header = TRUE,
                    sep = ";")
# Load truncation bounds and identifiers of individuals
# used in the paper - these may differ from metadata
# as new data is added
url_root <- "https://raw.githubusercontent.com/lbelzile/supercentenarian/master/auxiliary/"
truncation_bounds_tab <- paste0(url_root, "IDL2021_truncation_bounds.txt")
id_IDL2021_tab <- paste0(url_root, "IDL2021_id.txt")
id_France_tab <- paste0(url_root, "France_id.txt")
id_IDL2021 <- unlist(read.table(id_IDL2021_tab,
                          header = FALSE))
id_France <- unlist(read.table(id_France_tab,
                         header = FALSE))
# Create table with countries and
# dates of the data collection periods
sampling_frame <-
  read.table(file = truncation_bounds_tab,
             header = TRUE) %>%
  as_tibble %>%
  mutate(country = factor(country),
         c1 = ymd(c1),
         c2 = ymd(c2))
# c1: earliest day for collection period (110+ years)
# c2: latest day for collection period (110+ years)

# Tidy IDL database
idl_full <- idl_raw %>% # convert to tibble
  as_tibble %>% # select specific columns
  select(IDNUMBER,
         AGEDAYS,
         AGEYEARS,
         SEX,
         BDATE,
         DDATE,
         UPDATE,
         VALIDATION) %>% #match names in script
  rename(id = IDNUMBER,
         ndays = AGEDAYS,
         gender = SEX,
         bdate = BDATE,
         ddate = DDATE,
         wave = UPDATE,
         ageyear = AGEYEARS,
         valid = VALIDATION) %>%
  mutate(country = factor(substr(id, 1, ifelse(substr(id, 1, 2) == "UK", 2, 3))),
         valid = factor(valid),
         gender = factor(gender)) %>%
  # Create country manually b/c one record is . (AUSTRIA)
    filter(country != "JPN") # exclude Japan

# Next, we split the analysis
# We first select only French and
# remove records of individuals who died
# before 1987 - this is done
# automatically by selecting the
# records in the list of IDs


# French semi-supercentenarians and
# supercentenarians
u_fr <- 38350L
idl_france <- idl_full %>%
  filter(country == "FRA") %>%
  mutate(bdate = dmy(bdate),
         ddate = dmy(ddate)) %>%
  filter(ddate >= dmy("01-01-1987")) %>%
  mutate(rtrunc = as.numeric(dmy("31-12-2017") - bdate),
         ltrunc = pmax(u_fr, dmy("01-01-1987") - bdate)
  ) %>%
  mutate(gender = fct_recode(gender,
                             "female"  = "F",
                             "male" = "M"))
write.csv(idl_france,
          file = "france.csv",
          row.names = FALSE)
francent <- idl_france %>%
  rename(birth = bdate,
         numdays = ndays) %>%
  select(birth, numdays, gender)
save(francent, file = "francent.rda", version = 2)


# For IDL2021, we have to account for the
# fact that USA doesn't record dates, only
# years of birth/death.
# After having computed the truncation
# bounds separately for this country, we
# can merge the two database and purge the IDs

idl2021_rest <-
  idl_full %>%
    filter(! country  %in% c("FRA","USA"),
           id %in% id_IDL2021) %>%
  mutate(bdate = dmy(bdate),
         ddate = dmy(ddate))

# We have the number of days survived, but not
# the dates of death/birth
# There is partial information about these
idl2021_USA <- idl_full %>%
  filter(country == "USA",
         ageyear >= 110) %>%
  rename(byear = bdate,
         dyear = ddate) %>%
  mutate(bdate = pmax(dmy(paste0("01-01-", byear)),
                      dmy(paste("01-01-", dyear)) - ndays),
         ddate = pmin(dmy(paste0("31-12-", dyear)),
                      dmy(paste("31-12-", byear)) + ndays)) %>%
    select(! c(dyear, byear))
# If a person died in 2010, the earliest they could
# have died is January 1st, which gives a lower bound
# for their birth date, etc.

idl2021 <- bind_rows(idl2021_rest,
                     idl2021_USA) %>%
  filter(ageyear >= 110,
         id %in% id_IDL2021) %>%
  select(! valid) %>% # all validated (flag)
  mutate(country = factor(country, exclude = NULL)) %>%
  # remove unused factor level
  left_join(y = sampling_frame, by = c("country")) %>%
  # # Add date at which people reached 110 years
  mutate(x110 = if_else(mday(bdate) == 29 & month(bdate) == 2,
                        bdate + days(1) + years(110), #leap year and people born on February 29th
                        bdate + years(110)),
        ) %>%
  arrange(country, ndays) %>%
  mutate(gender = fct_recode(gender,
                             "female"  = "F",
                             "male" = "M")) %>%
  mutate(country = fct_recode(country,
                              "OS" = "AUT",
                              "BE" = "BEL",
                              "QC" = "CAN",
                              "DE" = "DEU",
                              "DN" = "DNK",
                              "ES" = "ESP",
                              "FI" = "FIN",
                              "FR" = "FRA",
                              "NO" = "NOR",
                              "SV" = "SWE",
                              "EW" = "UK",
                              "US" = "USA")) %>%
  mutate(ltrunc = pmax(if_else(country == "US",
                               40175L,
                               as.integer(x110 - bdate)),
                               as.integer(c1 - bdate)),
         rtrunc = as.integer(c2 - bdate))
  # Fix this operation for the US, which
  # doesn't have known birth/death dates
 

# Sanity checks: right truncation must be less than number of days
isTRUE(all(with(idl2021, rtrunc >= ndays)))
isTRUE(all(with(idl2021, ltrunc <= ndays)))
colnames(idl2021)
idl2021_csv <- idl2021 %>% select(! c(c1, c2, x110))
write.csv(idl2021_csv, # removed unused information,
          file = "idl2021.csv",
          row.names = FALSE)
u110 <- with(idl2021, min(ndays) - 1L)
idlex <- idl2021 %>%   
  transmute(slow = as.numeric(pmax(0, ltrunc - u110))/365.25,
            supp = as.numeric(rtrunc - u110)/365.25,
            datu = (ndays - u110)/365.25,
            gender = gender)
save(idlex, file = "IDL2021.rda", version = 2)

