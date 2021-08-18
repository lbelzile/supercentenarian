## Script to recreate `italcent.rda` datafile
## based on the Istat spreadsheet
# TODO set working directory to location of data file
# setwd('')
filename <- "ISTAT_Italian_SSC_1896-1910.xlsx"

library(lubridate)
library(readxl)
library(tidyverse)
#Import dataset and split in two (since years for cohorts pre-1899 are saved as DD/MM/YYYY and after that as MM/DD/YYYY
# Type of column entries (for all but the first column)
ctype <- c("date","numeric","numeric","text","numeric")
# Rename the columns to remove spaces, etc.
column_names <- c("birth_date",
                  "death_date",
                  "age",
                  "days",
                  "gender",
                  "cohort")
istat1 <- readxl::read_xlsx(filename,
                            col_names = column_names,
                            range = "A4:F26", 
                            col_types = c("text", ctype)) %>%
  mutate(birth_date = dmy(birth_date),
         death_date = ymd(death_date))
istat2 <- readxl::read_xlsx(filename,
                            col_names = column_names,
                            range = "A27:F3839",
                            col_types = c("date", ctype)) %>%
  mutate(birth_date = ymd(birth_date),
         death_date = ymd(death_date))
istat <- rbind(istat1, istat2) %>%
  mutate(rightcens = is.na(death_date),
         cohort = as.integer(cohort),
         gender = factor(gender, 
                         levels = c("f","m"),
                         labels = c("female","male")),
         numdays = as.integer(if_else(rightcens, 
                                      dmy("01-01-2016"),
                                      death_date) - birth_date)) %>%
  rename(birth = birth_date) %>%
  select(birth, gender, rightcens, numdays)
# Save file to rda format
save(italcent, 
     file = "italcent.rda",
     version = 2)
