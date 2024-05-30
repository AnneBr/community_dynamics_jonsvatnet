library(tidyverse)
library(lubridate)
library(glmmTMB)

# Downloading data
download_url <- "https://api.gbif.org/v1/occurrence/download/request/0239745-200613084148143.zip"

temp <- tempfile()
download.file(url=download_url, destfile=temp,quiet=FALSE)
alldata <- read.table(unzip(temp, "occurrence.txt"),sep="\t",header=T, fileEncoding = "UTF-8")

write_rds(alldata, "alldata.rds")
alldata <- read_rds("alldata.rds")

# How many variables are empty, i.e. NA?
alldata_sum <- alldata %>% 
  # Check the total number of NAs in each column
  summarise_all(.funs = ~ sum(is.na(.))) %>% 
  # Pivot to long format
  pivot_longer(gbifID:iucnRedListCategory) %>% 
  # So we can keep only columns where the number of NAs is unequal to the
  # total number of records
  filter(value != nrow(alldata))

# Keep non-empty columns
alldata <- alldata %>% 
  dplyr::select(alldata_sum$name)

# How many variables have only one value? (Same procedure as above)
alldata_sum <- alldata %>% 
  summarise_all(.funs = ~ n_distinct(.)) %>% 
  pivot_longer(gbifID:level2Name) %>% 
  filter(value != 1)

alldata <- alldata %>% 
  dplyr::select(alldata_sum$name)

# Keep only "ACCEPTED" entries (sounds most reliable)
alldata <- alldata %>% 
  filter(taxonomicStatus == "ACCEPTED") %>% 
  dplyr::select(-taxonomicStatus)

# Non-interesting columns
alldata <- alldata %>% 
  dplyr::select(-issue, -fieldNumber)

# Take out modifie/identifie columns
alldata <- alldata %>% 
  dplyr::select(-colnames(alldata)[str_detect(colnames(alldata), "ifie")])

# Keep records with the same sampling protocol 
# THIS ONE SOMETIMES DOESN'T WORK, in that case here is a possible solution:
#1. Open R-studio
#2. In the R console, write file.edit(file.path("~", ".Rprofile"))
# the Rprofile file opens to the left. In it, you write
# Sys.setlocale("LC_ALL", "nb_NO.UTF8")
# save it, close everything and reopen R-studio. 
# Now, if you write "Sys.getlocale() it should be "nb_NO.UTF8" on all.
alldata <- alldata %>% 
  filter(samplingProtocol == unique(alldata$samplingProtocol)[1]) %>% 
  dplyr::select(-samplingProtocol)

# Keep records of the same quantity
#alldata <- alldata %>% 
#  filter(organismQuantity == "1/1") %>% 
#  dplyr::select(-organismQuantity)

# We do not need this
alldata <- alldata %>% 
  dplyr::select(-organismQuantityType)

# What are the different taxon ranks?
table(alldata$taxonRank)

# Keep those with "species" rank
alldata_sp <- alldata %>% 
  filter(taxonRank == "SPECIES")

# There is one locality and a few depths that do not contain much data, remove!
alldata <- alldata %>% 
  filter(locality != unique(alldata$locality)[4],
         !(depth %in% c(5, 10, 15)))

# Since we modified the "alldata" object, we remake the species filter
alldata_sp <- alldata %>% 
  filter(taxonRank == "SPECIES")

# Now things look better!
alldata_sp_sum <- alldata_sp %>% 
  group_by(locality, depth, species, eventDate) %>%
  summarise(count = n()) %>% 
  ungroup()

alldata_sp_tmp <- inner_join(filter(alldata_sp_sum,
                                    count > 1), alldata_sp) %>% 
  group_by(locality, depth, species, eventDate) %>% 
  summarise(individualCountdiff = abs(diff(individualCount))) %>% 
  ungroup()

summary(alldata_sp_tmp$individualCountdiff)

alldata_sp_sample <- alldata_sp %>% 
  group_by(locality, depth, species, eventDate) %>% 
  summarise(occurrenceID = sample(occurrenceID, 1)) %>% 
  ungroup()

alldata_sp_sample <- inner_join(alldata_sp, alldata_sp_sample)

# Add numFactor for day in sampling season (time), and nf_year
# glmmTMB needs numFactors

all_dall <- alldata_sp_sample %>% 
  mutate(new_date = date(paste("1900", month, day, sep = "-")),
         duration = as.duration(interval(min(new_date), new_date)),
         duration = day(as.period(duration, unit = "days")),
         time = numFactor(duration),
         nf_year = numFactor(year - min(year))) 
all_dall

# collated data and all species (so without depth)
data_agg <- all_dall %>% 
  filter(year > 1993) %>%  ## data from before 1994 is less consistent/more patchy
  mutate(new_date = date(paste("1900", month, day, sep = "-")),
         duration = as.duration(interval(min(new_date), new_date)),
         duration = day(as.period(duration, unit = "days")),
         time = numFactor(duration),
         nf_year = numFactor(year - min(year)))

data_agg <- data_agg %>% 
  group_by(locality, phylum, species, time, year, month, day, nf_year, depth) %>% 
  summarise(individualCount = sum(individualCount), .groups = "drop")

write_rds(data_agg, "data_agg_depth.rds")
