#####################################################
#   Code to extract observation level data from LOTUS
#   Requires approved NDA access 
#   This version uses Free TDS as the driver
#
######################################################
library(DBI)
library(odbc)
library(glue)
library(dplyr)

working_dir <- "set your working dir"
setwd(working_dir)

driver.name <- "/opt/homebrew/opt/freetds/lib/libtdsodbc.so"
db.name<-"LOTUS"
port<-"11433"
user.name<-Sys.getenv("LODS_USER")
pwd<-Sys.getenv("LODS_PASS") 
host.name<-Sys.getenv("LODS_LOC")

con <- dbConnect(
  odbc::odbc(),
  Driver = driver.name,
  Server = host.name,
  Port = port,
  Database = db.name,
  UID = user.name,
  PWD = pwd,
  TDS_Version = "auto",
  Protocol = "TCPIP",
  Encrypt = "yes",
  TrustServerCertificate = "no",
  ConnectionTimeout = 30
)

as_catch <- dbGetQuery(
  con,
  "
  SELECT *
  FROM LOTUS.NEWOBS.LDS_CATCH_AS_V
  "
)

as_gear <- dbGetQuery(
  con,
  "
  SELECT *
  FROM LOTUS.NEWOBS.LDS_GEAR_CFG_AS_V
  "
)

DBI::dbDisconnect(con)

tmp <- as_catch %>%
  left_join(
    as_gear,
    by = c("TRIP_NUM", "SET_NUM")
  )

keep_cols <- c("TRIP_NUM", "SET_NUM", "PERMIT_NUM.x", "HAUL_BEGIN_DTM.x", 
        "MAINLN_LEN_RPTD", "BRNCHLN_LEN", "NUM_FLTS", "HKS_PER_FLT",
        "NUM_HKS_SET", "FLT_NUM", "HK_NUM", "ENGLISH_NAME.x")

new_names <- c("TRIP_NUM", "SET_NUM", "PERMIT_NUM", "HAUL_BEGIN_DTM", 
        "MAINLN_LEN_RPTD", "BRNCHLN_LEN", "NUM_FLTS", "HKS_PER_FLT",
        "NUM_HKS_SET", "FLT_NUM", "HK_NUM", "ENGLISH_NAME")       

species <- c("ALBACORE", "OCEANIC WHITETIP SHARK", "GREEN SEA TURTLE", "WAHOO", "MAHIMAHI")   

as_data <- tmp %>%
  filter(ENGLISH_NAME.x %in% species) %>%
  select(all_of(setNames(keep_cols, new_names)))

####Filter for missing data and add flags and fixes 
reg_date <- as.POSIXct("2011-09-23 00:00:00")

as_data <- as_data %>%
  mutate(
    # Determine if 70 m had to be added
    REG_TURTLE = if_else(
      HAUL_BEGIN_DTM > reg_date,
      70L,
      0L
    ),
    # convert miles to meters
    MAINLN_LEN_RPTD_M = MAINLN_LEN_RPTD * 1852,
    # Fill with grand mean 
    MAINLN_LEN_RPTD_M = coalesce(
      MAINLN_LEN_RPTD_M,
      mean(MAINLN_LEN_RPTD_M, na.rm = TRUE)
    ),
    #fill with grand mean
    NUM_FLTS = coalesce(
       NUM_FLTS,
      mean(NUM_FLTS, na.rm = TRUE)  
    ),
    #fill with grand mean
    HKS_PER_FLT = coalesce(
       HKS_PER_FLT,
      mean(HKS_PER_FLT, na.rm = TRUE)  
    ),
    HK_NUM = if_else(
      HK_NUM <= HKS_PER_FLT,
      HK_NUM,                                    # True: Do nothing
      ((HK_NUM - 1) %% HKS_PER_FLT) + 1      # False: Use the MOD function
    ),
    #estimate meters of main line per basket
    CALC_BSK_LEN = MAINLN_LEN_RPTD_M / (NUM_FLTS  - 1)
  )


save(as_data, file = "as_data.Rdata" )
