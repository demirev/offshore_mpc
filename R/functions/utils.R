library(magrittr)
library(stringr)
library(dplyr)
library(readr)

# IO ----------------------------------------------------------------------
listFilesType <- function(dir = "data/HFCS_UDB_1_2_ASCII/", type = "H") {
  
  # lists all files in `dir`` that match the pattertn type+number.csv
  # to be uesed in lapply loops when keeping all data in memory all the time
  # is unwarranted
  
  datasets <- dir %>%
    list.files() %>%
    str_extract(paste0("^",type,"\\d\\.csv")) %>%
    na.omit %>%
    as.list
  
  names(datasets) <- dir %>%
    list.files() %>%
    str_extract(paste0("^",type,"\\d")) %>%
    na.omit
  
  return(datasets)
}

importType <- function(dir = "data/HFCS_UDB_1_2_ASCII/", type = "H") {
  
  # imports all files in `dir`` that match the pattertn type+number.csv
  # outputs a list of dataframes
  
  files <- listFilesType(dir, type)
  
  datasets <- files %>%
    lapply(function(fl){
      fl <- paste0(dir,fl)
      read_csv(fl)
    })
  
  names(datasets) <- names(files)
  
  return(datasets)
}

importTypes <- function(dir = "data/HFCS_UDB_1_2_ASCII/", 
                        types = c("H","PN")) {
  
  # wrapper around importType that imports multiple types at once
  # outputs a list of lists of the type of `importType`
  
  datasets <- types %>% 
    lapply(function(t){
      importType(dir, t)
    })
  
  names(datasets) <- types
  
  return(datasets)
}

importSurveys <- function(dir = "data/", folders = list.files("data/"), 
                          types = c("H", "PN")) {
  
  datasets <- folders %>%
    lapply(function(d){
      importTypes(paste0(dir,d,"/"), types)
    })
  
  names(datasets) <- folders
  
  return(datasets)
}


# MI and bootstrap --------------------------------------------------------
mi_point <- function(datasets, FUN) {
  
}

