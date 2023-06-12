
#' @title Filter and summarize MC-IF cell-level data into patient-level
#' @description This function takes a data frame of cell-level spatial scores and filters out low quality cores as determined by patient information
#' in a specified xlsx file. The remaining cores are summarized into patient-level data for the markers CD30, CD4, CD68, CD20, and PD1.
#' The column names of a given data are hard coded, further changes are needed if data format is changed
#' @param data A data frame with spatial scores, which are returned by the spatial_score function.
#' @param patient The path to an xlsx file containing patient information, including related slide ID and low quality cores.
#' @param core The core used to calculate patient-level data.
#' @return A data frame with columns for patient ID, and marker signal for CD30, CD4, CD68, CD20, and PD1.
#' @importFrom magrittr %>%
#' @author Yifan Yin, Aixiang Jiang
#' @export

convert_patient <- function(data, patient, core = "A"){

  lib.info <- readxl::read_xlsx(path = patient, sheet = core) ## need to give sheet name

  # get low quality ID
  lib.info <- lib.info %>%
    tidyr::separate(Low, into = c('Low1','Low2'), sep = ',', remove = FALSE)
  Low.ID <-  c(lib.info$Low1, lib.info$Low2)
  Low.ID <- Low.ID[!is.na(Low.ID)]

  # transform the dataframe to make the Low IDs into one col
  lib.info <- lib.info %>%
    tidyr::pivot_longer(values_to = 'SlideID', cols = c('ID1','ID2'))

  # remove low quality IDs
  lib.info <- lib.info %>%
    dplyr::filter(is.na(Low)) %>%
    dplyr::select(-c(Low, Low1, Low2,name))

  # keep data by cores
  data <- data %>%
    dplyr::filter(grepl(paste0(core,'1'),Sample_Name))

  # Extract slideID from SampleName
  data <- data %>%
    tidyr::separate(Sample_Name, into = 'Sample_Name1', sep = ']', remove = FALSE) %>%
    tidyr::separate(Sample_Name1, into = c('Name1', 'Name2','Name3'), sep = ',') %>%
    tidyr::unite('SlideID',Name2:Name3, sep = '.') %>%
    dplyr::select(-Name1)

  # remove low quality slides
  data <- data %>%
    dplyr::filter(!SlideID %in% Low.ID)

  # add patient information
  data <- dplyr::left_join(data, lib.info, by = 'SlideID')

  # keep only HRS cells and relapse samples
  #### AJ note: for some reason, the following two rows only return 49 rows when I set core = 'B", which is not right
  # HRSCenter.data <- data %>%
  #   dplyr::filter(data$CD30pos_HRS_pos == 1)

  ## AJ: change to a different way
  HRSCenter.data = subset(data, data$CD30pos_HRS.dist.1 == 1)  ## now, for core = "B", it is 4247 rows, looks correct

  ## AJ: the following 5 chunk of code might be better replaced with a simple in-line or independent function
  ##      however, since they work, I keep as they are for now
  CD30pos_HRS_spatial_score_patient <- HRSCenter.data %>%
    dplyr::group_by(BCCA_ID) %>%
    dplyr::top_frac(0.1, wt = CD30pos_HRS_spatial_score) %>%
    dplyr::summarise(CD30pos_HRS_spatial_score = mean(CD30pos_HRS_spatial_score))

  CXCR5pos_HRS_spatial_score_patient <- HRSCenter.data%>%
    dplyr::group_by(BCCA_ID) %>%
    dplyr::top_frac(0.1, wt = CXCR5pos_HRS_spatial_score) %>%
    dplyr::summarise(CXCR5pos_HRS_spatial_score = mean(CXCR5pos_HRS_spatial_score)) %>%
    dplyr::select(-BCCA_ID)

  CD68pos_Mac_spatial_score_patient <- HRSCenter.data %>%
    dplyr::group_by(BCCA_ID) %>%
    dplyr::top_frac(0.1, wt = CD68pos_Mac_spatial_score) %>%
    dplyr::summarise(CD68pos_Mac_spatial_score = mean(CD68pos_Mac_spatial_score))%>%
    dplyr::select(-BCCA_ID)

  PD1pos_T_spatial_score_patient <- HRSCenter.data %>%
    dplyr::group_by(BCCA_ID) %>%
    dplyr::top_frac(0.1, wt = PD1pos_T_spatial_score) %>%
    dplyr::summarise(PD1pos_T_spatial_score = mean(PD1pos_T_spatial_score))%>%
    dplyr::select(-BCCA_ID)

  CXCR5pos_B_spatial_score_patient <- HRSCenter.data%>%
    dplyr::group_by(BCCA_ID) %>%
    dplyr::top_frac(0.1, wt = CXCR5pos_B_spatial_score) %>%
    dplyr::summarise(CXCR5pos_B_spatial_score = mean(CXCR5pos_B_spatial_score))%>%
    dplyr::select(-BCCA_ID)

  MC_IF.patient <- cbind(CD30pos_HRS_spatial_score_patient,
                       CXCR5pos_HRS_spatial_score_patient,
                       CD68pos_Mac_spatial_score_patient,
                       PD1pos_T_spatial_score_patient,
                       CXCR5pos_B_spatial_score_patient)

  return(MC_IF.patient)
}
