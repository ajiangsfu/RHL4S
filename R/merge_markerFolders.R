
#' @title Merge IHC data into a single data frame
#' @description
#' Please notice that the subfolders' names and content are hard coded, further changes are needed if data names and formats are changed
#'
#' @param path the location where the IHC data is stored. Under the path, there should be 6 separate folders which contain the IHC output
#'        for each markers. The name for each subfolder is the marker name e.g. CD30, CD20 and etc.
#' @param markerFolders the IHC raw data folder needed for RHL4S, you could have other subfolders, but at least these subfolders should be included
#' @return
#'        a dataframe which contains the marker signal for CD30, CD4, CD68, CD20 and PD1
#' @author Yifan Yin, Aixiang Jiang
#'
#' @export
#'
merge_markerFolders<- function(path, markerFolders = c("CD30", "CD4", "CD68", "CD20","PD1","CXCR5")){

  ## AJ: check to see if I have all subfolders I need, if not, return error message and stop
  markerFolders = intersect(markerFolders, c("CD30", "CD4", "CD68", "CD20","PD1","CXCR5"))
  if(length(markerFolders) < 6){
    return("Error: please provide a valid path and 6 subfolders named as CD30, CD4, CD68, CD20, PD1, CXCR5")
  }

  # get the txt file names in each marker(the file names are identical in all single marker folders)
  file.in.dir <- list.files(path = paste0(path, '/', markerFolders[1]))

  ## AJ: this function contains a lot of hard-coding, change is needed if the data file name and content format are changed

  # extract the segmentation files in CD30
  seg.data <- file.in.dir[grep(']_cell_seg_data.txt',file.in.dir)]
  # extract slide ID
  sample.data <- data.frame(IHC_ID = seg.data)

  # merge the data into a single dataframe
  all.data <- data.frame()
  for (i in 1:nrow(sample.data)) {
    # read summary table
    print(paste0("Processing sample: ", sample.data$IHC_ID[i]))
    s.data <- data.table::data.table()
    for (m in markerFolders) {
      ## AJ: the following "\\" did not work, change to "/"
      #m.cellTablePath <- paste0(path,'\\',  m, "\\", sample.data$IHC_ID[i]). ## this line caused error
      m.cellTablePath <- paste0(path,'/',  m, "/", sample.data$IHC_ID[i])
      m.data <- data.table::fread(m.cellTablePath)

      m.data <- cbind(m.data,
                      Marker_pos = ifelse(m.data$Phenotype == paste0(m, "+"), "POS", "Other"))
      colnames(m.data)[ncol(m.data)] <- m
      m.data <- m.data[,c("Sample Name", "Tissue Category", "Cell X Position", "Cell Y Position", "Entire Cell Area (square microns)", m,
                          grep(paste0(m, ".*Mean"), colnames(m.data), value = TRUE)), with = FALSE]
      colnames(m.data) <- gsub(" ", "_", gsub(" \\(Normalized Counts, Total Weighting\\)", "", gsub(" \\(Opal \\d+\\)", "", colnames(m.data))))

      # Remove rows with duplicated coordinates, since they can't be distinguished as individual cells
      m.comb <- paste(m.data$Cell_X_Position, m.data$Cell_Y_Position, sep = "_")
      m.dup <- m.comb[which(duplicated(m.comb))]
      if (length(m.dup) > 0) {
        print(paste0("Note: ", length(m.dup), " cells with duplicated coordinates removed"))
        m.data <- m.data[-which(m.comb %in% m.dup),]
      }

      if (nrow(s.data) < 1) {
        s.data <- m.data
      } else {
        s.data <- merge(s.data, m.data)
        if (nrow(s.data) != nrow(m.data)) {
          print(paste0("Warning: some non-overlapping x/y coordinates identified, ", nrow(m.data) - nrow(s.data), " cells dropped"))
        }
      }

      rm(m.cellTablePath, m.data, m.comb, m.dup)
    }
    s.data <- cbind(s.data, "Cell_ID" = as.character(rownames(s.data)))
    if (nrow(all.data) < 1) { all.data <- s.data } else { all.data <- rbind(all.data, s.data) }
    rm(s.data)
  }
  return(all.data)
}

