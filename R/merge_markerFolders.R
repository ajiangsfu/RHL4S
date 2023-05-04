#' @title Merge MC-IF data into a single data frame
#' @description This function reads in MC-IF data stored in six separate subfolders, one for each marker (e.g. CD4, CD20, CD30, CD68, CXCR5, and PD1), 
#' and merges them into a single data frame. Note that the function assumes that the folder structure and file names are in a specific format and any deviations
#'  will require changes to the code.
#' @param path The path to the folder containing the MC-IF data subfolders. Under the path, there should be 6 separate folders which contain the MC-IF output
#'        for each marker, and their names should be CD4, CD20, CD30, CD68, CXCR5, and PD1.
#' @param markerFolders A character vector of the marker subfolder names within the main MC-IF data folder. This should include all six marker subfolders.
#' @return A data frame containing the merged MC-IF data for all six markers.
#' @author Yifan Yin, Aixiang Jiang
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
  sample.data <- data.frame(MC_IF_ID = seg.data)

  # merge the data into a single dataframe
  all.data <- data.frame()
  for (i in 1:nrow(sample.data)) {
    # read summary table
    print(paste0("Processing sample: ", sample.data$MC_IF_ID[i]))
    s.data <- data.table::data.table()
    for (m in markerFolders) {
      ## AJ: the following "\\" did not work, change to "/"
      #m.cellTablePath <- paste0(path,'\\',  m, "\\", sample.data$MC_IF_ID[i]). ## this line caused error
      m.cellTablePath <- paste0(path,'/',  m, "/", sample.data$MC_IF_ID[i])
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

