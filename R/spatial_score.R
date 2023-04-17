#' @title Calculate cell level spatial scores
#' @description The column names of a given data are hard coded, further changes are needed if data format is changed
#' @param data a dataframe which includes the merged results for CD30, CD4, CD68, CD20, CXCR5 and PD1
#' @import data.table
#' @return a dataframe which contains the spatial score at cell level
#' @author Alex Xu, Aixiang Jiang, Yifan Yin
#' @export

spatial_score <- function(data){

  ## AJ: found the problem: the output from merge_markerFolders is in data.frame
  ##     but here, data should be in data.table
  data = data.table::as.data.table(data)

  # add cell types based on marker positivity
  # CD30+ HRS cells
  data <- dplyr::mutate(data,CD30pos_HRS_pos = dplyr::if_else(CD30 == 'POS' & `Entire_Cell_Area_(square_microns)` > 40, 1,0))


  ## CXCR5+CD30+ HRS
  data <- dplyr::mutate(data,CXCR5pos_HRS_pos = dplyr::if_else(CD30 =='POS' & CXCR5 =='POS' & CD4 =='Other' & CD20 == 'Other' & CD68 == 'Other' & `Entire_Cell_Area_(square_microns)` > 40, 1,0))

  ## CD68+ CD4+ CD30- Macrophage
  data <- dplyr::mutate(data,CD68pos_Mac_pos = dplyr::if_else(CD68 == 'POS' & CD30 =='Other'& CD4 == 'POS', 1,0))


  ## CD4+PD1+ CD30 CD68-
  data <- dplyr::mutate(data, PD1pos_T_pos = dplyr::if_else(CD4 == 'POS' & PD1 == 'POS' & CD68 == 'Other', 1,0))


  ## CXCR5+CD20+CD30-CD68-CD4-
  data <- dplyr::mutate(data,CXCR5pos_B_pos = dplyr::if_else(CD20 == 'POS' & CXCR5 == 'POS' & CD4 =='Other' & CD30 =='Other' & CD68 == 'Other', 1,0))

  # number of NN
  NNtarget <- 5

  # define the cell types
  targetpositives <- c('CD30pos_HRS','CXCR5pos_HRS','CD68pos_Mac','PD1pos_T','CXCR5pos_B')
  # markers
  panel_proteins <- c('CD68','CD20','CD30','CXCR5','CD4','PD1')

  # keep the relevant data
  PointDatacols <- c("Cell_X_Position", "Cell_Y_Position", "Cell_ID",

                     panel_proteins, paste0(targetpositives,'_pos'))

  # keep the location for each cell (X Y and object number)
  PointDatalocation <- c("Cell_X_Position", "Cell_Y_Position", "Cell_ID")

  for (ROIImage in unique(data$Sample_Name)){
    print(ROIImage)
    ROIsubset <- data
    # Extract ROI data
    #ROIsubset <- ROIsubset[Sample_Name==ROIImage,..PointDatacols]
    ## AJ: the above line caused multiple problems within package, change to the following
    ROIsubset <- ROIsubset[which(data$Sample_Name == ROIImage), PointDatacols, with = FALSE]

    # Make dummy cells

    ROIdummy <- data.table::copy(ROIsubset[1:5])

    ROIdummy[,Cell_ID := 100000000:100000004]

    ROIdummy[, which(colnames(ROIdummy) %like% "_pos") := 1]

    ROIdummy[, which(colnames(ROIdummy) %like% "Position") := 100000000:100000004]

    ROIsubset <- rbind(ROIsubset, ROIdummy)


    # convert the dt into a ppp obj for each cell type i
    X <- spatstat.geom::ppp(x = ROIsubset[,Cell_X_Position],
             y = ROIsubset[,Cell_Y_Position],
             window = spatstat.geom::owin(c(ROIsubset[,min(Cell_X_Position)], ## AJ: spatstat.geom:: was missed
                             ROIsubset[,max(Cell_X_Position)]),
                           c(ROIsubset[,min(Cell_Y_Position)],
                             ROIsubset[,max(Cell_Y_Position)])),
             unitname = "micrometer")

    # Prepare dataset for cell subset
    NNsubi <- ROIsubset[,..PointDatalocation]

    # Compare to each cell type j

    for(j in targetpositives){

      # Generate points for targeted cell types
      index_j <- ROIsubset[[paste0(j,"_pos")]]==1

      # convert targeted cell types into ppp obj
      Y <- spatstat.geom::ppp(x = ROIsubset[index_j,Cell_X_Position],
               y = ROIsubset[index_j,Cell_Y_Position],
               window = spatstat.geom::owin(c(ROIsubset[,min(Cell_X_Position)], ## AJ: spatstat.geom:: was missed
                               ROIsubset[,max(Cell_X_Position)]),
                             c(ROIsubset[,min(Cell_Y_Position)],
                               ROIsubset[,max(Cell_Y_Position)])),
               unitname = "micrometer")

      # Find spatial nearest neighbors of i in j for the targeted cell types
      NNdataj <- as.data.table(spatstat.geom::nncross(X,Y,k=1:NNtarget, # AJ: spatstat.geom:: was missed
                                       iX=as.integer(ROIsubset$Cell_ID),
                                       iY=as.integer(ROIsubset[index_j]$Cell_ID)))

      # Find cells by object number instead of index
      NNdataj[, which(names(NNdataj) %like% "which") := lapply(.SD, function(x) ROIsubset[index_j][x,Cell_ID]), .SDcols = names(NNdataj) %like% "which"]

      # Name the columns for cell type j
      colnames(NNdataj) <- paste0(j, ".", colnames(NNdataj))

      # Store data
      NNsubi <- cbind(NNsubi, NNdataj)
    }

    # calcualte the centroids for each cell
    Xcentroids <- vector("list",length(targetpositives))
    Ycentroids <- vector("list",length(targetpositives))
    NNcentroiddistances <- vector("list", length(targetpositives))

    # For every cell type
    for (cluster in 1:length(targetpositives))
    {
      # Identify the columns containing spatial data for type
      coltargets <- names(NNsubi)[names(NNsubi) %like% paste0("^",targetpositives[cluster],".which")]
      NNminiX <- NNsubi[,..coltargets]
      NNminiY <- NNsubi[,..coltargets]

      # Replace ObjectNumber with X/Y coordinates
      NNminiX[,(coltargets) := lapply(.SD, function(x) ROIsubset[match(x, Cell_ID),Cell_X_Position]), .SDcols=coltargets]
      NNminiY[,(coltargets) := lapply(.SD, function(x) ROIsubset[match(x, Cell_ID),Cell_Y_Position]), .SDcols=coltargets]

      # calculate the centroids and save it into the list
      Xcentroids[[cluster]] <- NNminiX[,rowMeans(NNminiX, na.rm=T)]

      Ycentroids[[cluster]] <- NNminiY[,rowMeans(NNminiY, na.rm=T)]

      # Calculate distance from each cell to each centroid for each cell type

      NNcentroiddistances[[cluster]] <- sqrt((Xcentroids[[cluster]]-NNsubi$Cell_X_Position)^2 + (Ycentroids[[cluster]]-NNsubi$Cell_Y_Position)^2)
    }

    # Store new distance data for centroids

    newnncentroidcols <- apply(expand.grid(c("NN.centroid.X.","NN.centroid.Y.","NN.centroiddist."), targetpositives), 1, paste, collapse="")

    idx <- order(c(seq_along(Xcentroids), seq_along(Ycentroids), seq_along(NNcentroiddistances)))

    NNsubi[, (newnncentroidcols) := (c(Xcentroids, Ycentroids, NNcentroiddistances))[idx]]

    # Substitute spatial calculations back in to ROI subset data

    spatialcolumns <- names(NNsubi)[(names(NNsubi) %like% "dist") | (names(NNsubi) %like% "NN")]

    spatialdata <- NNsubi[,..spatialcolumns]

    spatialdata <- spatialdata[1:(nrow(spatialdata)-NNtarget)]

    data[Sample_Name==ROIImage, (spatialcolumns) := spatialdata]

  }


  cap=50

  for (k in targetpositives){

    print(k)

    regexin <- paste0("centroiddist[.]",k,"$|","^",k,"[.]")

    distancecols <- colnames(data)[grep(regexin, colnames(data))]

    censoreddistance <- censor(data[,..distancecols],cap=cap)

    # scaledcols <- paste0(names(censoreddistance),".scaled",cap)

    # data[,(distancecols) := (censoreddistance)/cap]

    data[,(distancecols) := (censoreddistance)/(cap)]


  }


  for (target in targetpositives){

    regexin <- paste0("^",target,"[.]")

    distancecols <- colnames(data)[grep(regexin, colnames(data))]

    data[, paste0(target, "_spatial_score") := (NNtarget-rowSums(.SD))/NNtarget, .SDcols=distancecols]
    # AJ: to be consistent to current naming system, change "summary_metric" to "spatial_score:"

  }

  return(data)
}
