#' @title Generate pSILAC class object
#'
#' @description Creates a pSILAC object. Needs a design table and precursor-level report from an DIA-MS data processing software such as
#' Spectronaut, DIA-NN, and Fragpipe. 
#' These can be provided as path or as data loaded in R environment. 
#'
#' @param dataset A tab-separated file or `data.frame` in R containing precursor-level data.
#' 
#' If provided as a file, `dataset` must include:
#' - Columns for light and heavy precursor intensities
#' - A unique precursor ID column
#' - A protein ID column
#' 
#' If provided as a `data.frame` in R, `dataset` must:
#' - Use precursor IDs as row names
#' - Include columns for light and heavy intensities
#' 
#' @param design A tab-separated file or `data.frame` in R defining the experimental design.
#' 
#' - **File format**: The first column must be the raw file name. The file should also contain:
#'   - A `sample` column for sample identifiers
#'   - A `time` column in numeric format (hours)
#'   - Optional columns:
#'     - `replicate` column (numeric), for replicate IDs if a replicate design analysis is desired
#'     - `color` column to customize data plotting
#' - For replicate samples to be averaged, they must have the same value in the `sample` column.
#' 
#' - **Data.frame format**: If loaded as a `data.frame` in R, `design` must use raw file names as row names.
#' - **Template**: Use `generate_design_template()` to create a customizable design table template.
#' 
#' @param inputDataType A character string specifying the type of input data. Options include `spectronaut`, `diann`, `fragpipe`, `maxquant`, or `openswath`.
#' 
#' @param aggregate.replicates Method for aggregating replicate values:
#' - Options: `median` (recommended) or `mean`. Set to `NA` if no replicates are present or should not be aggregated.
#' - For replicate averaging, all replicates should have the same value in the `sample` column of `design`.
#' 
#' @param filterPeptides Logical; if `TRUE`, removes peptides without any lysine (K) or arginine (R) residues from the analysis (recommended).
#' 
#' @param noiseCutoff Numeric; specifies a cutoff for light and heavy intensities. If `NULL`, no cutoff is applied.
#' 
#' @param ncores Integer; the number of cores to use for parallel computation. Defaults to `1`.
#' 
#' @param requant Relevant for OpenSwath data only. Specifies handling of requantified values with a score greater than `0.05`.
#' - Options:
#'   - `keep`: Retain requantified values
#'   - `remove`: Exclude requantified values (default)
#'   - `impute`: Replace requantified values with imputed values
#' - **Data requirements**: Assumes columns are prepended with `Intensity_` and `score_` for each sample/timepoint in `dataset`.
#' 
#' @import dplyr
#' @importFrom parallel makeCluster stopCluster detectCores
#' @importFrom stringr str_detect
#' 
#' @return a pSILAC object.
#'
#' @export
generatepSILACObject <- function(dataset, design, inputDataType = "spectronaut", requant="remove", aggregate.replicates=NA, filterPeptides=T, ncores=1, imputeMethod="normD1", noiseCutoff = 8){
  
  if (is.character(dataset) && length(dataset) == 1) {
    
    if (inputDataType %in% c("spectronaut", "diann", "fragpipe")) {
      
      message(paste(Sys.time(), "Reading input data...", sep = " "))
      
      # For Spectronaut, DIA-NN, and FragPipe formats, do not use row names
      dataset <- read.delim(dataset, header = TRUE, row.names = NULL, stringsAsFactors = FALSE,
                            na.strings = c("", " ", "NA", "N.A.", "NaN", "n.a.", "Filtered", "na", "#N/A", "Fil", "#DIV/0!"))
      
    } else if (inputDataType %in% c("maxquant", "openswath")) {
      
      message(paste(Sys.time(), "Reading input data...", sep = " "))
      
      # For MaxQuant and OpenSWATH formats, treat the first column as row names
      dataset <- read.delim(dataset, header = TRUE, row.names = 1, stringsAsFactors = FALSE,
                            na.strings = c("", " ", "NA", "N.A.", "NaN", "n.a.", "Filtered", "na", "#DIV/0!"))
    }
  } else {
    message(paste(Sys.time(), "Using existing data frame as dataset...", sep = " "))
  }
  
  # Load design data if a file path is provided
  if(is.character(design) & length(design)==1)
    design <- read.delim(design, header=T, row.names=1, stringsAsFactors=F)
  
  # Check if required columns exist in design
  if(!all(c("time", "sample") %in% colnames(design)))
    stop("The design data.frame should (at least) have 'sample' and 'time' columns.")
  
  # Ensure the 'time' column is numeric
  if(!is.numeric(design$time))
    stop("The 'time' column of the design data.frame should be numeric.")
  
  # Aggregate replicates if specified
  if(!is.na(aggregate.replicates))
    aggregate.replicates <- match.arg(aggregate.replicates, c("median", "mean"))
  
  # Error if replicates exist without proper specification
  if(is.na(aggregate.replicates) && (length(unique(table(design$sample))) != 1 || length(unique(table(design$time))) != 1 ||
                                     !all(table(paste(design$sample, design$time)) == 1))) {
    stop("If working with replicates, set 'aggregate.replicates'. Otherwise, ensure each sample has one entry per time point.")
  }
  
  # Sort design data by sample and time
  design <- design[order(design$sample, design$time),]
  
  if(inputDataType == "spectronaut"){
    
    # Define the required columns for Spectronaut data
    required_columns <- c("PG.ProteinGroups", "EG.PrecursorId")
    
    # Check if all required columns are present in the dataset
    if(!all(required_columns %in% colnames(dataset))) {
      missing_columns <- required_columns[!required_columns %in% colnames(dataset)]
      stop(paste("The following required columns are missing from the dataset based on selected input data type:", paste(missing_columns, collapse = ", ")))
    }
    
    message(paste(Sys.time(), "Processing Spectronaut data...", sep = " "))
    
    dataset <- dataset %>% 
      dplyr::rename("Proteins" = PG.ProteinGroups) %>%
      dplyr::rename("id" = EG.PrecursorId) %>%
      dplyr::rename_with(~gsub("(.*).EG.Channel1Quantity", "Intensity.L.\\1", .), .cols = ends_with("Channel1Quantity")) %>%
      dplyr::rename_with( ~gsub("(.*).EG.Channel2Quantity", "Intensity.H.\\1", .), .cols = ends_with("Channel2Quantity")) %>%
      dplyr::select(id, Proteins, dplyr::starts_with("Intensity"))
    
    row.names(dataset) <- dataset$id
    
  } else if (inputDataType == "diann" | inputDataType == "fragpipe"){
    
    
    # Define the required columns for Spectronaut data
    required_columns <- c("Protein.Group", "Precursor.Id")
    
    # Check if all required columns are present in the dataset
    if(!all(required_columns %in% colnames(dataset))) {
      missing_columns <- required_columns[!required_columns %in% colnames(dataset)]
      stop(paste("The following required columns are missing from the dataset based on selected input data type:", paste(missing_columns, collapse = ", ")))
    }
    
    message(paste(Sys.time(), "Processing DIA-NN data...", sep = " "))
    
    dataset <- dataset %>% 
      dplyr::rename("Proteins" = Protein.Group) %>%
      dplyr::rename("id" = Precursor.Id) %>%
      dplyr::rename_with(.cols = ends_with(".L"), ~gsub("(.*).L", "Intensity.L.\\1", .)) %>%
      dplyr::rename_with(.cols = ends_with(".H"), ~gsub("(.*).H", "Intensity.H.\\1", .)) %>%
      dplyr::select(id, Proteins, dplyr::starts_with("Intensity."))
    
    row.names(dataset) <- dataset$id
  }else {
    message(paste(Sys.time(), "Using dataset as provided without additional processing...", sep = " "))
  }
  
  # Identify if dataset is in the expected format
  isStandardFormat <- c("Proteins" %in% colnames(dataset), length(grep("heavy", row.names(dataset)))==0,
                        length(grep("H.", colnames(dataset), fixed=T))>0 & length(grep("L.", colnames(dataset), fixed=T))>0)
  
  if(!all(isStandardFormat) & any(isStandardFormat)) stop("Could not identify the type of data.")
  isStandardFormat <- all(isStandardFormat)
  
  # Process standard formatted data
  if(isStandardFormat) {
    
    # Attempt to locate the columns for heavy and light intensities using a standard naming convention
    snh <- paste("Intensity.H.", row.names(design), sep="")
    snl <- paste("Intensity.L.", row.names(design), sep="")
    
    # If the expected naming is not found in dataset columns, try an alternative naming convention
    if(!all(c(snh %in% colnames(dataset), snl %in% colnames(dataset)))){
      
      # Alternative naming convention for heavy and light intensities
      snh <- paste("H.", row.names(design), sep="")
      snl <- paste("L.", row.names(design), sep="")
      
      # If neither naming convention is found, stop with an error indicating missing intensities
      if(!all(c(snh %in% colnames(dataset), snl %in% colnames(dataset))))
        stop("Could not find the light and heavy intensities for all samples of the design data.frame.")
    }
    
    # subset data to select columns with light intensities
    channel_1 <- dataset[,snl]
    
    # subset data to select columns with heavy intensities
    channel_2 <- dataset[,snh]
    
    # add "heavy" to the row names of the channel 2 data
    row.names(channel_2) <- paste(row.names(channel_2), "heavy", sep="")
    
    # rename the light and heavy data using the raw file names
    colnames(channel_1) <- row.names(design)
    colnames(channel_2) <- row.names(design)
    
    # combine light and heavy intensities
    e <- rbind(channel_1, channel_2)
    
    rm(channel_1, channel_2)
    
  } else {
    # Process Openswath format
    if(!all(row.names(design) %in% colnames(dataset))){
      if(!all(paste("Intensity", row.names(design), sep="_") %in% colnames(dataset))){
        stop(paste("Could not find the intensity for all samples of the design data.frame. Missing samples:",
                   paste(head(row.names(design)[which(!(paste("Intensity", row.names(design), sep="_") %in% colnames(dataset)))]), "...", collapse=", ")
        ))
      }
      e <- dataset[,paste("Intensity", row.names(design), sep="_")]
      colnames(e) <- gsub("^Intensity_", "", colnames(e))
      requant <- match.arg(requant, c("remove", "keep", "impute"))
      if(requant == "impute") stop("Imputation not yet implemented")  
      if(requant != "keep" & !all(paste("score", row.names(design), sep="_") %in% colnames(dataset))){
        warning("Missing scores for samples; using all intensities without filtering.")
        requant <- "keep"
      }
      if(requant == "remove") {
        s <- dataset[, paste("score", row.names(design), sep="_")]
        colnames(s) <- gsub("^score_", "", colnames(s))
        s <- s[colnames(e)]
        e[s > 0.05] <- NA
        message(paste("Replaced", sum(s > 0.05, na.rm=T), "(", round(100 * sum(s > 0.05, na.rm=T) / (ncol(s) * nrow(s))), "% ) badly requantified datapoints with NAs..."))
        rm(s)
      }
      e <- e[, row.names(design)]
    } else {
      e <- dataset[, row.names(design)]
    }
  }
  
  # Impute missing values if selected
  if(requant == "impute") e <- imputeAll(e, paste(design$condition, design$time), imputeMethod)
  
  # Filter out peptides without "K" or "R" if specified
  if(filterPeptides){
    w <- unique(c(grep("K", row.names(e)), grep("R", row.names(e))))
    if(length(w) < nrow(e)) {
      message(paste(nrow(e) - length(w), "peptides without K or R were discarded."))
      e <- e[w,]
    }
  }
  
  # Handle replicate aggregation if enabled
  if(!is.na(aggregate.replicates)){
    nnames <- paste(design$sample, design$time, sep=".")
    if(length(unique(nnames)) == length(nnames)){
      message("Replicate aggregation activated, but no replicates found.")
    } else {
      message(paste("Aggregating replicates using", aggregate.replicates, "..."))
      design <- aggregate(design, by=list(name=nnames), FUN=function(x){ x[[1]] })
      o <- order(design$sample, design$time)
      design <- design[o,]
      row.names(design) <- design$name
      e2 <- matrix(0, nrow=nrow(e), ncol=length(unique(nnames)))
      if(is.null(ncores)){
        library(parallel)
        ncores <- detectCores() - 1
      } else {
        if(ncores > 1) library(parallel)
      }
      if(ncores > 1){
        library(parallel)
        cl <- makeCluster(ncores)
        clusterExport(cl, c("nnames", "e", "aggregate.replicates"), environment())
        e <- parSapply(cl, unique(nnames), FUN=function(x){ apply(e[, which(nnames == x)], 1, na.rm=T, FUN=aggregate.replicates) })
        stopCluster(cl)
      } else {
        e <- sapply(unique(nnames), FUN=function(x){ apply(e[, which(nnames == x)], 1, na.rm=T, FUN=aggregate.replicates) })
      }
      e <- e[, o]
    }
  }
  
  # Finalize design columns
  if(!("name" %in% colnames(design))) design$name <- paste(design$sample, design$time, sep=".")
  
  # rename the columns using the sample name in the deign table
  colnames(e) <- design$name
  
  # set a default color if none is specified
  if(!("color" %in% colnames(design))) design$color <- "black"
  
  # Create output object with data and settings
  object <- list(
    design = design,
    info = list(creationDate = date(), requant = requant, aggregate.replicates = aggregate.replicates, protk.method = "", pSILAC.call = match.call(), protk.call = NULL),
    light = e[grep("heavy", row.names(e), invert=T),],
    heavy = e[grep("heavy", row.names(e), invert=F),],
    peptides = .getPeptideDefinition(dataset, isStandardFormat),
    RIA = NULL,
    hol = NULL,
    NLI = NULL,
    RIA.kloss = NULL,
    hol.kloss = NULL,
    NLI.kloss = NULL, 
    protein.kloss = NULL, 
    protein.kdeg = NULL,
    protein.halflife = NULL
  )
  
  # Additional cleanup for standard data types
  if(isStandardFormat){
    
    # Apply filtering to the data based on noiseCutoff if it's not NULL
    if(!is.null(noiseCutoff)) {
      
      object$light[object$light < noiseCutoff] <- NA
      object$heavy[object$heavy < noiseCutoff] <- NA
    }
    
    # Identify rows in 'light' and 'heavy' matrices where not all values are NA
    allNA <- c(row.names(object$light)[apply(object$light, 1, FUN=function(x){ !all(is.na(x)) })],
               row.names(object$heavy)[apply(object$heavy, 1, FUN=function(x){ !all(is.na(x)) })])
    
    # Filter 'light' and 'heavy' matrices to retain only rows with non-NA values found in 'rr'
    object$light <- object$light[row.names(object$light) %in% allNA, ]
    object$heavy <- object$heavy[row.names(object$heavy) %in% allNA, ]
    
    # Filter peptides to retain only those corresponding to the filtered rows in 'rr'
    object$peptides <- object$peptides[allNA, ]
  }
  
  # Update row names and calculate RIA and hol metrics
  row.names(object$heavy) <- gsub("heavy", "", row.names(object$heavy), fixed=T)
  g <- row.names(object$heavy)
  g <- g[g %in% row.names(object$light)]
  
  # generate the RIA matrix for nls fitting
  object$RIA <- as.data.frame(object$light[g, ] / (object$heavy[g, ] + object$light[g, ]))
  
  # generate the log-transformed H/L ratio for linear fitting 
  object$hol <- as.data.frame(log(object$heavy[g, ] / object$light[g, ] + 1))
  
  # Define class and return object
  class(object) <- "pSILAC"
  
  # normalize light channels to generate the NLI matrix for nls fitting
  object <- normalizeLightChannel(object)
  
  return(object)
}
