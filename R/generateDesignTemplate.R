#' Generate design table template
#'
#' @param dataset path to the data exported from MS data processing file
#' @param inputDataType can be "spectronaut", "diann", or "fragpipe", output from other software are not supported (for this function)
#' 
#' @export
generateDesignTemplate <- function(dataset, inputDataType = "spectronaut"){
  
  # Load dataset if a file path is provided
  if(is.character(dataset) & length(dataset)==1 & inputDataType %in% c("spectronaut", "diann", "fragpipe"))
    dataset <- read.delim(dataset, header=T, row.names=NULL, stringsAsFactors=F, na.strings=c(""," ","NA","N.A.","NaN","n.a.", "Filtered", "na", "#N/A", "Fil", "#DIV/0!"))
  
  # Load dataset if a file path is provided
  if(is.character(dataset) & length(dataset)==1 & inputDataType %in% c("maxquant", "openswath"))
    dataset <- read.delim(dataset, header=T, row.names=1, stringsAsFactors=F, na.strings=c(""," ","NA","N.A.","NaN","n.a.", "Filtered", "na", "#N/A", "Fil", "#DIV/0!"))
  
  if(inputDataType == "spectronaut"){
    
    # Define the required columns for Spectronaut data
    required_columns <- c("PG.ProteinGroups", "EG.PrecursorId")
    
    # Check if all required columns are present in the dataset
    if(!all(required_columns %in% colnames(dataset))) {
      missing_columns <- required_columns[!required_columns %in% colnames(dataset)]
      stop(paste("The following required columns are missing from the dataset based on selected input data type:", paste(missing_columns, collapse = ", ")))
    }
    
    raw_file_names <-  dataset %>% 
      dplyr::select(ends_with(".EG.Channel1Quantity")) %>%
      dplyr::rename_with( ~gsub("(.*).EG.Channel1Quantity", "Intensity.L.\\1", .), .cols = ends_with("Channel1Quantity")) %>%
      dplyr::rename_with(~gsub("Intensity.L.", "", .), .cols = everything()) %>%
      colnames()

  }else if (inputDataType == "diann" | inputDataType == "fragpipe"){
    
    # Define the required columns for Spectronaut data
    required_columns <- c("Protein.Group", "Precursor.Id")
    
    # Check if all required columns are present in the dataset
    if(!all(required_columns %in% colnames(dataset))) {
      missing_columns <- required_columns[!required_columns %in% colnames(dataset)]
      stop(paste("The following required columns are missing from the dataset based on selected input data type:", paste(missing_columns, collapse = ", ")))
    }
    
    
    raw_file_names <- dataset %>%
      dplyr::select(ends_with(".L", ignore.case = FALSE)) %>%
      dplyr::rename_with(.cols = ends_with(".L"), ~gsub("(.*).L", "Intensity.L.\\1", .)) %>%
      dplyr::rename_with(~gsub("Intensity.L.", "", .), .cols = everything()) %>%
      colnames()
  }
  
  n_samples <- length(raw_file_names)
  
  # make a design table to be exported and modified in excel
  design <- data.frame(
    raw_file = raw_file_names, 
    condition = vector(mode = "character", length = n_samples), 
    sample = vector(mode = "character", length = n_samples), 
    time = vector(mode = "numeric", length = n_samples), 
    replicate = vector(mode = "numeric", length = n_samples), 
    color = vector(mode = "character", length = n_samples)) 
  
  write.table(
    design, 
    file = "design_table_template.tsv",       
    sep = "\t",                      
    row.names = FALSE,                
    quote = FALSE                    
  )
  
  message(paste(Sys.time(), "The design template has been generated and save in the working folder.", sep = " "))
  
}