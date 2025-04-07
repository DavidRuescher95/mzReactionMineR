#' mzmine_to_se
#'
#' Convert mzMine feature table to a SummarizedExperiment object. Assumes the
#'     "Export to CSV (modular)" was used
#'
#' @importFrom utils read.csv
#' @importFrom dplyr %>% mutate group_by
#' @importFrom SummarizedExperiment SummarizedExperiment assays rowData colData
#' @param path_to_file path to the mzmine feature table
#' @param sample_meta_data data.frame. sample meta data to become colData
#'     the first columnn has to be the sample names
#' @param intensity_type character. either "area" or "height"
#' @param rowData_cols character vector. column names to be used as rowData
#' @param remove_adducts logical. if TRUE, remove adducts from the feature
#'     table. Requires that the ion_identities.iin_id column is present in the
#'     data.
#'
#' @returns a summarizedExperiment Object
#' @export
#'

mzmine_to_se <- function(
    path_to_file,
    sample_meta_data,
    intensity_type = c("area","height"),
    rowData_cols = c("id","rt","mz","ion_mobility"),
    remove_adducts = FALSE
) {

  # load feature data

  print("Loading mzMine feature table")

  features <- read.csv(
    path_to_file
  ) %>%
    mutate(
      id = as.character(id) # convert id from int to character
    )

  # extract columns with height or area

  print("Processing feature table")

  if(remove_adducts) {
    features <- features[is.na(features$ion_identities.iin_id), ]
  }

  feature_area <- features[,names(features)[grepl(paste0("[.]",intensity_type),names(features))]]

  # rename to fit the file names in meta_data

  names(feature_area) <- gsub(paste0("[.]",intensity_type), "", names(feature_area))
  names(feature_area) <- gsub(paste0("datafile[.]"), "", names(feature_area))

  sample_meta_data[,1] <- make.names(sample_meta_data[,1])

  samples <- intersect(
    names(feature_area),
    sample_meta_data[,1]
  )


  # generate summarized experiment

  print("Creating SummarizedExperiment object")

  se <- SummarizedExperiment(
    rowData = features[,rowData_cols],
    assays = list(
      raw = as.matrix(feature_area[,samples])
    ),
    colData = sample_meta_data[
      match(samples, sample_meta_data[,1]),
    ]
  )

  return(se)

}


