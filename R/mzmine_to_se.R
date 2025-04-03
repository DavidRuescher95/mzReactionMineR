library("dplyr",include.only = ("%>%"))

#' Title
#' @importFrom utils read.csv
#' @importFrom dplyr %>% mutate group_by
#' @importFrom SummarizedExperiment SummarizedExperiment assays rowData colData
#'     assays<-
#' @importFrom limma normalizeVSN
#' @param path_to_file path to the mzmine feature table
#' @param path_to__meta_data data.frame with sample meta data
#'     the first columnn has to be the sample names
#' @param intensity_type either "area" or "height"
#' @param normalize logical. Whether to normalize the data
#' @param transform logical. Whether to transform the normalized data
#' @param vsn logical. Whether to perform variance stabilizing normalization
#'
#' @returns a summarizedExperiment Object
#' @export
#'

mzmine_to_se <- function(
    path_to_file,
    sample_meta_data,
    intensity_type = c("area","height"),
    normalize = TRUE,
    transform = TRUE,
    vsn = TRUE
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

  feature_area <- features[,names(features)[grepl(paste0("[.]",intensity_type),names(features))]]

  # rename to fit the file names in meta_data

  names(feature_area) <- gsub(paste0("[.]",intensity_type), "", names(feature_area))
  names(feature_area) <- gsub(paste0("datafile[.]"), "", names(feature_area))

  # generate summarized experiment

  print("Creating SummarizedExperiment object")

  se <- SummarizedExperiment(
    rowData = features[,c("id","rt","mz","ion_mobility")],
    assays = list(
      raw = as.matrix(feature_area)
    ),
    colData = sample_meta_data[
      match(names(feature_area),make.names(sample_meta_data[,1])),
    ]
  )

  # processing

  if (normalize) {

    print("Normalizing data")

    assays(se)$normalized <- normalizePQN(assays(se)$raw,"median")
  }


  if (transform) {

    print("log2(x+1) transformation of normalized data")

    assays(se)$transformed <- log2(assays(se)$normalized + 1)
  }

  if(vsn) {

    print("Variance stabilizing normalization of raw data")

    assays(se)$vsn <- normalizeVSN(assays(se)$transformed)

  }

  return(se)

}


