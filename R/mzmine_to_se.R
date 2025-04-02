library("dplyr",include.only = ("%>%"))

#' Title
#' @importFrom utils read.csv
#' @importFrom dplyr %>% mutate group_by
#' @importFrom SummarizedExperiment SummarizedExperiment assays rowData colData
#'     assays<-
#' @param path_to_file path to the mzmine feature table
#' @param sample_meta_data data.frame containing sample metadata
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
    normalize = c(TRUE, FALSE),
    transform = c(TRUE, FALSE),
    vsn = c(TRUE,FALSE)
) {

  features <- read.csv(
    path_to_file
  ) %>%
    mutate(
      id = as.character(id) # convert id from int to character
    )

  feature_area <- features[,names(features)[grepl(paste0("[.]",intensity_type),names(features))]]

  names(feature_area) <- gsub(paste0("[.]",intensity_type), "", names(feature_area))
  names(feature_area) <- gsub(paste0("datafile[.]"), "", names(feature_area))

  se <- SummarizedExperiment(
    rowData = features[,c("id","rt","mz","ion_mobility")],
    assays = list(
      raw = feature_area
    ),
    colData = sample_meta_data[
      match(names(feature_area),make.names(sample_meta_data[,1])),
    ]
  )

  if (normalize) {
    assays(se)$normalized <- normalizePQN(assays(se)$raw)
  }


  if (transform) {
    assays(se)$transformed <- log2(assays(se)$normalized + 1)
  }

  if(vsn) {
    assays(se)$vsn <- limma::normalizeVSN(assays(se)$transformed)
  }

}


