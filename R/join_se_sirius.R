#' join_se_sirius
#'
#' A convenience function that joins a SIRIUS output file to a
#' SummarizedExperiment object.
#'
#' @importFrom SummarizedExperiment rowData SummarizedExperiment
#' @importFrom dplyr left_join
#' @importFrom utils read.delim
#'
#' @param object a SummarizedExperiment object
#' @param path_to_sirius path to the SIRIUS output file that should be joined
#'
#' @returns a SummarizedExperiment object with the SIRIUS data joined to the
#' rowData
#' @export
#'
join_se_sirius <- function(
    object = NULL,
    path_to_sirius = NULL
) {
  canopus_structure_summary <- utils::read.delim(
    path_to_sirius
  ) %>%
    dplyr::rename(
      "id" = .data$mappingFeatureId
    ) %>%
    dplyr::mutate(
      id = as.character(.data$id)
    )
  new_rowData <- dplyr::left_join(
    base::as.data.frame(
      SummarizedExperiment::rowData(object)
    ),
    canopus_structure_summary
    )
  SummarizedExperiment::rowData(object) <- new_rowData
  return(object)
}
utils::globalVariables(c(".",".data"))
