#' get_colData
#'
#' Convenience function that extracts the colData from a SummarizedExperiment
#'     object as a data.frame
#'
#' @importFrom SummarizedExperiment colData SummarizedExperiment
#' @param object a SummarizedExperiment object
#'
#' @returns colData as data.frame
#' @export
#'
get_colData <- function(
    object = NULL
) {
  data <- base::as.data.frame(
    SummarizedExperiment::colData(object)
  )
  return(data)
}
