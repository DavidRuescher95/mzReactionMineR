#' se_to_long
#'
#' Generates a long formdata.frame from a summarized experiment produced by
#'     mzmine_to_se
#'
#' @importFrom dplyr %>% inner_join
#' @importFrom tidyr pivot_longer
#' @importFrom mzReactionMineR get_colData get_rowData
#' @importFrom tidyselect all_of
#' @importFrom SummarizedExperiment rowData SummarizedExperiment
#' @param object a SummarizedExperiment object
#' @param assay a character string indicating the assay to be used
#'
#' @returns a data.frame in long format
#' @export
#'
se_to_long <- function(
  object,
  assay
) {
  df <- cbind(
    mzReactionMineR::get_rowData(object), assays(object)[[assay]]
  ) %>%
    tidyr::pivot_longer(
      cols = -tidyselect::all_of(colnames(rowData(object))),
      names_to = "filename",
      values_to = "value"
    ) %>%
    dplyr::inner_join(
      mzReactionMineR::get_colData(object)
    )
  return(df)
}
