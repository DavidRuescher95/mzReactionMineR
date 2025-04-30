#' filterSE
#'
#' A function that filters a SummarizedExperiment object based on various
#'     measures.
#'
#' @importFrom dplyr %>% mutate filter group_by summarize ungroup pull sym n
#'     inner_join
#' @importFrom tidyr pivot_longer
#' @importFrom SummarizedExperiment SummarizedExperiment rowData assays colData
#' @importFrom tidyselect all_of
#'
#' @param object A SummarizedExperiment object
#' @param assay Assay name to filter on
#' @param sample_col Character. Name of the column in the colData that
#'     that corresponds to the assay columns
#' @param group_col Character. Either "none" or a columnname in the colData
#' @param not_in Character. Either "none" or a value in the group_col.
#'     Specify, if an id should not be present in a group.
#' @param min_abundance Numeric. Minimum abundance for a feature to be present
#'     in a sample
#' @param min_pct Numeric. Minimum percentage of samples a feature must be
#'     present in. If group_col is specified, the percentage is calculated for
#'     for each group.
#' @param mz_range Numeric. A vector of length 2. The range of kept m/z values.
#' @param rt_range Numeric. A vector of length 2. The range of kept rt values.
#' @param mobility_range Numeric. A vector of length 2.
#'     The range of kept ion mobility values.
#' @param id_col Character. The respective column name in the rowData.
#' @param rt_col Character. The respective column name in the rowData.
#' @param mz_col Character. The respective column name in the rowData.
#' @param ion_mobility_col Character. The respective column name in the rowData.
#'
#' @returns A filtered summarizeExperiment object
#' @export
filterSe_ims <- function(
    object = NULL,
    assay = NULL,
    sample_col = "filename",
    group_col = "none",
    not_in = "none",
    min_abundance = 0,
    min_pct = 0.8,
    mz_range = c(0,Inf),
    rt_range = c(0,Inf),
    mobility_range = c(0,Inf),
    id_col = "id",
    rt_col = "rt",
    mz_col = "mz",
    ion_mobility_col = "ion_mobility"
) {

  # get column names of rowData

  id_cols <- names(as.data.frame(rowData(object)))

  # make long data frame of assay for easy filtering

  data_long <- cbind(
    as.data.frame(rowData(object)),
    as.data.frame(assays(object)[[assay]]) %>%
      replace(is.na(.data),0) %>%
      replace(.data > min_abundance, 1)
  ) %>%
    pivot_longer(
      cols = -all_of(id_cols),
      names_to = sample_col,
      values_to = "Value"
    ) %>%
    inner_join(
      as.data.frame(colData(object))
    ) %>%
    filter(
      !!sym(mz_col) > mz_range[1] & !!sym(mz_col) < mz_range[2],
      !!sym(rt_col) > rt_range[1] & !!sym(rt_col) < rt_range[2],
      !!sym(ion_mobility_col) > mobility_range[1] &
        !!sym(ion_mobility_col) < mobility_range[2]
    )

  # remove features with values < min_pct

  if(group_col != "none") {

    if(all(not_in != "none")) {

      ids <- data_long %>%
        filter(
          !!sym(group_col) %in% not_in
        ) %>%
        group_by(!!sym(id_col)) %>%
        summarize(
          Value = sum(.data$Value),
          na.rm = TRUE
        ) %>%
        filter(
          Value == 0
        ) %>%
        ungroup() %>%
        pull(!!sym(id_col)) %>%
        unique()

      data_long_grouped <- data_long %>%
        filter(
          !!sym(id_col) %in% ids
        ) %>%
        group_by(!!sym(id_col), !!sym(group_col))

    } else {

      data_long_grouped <- data_long %>%
        group_by(!!sym(id_col), !!sym(group_col))

    }

  } else {

    data_long_grouped <- data_long %>%
      group_by(!!sym(id_col))

  }

  ids <- data_long_grouped %>%
    summarize(
      Value = sum(.data$Value),
      cutoff = floor(n()*min_pct),
      na.rm = TRUE
    ) %>%
    filter(
      Value >= .data$cutoff
    ) %>%
    ungroup() %>%
    pull(!!sym(id_col)) %>%
    unique()

  result <- object[rowData(object)[[id_col]] %in% ids,]

  return(result)

}

utils::globalVariables(c(".",".data","Value"))
