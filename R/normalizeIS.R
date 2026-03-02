#' normalizeIS
#'
#' A function to normalize the intensities of an assay in a SummarizedExperiment
#'  object based on the intensities of an internal standard. The internal
#'  standard can be specified by its id, or by its rt and mz values. The
#'  function will find the closest feature to the provided rt and mz values, and
#'  use its intensities for normalization.
#'
#' @importFrom dplyr %>% select mutate filter slice_min pull between across
#' @importFrom SummarizedExperiment SummarizedExperiment assays rowData
#'
#' @param object A SummarizedExperiment object
#' @param assay Character. The name of the assay to be normalized.
#' @param rt Numeric. The theoretical retention time of the internal standard.
#' @param mz Numeric. The theoretical m/z of the internal standard.
#' @param id Character. The id of the internal standard. If "none", the function
#'        will find the closest feature to the provided rt and mz values.
#'        Default is "none".
#' @param mz_tolerance Numeric vector of length 2. Absolute and relative (ppm)
#'        m/z tolerance for finding the internal standard.
#' @param rt_tolerance Numeric. Absolute retention time tolerance for finding
#'        the internal standard.
#' @param new_assay_name Character. The name of the new assay that will contain
#'        the normalized intensities. Default is "is_normalized".
#' @param remove Logical. Whether to remove the internal standard from the
#'        object after normalization. Default is TRUE.
#' @param id_col Character. The name of the column in rowData that contains the feature ids.
#' @param rt_col Character. The name of the column in rowData that contains the retention times.
#' @param mz_col Character. The name of the column in rowData that contains the m/z values.
#'
#' @returns A SummarizedExperiment object with a new assay containing the
#'        normalized intensities.
#' @export
normalizeIS <- function(
    object,
    assay = NULL,
    rt = NULL,
    mz = NULL,
    id = "none",
    mz_tolerance = c(0.005, 10),
    rt_tolerance = 0.1,
    new_assay_name = "is_normalized",
    remove = TRUE,
    id_col = "id",
    rt_col = "rt",
    mz_col = "mz"
) {

  is_id <- get_is_id(
    object = object,
    rt = rt,
    mz = mz,
    id = id,
    mz_tolerance = mz_tolerance,
    rt_tolerance = rt_tolerance,
    id_col = id_col,
    rt_col = rt_col,
    mz_col = mz_col
  )

  is_intensities <- get_is_intensities(
    object = object,
    assay = assay,
    id = is_id
  )

  tmp <- object
  assays(tmp, withDimnames = FALSE)[[new_assay_name]] <- t(apply(
    t(assays(tmp)[[assay]]),2,function(col) {
      col/is_intensities
    }))

  if(remove) {
    tmp <- tmp[-which(rowData(tmp)[[id_col]] == is_id), ]
  }

  return(tmp)

}

# function to find the id of the internal standard

get_is_id <- function(
    object,
    rt = NULL,
    mz = NULL,
    id = "none",
    mz_tolerance = c(0.005, 10),
    rt_tolerance = 0.1,
    id_col = "id",
    rt_col = "rt",
    mz_col = "mz"
) {

  if(id != "none") {

    if(id %in% rowData(object)[[id_col]]) {
      is_id <- id
    } else {
      stop(paste0("The provided id (", id, ") is not present in the rowData of the object."))
    }

  } else {

    mz_range <- calc_mz_range(mz, mz_tolerance)
    rt_range <- calc_rt_range(rt, rt_tolerance)

    # find is_id
    is_id <- rowData(object) %>%
      as.data.frame() %>%
      filter(
        between(.data[[mz_col]], mz_range[1], mz_range[2]),
        between(.data[[rt_col]], rt_range[1], rt_range[2])
      ) %>%
      mutate(
        mz_diff = abs(.data[[mz_col]]-mz),
        rt_diff = abs(.data[[rt_col]]-rt)
      ) %>%
      slice_min(
        across(c(rt_diff,mz_diff)), n=1
      ) %>%
      pull(.data[[id_col]])

  }

  return(is_id)

}

# function to extract intensities of the internal standatd

get_is_intensities <- function(
    object,
    assay = NULL,
    id,
    id_col = "id"
) {

  values <- assays(object[which(rowData(object)[[id_col]] == id), ])[[assay]]
  return(values)

}

# function to calculate mz range based on absolute tolerance and ppm

calc_mz_range <- function(
    mz,
    tolerance = c(0.005, 10)
) {

  mz_tolerance <- max(mz_tolerance[1], (mz_tolerance[2] * (mz/1e6)))
  mz_range <- c(mz - mz_tolerance, mz + mz_tolerance)
  return(mz_range)

}

# function to calculate rt range based on absolute tolerance and ppm

calc_rt_range <- function(
    rt,
    tolerance = 0.1
) {

  rt_range <- c(rt - mz_tolerance, rt + mz_tolerance)
  return(rt_range)

}


utils::globalVariables(c(".data", ".", "mz_diff", "rt_diff"))
