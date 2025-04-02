library("dplyr",include.only = ("%>%"))

#' normalizePQN
#'
#' Probabilistsic quotient normalization (PQN) is a method to normalize data by
#'    dividing each value within a sample by a normalization factor.
#'
#'    \deqn{normalization\_factor = median \left( \frac{value_{i,j}}{median(feature_i)} \right)}
#'
#'    Where \eqn{value_{i,j}} is the value of the i-th feature in the j-th sample and
#'    \eqn{median(feature_i)} is the median of the i-th feature across all samples.
#'
#'    If measure = "mean", the normalization factor is calculated as:
#'
#'    \deqn{normalization\_factor = median \left( \frac{value_{i,j}}{mean(feature_i)} \right)}
#'
#'    The normalized matrix is then simply the quotient of a \eqn{value_{i,j}} by
#'    the \eqn{normalization_factor_{j}} of a given sample j.
#'
#'    \deqn{normalized\_data = \frac{value_{i,j}}{normalization\_factor_j}}
#'
#' @importFrom dplyr %>%
#' @importFrom MatrixGenerics rowMedians
#' @importFrom stats median
#'
#' @param data Input data to be normalized. Usually a numeric matrix.
#' @param measure What measure to be used for the normalization factor.
#'     Can be "median" or "mean".
#'
#' @returns A matrix containing normalized values.
#' @export
#

normalizePQN <- function(data, measure = c("median", "mean")) {

  raw_data <- as.matrix(data)

  # normalization_factor = median ( value_i_j/ median(feature_i) )
    # where value_i_n is the value of the i-th feature in the j-th sample
    # and median(feature_i) is the median of the i-th feature across all samples

  switch(
    measure,
    median = {
      normalization_factor <- apply(
        raw_data /rowMedians(raw_data, na.rm = TRUE), 2, median,
        na.rm = TRUE
      )
    },
    mean = {
      normalization_factor <- apply(
        raw_data / rowMeans(raw_data, na.rm = TRUE), 2, median,
        na.rm = TRUE
      )
    }
  )

  normalized_data <- apply(raw_data, 1, function(x) x/normalization_factor) %>% t()

  rownames(normalized_data) <- rownames(raw_data)

  return(normalized_data)

}
