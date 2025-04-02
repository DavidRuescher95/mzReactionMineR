library("dplyr",include.only = ("%>%"))

#' Align peaks across multiple samples
#'
#' This function aligns peaks across multiple samples based on m/z and retention time (rt) tolerances.
#' It is a simplified implementation of mzMine's Join Aligner algorithm.
#' mzmine.github.io/mzmine_documentation/module_docs/align_join_aligner/join_aligner.html
#'
#' @importFrom dplyr %>% select mutate bind_rows sym row_number if_else filter
#'     ungroup group_by arrange mutate_if
#' @importFrom rlang :=
#' @param input A named list of data frames, each containing peak information with columns for id, rt, and mz.
#' @param mz_tolerance Numeric,
#' the tolerance for m/z alignment window.
#' @param rt_tolerance Numeric,
#' the tolerance for retention time alignment window.
#' @param mz_weight Numeric,
#' the weight for m/z in the alignment score calculation.
#' @param rt_weight Numeric,
#' the weight for retention time in the alignment score calculation.
#' @param id_col Character,
#' the name of the column containing peak IDs.
#' @param rt_col Character,
#' the name of the column containing retention times.
#' @param mz_col Character,
#' the name of the column containing m/z values.
#' @return A data frame with aligned peaks across all samples.
#' The data frame will have columns for id, rt, mz.
#' The original id for each sample is stored in separate columns.
#'
#'

join_aligner <- function(
  input,
  mz_tolerance = 0.05,
  rt_tolerance = 0.1,
  mz_weight = 3,
  rt_weight = 1,
  id_col = "id",
  rt_col = "rt",
  mz_col = "mz"
) {
  # Check if input is a list
  if (!is.list(input)) {
    stop("Input must be a list")
  }

  # Ensure input is a named list
  if (is.null(names(input))) {
    warning("Input must be a named list\n renaming...")
    names(input) <- paste("Sample", seq_along(input), sep = "_")
  }
  #

  # Convert only columns that can be converted to numeric
  input <- lapply(input, function(df) {
    df[, c(id_col, rt_col, mz_col)] %>%
      mutate_if(is.character, as.numeric)
  })

  # Order input by the number of peaks in each sample (descending)
    # technically unncecessary
  new_order <- names(sort(unlist(lapply(input, nrow)), decreasing = TRUE))
  input <- input[new_order]

  # Initialize master list with the first sample
  master_list <- input[[1]]

  # Iterate over each sample
  for (sample_id in seq_along(input)) {
    sample <- input[[sample_id]]
    sample_name <- names(input)[sample_id]

    # Iterate over each peak in the master list
    for (i in seq_len(nrow(master_list))) {
      # Define the alignment window (AW)

      master_mz <- master_list$mz[i]
      master_rt <- master_list$rt[i]

      mz_range <- c(master_mz - mz_tolerance, master_mz + mz_tolerance)
      rt_range <- c(master_rt - rt_tolerance, master_rt + rt_tolerance)

      # Find peaks within the alignment window
      candidate_peaks <- sample %>%
        filter(mz >= mz_range[1] & mz <= mz_range[2] & rt >= rt_range[1] & rt <= rt_range[2])

      if (nrow(candidate_peaks) > 0) {
        # Compute scores for each candidate peak
        scores <- rowSums(cbind(
          (1 - abs(candidate_peaks$mz - master_mz) / mz_tolerance) * mz_weight,
          (1 - abs(candidate_peaks$rt - master_rt) / rt_tolerance) * rt_weight
        ))

        # Find the peak with the best score
        best_match <- candidate_peaks[which.max(scores), ]
        best_score <- max(scores)

        # Align the best match peak
        master_list[i, "mz"] <- mean(c(master_mz, best_match$mz))
        master_list[i, "rt"] <- mean(c(master_rt, best_match$rt))
        master_list[i, sample_name] <- best_match[[id_col]]
        master_list[i, "score"] <- best_score


      }
    }

    master_list <- master_list %>%
      group_by(!!sym(sample_name)) %>%
      mutate(
        score = replace(.data$score, is.na(.data$score), 0),
        !!sym(sample_name) := if_else(.data$score == max(.data$score, na.rm = TRUE), !!sym(sample_name), NA_real_)
      ) %>%
      ungroup() %>%
      select(-.data$score)  # Remove the score column

    # Find missing features in the master list
    missing_features <- sample %>%
      filter(!id %in% master_list[[sample_name]])

    # Add missing features to the master list
    if (nrow(missing_features) > 0) {
      missing_features[[sample_name]] <- missing_features[[id_col]]
      master_list <- bind_rows(master_list, missing_features) %>%
        arrange(rt, mz) %>%
        mutate(id = row_number())
    }


  }

  return(master_list)
}

utils::globalVariables(c(".data", "."))
