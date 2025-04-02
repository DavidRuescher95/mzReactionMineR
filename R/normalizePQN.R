library(MatrixGenerics)

normalizePQN <- function(data, measure = c("median", "mean")) {
  
  raw_data <- as.matrix(data)
  
  # normalization_factor = median ( value_i_j/ median(feature_i) )
    # where value_i_n is the value of the i-th feature in the j-th sample
    # and median(feature_i) is the median of the i-th feature across all samples
  
  switch(
    measure,
    median = {
      normalization_factor <- apply(
        raw_data / rowMedians(raw_data, na.rm = TRUE), 2, median,
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