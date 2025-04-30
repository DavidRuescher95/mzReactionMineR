#' mass_difference_network
#'
#' Generates a mass difference network from a vector of m/z values.
#'
#' @param x a numeric vector of mz values. Required.
#' @param y a numeric vector of mz values. Optional.
#'     If not provided, the function will use x for both sets,
#'     generating a bipartite graph.
#'
#' @returns An adjacency matrix representing the mass difference network.
#' @export
#'
#'
mass_difference_network <- function(
    x,
    y = NULL
) {

  if(is.null(y)) {
    y <- x
  } else {
    print("Two sets provided. Generating a bipartite graph.")
  }

  # generate adjacency matrix

  adjacency_matrix <- matrix(0, nrow = length(y), ncol = length(x))

  # change colnames
  if(!is.null(names(x))) {
    colnames(adjacency_matrix) <- names(x)
  }

  # change colnames
  if(!is.null(names(y))) {
    rownames(adjacency_matrix) <- names(y)
  }

  # calculate the mass difference between each pair of points

  for(i in seq_along(y)) {

    adjacency_matrix[i,] <- abs(y[i] - x)

  }

  return(adjacency_matrix)

}
