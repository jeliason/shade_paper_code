#' Plot radial basis functions
#'
#' @param x vector of spatial distances at which to evaluate basis functions
#' @param rbfs list of radial basis functions
#'
#' @return ggplot2 object
#' @export
plot_rbfs <- function(x,rbfs) {
  d <- data.frame(x = x)
  
  # Calculate the values for each basis function at each point
  for (i in 1:length(rbfs)) {
    d[, paste0("rbf_", i)] <- rbfs[[i]](x)
  }
  
  # Reshape data for ggplot
  d <- tidyr::pivot_longer(d, cols = starts_with("rbf"), names_to = "rbf", values_to = "y")
  
  # Plot
  ggplot2::ggplot(d, ggplot2::aes(x = x, y = y, color = rbf)) +
    ggplot2::geom_line()
  
}