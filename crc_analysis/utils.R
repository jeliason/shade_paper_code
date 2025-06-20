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

convert_tma_metadata_final <- function(pt_data) {
  pt_data %>%
    rename(
      LA_A = TMA_A_LA,
      LA_B = TMA_B_LA,
      Diffuse_A = TMA_A_Diffuse,
      Diffuse_B = TMA_B_Diffuse
    ) %>%
    separate(`TMA spot`, into = c("spot1", "spot2"), sep = ",", convert = TRUE) %>%
    pivot_longer(
      cols = c("spot1", "spot2"),
      names_to = "which_spot",
      values_to = "spot_number"
    ) %>%
    mutate(
      spot_number = as.integer(spot_number),
      spot_A = paste0(spot_number, "_A"),
      spot_B = paste0(spot_number, "_B"),
      Category_A = case_when(
        LA_A == 2 ~ 1,
        Diffuse_A == 2 ~ 2,
        LA_A == 1 & which_spot == "spot1" ~ 1,
        Diffuse_A == 1 & which_spot == "spot2" ~ 2,
        TRUE ~ NA_real_
      ),
      Category_B = case_when(
        LA_B == 2 ~ 1,
        Diffuse_B == 2 ~ 2,
        LA_B == 1 & which_spot == "spot1" ~ 1,
        Diffuse_B == 1 & which_spot == "spot2" ~ 2,
        TRUE ~ NA_real_
      )
    ) %>%
    pivot_longer(
      cols = c(spot_A, spot_B, Category_A, Category_B),
      names_to = c("what", "copy"),
      names_pattern = "(spot|Category)_(.)",
      values_to = "value",
      values_transform = list(value = as.character)
    ) %>%
    pivot_wider(names_from = what, values_from = value) %>%
    mutate(
      Category = as.integer(Category)
    ) %>%
    filter(!is.na(Category)) %>%
    rename(Spot = spot) %>%
    mutate(
      Patient = as.integer(factor(Patient)),
      Group = as.integer(factor(Group)),
      Category = as.integer(Category)
    ) %>%
    select(Spot, Patient, Group, Category)
}