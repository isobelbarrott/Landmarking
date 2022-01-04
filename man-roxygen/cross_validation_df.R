#' @param cross_validation_df List of data frames containing the cross-validation fold each individual is assigned to. Each data frame in the list should be
#' named according to the landmark time `x_L` that they correspond. Each data frame should contain the columns \code{individual_id} and a column \code{cross_validation_number}
#' which contains the cross-validation fold of the individual. An alternative to setting parameter `k` for performing cross-validation;
#' if both are missing no cross-validation is used.
