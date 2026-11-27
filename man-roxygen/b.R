#' @param b  Integer specifying the number of bootstrap samples to take when calculating standard error of c-index and Brier score.
#' The bootstrap resamples the individuals used to assess performance while holding the fitted models and their predictions fixed,
#' so the standard errors do not include the uncertainty from fitting the models. Use \code{set.seed()} beforehand for reproducible standard errors.
