#' Assign a k-fold cross-validation number
#'
#' Randomly assigns a k-fold cross-validation number to each individual in a dataset.
#'
#' @param data_long Data frame in long format i.e. there may be more than one row per individual
#' @template individual_id
#' @template k
#' @param seed The value of the seed (default is 1)
#' @return Data frame `data_long` updated to contain a new column `cross_validation_number`
#' indicating the fold to which the individual has been assigned.
#' @author Isobel Barrott \email{isobel.barrott@@gmail.com}
#' @examples
#' data(data_repeat)
#' data_repeat <- add_cv_number(data_long = data_repeat,
#'                                       individual_id = "id",
#'                                       k = 10)
#' @details This function randomly divides the \code{n} individual IDs
#' into \code{k} groups, each with \code{n/k} members (or as close to this number as possible).
#' @export



add_cv_number <- function(data_long, individual_id, k, seed = 1) {
  set.seed(seed)
  if (!(is.data.frame(data_long))) {
    stop("data_long should be a data frame")
  }
  if (!(individual_id %in% colnames(data_long))) {
    stop("individual_id should be a column name in data_long")
  }
  ids <- data_long[[individual_id]]
  unique_ids <- unique(ids)
  n <- length(unique_ids)
  sample_size <- rep(round(n / k, 0), k)
  diff <- n - round(n / k, 0) * k
  if (diff != 0) {
    sample_size[1:abs(diff)] <- sample_size[1:abs(diff)] + sign(diff)
  }
  sample_size_end <- cumsum(sample_size)
  sample_size_start <- c(1, cumsum(sample_size)[-k] + 1)
  samples <- sample(1:n, size = n)
  df <- data.frame(cross_validation_number = rep(NA, n))
  df[[individual_id]] <- methods::as(unique_ids, class(ids))
  for (i in 1:k) {
    df[samples[sample_size_start[i]:sample_size_end[i]], "cross_validation_number"] <-
      i
  }
  df <- dplyr::left_join(data_long, df, by = individual_id)
  df
}
