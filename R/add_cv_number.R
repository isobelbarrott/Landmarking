#' Assign a k-fold cross-validation number
#'
#' Randomly assigns a k-fold cross-validation number to each individual in a dataset.
#'
#' @param data Data frame; there may be more than one row per individual
#' @template patient_id
#' @param k Integer specifying the number of folds
#' @return Data frame `data` updated to contain a new column `cross_validation_number`
#' indicating the fold to which the individual has been assigned.
#' @author Isobel Barrott \email{isobel.barrott@@gmail.com}
#' @examples
#' data(data_repeat_outcomes)
#' data_landmark_cv<-add_cv_number(data=data_landmark, patient_id="id", k=10)
#' @export

add_cv_number <- function(data, patient_id, k) {
  if (!(is.data.frame(data))) {
    stop("data should be a data frame")
  }
  if (!(patient_id %in% colnames(data))) {
    stop("patient_id should be a column name in data")
  }
  ids <- data[[patient_id]]
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
  df[[patient_id]] <- methods::as(unique_ids, class(ids))
  for (i in 1:k) {
    df[samples[sample_size_start[i]:sample_size_end[i]], "cross_validation_number"] <-
      i
  }
  df <- dplyr::left_join(data, df, by = patient_id)
  df
}
