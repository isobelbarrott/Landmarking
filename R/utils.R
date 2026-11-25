return_LOCF_by_variable <- function(data_long,
                                    i,
                                    covariates,
                                    covariates_time,
                                    individual_id,
                                    x_L) {
  var <- covariates[i]
  time <- covariates_time[i]
  data_var <-
    data_long[data_long[[time]] <= x_L, ][, c(individual_id, var, time)]
  data_var <-
    data_var[!is.na(data_var[[var]]), ]#removes NA values
  data_var <-
    data_var[order(data_var[[time]], decreasing = TRUE),]
  data_var <-
    data_var[!duplicated(data_var[[individual_id]]),]
  return(data_var[, c(individual_id, var)])
}

# Return an argument used to fit a landmark model. Models fitted with
# Landmarking >= 1.0.3 store their resolved arguments in `args`; for older
# objects fall back to evaluating the stored call.
get_landmark_arg <- function(landmark_model, name, default = NULL) {
  if (!is.null(landmark_model$args) && name %in% names(landmark_model$args)) {
    return(landmark_model$args[[name]])
  }
  value <- landmark_model$call[[name]]
  if (is.null(value)) {
    return(default)
  }
  eval(value)
}
