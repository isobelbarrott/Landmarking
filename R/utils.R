return_LOCF_by_variable <- function(data_long,
                                    i,
                                    covariates,
                                    covariates_time,
                                    individual_id,
                                    x_L) {
  var <- covariates[i]
  time <- covariates_time[i]
  data_var <-
    data_long[data_long[[time]] < x_L,][, c(individual_id, var, time)]
  data_var <-
    data_var[!is.na(data_var[[var]]),]#removes NA values
  data_var <-
    data_var[order(data_var[[time]], decreasing = TRUE), ]
  data_var <-
    data_var[!duplicated(data_var[[individual_id]]), ]
  return(data_var[, c(individual_id, var)])
}
