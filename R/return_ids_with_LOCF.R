#' Select individuals in a dataset with a last observation carried forward (LOCF) at a landmark time
#'
#' To fit the LOCF model, all individuals must have at least one
#' non-`NA` entry by landmark time `x_L` for all covariates.
#' This function selects these individuals and removes the other rows.
#'
#' @details Individuals have a LOCF if there is a non-`NA` entry for each of the covariates in
#' `covariates` up until (not including) time `x_L`.
#'
#' @param data_long Data frame with repeated measurements data in long format
#' @template individual_id
#' @template x_L
#' @template covariates
#' @template covariates_time
#' @return List of data frames which correspond to each landmark time `x_L`.
#' Each data frame is an updated version of `data_long` which contains only rows
#' of individuals with a LOCF at age `x_L`, other rows are removed.
#' @author Isobel Barrott \email{isobel.barrott@@gmail.com}
#' @examples
#' library(Landmarking)
#' data(data_repeat_outcomes)
#' data_repeat_outcomes <-
#'   return_ids_with_LOCF(
#'     data_long = data_repeat_outcomes,
#'     individual_id = "id",
#'     covariates =
#'       c("ethnicity", "smoking", "diabetes", "sbp_stnd", "tchdl_stnd"),
#'     covariates_time =
#'       c(rep("response_time_sbp_stnd", 4), "response_time_tchdl_stnd"),
#'     x_L = c(60,61)
#'   )
#' @export



return_ids_with_LOCF <-
  function(data_long,
           individual_id,
           x_L,
           covariates,
           covariates_time) {
    if (!(inherits(data_long,"data.frame"))) {
      stop("data_long should be a data frame")
    }
    if (!(inherits(individual_id,"character"))) {
      stop("individual_id should have class character")
    }
    if (!(inherits(covariates,"character"))) {
      stop("covariates should have class character")
    }
    if (!(inherits(covariates_time,"character"))) {
      stop("covariates_time should have class character")
    }

    for (col in c(covariates,
                  individual_id,
                  covariates_time)) {
      if (!(col %in% names(data_long))) {
        stop(col, " is not a column name in data_long")
      }
      # if(any(is.na(data_long[[col]]))){
      #   stop(col, " contains NA values")
      # }
    }

    if (!(length(covariates_time) %in% c(length(covariates), 1))) {
      stop("Length of covariates_time should be equal to length of covariates or 1")
    }

    if (length(covariates_time) == 1) {
      covariates_time <- rep(covariates_time, times = length(covariates))
    }
    if (!(inherits(x_L,"numeric"))) {
      stop("x_L should have class numeric")
    }

    out_list <- lapply(x_L, function(x_l) {
      #Pick out IDs with one LOCF before Ln for each exposure
      LOCF_IDs_by_variable <-
        lapply(1:length(c(covariates)), function(i) {
          var <- c(covariates)[i]
          time <- c(covariates_time)[i]
          data_var <- data_long[data_long[[time]] <= x_l, ]
          data_var <-
            data_var[!is.na(data_var[[var]]), ]
          data_var <-
            data_var[order(data_var[[time]], decreasing = TRUE), ]
          data_var <-
            data_var[!duplicated(data_var[[individual_id]]), ]
          return(data_var[[individual_id]])
        })
      LOCF_IDs_by_variable <-
        Reduce(intersect, LOCF_IDs_by_variable)
      data_long <-
        data_long[data_long[[individual_id]] %in% LOCF_IDs_by_variable, ]
      data_long
    })
    if (length(x_L) == 1) {
      out_list[[1]]
    } else{
      names(out_list) <- x_L
      out_list
    }
  }
