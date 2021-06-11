#' Select individuals in a dataset with a last observation carried forward (LOCF) at a particular time
#'
#' To fit the LOCF model, all individuals must have at least one
#' non-`NA` entry by age `x_L` for all covariates. This function selects these individuals and removes the other rows.
#'
#' @details Individuals have a LOCF if there is a non-`NA` entry for each of the covariates in
#' `covariates` up until (and including) time `x_L`.
#'
#' @param data Data frame with repeated measurements in long format; each row corresponds to an assessment entry
#' @template patient_id
#' @template x_L
#' @template covariates
#' @template covariates_time
#' @return Data frame `data` updated to contain only rows of individuals with a LOCF at age `x_L`, other rows are removed
#' @author Isobel Barrott \email{isobel.barrott@@gmail.com}
#' @examples data(data_repeat_outcomes)
#' data_repeat_outcomes<-return_ids_with_LOCF(data=data_repeat_outcomes,
#'   patient_id="id",
#'   covariates=c("ethnicity","smoking","diabetes","deprivation",
#'   "atrial_fibrillation","sbp_stnd","tchdl_stnd"),
#'   covariates_time=c(rep("response_time_sbp_stnd",6),"response_time_tchdl_stnd"),
#'   x_L=60)
#' @export


return_ids_with_LOCF <-
  function(data,
           patient_id,
           x_L,
           covariates,
           covariates_time) {
    for (col in c(covariates,
                  patient_id,
                  covariates_time)) {
      if (!(col %in% names(data))) {
        stop(col, " is not a column name in data")
      }
    }

    if (!(length(covariates_time) %in% c(length(covariates),1))){
      stop("Length of covariates_time should be equal to length of covariates or 1")}

    if (length(covariates_time)==1){covariates_time<-rep(covariates_time,times=length(covariates))}

    #Pick out IDs with one LOCF before Ln for each exposure
    LOCF_IDs_by_variable <-
      lapply(1:length(c(covariates)), function(i) {
        var <- c(covariates)[i]
        time <- c(covariates_time)[i]
        data_var <- data[data[[time]] <= x_L,]
        data_var <-
          data_var[!is.na(data_var[[var]]),]
        data_var <-
          data_var[order(data_var[[time]], decreasing = TRUE),]
        data_var <-
          data_var[!duplicated(data_var[[patient_id]]),]
        return(data_var[[patient_id]])
      })
    LOCF_IDs_by_variable <- Reduce(intersect, LOCF_IDs_by_variable)
    data <- data[data[[patient_id]] %in% LOCF_IDs_by_variable,]
    data
  }
