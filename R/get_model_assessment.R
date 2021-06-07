#' Compute C-index and Brier score
#'
#' Performs model assessment by computing the C-index and Brier score at time `x_hor`. There is the option
#' to calculate their standard errors using bootstraping.
#'
#' @param data Data frame containing survival outcomes and the event predictions from the model,
#' there should be one row for each individual
#' @template event_status
#' @template event_time
#' @template x_hor
#' @template patient_id
#' @param event_prediction Character string specifying the column name in `data` containing the prediction of the event of interest
#' @param return_c_index Boolean indicating whether the C-index should be calculated
#' @param return_brier_score Boolean indicating whether the Brier score should be calculated
#' @param b Integer specifying the number of bootstrap samples to take
#' @param standard_error Boolean indicating whether to calculate the standard error
#' @return List containing C-index, Brier score and their standard errors
#' @details There are two factors in assessing the performance of a prediction model; its
#' discrimination and its calibration. The c-index is one method to assess
#' discrimination, this refers to the ability of the model to separate individuals into
#' those that will have an event and those that will not. The c-index at a horizon time `x_hor`
#' looks at the pairs of individuals where one individual has the event at a time T and the other has not had the event at time T.
#' It is calculated as the proportion of these pairs where their relative risk prediction agrees with the
#' actual outcomes for the two individuals. This is extended to the competing risks case
#' by comparing individuals where one had the event of interest at time T and the other individual either
#' did not experience the event before this time T or experienced a competing event.
#'
#' The Brier score is one method to assess calibration, this refers to the agreement between the risk prediction and
#' the outcome. The Brier score is calculated as the average mean squared error of the predicted risk and the event outcome (where
#' an event is 1 and not experiencing the event is 0). This is extended to the competing risks case by including the competing risk events as
#' not experiencing the event.
#'
#' For both the c-index and Brier score calculations, inverse probability censoring weighting (IPCW) is used to create weights
#' which account for the occurence of censoring. The censoring model assumes for this function is the Kaplan Meier model, i.e. censoring occurs
#' independently of covariates.
#'
#' The c-index is calculated using the `cindex` function in package `pec`. The Brier score is calculated using
#' `pec` function in package `pec`.
#'
#'
#' @author Isobel Barrott \email{isobel.barrott@@gmail.com}
#' @examples \dontrun{data(data_landmark_LOCF)
#' get_model_assessment(data=data_landmark_LOCF,
#'   patient_id="id",
#'   event_prediction="event_prediction",
#'   event_status="event_status",
#'   event_time="event_time",
#'   x_hor=65,
#'   return_c_index = TRUE,
#'   return_brier_score = TRUE,
#'   b=10,
#'   standard_error = TRUE)}
#' @export

get_model_assessment <-
  function(data,
           patient_id,
           event_prediction,
           event_status,
           event_time,
           x_hor,
           return_c_index = TRUE,
           return_brier_score = TRUE,
           b,
           standard_error = FALSE) {
    # for (col in c(event_prediction,
    #               event_status,
    #               event_time)) {
    #   if (!(col %in% names(data))) {
    #     stop(col, " is not a column name in data")
    #   }
    # }
    # if (!(is.numeric(x_hor))) {
    #   stop("x_hor should be numeric")
    # }
    # if (!(is.numeric(b))) {
    #   stop("b should be numeric")
    # }
    # if (length(unique(data[[patient_id]]))!=dim(data)[1]) {
    #   stop("There should be one row for each individual in data")
    # }
    # n <- dim(data)[1]
    # Hist<-prodlim::Hist
    # Surv<-survival::Surv

    # if (return_c_index == TRUE) {
    #
    #   c_index <-
    #     pec::cindex(
    #       object = matrix(data$event_prediction),
    #       formula = Hist(event_time, event_status) ~ 1,
    #       cause = 1,
    #       data = data,
    #       cens.model = "marginal",
    #       splitMethod = "none",
    #       B = 0,
    #       exact = FALSE,
    #       verbose = FALSE,
    #       eval.times = x_hor
    #     )$AppCindex$matrix
    #
    #   if (standard_error == TRUE) {
    #     rho_hat_star <-
    #       unlist(lapply(1:b, function(i) {
    #         c_index_bootstrap(
    #           n = n,
    #           data = data,
    #           event_prediction =
    #             event_prediction,
    #           event_time =
    #             event_time,
    #           event_status =
    #             event_status,
    #           x_hor =
    #             x_hor
    #         )
    #       }))
    #     c_index_error <-
    #       sqrt((1 / (b - 1)) * (sum((
    #         rho_hat_star - mean(rho_hat_star)
    #       ) ^ 2)))
    #   }
    # }
    #
    # if (return_brier_score == TRUE) {
    #   brier_score <-
    #     pec::pec(
    #       object = matrix(cbind(0, data[[event_prediction]]), ncol =
    #                         2),
    #       formula = Hist(event_time, event_status) ~ 1,
    #       cause = 1,
    #       data = data,
    #       cens.model = "marginal",
    #       splitMethod = "none",
    #       B = 0,
    #       exact = FALSE,
    #       verbose = FALSE,
    #       reference = FALSE,
    #       times = x_hor
    #     )$AppErr$matrix[2]
    #
    #   if (standard_error == TRUE) {
    #     rho_hat_star <-
    #       unlist(lapply(1:b, function(i) {
    #         brier_bootstrap(
    #           n = n,
    #           data = data,
    #           event_prediction =
    #             event_prediction,
    #           event_time =
    #             event_time,
    #           event_status =
    #             event_status,
    #           x_hor =
    #             x_hor
    #         )
    #       }))
    #     brier_score_error <-
    #       sqrt((1 / (b - 1)) * (sum((
    #         rho_hat_star - mean(rho_hat_star)
    #       ) ^ 2)))
    #
    #
    #   }
    # }
    #
    # if (standard_error == FALSE) {
    #   return(list(c_index = c_index, brier_score = brier_score))
    # }
    # if (standard_error == TRUE) {
    #   return(list(
    #     c_index = c_index,
    #     c_index_standard_error = c_index_error,
    #     brier_score = brier_score,
    #     brier_score_standard_error = brier_score_error
    #   ))
    # }
    1
  }
#
# c_index_bootstrap <-
#   function(n,
#            data,
#            event_prediction,
#            event_time,
#            event_status,
#            x_hor) {
#     bootstrap_sample <- sample(n, replace = TRUE)
#     Hist<-prodlim::Hist
#     Surv<-survival::Surv
#     pec::cindex(
#       object = matrix(data[[event_prediction]][bootstrap_sample]),
#       formula = Hist(event_time, event_status) ~ 1,
#       cause = 1,
#       data = data[bootstrap_sample,],
#       cens.model = "marginal",
#       splitMethod = "none",
#       B = 0,
#       exact = FALSE,
#       verbose = FALSE,
#       eval.times = x_hor
#     )$AppCindex$matrix
#   }
#
# brier_bootstrap <-
#   function(n,
#            data,
#            event_prediction,
#            event_time,
#            event_status,
#            x_hor) {
#     bootstrap_sample <- sample(n, replace = TRUE)
#     Hist<-prodlim::Hist
#     Surv<-survival::Surv
#
#     pec::pec(
#       object = matrix(cbind(0, data[[event_prediction]][bootstrap_sample]), ncol =
#                         2),
#       formula = Hist(event_time, event_status) ~ 1,
#       cause = 1,
#       data = data[bootstrap_sample,],
#       cens.model = "marginal",
#       splitMethod = "none",
#       B = 0,
#       exact = FALSE,
#       verbose = FALSE,
#       reference = FALSE,
#       times = x_hor
#     )$AppErr$matrix[2]
#   }
