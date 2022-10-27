#' Compute C-index and Brier score
#'
#' Performs model assessment by computing the C-index and Brier score at time `x_hor`. There is the option
#' to calculate their standard errors using bootstraping.
#'
#' @param data Data frame containing survival outcomes and the event predictions from the model,
#' there should be one row for each individual
#' @param event_status Character string specifying the column name in `data` which contains the event status (where 0=censoring, 1=event of interest, if there are competing events these are labelled 2 or above). Events at time `x_hor` should be labelled censored.
#' @param event_time Character string specifying the column name in `data` which contains the event time.
#' @template x_hor
#' @param individual_id  Character string specifying the column name in `data` which contains the individual identifiers
#' @template event_prediction
#' @template b
#' @return List containing C-index, Brier score and their standard errors
#' @details There are two factors in assessing the performance of a prediction model; its
#' discrimination and its calibration. The c-index is a commonly used metric which assesses
#' discrimination, this refers to the ability of the model to separate individuals into
#' those that will have an event and those that will not. The c-index at a horizon time `x_hor`
#' looks at the pairs of individuals where one individual has the event at a time T and the other has not had the event at time T.
#' It is calculated as the proportion of these pairs where their relative risk prediction agrees with the
#' actual outcomes for the two individuals. This is extended to the competing risks case
#' by comparing individuals where one had the event of interest at time T and the other individual either
#' did not experience the event before this time T or experienced a competing event.
#'
#' The Brier score gives an indication of the calibration of a model (and its discrimination to an extent), this refers to the agreement between the risk prediction and
#' the outcome. The Brier score is calculated as the average mean squared error of the predicted risk and the event outcome (where
#' an event is 1 and not experiencing the event is 0). This is extended to the competing risks case by including the competing risk events as
#' not experiencing the event.
#'
#' For both the c-index and Brier score calculations, inverse probability censoring weighting (IPCW) is used to create weights
#' which account for the occurrence of censoring. The censoring model assumes for this function is the Kaplan Meier model, i.e. censoring occurs
#' independently of covariates.
#'
#' The c-index is calculated using the `cindex` function in package `pec`. The Brier score is calculated using
#' `pec` function in package `pec`.
#'
#'
#' @author Isobel Barrott \email{isobel.barrott@@gmail.com}
#' @examples \dontrun{
#' library(Landmarking)
#' data(data_repeat_outcomes)
#' data_model_landmark_LOCF <-
#'   fit_LOCF_landmark(
#'     data_long = data_repeat_outcomes,
#'     x_L = c(60, 61),
#'     x_hor = c(65, 66),
#'     covariates =
#'       c("ethnicity", "smoking", "diabetes", "sbp_stnd", "tchdl_stnd"),
#'     covariates_time =
#'       c(rep("response_time_sbp_stnd", 4), "response_time_tchdl_stnd"),
#'     k = 10,
#'     individual_id = "id",
#'     event_time = "event_time",
#'     event_status = "event_status",
#'     survival_submodel = "cause_specific"
#'   )
#' get_model_assessment(data = data_model_landmark_LOCF[["60"]]$data,
#'   individual_id = "id",
#'   event_prediction = "event_prediction",
#'   event_status = "event_status",
#'   event_time = "event_time",
#'   x_hor = 65,
#'   b = 100)}
#' @export


get_model_assessment <-
  function(data,
           individual_id,
           event_prediction,
           event_status,
           event_time,
           x_hor,
           b) {
    if (!(inherits(data,"data.frame"))) {
      stop("data should be a data frame")
    }
    if (!(inherits(individual_id,"character"))) {
      stop("individual_id should have class character")
    }
    if (!(inherits(event_prediction,"character"))) {
      stop("event_prediction should have class character")
    }
    if (!(inherits(event_status,"character"))) {
      stop("event_status should have class character")
    }
    if (!(inherits(event_time,"character"))) {
      stop("event_time should have class character")
    }
    for (col in c(individual_id,
                  event_prediction,
                  event_status,
                  event_time)) {
      if (!(col %in% names(data))) {
        stop(col, " is not a column name in data")
      }
      if(any(is.na(data[[col]]))){
        stop(col, " contains NA values")
      }
    }

    if (!(inherits(x_hor,"numeric"))) {
      stop("x_hor should be numeric")
    }
    if (!is.na(b)) {
      if (!(inherits(b,"numeric"))) {
        stop("b should be numeric")
      }
      standard_error <- TRUE
    }
    else{
      standard_error <- FALSE
    }
    if (length(unique(data[[individual_id]])) != dim(data)[1]) {
      stop("There should be one row for each individual in data")
    }

    n <- dim(data)[1]
    Hist <- prodlim::Hist
    Surv <- survival::Surv

    data[["event_time"]] <- data[[event_time]]
    data[["event_status"]] <- data[[event_status]]
    data[["event_prediction"]] <- data[[event_prediction]]

    if (setequal(data[[event_status]],0:1)){ # survival probabilities are needed for the non-competing risks model
      data[["event_prediction"]] <- 1 - data[["event_prediction"]]
    }

    c_index <-
      pec::cindex(
        object = matrix(data[["event_prediction"]]),
        formula = Hist(event_time, event_status) ~ 1,
        cause = 1,
        data = data,
        cens.model = "marginal",
        splitMethod = "none",
        B = 0,
        exact = FALSE,
        verbose = FALSE,
        eval.times = x_hor
      )$AppCindex$matrix

    if (standard_error == TRUE) {
      rho_hat_star <-
        unlist(lapply(1:b, function(i) {
          c_index_bootstrap(
            n = n,
            data = data,
            event_prediction =
              event_prediction,
            event_time =
              event_time,
            event_status =
              event_status,
            x_hor =
              x_hor
          )
        }))
      c_index_error <-
        sqrt((1 / (b - 1)) * (sum((
          rho_hat_star - mean(rho_hat_star)
        ) ^ 2)))
    }

    brier_score <-
      pec::pec(
        object = matrix(cbind(0, data[["event_prediction"]]), ncol =
                          2),
        formula = Hist(event_time, event_status) ~ 1,
        cause = 1,
        data = data,
        cens.model = "marginal",
        splitMethod = "none",
        B = 0,
        exact = FALSE,
        verbose = FALSE,
        reference = FALSE,
        times = x_hor
      )$AppErr$matrix[2]

    if (standard_error == TRUE) {
      rho_hat_star <-
        unlist(lapply(1:b, function(i) {
          brier_bootstrap(
            n = n,
            data = data,
            event_prediction =
              event_prediction,
            event_time =
              event_time,
            event_status =
              event_status,
            x_hor =
              x_hor
          )
        }))
      brier_score_error <-
        sqrt((1 / (b - 1)) * (sum((
          rho_hat_star - mean(rho_hat_star)
        ) ^ 2)))


    }


    if (standard_error == FALSE) {
      c_index_error <- NA
      brier_score_error <- NA
    }

    return(
      list(
        c_index = c_index,
        c_index_standard_error = c_index_error,
        brier_score = brier_score,
        brier_score_standard_error = brier_score_error
      )
    )
  }

c_index_bootstrap <-
  function(n,
           data,
           event_prediction,
           event_time,
           event_status,
           x_hor) {
    bootstrap_sample <- sample(n, replace = TRUE)
    Hist <- prodlim::Hist
    Surv <- survival::Surv

    data[["event_time"]] <- data[[event_time]]
    data[["event_status"]] <- data[[event_status]]
    data[["event_prediction"]] <- data[[event_prediction]]

    pec::cindex(
      object = matrix(data[["event_prediction"]][bootstrap_sample]),
      formula = Hist(event_time, event_status) ~ 1,
      cause = 1,
      data = data[bootstrap_sample, ],
      cens.model = "marginal",
      splitMethod = "none",
      B = 0,
      exact = FALSE,
      verbose = FALSE,
      eval.times = x_hor
    )$AppCindex$matrix
  }

brier_bootstrap <-
  function(n,
           data,
           event_prediction,
           event_time,
           event_status,
           x_hor) {
    bootstrap_sample <- sample(n, replace = TRUE)

    data[["event_time"]] <- data[[event_time]]
    data[["event_status"]] <- data[[event_status]]
    data[["event_prediction"]] <- data[[event_prediction]]

    pec::pec(
      object = matrix(cbind(0, data[["event_prediction"]][bootstrap_sample]), ncol =
                        2),
      formula = Hist(event_time, event_status) ~ 1,
      cause = 1,
      data = data[bootstrap_sample, ],
      cens.model = "marginal",
      splitMethod = "none",
      B = 0,
      exact = FALSE,
      verbose = FALSE,
      reference = FALSE,
      times = x_hor
    )$AppErr$matrix[2]
  }
