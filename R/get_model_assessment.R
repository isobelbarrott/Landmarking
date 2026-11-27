#' Compute the time-dependent AUC and Brier score
#'
#' Performs model assessment by computing the time-dependent area under the ROC curve (AUC)
#' and the Brier score at time `x_hor`, with their standard errors.
#'
#' @param data Data frame containing survival outcomes and the event predictions from the model,
#' there should be one row for each individual
#' @param event_status Character string specifying the column name in `data` which contains the event status (where 0=censoring, 1=event of interest, if there are competing events these are labelled 2 or above).
#' @param event_time Character string specifying the column name in `data` which contains the event time. Follow-up
#' may be censored at `x_hor`: individuals without an event on or before `x_hor` are treated as event-free at `x_hor`.
#' @template x_hor
#' @param individual_id  Character string specifying the column name in `data` which contains the individual identifiers
#' @template event_prediction
#' @template b
#' @return List containing the AUC (`auc`) and its standard error (`auc_standard_error`), the Brier score (`brier_score`)
#' and its standard error (`brier_score_standard_error`).
#' @details There are two factors in assessing the performance of a prediction model; its
#' discrimination and its calibration. Discrimination refers to the ability of the model to separate individuals into
#' those that will have the event of interest by time `x_hor` and those that will not. It is assessed with the
#' time-dependent AUC at `x_hor`: the probability that an individual who has the event of interest by `x_hor` (a case)
#' has a higher predicted risk than an individual who does not (a control). With competing risks, individuals who have a competing
#' event by `x_hor` are controls. The AUC is used rather than the C-index because the C-index is not a proper measure
#' for evaluating predicted risks at a fixed time (Blanche, Kattan and Gerds, 2019).
#'
#' The Brier score is the mean squared difference between the predicted risk and the event outcome by `x_hor` (where
#' an event of interest is 1 and not experiencing the event, including having a competing event, is 0). It reflects both calibration and
#' discrimination. Lower values are better.
#'
#' Inverse probability of censoring weighting (IPCW) is used to account for censoring, with a Kaplan-Meier model for censoring,
#' i.e. censoring is assumed to occur independently of covariates. The standard errors are based on the influence function
#' (Blanche et al., 2015) and treat the predicted risks as fixed, so they do not include the uncertainty from fitting the model.
#'
#' The metrics are calculated using the `Score` function in package `riskRegression`.
#'
#' In versions of Landmarking before 1.1.0, this function reported the C-index and Brier score, calculated with package `pec`
#' on follow-up censored at `x_hor`. Individuals who were event-free at `x_hor` were then treated as censored before `x_hor`,
#' which made the Brier score too low.
#'
#' @references
#' Blanche P, Kattan MW, Gerds TA (2019). The c-index is not proper for the evaluation of t-year predicted risks.
#' Biostatistics, 20(2), 347-357.
#'
#' Blanche P, Proust-Lima C, Loubere L, Berr C, Dartigues JF, Jacqmin-Gadda H (2015). Quantifying and comparing dynamic
#' predictive accuracy of joint models for longitudinal marker and time-to-event in presence of censoring and competing risks.
#' Biometrics, 71(1), 102-113.
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
#'   x_hor = 65)}
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

    if (!(is.numeric(x_hor))) {
      stop("x_hor should be numeric")
    }
    if (!missing(b)) {
      warn_b_deprecated(b)
    }
    if (length(unique(data[[individual_id]])) != dim(data)[1]) {
      stop("There should be one row for each individual in data")
    }

    status <- data[[event_status]]
    if (!is.null(levels(status))) {
      status <- as.numeric(levels(status))[status]
    }
    data_score <- data.frame(event_time = data[[event_time]],
                             event_status = status)

    #Individuals without an event on or before x_hor are event-free at x_hor.
    #When follow-up has been censored at x_hor their event time equals x_hor,
    #which Score() would treat as censored before x_hor, so move them just
    #past the horizon. This does not change the censoring weights before x_hor.
    event_free_at_x_hor <-
      data_score$event_time >= x_hor & data_score$event_status == 0
    data_score$event_time[event_free_at_x_hor] <- x_hor + 1

    Hist <- prodlim::Hist
    score <- riskRegression::Score(
      object = list(landmark_model = matrix(data[[event_prediction]], ncol = 1)),
      formula = Hist(event_time, event_status) ~ 1,
      data = data_score,
      times = x_hor,
      cause = 1,
      metrics = c("auc", "brier"),
      null.model = FALSE,
      se.fit = TRUE,
      cens.model = "km",
      progress.bar = NULL
    )
    auc <- as.data.frame(score$AUC$score)
    auc <- auc[auc$model == "landmark_model", ]
    brier <- as.data.frame(score$Brier$score)
    brier <- brier[brier$model == "landmark_model", ]

    return(
      list(
        auc = auc$AUC,
        auc_standard_error = auc$se,
        brier_score = brier$Brier,
        brier_score_standard_error = brier$se
      )
    )
  }

# Warn that the bootstrap argument b is no longer used (Landmarking >= 1.1.0)
warn_b_deprecated <- function(b) {
  if (!is.null(b) && !all(is.na(b))) {
    warning("b is no longer used: standard errors are now calculated with ",
            "riskRegression::Score() (see help(get_model_assessment))",
            call. = FALSE)
  }
}
