#' Smooth calibration curve and calibration metrics for a landmark model
#'
#' Estimates a smooth calibration curve for the landmark model at landmark time `x_L`, and summarises
#' it with the integrated calibration index (ICI), E50 and E90.
#'
#' @param x Object inheriting the class `landmark`, this should be the output from either `fit_LME_landmark` or `fit_LOCF_landmark`.
#' @param x_L Numeric specifying the landmark time. This indicates which landmark model in `x` to use.
#' @param span Numeric specifying the span of the `stats::loess` smoother (default 0.75).
#' @return List containing `curve`, a data frame of predicted risks (`predicted`) and the smoothed observed risks (`observed`)
#' at 100 points over the range of the predicted risks, and the calibration metrics `ICI`, `E50` and `E90`.
#' @details The observed risk of the event of interest (the event with `event_status=1`) by the horizon time is estimated for each individual
#' with a pseudo-observation: the jackknife pseudo-value of the Aalen-Johansen estimate of the cumulative incidence of the event of interest
#' (or of one minus the Kaplan-Meier estimate when there are no competing events) at the horizon time. The pseudo-observations are
#' regressed on the predicted risks with a loess smoother to give the calibration curve, which accounts for censoring and competing risks.
#' For a perfectly calibrated model the curve lies on the line y=x.
#'
#' For each individual, the absolute difference between their predicted risk and the smoothed observed risk at that prediction is
#' calculated. The ICI is the mean of these differences, and E50 and E90 are their median and 90th percentile
#' (Austin et al., 2020, 2022). Smaller values indicate better calibration.
#'
#' If the model was fitted with cross-validation, the predicted risks are the out-of-fold predictions, so the calibration is
#' assessed on individuals who were not used to fit the models that predicted them.
#' @references
#' Austin PC, Harrell FE, van Klaveren D (2020). Graphical calibration curves and the integrated calibration index (ICI)
#' for survival models. Statistics in Medicine, 39(21), 2714-2742.
#'
#' Austin PC, Putter H, Giardiello D, van Klaveren D (2022). Graphical calibration curves and the integrated calibration index (ICI)
#' for competing risk models. Diagnostic and Prognostic Research, 6, 2.
#' @examples
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
#' calibration <- get_calibration(x = data_model_landmark_LOCF, x_L = 60)
#' calibration[c("ICI", "E50", "E90")]
#' @export
get_calibration <- function(x, x_L, span = 0.75) {
  if (!inherits(x, "landmark")) {
    stop("x must have class 'landmark'")
  }
  if (!is.numeric(x_L)) {
    stop("x_L should be numeric")
  }
  if (!(as.character(x_L) %in% names(x))) {
    stop("x_L must be the name of an element in 'x'")
  }
  if (!is.numeric(span) || length(span) != 1 || span <= 0) {
    stop("span should be a positive number")
  }
  x <- x[[as.character(x_L)]]
  event_status <- get_landmark_arg(x, "event_status")
  event_time <- get_landmark_arg(x, "event_time")
  x_hor <- x$x_hor
  if (is.null(x_hor)) {
    #Models fitted with Landmarking < 1.1.0 do not store the horizon, but
    #their follow-up is censored at the horizon
    x_hor <- max(x$data[[event_time]])
  }
  get_smooth_calibration(
    event_prediction = x$data$event_prediction,
    event_time = x$data[[event_time]],
    event_status = x$data[[event_status]],
    x_hor = x_hor,
    span = span
  )
}

# Pseudo-observation based smooth calibration curve and ICI/E50/E90
get_smooth_calibration <- function(event_prediction,
                                   event_time,
                                   event_status,
                                   x_hor,
                                   span = 0.75) {
  if (!is.null(levels(event_status))) {
    event_status <- as.numeric(levels(event_status))[event_status]
  }
  data_pseudo <- data.frame(
    event_time = pmin(event_time, x_hor),
    event_status = ifelse(event_time > x_hor, 0, event_status)
  )
  if (!any(data_pseudo$event_status == 1)) {
    stop("There are no events of interest (event_status = 1) by the horizon time")
  }
  Hist <- prodlim::Hist
  fit <- prodlim::prodlim(Hist(event_time, event_status) ~ 1, data = data_pseudo)
  if (all(data_pseudo$event_status %in% c(0, 1))) {
    #Pseudo-values of the survival probability; the risk is one minus these
    pseudo <- 1 - as.numeric(prodlim::jackknife(fit, times = x_hor))
  } else {
    pseudo <- as.numeric(prodlim::jackknife(fit, times = x_hor, cause = 1))
  }

  data_smooth <- data.frame(pseudo = pseudo, event_prediction = event_prediction)
  smooth <- stats::loess(pseudo ~ event_prediction, data = data_smooth, span = span, degree = 2)
  observed_risk <- function(prediction) {
    observed <- stats::predict(smooth, newdata = data.frame(event_prediction = prediction))
    pmin(pmax(observed, 0), 1)
  }

  difference <- abs(observed_risk(event_prediction) - event_prediction)
  predicted <- seq(min(event_prediction), max(event_prediction), length.out = 100)
  list(
    curve = data.frame(predicted = predicted, observed = observed_risk(predicted)),
    ICI = mean(difference),
    E50 = stats::median(difference),
    E90 = unname(stats::quantile(difference, 0.9))
  )
}
