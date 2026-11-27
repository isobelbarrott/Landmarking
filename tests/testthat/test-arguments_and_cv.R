# Targeted tests for argument handling, cross-validation folds, alternative
# survival submodels and predict()/plot() on fitted landmark models.

data(data_repeat_outcomes)

covs <- c("ethnicity", "smoking", "diabetes", "sbp_stnd", "tchdl_stnd")
cov_times <- c(rep("response_time_sbp_stnd", 4), "response_time_tchdl_stnd")

quietly <- function(expr) suppressWarnings(suppressMessages(expr))

fit_locf <- function(...) {
  quietly(fit_LOCF_landmark(
    individual_id = "id",
    event_time = "event_time",
    event_status = "event_status",
    ...
  ))
}

fit_lme <- function(data_long, individual_id = "id",
                    responses_LME = c("sbp_stnd", "tchdl_stnd"),
                    responses_LME_time = c("response_time_sbp_stnd", "response_time_tchdl_stnd"),
                    ...) {
  quietly(fit_LME_landmark(
    data_long = data_long,
    x_L = 60,
    x_hor = 65,
    predictors_LME = c("ethnicity", "smoking", "diabetes"),
    predictors_LME_time = "response_time_sbp_stnd",
    responses_LME = responses_LME,
    responses_LME_time = responses_LME_time,
    individual_id = individual_id,
    standardise_time = TRUE,
    lme_control = nlme::lmeControl(maxIter = 100, msMaxIter = 100),
    event_time = "event_time",
    event_status = "event_status",
    survival_submodel = "cause_specific",
    ...
  ))
}

# Long-format data restricted to the risk set at landmark age 60
risk_set_60 <- quietly(find_LOCF_risk_set(
  data_long = data_repeat_outcomes,
  x_L = 60,
  x_hor = 65,
  covariates = covs,
  covariates_time = cov_times,
  individual_id = "id",
  event_time = "event_time",
  event_status = "event_status"
))[["60"]]

# Compare predictions for the same individuals, regardless of row order
expect_same_predictions <- function(pred, fitted, individual_id = "id") {
  idx <- match(as.character(fitted[[individual_id]]), as.character(pred[[individual_id]]))
  expect_false(anyNA(idx))
  expect_equal(pred$event_prediction[idx], fitted$event_prediction, tolerance = 1e-8)
}

test_that("integer x_L, x_hor, k and b are accepted and match doubles", {
  m_int <- fit_locf(data_long = data_repeat_outcomes, x_L = 60L, x_hor = 65L,
                    covariates = covs, covariates_time = cov_times, k = 5L, b = 3L,
                    survival_submodel = "cause_specific")
  m_dbl <- fit_locf(data_long = data_repeat_outcomes, x_L = 60, x_hor = 65,
                    covariates = covs, covariates_time = cov_times, k = 5,
                    survival_submodel = "cause_specific")
  expect_equal(m_int[["60"]]$data$event_prediction, m_dbl[["60"]]$data$event_prediction)
  expect_true(is.numeric(m_int[["60"]]$prediction_error$c_index_standard_error))
})

test_that("an individual_id other than 'id' works with random slopes (two responses)", {
  skip_on_cran() # slow: fits LME models
  data_renamed <- data_repeat_outcomes
  names(data_renamed)[names(data_renamed) == "id"] <- "patient_id"

  m_id <- fit_lme(data_repeat_outcomes, individual_id = "id", k = 3)
  m_pid <- fit_lme(data_renamed, individual_id = "patient_id", k = 3)

  expect_true(all(c("sbp_stnd_slope", "tchdl_stnd_slope") %in% names(m_pid[["60"]]$data)))
  expect_equal(m_pid[["60"]]$data$event_prediction, m_id[["60"]]$data$event_prediction)
})

test_that("an individual_id other than 'id' works with random slopes (one response)", {
  data_renamed <- data_repeat_outcomes
  names(data_renamed)[names(data_renamed) == "id"] <- "patient_id"

  m <- fit_lme(data_renamed, individual_id = "patient_id",
               responses_LME = "sbp_stnd",
               responses_LME_time = "response_time_sbp_stnd")
  expect_true("sbp_stnd_slope" %in% names(m[["60"]]$data))
  expect_true(all(m[["60"]]$data$event_prediction >= 0 & m[["60"]]$data$event_prediction <= 1))
})

test_that("predict() works when arguments were passed as variables (LOCF)", {
  fit_in_function <- function() {
    my_covs <- covs
    my_cov_times <- cov_times
    fit_locf(data_long = data_repeat_outcomes, x_L = 60, x_hor = 65,
             covariates = my_covs, covariates_time = my_cov_times,
             survival_submodel = "cause_specific")
  }
  m <- fit_in_function()
  expect_false(exists("my_covs"))

  pred <- quietly(predict(m, x_L = 60, x_hor = 65, newdata = risk_set_60))
  # Without cross-validation, in-sample predictions should reproduce the fit
  expect_same_predictions(pred, m[["60"]]$data)
})

test_that("predict() works when arguments were passed as variables (LME, renamed id)", {
  skip_on_cran() # slow: fits LME models
  data_renamed <- data_repeat_outcomes
  names(data_renamed)[names(data_renamed) == "id"] <- "patient_id"
  newdata <- risk_set_60
  names(newdata)[names(newdata) == "id"] <- "patient_id"

  fit_in_function <- function() {
    my_responses <- c("sbp_stnd", "tchdl_stnd")
    fit_lme(data_renamed, individual_id = "patient_id", responses_LME = my_responses)
  }
  m <- fit_in_function()

  pred <- quietly(predict(m, x_L = 60, x_hor = 65, newdata = newdata))
  expect_same_predictions(pred, m[["60"]]$data, individual_id = "patient_id")
})

test_that("a single cross_validation_df can be used for several landmark ages", {
  # Long format: one row per measurement, so several rows per individual
  cv_df <- add_cv_number(data_long = data_repeat_outcomes, individual_id = "id", k = 5)

  m <- fit_locf(data_long = data_repeat_outcomes, x_L = c(60, 61), x_hor = c(65, 66),
                covariates = covs, covariates_time = cov_times,
                cross_validation_df = cv_df, survival_submodel = "cause_specific")

  expect_named(m, c("60", "61"))
  # One row per individual (long-format folds must not duplicate rows)
  expect_equal(nrow(m[["60"]]$data), length(unique(risk_set_60$id)))
  # The folds used are those supplied
  folds <- dplyr::distinct(cv_df[, c("id", "cross_validation_number")])
  used <- m[["60"]]$data
  expect_equal(used$cross_validation_number,
               folds$cross_validation_number[match(as.character(used$id), folds$id)])
})

test_that("missing folds are detected at every landmark age, not only the last", {
  cv_df <- dplyr::distinct(
    add_cv_number(data_long = data_repeat_outcomes, individual_id = "id", k = 5)[, c("id", "cross_validation_number")]
  )
  dropped_id <- as.character(risk_set_60$id[1])
  cv_list <- list("60" = cv_df[cv_df$id != dropped_id, ], "61" = cv_df)

  expect_error(
    fit_locf(data_long = data_repeat_outcomes, x_L = c(60, 61), x_hor = c(65, 66),
             covariates = covs, covariates_time = cov_times,
             cross_validation_df = cv_list, survival_submodel = "cause_specific"),
    "landmark age 60"
  )
})

test_that("standard_cox submodel fits and predict() reproduces the fit", {
  data_cox <- data_repeat_outcomes
  data_cox$event_status <- factor(ifelse(data_cox$event_status == "1", 1, 0))

  m <- fit_locf(data_long = data_cox, x_L = 60, x_hor = 65,
                covariates = covs, covariates_time = cov_times,
                survival_submodel = "standard_cox")
  p <- m[["60"]]$data$event_prediction
  expect_true(all(p >= 0 & p <= 1))
  expect_true(m[["60"]]$prediction_error$c_index > 0.5)

  newdata <- data_cox[data_cox$id %in% risk_set_60$id, ]
  pred <- quietly(predict(m, x_L = 60, x_hor = 65, newdata = newdata))
  expect_same_predictions(pred, m[["60"]]$data)
})

test_that("fine_gray submodel fits", {
  m <- fit_locf(data_long = data_repeat_outcomes, x_L = 60, x_hor = 65,
                covariates = covs, covariates_time = cov_times,
                survival_submodel = "fine_gray")
  p <- m[["60"]]$data$event_prediction
  expect_true(all(p >= 0 & p <= 1))
  expect_s3_class(m[["60"]]$model_survival, "FGR")
})

test_that("plot() works when arguments were passed as variables", {
  skip_if_not_installed("mstate")
  fit_in_function <- function() {
    my_covs <- covs
    my_cov_times <- cov_times
    fit_locf(data_long = data_repeat_outcomes, x_L = 60, x_hor = 65,
             covariates = my_covs, covariates_time = my_cov_times,
             survival_submodel = "cause_specific")
  }
  m <- fit_in_function()
  expect_s3_class(plot(m, x_L = 60, n = 5), "ggplot")
})

# A new individual, as in the how_to_use vignette: character/numeric
# covariates rather than factors, and only one ethnicity present
new_individual <- data.frame(
  id = c("3001", "3001", "3001"),
  response_time_sbp_stnd = c(57, 58, 59),
  smoking = c(0, 0, 0),
  diabetes = c(0, 0, 0),
  ethnicity = c("Indian", "Indian", "Indian"),
  sbp_stnd = c(0.45, 0.87, 0.85),
  tchdl_stnd = c(-0.7, 0.24, 0.3),
  response_time_tchdl_stnd = c(57, 58, 59)
)

# The same individual alongside 50 individuals from the risk set
with_others <- function(new_individual) {
  others <- risk_set_60[risk_set_60$id %in% unique(risk_set_60$id)[1:50], names(new_individual)]
  others[] <- lapply(others, function(x) if (is.factor(x)) as.character(x) else x)
  others$smoking <- as.numeric(others$smoking)
  others$diabetes <- as.numeric(others$diabetes)
  rbind(new_individual, others)
}

expect_same_alone_and_with_others <- function(model) {
  pred_alone <- quietly(predict(model, x_L = 60, x_hor = 64, newdata = new_individual))
  pred_group <- quietly(predict(model, x_L = 60, x_hor = 64, newdata = with_others(new_individual)))
  expect_equal(nrow(pred_alone), 1)
  p <- pred_alone$event_prediction
  expect_true(p > 0 && p < 1)
  expect_equal(p, pred_group$event_prediction[as.character(pred_group$id) == "3001"])
}

test_that("predict() works for a single new individual (LOCF)", {
  m <- fit_locf(data_long = data_repeat_outcomes, x_L = 60, x_hor = 65,
                covariates = covs, covariates_time = cov_times,
                survival_submodel = "cause_specific")
  expect_same_alone_and_with_others(m)
})

test_that("predict() works for a single new individual (LME)", {
  skip_on_cran() # slow: fits LME models
  m <- fit_lme(data_repeat_outcomes)
  expect_same_alone_and_with_others(m)
})

test_that("add_cv_number() does not change the global random number state", {
  set.seed(123)
  expected <- stats::runif(1)
  set.seed(123)
  folds_1 <- add_cv_number(data_long = data_repeat_outcomes, individual_id = "id", k = 5)
  expect_equal(stats::runif(1), expected)
  # Fold assignment is still determined by the seed argument
  folds_2 <- add_cv_number(data_long = data_repeat_outcomes, individual_id = "id", k = 5)
  expect_equal(folds_1$cross_validation_number, folds_2$cross_validation_number)
})

test_that("x_hor must be greater than x_L", {
  expect_error(
    fit_locf(data_long = data_repeat_outcomes, x_L = c(60, 61), x_hor = c(65, 61),
             covariates = covs, covariates_time = cov_times,
             survival_submodel = "cause_specific"),
    "x_hor should be greater than x_L"
  )
  expect_error(
    fit_LME_landmark(
      data_long = data_repeat_outcomes, x_L = 60, x_hor = 60,
      predictors_LME = c("ethnicity", "smoking", "diabetes"),
      predictors_LME_time = "response_time_sbp_stnd",
      responses_LME = c("sbp_stnd", "tchdl_stnd"),
      responses_LME_time = c("response_time_sbp_stnd", "response_time_tchdl_stnd"),
      individual_id = "id", event_time = "event_time", event_status = "event_status",
      survival_submodel = "cause_specific"
    ),
    "x_hor should be greater than x_L"
  )
})

test_that("LOCF model without cross-validation: in-sample message, predict() and plot() checks", {
  suppressWarnings(suppressMessages(expect_message(
    m <- fit_LOCF_landmark(
      data_long = data_repeat_outcomes, x_L = 60, x_hor = 65,
      covariates = covs, covariates_time = cov_times,
      individual_id = "id", event_time = "event_time", event_status = "event_status",
      survival_submodel = "cause_specific"
    ),
    "in-sample"
  )))

  # ...but not when cross-validation is used
  messages <- character()
  withCallingHandlers(
    suppressWarnings(fit_LOCF_landmark(
      data_long = data_repeat_outcomes, x_L = 60, x_hor = 65,
      covariates = covs, covariates_time = cov_times, k = 5,
      individual_id = "id", event_time = "event_time", event_status = "event_status",
      survival_submodel = "cause_specific"
    )),
    message = function(cnd) {
      messages <<- c(messages, conditionMessage(cnd))
      invokeRestart("muffleMessage")
    }
  )
  expect_false(any(grepl("in-sample", messages)))

  expect_error(predict(m, x_L = 60, x_hor = 60, newdata = new_individual),
               "x_hor should be greater than x_L")

  # Too many tied predictions to form the requested number of groups
  m_tied <- m
  m_tied[["60"]]$data$event_prediction <-
    rep(c(0.1, 0.2), length.out = nrow(m_tied[["60"]]$data))
  expect_error(plot(m_tied, x_L = 60, n = 5), "tied values of event_prediction")

  # An event type with no events
  m_no_event_1 <- m
  status <- as.numeric(as.character(m_no_event_1[["60"]]$data$event_status))
  m_no_event_1[["60"]]$data$event_status <- ifelse(status == 1, 0, status)
  expect_error(plot(m_no_event_1, x_L = 60, n = 5), "Not enough events \\(for event 1\\)")
})

test_that("LME model works when a fold has no individuals in a category of a predictor", {
  skip_on_cran() # slow: fits LME models
  # All individuals of the rarest ethnicity are in folds 2 and 3, so fold 1 has
  # none, but every training set (the other folds) still contains some
  individuals <- data_repeat_outcomes[!duplicated(data_repeat_outcomes$id), c("id", "ethnicity")]
  set.seed(42)
  individuals$cross_validation_number <- sample(1:3, nrow(individuals), replace = TRUE)
  indian <- individuals$ethnicity == "Indian"
  individuals$cross_validation_number[indian] <- rep(2:3, length.out = sum(indian))
  cv_df <- individuals[, c("id", "cross_validation_number")]

  m <- fit_lme(data_repeat_outcomes, cross_validation_df = cv_df)
  fitted <- m[["60"]]$data
  expect_false(any(fitted$ethnicity[fitted$cross_validation_number == 1] == "Indian"))
  expect_true(any(fitted$cross_validation_number == 1))
  expect_true(all(fitted$event_prediction >= 0 & fitted$event_prediction <= 1))
})

test_that("LME model warns about individuals without all responses by the landmark age", {
  skip_on_cran() # slow: fits LME models
  # Individuals in the LME risk set at age 60 (predictors measured by age 60)
  # who have no total/HDL cholesterol measurement on or before age 60
  d <- data_repeat_outcomes
  risk_set_ids <- unique(d$id[d$event_time > 60 & d$response_time_sbp_stnd <= 60])
  n_expected <- sum(!(risk_set_ids %in% d$id[d$response_time_tchdl_stnd <= 60]))
  expect_gt(n_expected, 0)

  suppressMessages(suppressWarnings(expect_warning(
    m <- fit_LME_landmark(
      data_long = data_repeat_outcomes,
      x_L = 60,
      x_hor = 65,
      predictors_LME = c("ethnicity", "smoking", "diabetes"),
      predictors_LME_time = "response_time_sbp_stnd",
      responses_LME = c("sbp_stnd", "tchdl_stnd"),
      responses_LME_time = c("response_time_sbp_stnd", "response_time_tchdl_stnd"),
      individual_id = "id",
      standardise_time = TRUE,
      lme_control = nlme::lmeControl(maxIter = 100, msMaxIter = 100),
      event_time = "event_time",
      event_status = "event_status",
      survival_submodel = "cause_specific"
    ),
    paste0("^", n_expected, " individuals have been removed from the survival submodel")
  )))
  expect_equal(nrow(m[["60"]]$data), length(risk_set_ids) - n_expected)

  # Character predictors give the same model as factor predictors, and
  # predict() works for a single new individual
  data_chr <- data_repeat_outcomes
  data_chr$ethnicity <- as.character(data_chr$ethnicity)
  m_chr <- fit_lme(data_chr)
  expect_equal(m_chr[["60"]]$data$event_prediction, m[["60"]]$data$event_prediction,
               tolerance = 1e-6)
  expect_same_alone_and_with_others(m_chr)
})
