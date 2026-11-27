# Fitting these models takes several minutes, so skip this file on CRAN
skip_on_cran()

set.seed(1)
data(data_repeat_outcomes)
data_repeat_outcomes <-
  find_LOCF_risk_set(
    data_long = data_repeat_outcomes,
    individual_id = "id",
    event_time = "event_time",
    event_status = "event_status",
    covariates = c("ethnicity", "smoking", "diabetes", "sbp_stnd", "tchdl_stnd"),
    covariates_time = c(rep("response_time_sbp_stnd", 4), "response_time_tchdl_stnd"),
    x_L = c(60,61),
    x_hor= c(65,66)
  )
data_model_landmark_LOCF <-
  fit_LOCF_landmark(
    data_long = data_repeat_outcomes,
    x_L = c(60, 61),
    x_hor = c(65, 66),
    covariates = c("ethnicity", "smoking", "diabetes", "sbp_stnd", "tchdl_stnd"),
    covariates_time = c(rep("response_time_sbp_stnd", 4), "response_time_tchdl_stnd"),
    k = 10,
    individual_id = "id",
    event_time = "event_time",
    event_status = "event_status",
    survival_submodel = "cause_specific"
  )
cross_validation_list <- lapply(data_model_landmark_LOCF, "[[", i = 1)
data_model_landmark_LME <-
  fit_LME_landmark(
    data_long = data_repeat_outcomes[["60"]],
    x_L = c(60),
    x_hor = c(65),
    cross_validation_df =
      cross_validation_list,
    predictors_LME = c("ethnicity", "smoking", "diabetes"),
    predictors_LME_time =
      "response_time_sbp_stnd",
    responses_LME = c("sbp_stnd", "tchdl_stnd"),
    responses_LME_time = c("response_time_sbp_stnd", "response_time_tchdl_stnd"),
    individual_id = "id",
    standardise_time = TRUE,
    lme_control = nlme::lmeControl(maxIter = 100, msMaxIter = 100),
    event_time = "event_time",
    event_status = "event_status",
    survival_submodel = "cause_specific"
  )
# Compare predictions and performance metrics against reference values from
# Landmarking 1.1.0 (nested cross-validation for the LME model; AUC, Brier
# score and IPA from riskRegression::Score()). The LOCF predictions are the
# same as in the original release. Whole-object RDS snapshots were not
# reproducible: the riskRegression::CSC fit serialises differently on every
# run, even with the same data and seed. A small tolerance allows for
# numerical differences between nlme/riskRegression versions.
reference <- readRDS(test_path("fixtures", "main_test_reference.rds"))

expect_matches_reference <- function(model, ref) {
  expect_setequal(names(model), names(ref))
  for (x_l in names(ref)) {
    fitted <- model[[x_l]]$data
    idx <- match(ref[[x_l]]$id, as.character(fitted$id))
    expect_false(anyNA(idx))
    expect_equal(fitted$event_prediction[idx], ref[[x_l]]$event_prediction, tolerance = 1e-3)
    expect_equal(unlist(model[[x_l]]$prediction_error), ref[[x_l]]$prediction_error, tolerance = 1e-3)
  }
}

test_that("fitting landmark (LOCF)", {
  skip_on_cran()
  expect_matches_reference(data_model_landmark_LOCF, reference$LOCF)
})

test_that("fitting landmark (LME)", {
  skip_on_cran()
  expect_matches_reference(data_model_landmark_LME, reference$LME)
})
