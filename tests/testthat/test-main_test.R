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
    fixed_effects = c("ethnicity", "smoking", "diabetes"),
    fixed_effects_time =
      "response_time_sbp_stnd",
    random_effects = c("sbp_stnd", "tchdl_stnd"),
    random_effects_time = c("response_time_sbp_stnd", "response_time_tchdl_stnd"),
    individual_id = "id",
    standardise_time = TRUE,
    lme_control = nlme::lmeControl(maxIter = 100, msMaxIter = 100),
    event_time = "event_time",
    event_status = "event_status",
    survival_submodel = "cause_specific"
  )
save_file <- function(code) {
  path <- tempfile(fileext = ".RDS")
  saveRDS(code,file = path)
  path
}

test_that("fitting landmark (LOCF)", {
  expect_snapshot_file(save_file(
    data_model_landmark_LOCF), "out_data_model_landmark_LOCF.RDS",cran = FALSE)
})

test_that("fitting landmark (LME)", {
  expect_snapshot_file(save_file(
    data_model_landmark_LME), "out_data_model_landmark_LME.RDS",cran = FALSE)
})
