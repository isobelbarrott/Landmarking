context("Test data")
library(Landmarking)

test_that("fit_LOCF_landmark_model works", {
  library(Landmarking)
  set.seed(1)
  data(data_repeat_outcomes)
  data_repeat_outcomes <-
    return_ids_with_LOCF(
      data = data_repeat_outcomes,
      patient_id = "id",
      covariates = c("ethnicity", "smoking", "diabetes", "sbp_stnd", "tchdl_stnd"),
      covariates_time = c(rep("response_time_sbp_stnd", 4), "response_time_tchdl_stnd"),
      x_L = c(60, 61)
    )
  data_model_landmark_LOCF <-
    fit_LOCF_landmark_model(
      data_long = data_repeat_outcomes,
      x_L = c(60, 61),
      x_hor = c(65, 66),
      covariates = c("ethnicity", "smoking", "diabetes", "sbp_stnd", "tchdl_stnd"),
      covariates_time = c(rep("response_time_sbp_stnd", 4), "response_time_tchdl_stnd"),
      k = 10,
      start_study_time = "start_time",
      end_study_time = "event_time",
      patient_id = "id",
      event_time = "event_time",
      event_status = "event_status",
      survival_submodel = "cause_specific"
    )
  expect_true(all(0<=(data_model_landmark_LOCF$`60`$data$event_prediction)&
                    1>=(data_model_landmark_LOCF$`60`$data$event_prediction)), "Event predictions are not between 0 and 1")
  expect_true(all(0<=(data_model_landmark_LOCF$`61`$data$event_prediction)&
                    1>=(data_model_landmark_LOCF$`61`$data$event_prediction)), "Event predictions are not between 0 and 1")
  expect_equal(class(data_model_landmark_LOCF), "landmark")
  expect_true(all(0<=(data_model_landmark_LOCF$`60`$prediction_error$brier_score)&
                    1>=(data_model_landmark_LOCF$`60`$prediction_error$brier_score)), "Brier score is not between 0 and 1")
  expect_true(all(0<=(data_model_landmark_LOCF$`60`$prediction_error$c_index)&
                    1>=(data_model_landmark_LOCF$`60`$prediction_error$c_index)), "C-index is not between 0 and 1")

  newdata <-
    data.frame(
      id = c(3001, 3001, 3001),
      response_time_sbp_stnd = c(57, 58, 59),
      smoking = c(0, 0, 0),
      diabetes = c(0, 0, 0),
      ethnicity = c("Indian", "Indian", "Indian"),
      sbp_stnd = c(0.45, 0.87, 0.85),
      tchdl_stnd = c(-0.7, 0.24, 0.3),
      response_time_tchdl_stnd = c(57, 58, 59)
    )
  newdata_predict <-predict(
    object = data_model_landmark_LOCF,
    x_L = 60,
    x_hor = 62,
    newdata = newdata,
    cv_fold = 1
  )
  expect_true(all(0<=(newdata_predict$event_prediction)&
                    1>=(newdata_predict$event_prediction)), "Event predictions are not between 0 and 1")

  data_model_landmark_LOCF_plot<-plot(data_model_landmark_LOCF,x_L=60,n=5)
  expect_true(is.numeric(test$data$actual)&is.numeric(test$data$predicted),
              "Values are numeric")
})

test_that("fit_LME_landmark_model works", {
  library(Landmarking)
  set.seed(1)
  data(data_repeat_outcomes)
  data_repeat_outcomes <-
    return_ids_with_LOCF(
      data = data_repeat_outcomes,
      patient_id = "id",
      covariates = c("ethnicity", "smoking","sbp_stnd"),
      covariates_time = c(rep("response_time_sbp_stnd", 3)),
      x_L = c(60)
    )
  data_repeat_outcomes<-data_repeat_outcomes[data_repeat_outcomes$id %in% 1:100,]
  data_model_landmark_LME <-
    fit_LME_landmark_model(
      data_long = data_repeat_outcomes,
      x_L = c(60),
      x_hor = c(65),
      start_study_time = "start_time",
      end_study_time = "event_time",
      fixed_effects = c("ethnicity", "smoking"),
      fixed_effects_time =
        "response_time_sbp_stnd",
      random_effects = c("sbp_stnd"),
      random_effects_time = "response_time_sbp_stnd",
      patient_id = "id",
      standardise_time = TRUE,
      lme_control = nlme::lmeControl(maxIter =
                                       100, msMaxIter = 100),
      event_time = "event_time",
      event_status = "event_status",
      survival_submodel = "standard_cox"
    )
  expect_true(all(0<=(data_model_landmark_LME$`60`$data$event_prediction)&
                    1>=(data_model_landmark_LME$`60`$data$event_prediction)), "Event predictions are not between 0 and 1")
  expect_true(all(0<=(data_model_landmark_LME$`61`$data$event_prediction)&
                    1>=(data_model_landmark_LME$`61`$data$event_prediction)), "Event predictions are not between 0 and 1")
  expect_equal(class(data_model_landmark_LME), "landmark")
  expect_true(all(0<=(data_model_landmark_LME$`60`$prediction_error$brier_score)&
                    1>=(data_model_landmark_LME$`60`$prediction_error$brier_score)), "Brier score is not between 0 and 1")
  expect_true(all(0<=(data_model_landmark_LME$`60`$prediction_error$c_index)&
                    1>=(data_model_landmark_LME$`60`$prediction_error$c_index)), "C-index is not between 0 and 1")

  newdata <-
    data.frame(
      id = c(3001, 3001, 3001),
      response_time_sbp_stnd = c(57, 58, 59),
      smoking = c(0, 0, 0),
      diabetes = c(0, 0, 0),
      ethnicity = c("Indian", "Indian", "Indian"),
      sbp_stnd = c(0.45, 0.87, 0.85),
      tchdl_stnd = c(-0.7, 0.24, 0.3),
      response_time_tchdl_stnd = c(57, 58, 59)
    )
  newdata_predict <-predict(
    object = data_model_landmark_LME,
    x_L = 60,
    x_hor = 62,
    newdata = newdata
  )
  expect_true(all(0<=(newdata_predict$event_prediction)&
                    1>=(newdata_predict$event_prediction)), "Event predictions are not between 0 and 1")
})
