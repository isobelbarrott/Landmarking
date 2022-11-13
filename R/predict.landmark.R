#' Predict the risk of an event for a new individual using the landmark model
#'
#' This function predicts the risk of an event for new data using the landmark model fitted by `fit_LME_landmark` or `fit_LOCF_landmark`.
#' The 'event' is defined as event for which `event_status` is 1.
#'
#' @param object Object inheriting the class `landmark`, this should be the output from either `fit_LME_landmark` or `fit_LOCF_landmark`. It should contain a list
#' of landmark models corresponding to different landmark times `x_L`.
#' @param x_L 	Numeric specifying the landmark time. This indicates which landmark model in `object` to use.
#' @param x_hor Numeric specifying the horizon time. The function assesses the risk of event before this time.
#' @param newdata Data frame containing new data to return the risk prediction of the event of interest. The data should be in in long format
#' and the columns must contain the covariates and time variables that are used to fit the model.
#' For the LME model this the variables `predictors_LME`, `responses_LME`, `predictors_LME_time`, and
#' `responses_LME_time`. For the LOCF model this is `covariates` and `covariates_time`.
#' @param cv_fold If cross validation is used to fit `fit_LME_landmark` or `fit_LOCF_landmark`, then the cross validation fold to use when making risk predictions needs to be specified.
#' @param \dots Arguments passed on to `riskRegression::predictRisk`
#' @return Data frame `newdata` updated to contained a new column `event_prediction`
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
#'  newdata <-
#'    rbind(
#'      data.frame(
#'        id = c(3001, 3001, 3001),
#'        response_time_sbp_stnd = c(57, 58, 59),
#'        smoking = c(0, 0, 0),
#'        diabetes = c(0, 0, 0),
#'        ethnicity = c("Indian", "Indian", "Indian"),
#'        sbp_stnd = c(0.45, 0.87, 0.85),
#'        tchdl_stnd = c(-0.7, 0.24, 0.3),
#'        response_time_tchdl_stnd = c(57, 58, 59)
#'      )
#'    )
#'  predict(object=data_model_landmark_LOCF,x_L=60,x_hor=62,newdata=newdata,cv_fold=1)
#' @export
predict.landmark <- function(object, x_L, x_hor, newdata, cv_fold = NA, ...) {
  if (!inherits(object,"landmark")) {
    stop("object must have class 'landmark'")
  }
  if (!(as.character(x_L) %in% names(object))) {
    stop("x_L must be the name of an element in 'object'")
  }

  if (!(inherits(x_hor,"numeric"))) {
    stop("x_hor should have class numeric")
  }

  model_longitudinal <- object[[as.character(x_L)]]$model_longitudinal

  if (model_longitudinal == "LOCF") {
    call <- as.list(object[[as.character(x_L)]]$call)

    covariates <- eval(call$covariates)
    covariates_time <- eval(call$covariates_time)
    individual_id <- eval(call$individual_id)
    for (covariate in covariates) {
      if (is.factor(object[[as.character(x_L)]]$data[[covariate]])) {
        newdata[[covariate]] <-
          factor(newdata[[covariate]], levels = levels(object[[as.character(x_L)]]$data[[covariate]]))
      }
    }

    data_model_longitudinal <-
      fit_LOCF_longitudinal(
        data_long = newdata,
        x_L = x_L,
        covariates = covariates,
        covariates_time =
          covariates_time,
        individual_id =
          individual_id
      )
    data_longitudinal <- data_model_longitudinal$data_longitudinal
  }

  if (model_longitudinal == "LME") {
    call <- as.list(object[[as.character(x_L)]]$call)

    responses_LME <- eval(call$responses_LME)
    predictors_LME <- eval(call$predictors_LME)
    predictors_LME_time <- eval(call$predictors_LME_time)
    responses_LME_time <- eval(call$responses_LME_time)
    individual_id <- eval(call$individual_id)
    random_slope_survival <-
      ifelse(
        is.null(eval(call$random_slope_survival)),
        TRUE,
        eval(call$random_slope_survival)
      )

    for (predictor_LME in predictors_LME) {
      if (is.factor(object[[as.character(x_L)]]$data[[predictor_LME]])) {
        newdata[[predictor_LME]] <-
          factor(newdata[[predictor_LME]], levels = levels(object[[as.character(x_L)]]$data[[predictor_LME]]))
      }
    }

    if (!is.na(cv_fold)) {
      model_LME <-
        object[[as.character(x_L)]]$model_LME[[as.character(cv_fold)]]
      model_LME_standardise_time <-
        object[[as.character(x_L)]]$model_LME_standardise_time[[as.character(cv_fold)]]
    }else{
      model_LME <- object[[as.character(x_L)]]$model_LME
      model_LME_standardise_time <-
        object[[as.character(x_L)]]$model_LME_standardise_time
    }

    data_model_longitudinal <-
      fit_LOCF_longitudinal(
        data_long = newdata,
        x_L = x_L,
        covariates = predictors_LME,
        covariates_time =
          predictors_LME_time,
        individual_id =
          individual_id
      )
    data_LOCF <- data_model_longitudinal$data_longitudinal

    response_type <-
      Reduce(c, lapply(responses_LME, function(i) {
        rep(i, dim(data_LOCF)[1])
      }))
    response <- as.numeric(rep(NA, length(response_type)))
    response_time <- as.numeric(rep(x_L, length(response_type)))
    data_predictors_LME <-
      do.call("rbind", replicate(
        n = length(responses_LME),
        data_LOCF[, c(individual_id, predictors_LME)],
        simplify = FALSE
      ))
    data_LOCF <-
      data.frame(data_predictors_LME,
                 response_type,
                 response,
                 response_time,
                 predict = 1)

    #Create validation and development dataset
    #####
    response_type <-
      Reduce(c, lapply(responses_LME, function(i) {
        rep(i, dim(newdata)[1])
      }))
    response <-
      as.numeric(Reduce(c, lapply(1:length(responses_LME), function(i) {
        newdata[, responses_LME[i]]
      })))
    response_time <-
      as.numeric(Reduce(c, lapply(1:length(responses_LME), function(i) {
        newdata[, responses_LME_time[i]]
      })))
    data_predictors_LME <-
      do.call("rbind", replicate(
        n = length(responses_LME),
        newdata[, c(individual_id, predictors_LME)],
        simplify = FALSE
      ))
    data_predictors_LME[[individual_id]] <-
      as.factor(data_predictors_LME[[individual_id]])
    data_LME <-
      data.frame(data_predictors_LME,
                 response_type,
                 response,
                 response_time)
    data_LME<-data_LME[data_LME$response_time<=x_L,]
    data_LME$predict <- 0
    data_longitudinal <- dplyr::bind_rows(data_LOCF, data_LME)

    mean_response_time <-
      object[[as.character(x_L)]]$model_LME_standardise_time$mean_response_time
    sd_response_time <-
      object[[as.character(x_L)]]$model_LME_standardise_time$sd_response_time
    data_longitudinal$response_time <-
      (data_longitudinal$response_time - mean_response_time) / sd_response_time

    response_predictions <-
      which(data_longitudinal$predict == 1)

    mixoutsamp_longitudinal <-
      mixoutsamp(model = model_LME,
                 newdata = data_longitudinal)


    data_longitudinal <-
      mixoutsamp_longitudinal$preddata[response_predictions,][, c(individual_id,
                                                                  predictors_LME,
                                                                  "response_type",
                                                                  "fitted")]
    data_longitudinal <-
      stats::reshape(
        data_longitudinal,
        timevar = "response_type",
        idvar = c(individual_id, predictors_LME),
        direction = "wide"
      )
    for (name in responses_LME) {
      names(data_longitudinal)[grep(paste0("fitted.", name), names(data_longitudinal))] <-
        name
    }
    if (random_slope_survival == TRUE) {

      if (random_slope_survival == TRUE) {
        if (length(responses_LME)==1){
          slopes_df<-mixoutsamp_longitudinal$random[,c("id","reffresponse_time")]
          slopes_df["reffresponse_time"]<-slopes_df["reffresponse_time"]+model_LME$coefficients$fixed["response_time"]

        }
        if (length(responses_LME)>1){
          slopes_df<-mixoutsamp_longitudinal$random[,c("id",paste0("reffresponse_type",responses_LME,":response_time"))]
          responses_LME_dummy<-paste0("response_type",responses_LME,":response_time")
          responses_LME_dummy[1]<-"response_time"
          for (i in 1:length(responses_LME)){
            slopes_df[,paste0("reffresponse_type",responses_LME[i],":response_time")]<-
              slopes_df[,paste0("reffresponse_type",responses_LME[i],":response_time")]+model_LME$coefficients$fixed[responses_LME_dummy[1]]
            if(i!=1){slopes_df[,paste0("reffresponse_type",responses_LME[i],":response_time")]<-
              slopes_df[,paste0("reffresponse_type",responses_LME[i],":response_time")]+model_LME$coefficients$fixed[responses_LME_dummy[i]]}
          }
        }
        names(slopes_df)<-c("id",paste0(responses_LME,"_slope"))
        data_longitudinal <- dplyr::left_join(data_longitudinal,
                                                  slopes_df,
                                                  by = individual_id)
      }
    }
  }

  if (!is.na(cv_fold)) {
    model_survival <-
      object[[as.character(x_L)]]$model_survival[[as.character(cv_fold)]]
  } else{
    model_survival <- object[[as.character(x_L)]]$model_survival
  }
  if (!inherits(model_survival,c("CauseSpecificCox", "FGR", "coxph"))){
    stop("Class of survival model should be 'CauseSpecificCox','FGR', or 'coxph'")
  }
  if (inherits(model_survival,c("CauseSpecificCox", "FGR"))) {
    data_longitudinal$event_prediction <- as.numeric(
      riskRegression::predictRisk(
        model_survival,
        cause = 1,
        newdata = data_longitudinal,
        times = x_hor,
        ...
      )
    )
  }
  if (inherits(model_survival,"coxph")) {
    data_longitudinal$event_prediction <- as.numeric(
      riskRegression::predictRisk(model_survival, times = x_hor, newdata = data_longitudinal, ...)
    )
  }
  return(data_longitudinal)
}


