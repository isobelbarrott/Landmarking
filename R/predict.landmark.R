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
#' For the LME model this the variables `fixed_effects`, `random_effects`, `fixed_effects_time`, and
#' `random_effects_time`. For the LOCF model this is `covariates` and `covariates_time`.
#' @param cv_fold If cross validation is used to fit `fit_LME_landmark` or `fit_LOCF_landmark`, then the cross validation fold to use when making risk predictions needs to be specified.
#' @param \dots Arguments passed on to `riskRegression::predictRisk`
#' @return Data frame `newdata` updated to contained a new column `event_prediction`
#' @examples
#' library(Landmarking)
#' data(data_repeat_outcomes)
#' data_repeat_outcomes <-
#'   return_ids_with_LOCF(
#'     data_long = data_repeat_outcomes,
#'     individual_id = "id",
#'     covariates =
#'       c("ethnicity", "smoking", "diabetes", "sbp_stnd", "tchdl_stnd"),
#'     covariates_time =
#'       c(rep("response_time_sbp_stnd", 4), "response_time_tchdl_stnd"),
#'     x_L = c(60,61)
#'   )
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
predict.landmark<-function(object,x_L,x_hor,newdata,cv_fold=NA,...){
  if(class(object)!="landmark"){stop("object must have class 'landmark'")}
  if(!(as.character(x_L) %in% names(object))){stop("x_L must be the name of an element in 'object'")}

  call<-lapply(as.list(object[[as.character(x_L)]]$call),eval)


  model_longitudinal<-object[[as.character(x_L)]]$model_longitudinal

  if (model_longitudinal=="LOCF"){
    covariates<-call$covariates
    covariates_time<-call$covariates_time
    individual_id<-call$individual_id
    for (covariate in covariates){
      if(is.factor(object[[as.character(x_L)]]$data[[covariate]])){
        newdata[[covariate]]<-factor(newdata[[covariate]],levels=levels(object[[as.character(x_L)]]$data[[covariate]]))
        }
      }


    data_model_longitudinal<-fit_LOCF_longitudinal(data_long=newdata,
                                                         x_L=x_L,
                                                         covariates=covariates,
                                                         covariates_time=covariates_time,
                                                         individual_id=individual_id)
    data_longitudinal<-data_model_longitudinal$data_longitudinal
  }

  if (model_longitudinal=="LME"){
    random_effects<-call$random_effects
    fixed_effects<-call$fixed_effects
    fixed_effects_time<-call$fixed_effects_time
    random_effects_time<-call$random_effects_time
    individual_id<-call$individual_id
    random_slope_as_covariate<-ifelse(is.null(call$random_slope_as_covariate),TRUE,call$random_slope_as_covariate)

    for (fixed_effect in fixed_effects){
      if(is.factor(object[[as.character(x_L)]]$data[[fixed_effect]])){
        newdata[[fixed_effect]]<-factor(newdata[[fixed_effect]],levels=levels(object[[as.character(x_L)]]$data[[fixed_effect]]))
      }
    }

    if (!is.na(cv_fold)){
      model_LME<-object[[as.character(x_L)]]$model_LME[[as.character(cv_fold)]]
      model_LME_standardise_time<-object[[as.character(x_L)]]$model_LME_standardise_time[[as.character(cv_fold)]]
    }else{
      model_LME<-object[[as.character(x_L)]]$model_LME
      model_LME_standardise_time<-object[[as.character(x_L)]]$model_LME_standardise_time
    }

    data_model_longitudinal<-fit_LOCF_longitudinal(data_long=newdata,
                                                         x_L=x_L,
                                                         covariates=fixed_effects,
                                                         covariates_time=fixed_effects_time,
                                                         individual_id=individual_id)
    data_LOCF<-data_model_longitudinal$data_longitudinal

    response_type <- Reduce(c,lapply(random_effects,function(i){rep(i,dim(data_LOCF)[1])}))
    response <- as.numeric(rep(NA,length(response_type)))
    response_time<-as.numeric(rep(x_L,length(response_type)))
    data_fixed_effects<-do.call("rbind",replicate(n=length(random_effects),data_LOCF[,c(individual_id,fixed_effects)],simplify=FALSE))
    data_LOCF<-data.frame(data_fixed_effects,response_type,response,response_time,predict=1)

    #Create validation and development dataset
    #####
    response_type <-
      Reduce(c, lapply(random_effects, function(i) {
        rep(i, dim(newdata)[1])
      }))
    response <-
      as.numeric(Reduce(c, lapply(1:length(random_effects), function(i) {
        newdata[, random_effects[i]]
      })))
    response_time <-
      as.numeric(Reduce(c, lapply(1:length(random_effects), function(i) {
        newdata[, random_effects_time[i]]
      })))
    data_fixed_effects<-
      do.call("rbind", replicate(
        n = length(random_effects),
        newdata[, c(individual_id, fixed_effects)],
        simplify = FALSE
      ))
    data_fixed_effects[[individual_id]]<-as.factor(data_fixed_effects[[individual_id]])
    data_LME <-
      data.frame(data_fixed_effects, response_type, response, response_time)

    data_LME$predict<-0
    data_longitudinal<-dplyr::bind_rows(data_LOCF,data_LME)

    mean_response_time<-object[[as.character(x_L)]]$model_LME_standardise_time$mean_response_time
    sd_response_time<-object[[as.character(x_L)]]$model_LME_standardise_time$sd_response_time
    data_longitudinal$response_time<-(data_longitudinal$response_time-mean_response_time)/sd_response_time

    response_predictions <-
      which(data_longitudinal$predict == 1)

    mixoutsamp_longitudinal <-
        mixoutsamp(model = model_LME,
                   newdata = data_longitudinal)


    data_longitudinal <-mixoutsamp_longitudinal$preddata[response_predictions, ][, c(individual_id,
                                                                                             fixed_effects,
                                                                                             "response_type",
                                                                                             "fitted")]
    data_longitudinal <-
      stats::reshape(
        data_longitudinal,
        timevar = "response_type",
        idvar = c(individual_id, fixed_effects),
        direction = "wide"
      )
    for (name in random_effects) {
      names(data_longitudinal)[grep(paste0("fitted.",name), names(data_longitudinal))] <-
        name
    }
    if(random_slope_as_covariate==TRUE){
      data_longitudinal<-dplyr::left_join(data_longitudinal,
                                                     mixoutsamp_longitudinal$random[,c(individual_id,paste0("reffresponse_type",random_effects,":response_time"))],by=individual_id)
      for (name in paste0(random_effects)){
        names(data_longitudinal)[grep(paste0("reffresponse_type",name,":response_time"),names(data_longitudinal))]<-paste0(name,"_slope")
      }
    }
  }

  if (!is.na(cv_fold)){
      model_survival<-object[[as.character(x_L)]]$model_survival[[as.character(cv_fold)]]
    }else{
      model_survival<-object[[as.character(x_L)]]$model_survival
    }
  if(!(class(model_survival) %in% c("CauseSpecificCox","FGR","coxph"))){stop("Class of survival model should be 'CauseSpecificCox','FGR', or 'coxph'")}
  if (class(model_survival) %in% c("CauseSpecificCox","FGR")){
    data_longitudinal$event_prediction <- as.numeric(
      riskRegression::predictRisk(
        model_survival,
        cause = 1,
        newdata = data_longitudinal,
        times = x_hor,...
      )
    )
  }
  if (class(model_survival)=="coxph") {
    data_longitudinal$event_prediction <- as.numeric(
      riskRegression::predictRisk(model_survival, times = x_hor, newdata = data_longitudinal,...)
    )
  }
  return(data_longitudinal)
}

