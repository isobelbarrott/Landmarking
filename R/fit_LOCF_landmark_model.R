
return_LOCF_by_variable <- function(data,
                                    i,
                                    covariates,
                                    covariates_time,
                                    patient_id, x_L) {
  var <- covariates[i]
  time <- covariates_time[i]
  data_var <-
    data[data[[time]] <= x_L,][, c(patient_id, var, time)]
  data_var <-
    data_var[!is.na(data_var[[var]]),]#removes NA values
  data_var <-
    data_var[order(data_var[[time]], decreasing = TRUE), ]
  data_var <-
    data_var[!duplicated(data_var[[patient_id]]), ]
  return(data_var[, c(patient_id, var)])
}

fit_LOCF_longitudinal_model <- function(data,
                                        x_L,
                                        covariates,
                                        covariates_time,
                                        cv_name=NA,
                                        patient_id) {
  call <- match.call()
  if (!(is.data.frame(data))) {
    stop("data should be a dataframe")
  }
  if (!(is.numeric(x_L))){
    stop("'x_L' should be numeric")}

  for (col in c(
    covariates,
    covariates_time,
    patient_id
  )) {
    if (!(col %in% names(data))) {
      stop(col, " is not a column name in data")
    }
  }

  if (!is.na(cv_name)) {
    if (!(cv_name %in% names(data))) {
      stop(cv_name, " is not a column name in data")
    }
    if (any(is.na(data[[cv_name]]))){
      stop("The column ",cv_name, " contains NA values")
    }
  }

  if (!(length(covariates_time) %in% c(length(covariates),1))){
    stop("Length of covariates_time should be equal to length of covariates or 1")}

  if (length(covariates_time)==1){covariates_time<-rep(covariates_time,times=length(covariates))}

  if(dim(return_ids_with_LOCF(data=data,
                              patient_id=patient_id,
                              x_L=x_L,
                              covariates=covariates,
                              covariates_time=covariates_time))[1]!=dim(data)[1]){
    stop(paste0("data contains individuals that do not have a LOCF for all covariates at landmark age ",x_L,".
                  Use function return_ids_with_LOCF to remove these individuals from the dataset data."))
  }

  data[[patient_id]] <- as.factor(data[[patient_id]])
  data_LOCF <- data

  #Pick out LOCF for each exposure#
  LOCF_values_by_variable <-
    lapply(1:length(covariates), function(x) {
      return_LOCF_by_variable(
        data = data_LOCF,
        i = x,
        covariates = covariates,
        covariates_time = covariates_time,
        patient_id = patient_id,
        x_L = x_L
      )
    })
  data_LOCF <- Reduce(merge, LOCF_values_by_variable)
  data_LOCF<-data_LOCF[match(unique(data[[patient_id]]),data_LOCF[[patient_id]]),]
  if (!is.na(cv_name)){
    data_LOCF <-
      dplyr::left_join(data_LOCF, unique(data[c(patient_id, cv_name)]), by =
                         patient_id)
  }
  data_LOCF<-data_LOCF[,order(match(names(data_LOCF),names(data)))]
  data_LOCF<-data_LOCF[order(match(data_LOCF[[patient_id]],data[[patient_id]])),]
  rownames(data_LOCF)<-NULL

  list(data_longitudinal = data_LOCF,
       model_longitudinal = "LOCF",
       call = call)
}

fit_survival_model <- function(data,
                               patient_id,
                               cv_name=NA,
                               covariates,
                               event_time,
                               event_status,
                               survival_submodel = c("standard_cox", "cause_specific", "fine_gray"),
                               x_hor) {
  #Checks
  #####
  if (!(is.data.frame(data))) {
    stop("data should be a dataframe")
  }
  if (!(is.numeric(x_hor))) {
    stop("x_hor should be numeric")
  }
  for (col in c(covariates,
                event_time,
                event_status,
                patient_id)) {
    if (!(col %in% names(data))) {
      stop(col, " is not a column name in data")
    }
  }



  if (!is.na(cv_name)) {
    if (!(cv_name %in% names(data))) {
      stop(cv_name, " is not a column name in data")
    }
    if (any(is.na(data[[cv_name]]))){
      stop("The column ",cv_name, " contains NA values")
    }
  }

  if(is.na(cv_name)){data[["cross_validation_number"]]<-1;cv_name<-"cross_validation_number"}

  survival_submodel <- match.arg(survival_submodel)

  if(survival_submodel %in% c("cause_specific", "fine_gray")){
    if(!(setequal(data[[event_status]],0:max(data[[event_status]])))){
      stop("event_status column should contain only values 0, 1, 2, ... for cause_specific or fine_gray survival submodel,
        or values 0 and 1 for standard_cox survival submodel")
    }
  }
  if(survival_submodel %in% c("standard_cox")){
    if(!(setequal(data[[event_status]],0:1))){
      stop("event_status column should contain only values 0, 1, 2, ... for cause_specific or fine_gray survival submodel,
        or values 0 and 1 for standard_cox survival submodel")
    }
  }

  data[[patient_id]] <- as.factor(data[[patient_id]])

  if(is.na(cv_name)){data[["cross_validation_number"]]<-1;cv_name<-"cross_validation_number"}

    cv_numbers <- unique(data[[cv_name]])
    model <- as.list(cv_numbers)
    names(model) <- cv_numbers

    data_cv <- lapply(cv_numbers, function(cv_number) {

      if (length(cv_numbers)==1) {
        data_test <- data
        data_train <- data
      }else{
        data_test <- data[data[[cv_name]] == cv_number, ]
        data_train <- data[data[[cv_name]] != cv_number, ]
      }

      model_survival <- c()

      if (survival_submodel == "standard_cox") {
        formula_survival <-
          as.formula(paste0("Surv(", event_time, ", ", event_status, "==1) ~",
                            paste0(covariates, collapse = "+")))
        model_survival <- coxph(formula_survival, data = data_train,x=TRUE)
        data_test$event_prediction <-riskRegression::predictRisk(model_survival, times = x_hor, newdata = data_test)
      }

      if (survival_submodel == "cause_specific") {
        if(length(unique(data_train[[event_status]]))!=length(unique(data[[event_status]]))){
          stop("Not enough competing risk events to train competing risks model")
        }

        formula_survival <-
          as.formula(paste0("Hist(", event_time, ", ", event_status, ") ~",
                            paste0(covariates, collapse = "+")))
        model_survival <-
          riskRegression::CSC(
            formula_survival,
            data = data_train,
            fitter = "coxph",
            cause = 1
          )
        data_test$event_prediction <- as.numeric(
          riskRegression::predictRisk(
            model_survival,
            cause = 1,
            newdata = data_test,
            times = x_hor
          )
        )
      }

      if (survival_submodel == "fine_gray") {
        if(length(unique(data_train[[event_status]]))!=length(unique(data[[event_status]]))){
          stop("Not enough competing risk events to train competing risks model")
        }
        formula_survival <-
          as.formula(paste0("Hist(", event_time, ", ", event_status, ") ~",
                            paste0(covariates, collapse = "+")))
        model_survival <-
          riskRegression::FGR(formula_survival, data = data_train, cause = 1)
        data_test$event_prediction <-
          as.numeric(
            riskRegression::predictRisk(
              model_survival,
              cause = 1,
              newdata = data_test,
              times = x_hor
            )
          )
      }

      return(list(data_test = data_test, model_survival = model_survival))
    })
    data_survival <-
      data.frame(Reduce(dplyr::bind_rows, lapply(data_cv, `[[`, 1)))
    model_survival <- lapply(data_cv, `[[`, 2)
    names(model_survival)<-cv_numbers
    if (length(cv_numbers)==1) {
      model_survival<-model_survival[[1]]
    }

  data_survival<-data_survival[order(match(data_survival[[patient_id]],data[[patient_id]])),]
  rownames(data_survival)<-NULL
  list(data_survival = data_survival, model_survival = model_survival)
}

#' Fit a landmarking model using a last observation carried forward (LOCF) model for the longitudinal submodel
#'
#' This function performs the two-stage landmarking analysis. In the first stage, the longitudinal submodel is fitted using the LOCF model and in the
#' second stage the survival submodel is fitted.
#'
#' @param data_long Data frame or list of data frames each corresponding to a landmark age x_L (each element of the list must be named the value of x_L it corresponds to)
#' Each data frame contains repeat measurements data and time-to-event data in long format.
#' @template x_L
#' @template x_hor
#' @template event_status
#' @template event_time
#' @template start_study_time
#' @template end_study_time
#' @template cross_validation_df
#' @template k
#' @template b
#' @template covariates
#' @template covariates_time
#' @template patient_id
#' @template survival_submodel
#' @return List containing `data`, `model_longitudinal`, `model_survival`, and `prediction_error`.
#'
#' `data` is a data frame that includes the `covariates` values
#' calculated using the LOCF model (see Details section) and the predicted
#' probability that the event of interest has occured by time \code{x_L}, labelled as \code{"event_prediction"}.
#' There is one row for each individual.
#'
#' `model_longitudinal` indicates that the longitudinal submodel is LOCF.
#'
#' `model_survival` contains the outputs from the function used to fit the survival submodel, including the estimated parameters of the model.
#' For a model using cross-validation, `model_survival` contains a list of outputs with each
#' element in the list corresponding to a different cross-validation fold. For more information on how the survival model is fitted
#' please see `?fit_survival_model` which is a function used within `fit_LOCF_landmark_model`.
#'
#' `prediction_error` contains a list indicating the c-index and Brier score at time `x_hor` and their standard errors if parameter `b` is used.
#' For more information on how the prediction error is calculated
#' please see `?get_model_assessment` which is the function used to do this within `fit_LOCF_landmark_model`.
#'
#' @details
#' There are two parts to fitting the landmark model: the longitudinal submodel and the survival submodel.
#'
#' For the longitudinal model, this function uses the most recent values of the covariates at the landmark
#' age \code{x_L}. This is the LOCF model.
#'
#' For the survival submodel, there are three choices of model:
#' * the standard Cox model, this is a wrapper function for \code{coxph} from the package \code{survival}
#' * the cause-specific model, this is a wrapper function for \code{CSC} from package \code{riskRegression}
#' * the Fine Gray model, this is a wrapper function for \code{FGR} from package \code{riskRegression}
#'
#' The latter two models estimate the probability of the event of interest in the presence of competing events.
#'
#' @author Isobel Barrott \email{isobel.barrott@@gmail.com}
#' @examples
#' library(Landmarking)
#'  data(data_repeat)
#'  data(data_outcomes)
#'  data_repeat$response_time_tchdl_stnd <-
#'    as.numeric((
#'      as.Date(data_repeat$response_date_tchdl_stnd, format = "yyyy-mm-dd") -
#'        as.Date(data_repeat$dob, format = "yyyy-mm-dd")
#'    ) / 365.25)
#'  data_repeat$response_time_sbp_stnd <-
#'    as.numeric((
#'      as.Date(data_repeat$response_date_sbp_stnd, format = "yyyy-mm-dd") -
#'        as.Date(data_repeat$dob, format = "yyyy-mm-dd")
#'    ) / 365.25)
#' start_time <-
#'   stats::aggregate(stats::as.formula(
#'   paste0("response_time_sbp_stnd", "~", "id")
#'   ), data_repeat, function(x) {
#'     min(x)
#'   })
#' names(start_time)[2] <- "start_time"
#' data_repeat <- dplyr::left_join(data_repeat, start_time, by = "id")
#'  data_repeat_outcomes <-
#'    dplyr::left_join(data_repeat, data_outcomes, by = "id")
#'  data_repeat_outcomes <-
#'    return_ids_with_LOCF(
#'      data = data_repeat_outcomes,
#'      patient_id = "id",
#'      covariates =
#'        c("ethnicity", "smoking", "diabetes", "sbp_stnd", "tchdl_stnd"),
#'      covariates_time =
#'        c(rep("response_time_sbp_stnd", 4), "response_time_tchdl_stnd"),
#'      x_L = c(60, 61)
#'    )
#'  data_model_landmark_LOCF <-
#'    fit_LOCF_landmark_model(
#'      data_long = data_repeat_outcomes,
#'      x_L = c(60, 61),
#'      x_hor = c(65, 66),
#'      covariates =
#'        c("ethnicity", "smoking", "diabetes", "sbp_stnd", "tchdl_stnd"),
#'      covariates_time =
#'        c(rep("response_time_sbp_stnd", 4), "response_time_tchdl_stnd"),
#'      k = 10,
#'      start_study_time = "start_time",
#'      end_study_time = "event_time",
#'      patient_id = "id",
#'      event_time = "event_time",
#'      event_status = "event_status",
#'      survival_submodel = "cause_specific"
#'    )
#' @importFrom stats as.formula
#' @importFrom survival Surv
#' @importFrom survival coxph
#' @importFrom prodlim Hist
#' @export

fit_LOCF_landmark_model<-function(data_long,
                                  x_L,
                                  x_hor,
                                  covariates,
                                  covariates_time,
                                  start_study_time,
                                  end_study_time,
                                  k,
                                  cross_validation_df,
                                  patient_id,
                                  event_time,
                                  event_status,
                                  survival_submodel=c("standard_cox", "cause_specific", "fine_gray"),
                                  b){
  call <- match.call()

  survival_submodel <- match.arg(survival_submodel)

  #Checks

  if (missing(k)){k_add<-FALSE}else{k_add<-TRUE}
  if (missing(cross_validation_df)){cross_validation_df_add<-FALSE}else{cross_validation_df_add<-TRUE}
  if (k_add==TRUE) {
    if (!(is.numeric(k))) {
      stop("k should be numeric")
    }
  }
  if (cross_validation_df_add==TRUE) {
    if(class(cross_validation_df)=="list") {
      if(!all(x_L %in% names(cross_validation_df))){stop("The names of elements in cross_validation_df list should be the landmark times in x_L")}
      if(any(Reduce("c",lapply(cross_validation_df,function(x){any(duplicated(dplyr::distinct(x[,c(patient_id,"cross_validation_number")])[,patient_id]))})))){
        stop("Cross validation folds should be the same for the same individual")}}else if(class(cross_validation_df)=="data.frame"){
        if(any(duplicated(dplyr::distinct(cross_validation_df[,c(patient_id,"cross_validation_number")])[,patient_id]))){
          stop("Cross validation folds should be the same for the same individual")
        }
    }else{stop("cross_validation_df should be either a data frame or a list")}
  }
  if (k_add==TRUE && cross_validation_df_add==TRUE){stop("Either use parameter k or cross_validation_df but not both")}
  if (k_add==FALSE && cross_validation_df_add==FALSE){cv_name<-NA}else{cv_name<-"cross_validation_number"}

  if (!(length(covariates_time) %in% c(length(covariates),1))){
    stop("Length of covariates_time should be equal to length of covariates or 1")}
  if (length(covariates_time)==1){covariates_time<-rep(covariates_time,times=length(covariates))}

  if (length(x_L)!=length(x_hor)){stop("Length of x_L should be the same as length of x_hor")}
  if (!(is.data.frame(data_long)||is.list(data_long))){stop("data_long should be a list or data.frame")}
  if (is.data.frame(data_long)){
    data_long<-lapply(x_L,function(x_l){data_long})
    names(data_long)<-x_L
  }
  if (is.list(data_long)) {
    if(!setequal(names(data_long),x_L)){stop("Names of elements in data_long should be landmark ages x_L")}
  }



  if (missing(b)){b<-NA}

  data_long_x_L<-lapply(1:length(x_L),function(i){
    x_l<-x_L[i]
    x_h<-x_hor[i]

    data_long_x_L<-data_long[[as.character(x_l)]]

    if(survival_submodel %in% c("cause_specific", "fine_gray")){
      if(!(setequal(data_long_x_L[[event_status]],0:max(data_long_x_L[[event_status]])))){
      stop("event_status column should contain only values 0, 1, 2, ... for cause_specific or fine_gray survival submodel,
          or values 0 and 1 for standard_cox survival submodel")
      }
    }
    if(survival_submodel %in% c("standard_cox")){
      if(!(setequal(data_long_x_L[[event_status]],0:1))){
        stop("event_status column should contain only values 0, 1, 2, ... for cause_specific or fine_gray survival submodel,
          or values 0 and 1 for standard_cox survival submodel")
      }
    }
    if(!is.numeric(data_long_x_L[[event_status]])){stop("Column event_status should have class numeric")}

    if (!all(is.numeric(x_l))){
      stop("'x_L' should be numeric")}
    if (!all(is.numeric(x_h))){
      stop("'x_hor' should be numeric")}
    for (col in c(
      covariates,
      covariates_time,
      patient_id,
      start_study_time,
      end_study_time,
      event_time,
      event_status
    )) {
      if (!(col %in% names(data_long_x_L))) {
        stop(col, " is not a column name in data_long")
      }
    }
    data_long_x_L[[patient_id]]<-as.factor(data_long_x_L[[patient_id]])

      if(dim(return_ids_with_LOCF(data=data_long_x_L,
                                  patient_id=patient_id,
                                  x_L=x_l,
                                  covariates=covariates,
                                  covariates_time=covariates_time))[1]!=dim(data_long_x_L)[1]){
        stop("data_long contains individuals that do not have a LOCF for all covariates.
             Use function return_ids_with_LOCF to remove these individuals from the dataset")
      }


    data_long_x_L<-data_long_x_L[data_long_x_L[[start_study_time]]<=x_l & data_long_x_L[[end_study_time]]>x_l,]
    data_long_x_L[[event_status]][data_long_x_L[[event_time]]>x_h]<-0
    data_long_x_L[[event_time]][data_long_x_L[[event_time]]>x_h]<-x_h

    if(cross_validation_df_add==TRUE){data_long_x_L<-
      dplyr::left_join(data_long_x_L,cross_validation_df[[as.character(x_l)]][,c(patient_id,"cross_validation_number")],by=patient_id)}
    return(data_long_x_L)
  })

  names(data_long_x_L)<-x_L

  if(k_add==TRUE){
    data_long_x_L_cv<-add_cv_number(data=Reduce("rbind",data_long_x_L), patient_id=patient_id, k=k)
    data_long_x_L<-lapply(data_long_x_L,function(x){dplyr::left_join(x,dplyr::distinct(data_long_x_L_cv[,c(patient_id,"cross_validation_number")]),by=patient_id)})
  }

  out<-lapply(1:length(x_L),function(i){
    x_l<-x_L[i]
    x_h<-x_hor[i]

    data_long<-data_long_x_L[[as.character(x_l)]]

    print(paste0("Fitting longitudinal submodel, landmark age ",x_l))
    data_model_longitudinal<-fit_LOCF_longitudinal_model(data=data_long,
                                                         x_L=x_l,
                                                         covariates=covariates,
                                                         covariates_time=covariates_time,
                                                         cv_name=cv_name,
                                                         patient_id=patient_id)

    print(paste0("Complete, landmark age ",x_l))

    data_events<-dplyr::distinct(data_long[,c(patient_id, event_status,event_time)])
    data_longitudinal<-dplyr::left_join(data_model_longitudinal$data_longitudinal,data_events,by=patient_id)

    print(paste0("Fitting survival submodel, landmark age ",x_l))
    data_model_survival<-fit_survival_model(data=data_longitudinal,
                                            patient_id=patient_id,
                                            cv_name=cv_name,
                                            covariates=covariates,
                                            event_time=event_time,
                                            event_status=event_status,
                                            survival_submodel = survival_submodel,
                                            x_hor=x_h)
    print(paste0("Complete, landmark age ",x_l))


    prediction_error<-get_model_assessment(data=data_model_survival$data_survival,
                                               patient_id=patient_id,
                                               event_prediction="event_prediction",
                                               event_status=event_status,
                                               event_time=event_time,
                                               x_hor=x_h,
                                               b=b)

    list(data=data_model_survival$data_survival,
         model_longitudinal=data_model_longitudinal$model_longitudinal,
         model_survival=data_model_survival$model_survival,
         prediction_error=prediction_error,
         call=call
    )
  })
  names(out)<-x_L
  class(out)<-"landmark"
  out
}

