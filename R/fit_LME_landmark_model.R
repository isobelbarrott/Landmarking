
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



fit_LME_longitudinal_model <- function(data,
                                       x_L,
                                       fixed_effects,
                                       random_effects,
                                       fixed_effects_time,
                                       random_effects_time,
                                       standardise_time=FALSE,
                                       random_slope = TRUE,
                                       cv_name=NA,
                                       patient_id,
                                       lme_control = nlme::lmeControl()) {
  call <- match.call()
  if (!(is.data.frame(data))) {
    stop("data should be a dataframe")
  }
  if (!(is.numeric(x_L))){
    stop("'x_L' should be numeric")}

  for (col in c(
    fixed_effects,
    random_effects,
    fixed_effects_time,
    random_effects_time,
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

  if(is.na(cv_name)){data[["cross_validation_number"]]<-1;cv_name<-"cross_validation_number"}

  if (!(length(fixed_effects_time) %in% c(length(fixed_effects),1))){
    stop("Length of fixed_effects_time should be equal to length of fixed_effects or 1")}

  if (length(fixed_effects_time)==1){fixed_effects_time<-rep(fixed_effects_time,times=length(fixed_effects))}

  if (!(length(random_effects_time) %in% c(length(random_effects),1))){
    stop("Length of random_effects_time should be equal to length of random_effects or 1")}

  if (length(random_effects_time)==1){random_effects_time<-rep(random_effects_time,times=length(random_effects))}

  if(dim(return_ids_with_LOCF(data=data,
                              patient_id=patient_id,
                              x_L=x_L,
                              covariates=fixed_effects,
                              covariates_time=fixed_effects_time))[1]!=dim(data)[1]){
    stop("data contains individuals that do not have a LOCF for all fixed_effects. Use function return_ids_with_LOCF to remove these individuals from the dataset data.")
  }


  data[[patient_id]] <- as.factor(data[[patient_id]])

  data_LOCF <- data
  data_LME <- data

  #Pick out LOCF for each exposure#
  LOCF_values_by_variable <-
    lapply(1:length(c(fixed_effects, random_effects)), function(x) {
      return_LOCF_by_variable(
        data = data_LOCF,
        i = x,
        covariates=c(fixed_effects, random_effects),
        covariates_time = c(fixed_effects_time, random_effects_time),
        patient_id = patient_id,
        x_L =
          x_L
      )
    })
  data_LOCF <- Reduce(merge, LOCF_values_by_variable)
  data_LOCF<-data_LOCF[match(unique(data[[patient_id]]),data_LOCF[[patient_id]]),]

  data_LOCF <-
    dplyr::left_join(data_LOCF, unique(data[c(patient_id, cv_name)]), by =
                         patient_id)


  #Create validation and development dataset
  #####
  response_type <-
    Reduce(c, lapply(random_effects, function(i) {
      rep(i, dim(data_LME)[1])
    }))

  response <-
    as.numeric(Reduce(c, lapply(1:length(random_effects), function(i) {
      data_LME[, random_effects[i]]
    })))
  response_time <-
    as.numeric(Reduce(c, lapply(1:length(random_effects), function(i) {
      data_LME[, random_effects_time[i]]
    })))
    data_fixed_effects<-
      do.call("rbind", replicate(
        n = length(random_effects),
        data_LME[, c(patient_id, fixed_effects, cv_name)],
        simplify = FALSE
      ))

  data_LME <-
    data.frame(data_fixed_effects, response_type, response, response_time)

  if(standardise_time==TRUE){
    mean_response_time<-mean(data_LME$response_time,na.rm=TRUE)
    sd_response_time<-stats::sd(data_LME$response_time,na.rm=TRUE)
  }else{
    mean_response_time<-0
    sd_response_time<-1
  }
  data_LME$response_time<-(data_LME$response_time-mean_response_time)/sd_response_time
  x_L<-(x_L-mean_response_time)/sd_response_time
  standardise_time<-list(mean_response_time=mean_response_time,sd_response_time=sd_response_time)

  data_LME_model_val <- data_LME[data_LME$response_time <= x_L,]
  data_LME_model_dev <- data_LME[!is.na(data_LME$response),]
  #####

  if (length(random_effects)==1){
    formula_weights <- NULL
    if (random_slope == FALSE) {formula_random<-as.formula(paste0(" ~ 1 |",patient_id))}else{formula_random<-as.formula(paste0(" ~ 1 + response_time |",patient_id))}
    formula_fixed<-as.formula(paste0(c(paste0("response~ 1 "), c("response_time", fixed_effects)), collapse = "+"))
  }
  if (length(random_effects)>1){
    formula_weights <- nlme::varIdent(form = ~ 1 | "response_type")
    if (random_slope == FALSE) {formula_random<-as.formula(paste0(" ~ -1 + response_type | ",patient_id))}else{
     formula_random<-as.formula(paste0(" ~-1 + response_type + response_type:response_time | ",patient_id))}
    formula_fixed<-as.formula(paste0(c(
    paste0(c(paste0("response~-1+ response_type"), c("response_time", fixed_effects)), collapse = "+"),
    paste0(paste0(paste0(
      c("response_time", fixed_effects), ":response_type")
    ), collapse = "+")
  ), collapse = "+"))
  }


  #cl<-parallel::makeCluster(parallel::detectCores()-1,type="SOCK")
  cv_numbers <- unique(data_LME_model_dev[[cv_name]])

  model_LME <- lapply(cv_numbers, function(cv_number) {
    if (length(cv_numbers)>1){
      data_dev_cv <- data_LME_model_dev[data_LME_model_dev[[cv_name]]!= cv_number,]}
    if (length(cv_numbers)==1){
      data_dev_cv<-data_LME_model_dev}
    model_LME_cv<-nlme::lme(fixed=formula_fixed,
                            random=formula_random,
                            data = data_dev_cv,
                            weights = formula_weights,
                            control = lme_control)
    model_LME_cv$call$fixed<-formula_fixed
    model_LME_cv
  })
  names(model_LME) <- cv_numbers


  response_type <- Reduce(c,lapply(random_effects,function(i){rep(i,dim(data_LOCF)[1])}))
  response <- as.numeric(rep(NA,length(response_type)))
  response_time<-as.numeric(rep(x_L,length(response_type)))
  data_fixed_effects<-do.call("rbind",replicate(n=length(random_effects),data_LOCF[,c(patient_id,fixed_effects,cv_name)],simplify=FALSE))
  data_LOCF_model_val<-data.frame(data_fixed_effects,response_type,response,response_time,predict=1)
  data_LME_model_val$predict<-0
  data_LME_model_val<-dplyr::bind_rows(data_LOCF_model_val,data_LME_model_val)

  data_LME <-
    lapply(cv_numbers, function(cv_number) {
      data_LME_model_val_cv <-
        data_LME_model_val[which(data_LME_model_val[[cv_name]] == cv_number), ]
      data_LME_model_val_cv<-droplevels(data_LME_model_val_cv)
      model_LME_cv <- model_LME[[cv_number]]
      response_predictions <-
        which(data_LME_model_val_cv$predict == 1)
      data_LME_model_val_cv <-
        mixoutsamp(model = model_LME_cv,
                   newdata = data_LME_model_val_cv)$preddata[response_predictions, ][, c(patient_id,
                                                                                         fixed_effects,
                                                                                         cv_name,
                                                                                         "response_type",
                                                                                         "fitted")]
      data_LME_model_val_cv <-
          stats::reshape(
            data_LME_model_val_cv,
            timevar = "response_type",
            idvar = c(patient_id, fixed_effects, cv_name),
            direction = "wide"
          )
      for (name in random_effects) {
        names(data_LME_model_val_cv)[grep(paste0("fitted.",name), names(data_LME_model_val_cv))] <-
          name
      }
      data_LME_model_val_cv[[cv_name]]<-cv_number
      data_LME_model_val_cv
  })
  data_LME <- do.call("rbind", data_LME)

  data_LME<-data_LME[match(unique(data[[patient_id]]),data_LME[[patient_id]]),]
  data_LME<-data_LME[,order(match(names(data_LME),names(data)))]
  rownames(data_LME)<-NULL

  if (length(cv_numbers)==1){model_LME<-model_LME[[1]]}
  if(length(unique(data[[cv_name]]))==1){data[[cv_name]]<-NULL}

  return(list(
      data_longitudinal = data_LME,
      model_longitudinal = "LME",
      call = call,
      model_LME = model_LME,
      model_LME_standardise_time=standardise_time
    ))
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
      stop("event_status column should contain only values 0, 1, and 2 for cause_specific or fine_gray survival submodel,
        or values 0 and 1 for standard_cox survival submodel")
    }
  }
  if(survival_submodel %in% c("standard_cox")){
    if(!(setequal(data[[event_status]],0:1))){
      stop("event_status column should contain only values 0, 1, and 2 for cause_specific or fine_gray survival submodel,
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

#' Fit a landmarking model using a linear mixed effects (LME) model for the longitudinal submodel
#'
#' This function performs the two-stage landmarking analysis. In the first stage the longitudinal submodel is fitted using the LME model and in the
#' second stage the survival submodel is fitted.
#' Integer specifying the number of bootstrap samples to take
#' @param data_long Data frame or list of data frames each corresponding to a landmark age x_L (each element of the list must be named the value of x_L it corresponds to)
#' Each data frame contains repeat measurements data and time-to-event data in long format.
#' @template x_L
#' @template x_hor
#' @param standardise_time Boolean indicating whether to standarise the time variable by subtracting the mean
#' and dividing by the standard deviation (see Details section for more information)
#' @template event_status
#' @template event_time
#' @template k
#' @template cross_validation_df
#' @template patient_id
#' @template start_study_time
#' @template end_study_time
#' @template fixed_effects
#' @template random_effects
#' @template fixed_effects_time
#' @template random_effects_time
#' @param random_slope Boolean indicating whether to include a random slope in the LME model
#' @param lme_control Object created using `nlme::lmeControl()`, which will be passed to the `control` argument of the `lme`
#' function
#' @template b
#' @template survival_submodel
#' @return List containing `data`, `model_longitudinal`, `model_LME`, `model_LME_standardise_time`, `model_survival`, and `prediction_error`.
#'
#' `data` is a data frame that includes the LOCF values of the `fixed_effects` and estimates
#' of the `random_effects` predicted from the LME model (see Details section). It also includes the predicted
#' probability that the event of interest has occured by time \code{x_L}, labelled as \code{"event_prediction"}.
#' There is one row for each individual.
#'
#' `model_longitudinal` indicates that the longitudinal submodel is LME.
#'
#' `model_LME` contains the output from
#' the `lme` function from package `nlme`. For a model using cross-validation,
#' `model_LME` contains a list of outputs with each
#' element in the list corresponds to a different cross-validation fold.
#' `prediction_error` contains a list indicating the c-index and Brier score at time `x_hor` and their standard errors if parameter `b` is used.
#' For more information on how the prediction error is calculated
#' please see `?get_model_assessment` which is the function used to do this within `fit_LME_landmark_model`.
#'
#' `model_LME_standardise_time` contains a list of two objects `mean_response_time` and `sd_response_time` if the parameter `standardise_time=TRUE` is used. This
#' is the mean and standard deviation use to normalise times when fitting the LME model.
#'
#' `model_survival` contains the outputs from
#' the survival submodel functions, including the estimated parameters of the model. For a model using cross-validation,
#' `model_survival` will contain a list of outputs with each
#' element in the list corresponding to a different cross-validation fold.
#' `model_survival` contains the outputs from the function used to fit the survival submodel, including the estimated parameters of the model.
#' For a model using cross-validation, `model_survival` contains a list of outputs with each
#' element in the list corresponding to a different cross-validation fold. For more information on how the survival model is fitted
#' please see `?fit_survival_model` which is the function used to do this within `fit_LME_landmark_model`.
#'
#' `prediction_error` contains a list indicating the c-index and Brier score at time `x_hor` and their standard errors if parameter `b` is used.
#' @details
#'
#' There are two parts to fitting the landmark model: the longitudinal submodel and the survival submodel.
#' This function id used to fit the LME model as the longitudinal submodel. This model will be described in more detail.
#'
#' For an individual \eqn{i}, the LME model can be written as
#'
#' \deqn{Y_i = X_i \beta + Z_i U_i + \epsilon_i}
#'
#' where
#' * \eqn{Y_i} is the vector of outcomes at different time points for the individual
#' * \eqn{X_i} is the matrix of covariates for the fixed effects at these time points
#' * \eqn{\beta} is the vector of coefficients for the fixed effects
#' * \eqn{Z_i} is the matrix of covariates for the random effects
#' * \eqn{U_i} is the matrix of coefficients for the random effects
#' * \eqn{\epsilon_i} is the error term, typically from N(0, \eqn{\sigma})
#'
#' By using an LME model to fit repeat measures data we can allow measurements from the same individuals to be
#' more similar than measurements from different individuals. This is done through the random intercept and/or
#' random slope.
#'
#' Extending this model to the case where there are multiple random effects, denoted \eqn{k}, we have
#'
#' \deqn{Y_{ik} = X_{ik} \beta_k + Z_{ik} U_{ik} + \epsilon_{ik}}
#'
#' Using this model we can allow a covariance structure within the random effects term \eqn{U_{ik}}, for example a sample from the
#' multivariate normal (MVN) distribution \eqn{MVN(0,\Sigma_u)}. This covariance structure means the value of one random effects variable informs about the
#' value of the other random effects variables, leading to more accurate predictions and allowing there to be missing data in the
#' random effects variables.
#'
#' The function \code{fit_LME_landmark_model} uses this covariance structure for the random effects when fitting the LME model.
#' To fit the LME model the function \code{lme} from the package \code{nlme} is used.
#' The fixed effects are calculated as the LOCF for the variables \code{fixed_effects} at the landmark age \code{x_L} and the random effects
#' are those stated in \code{random_effects} and at times \code{random_effects_time}. This model is used to predict the
#' values of the random effects at the landmark time \code{x_L}.
#'
#' It is important to distinguish between the validation set and the development set for fitting the LME model. The development set includes
#' all the repeat measurements (including those after the landmark age \code{x_L}). Conversely, the validation set only includes
#' the repeat measurements recorded up until and including the landmark age \code{x_L}.
#'
#'
#' There is an important consideration about fitting the linear mixed effects model. As the variable \code{random_effects_time}
#' gets further from 0, the random effects coefficients get closer to 0. This causes computational issues
#' as the elements in the covariance matrix of the random effects, \eqn{\Sigma_u}, are constrained to
#' be greater than 0. Using parameter \code{standard_time=TRUE} can prevent this issue by standardising the
#' time variables to ensure that the \code{random_effects_time} values are not too close to 0.
#'
#' The predictions of the random effects at the landmark age, in addition to the LOCF values for the fixed effects,
#' are used as the covariates for the survival submodel.
#' The values of these covariates can be viewed in the data frame `data` that are returned by this function.
#'
#' For the survival submodel, there are three choices of model:
#' * the standard Cox model, this is a wrapper function for \code{coxph} from the package \code{survival}
#' * the cause-specific model, this is a wrapper function for \code{CSC} from package \code{riskRegression}
#' * the Fine Gray model, this is a wrapper function for \code{FGR} from package \code{riskRegression}
#'
#' The latter two models estimate the probability of the event of interest in the presence of competing events.
#'
#'
#' @author Isobel Barrott \email{isobel.barrott@@gmail.com}
#' @examples \dontrun{
#' #' library(Landmarking)
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
#' data_model_landmark_LME <-
#'   fit_LME_landmark_model(
#'     data_long = data_repeat_outcomes,
#'     x_L = c(60, 61),
#'     x_hor = c(65, 66),
#'     k = 10,
#'     start_study_time = "start_time",
#'     end_study_time = "event_time",
#'     fixed_effects = c("ethnicity", "smoking", "diabetes"),
#'     fixed_effects_time = "response_time_sbp_stnd",
#'     random_effects = c("sbp_stnd", "tchdl_stnd"),
#'     random_effects_time = c("response_time_sbp_stnd", "response_time_tchdl_stnd"),
#'     patient_id = "id",
#'     standardise_time = TRUE,
#'     lme_control = nlme::lmeControl(maxIter = 100, msMaxIter = 100),
#'     event_time = "event_time",
#'     event_status = "event_status",
#'     survival_submodel = "cause_specific"
#'   )
#' }
#' @importFrom stats as.formula
#' @importFrom prodlim Hist
#' @export


fit_LME_landmark_model<-function(data_long,
                                 x_L,
                                 x_hor,
                                 fixed_effects,
                                 random_effects,
                                 fixed_effects_time,
                                 random_effects_time,
                                 random_slope = TRUE,
                                 standardise_time=FALSE,
                                 patient_id,
                                 start_study_time,
                                 end_study_time,
                                 k,
                                 cross_validation_df,
                                 lme_control = nlme::lmeControl(),
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

  if (!(length(fixed_effects_time) %in% c(length(fixed_effects),1))){
    stop("Length of fixed_effects_time should be equal to length of fixed_effects or 1")}
  if (length(fixed_effects_time)==1){fixed_effects_time<-rep(fixed_effects_time,times=length(fixed_effects))}
  if (!(length(random_effects_time) %in% c(length(random_effects),1))){
    stop("Length of random_effects_time should be equal to length of random_effects or 1")}
  if (length(random_effects_time)==1){random_effects_time<-rep(random_effects_time,times=length(random_effects))}

  if(length(x_L)!=length(x_hor)){stop("Length of x_L should be the same as length of x_hor")}

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


      data_long_x_l<-data_long[[as.character(x_l)]]
      if(survival_submodel %in% c("cause_specific", "fine_gray")){
        if(!(setequal(data_long_x_l[[event_status]],0:2))){
          stop("event_status column should contain only values 0, 1, and 2 for cause_specific or fine_gray survival submodel,
          or values 0 and 1 for standard_cox survival submodel")
        }
      }
      if(survival_submodel %in% c("standard_cox")){
        if(!(setequal(data_long_x_l[[event_status]],0:1))){
          stop("event_status column should contain only values 0, 1, and 2 for cause_specific or fine_gray survival submodel,
          or values 0 and 1 for standard_cox survival submodel")
        }
      }

      if(!is.numeric(data_long_x_l[[event_status]])){stop("Column event_status should have class numeric")}

      if (!all(is.numeric(x_l))){
        stop("'x_L' should be numeric")}
      if (!all(is.numeric(x_h))){
        stop("'x_hor' should be numeric")}
      for (col in c(
        fixed_effects,
        fixed_effects_time,
        random_effects,
        random_effects_time,
        patient_id,
        start_study_time,
        end_study_time,
        event_time,
        event_status
      )) {
        if (!(col %in% names(data_long_x_l))) {
          stop(col, " is not a column name in data_long")
        }
      }
      data_long_x_l[[patient_id]]<-as.factor(data_long_x_l[[patient_id]])

      if(dim(return_ids_with_LOCF(data=data_long_x_l,
                                  patient_id=patient_id,
                                  x_L=x_l,
                                  covariates=fixed_effects,
                                  covariates_time=fixed_effects_time))[1]!=dim(data_long_x_l)[1]){
        stop("data_long contains individuals that do not have a LOCF for all fixed_effects.
             Use function return_ids_with_LOCF to remove these individuals from the dataset")
      }

      data_long_x_l<-data_long_x_l[data_long_x_l[[start_study_time]]<=x_l & data_long_x_l[[end_study_time]]>x_l,]
      data_long_x_l[[event_status]][data_long_x_l[[event_time]]>x_h]<-0
      data_long_x_l[[event_time]][data_long_x_l[[event_time]]>x_h]<-x_h

      if(cross_validation_df_add==TRUE){data_long_x_l<-
        dplyr::left_join(data_long_x_l,cross_validation_df[[as.character(x_l)]][,c(patient_id,"cross_validation_number")],by=patient_id)}
      return(data_long_x_l)
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
    print(paste0("Fitting longitudinal submodel, landmark age ", x_l))
    data_model_longitudinal<-fit_LME_longitudinal_model(data=data_long,
                                                        x_L=x_l,
                                                        fixed_effects=fixed_effects,
                                                        random_effects=random_effects,
                                                        fixed_effects_time=fixed_effects_time,
                                                        random_effects_time=random_effects_time,
                                                        random_slope = random_slope,
                                                        standardise_time=standardise_time,
                                                        cv_name=cv_name,
                                                        patient_id=patient_id,
                                                        lme_control = lme_control)
    print(paste0("Complete, landmark age ",x_l))

    data_events<-dplyr::distinct(data_long[,c(patient_id, event_status,event_time)])
    data_longitudinal<-dplyr::left_join(data_model_longitudinal$data_longitudinal,data_events,by=patient_id)

    print(paste0("Fitting survival submodel, landmark age ",x_l))

    data_model_survival<-fit_survival_model(data=data_longitudinal,
                                            patient_id=patient_id,
                                            cv_name=cv_name,
                                            covariates=c(fixed_effects,random_effects),
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
         model_LME=data_model_longitudinal$model_LME,
         model_LME_standardise_time=data_model_longitudinal$model_LME_standardise_time,
         model_survival=data_model_survival$model_survival,
         prediction_error=prediction_error,
         call=call)
  })
  names(out)<-x_L
  class(out)<-"landmark"
  out
}
