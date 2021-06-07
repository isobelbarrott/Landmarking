

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
                                        x_hor,
                                        covariates,
                                        covariates_time,
                                        cv_name,
                                        patient_id) {
  call <- match.call()
  if (!(is.data.frame(data))) {
    stop("data should be a dataframe")
  }
  if (!(is.numeric(x_L))){
    stop("'x_L' should be numeric")}
  if (!(is.numeric(x_hor))){
    stop("'x_hor' should be numeric")}

  for (col in c(
    covariates,
    covariates_time,
    patient_id
  )) {
    if (!(col %in% names(data))) {
      stop(col, " is not a column name in data")
    }
  }

  if (!missing(cv_name)) {
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
    stop("data contains individuals that do not have a LOCF for all covariates. Use function return_ids_with_LOCF to remove these individuals from the dataset data.")
  }

  data[[patient_id]] <- as.character(data[[patient_id]])
  data_LOCF <- data
  data_LME <- data

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
  if (!missing(cv_name)){
    data_LOCF <-
      dplyr::left_join(data_LOCF, unique(data[c(patient_id, cv_name)]), by =
                         patient_id)
  }

  data_LOCF<-data_LOCF[,order(match(names(data_LOCF),names(data)))]

  list(data_longitudinal = data_LOCF,
       model_longitudinal="LOCF",
       call = call)
}




fit_LME_longitudinal_model <- function(data,
                                       x_L,
                                       x_hor,
                                       fixed_effects,
                                       random_effects,
                                       fixed_effects_time,
                                       random_effects_time,
                                       standardise_time=TRUE,
                                       random_intercept = TRUE,
                                       random_slope = TRUE,
                                       cv_name,
                                       patient_id,
                                       lme_control = nlme::lmeControl()) {
  call <- match.call()
  if (!(is.data.frame(data))) {
    stop("data should be a dataframe")
  }
  if (!(is.numeric(x_L))){
    stop("'x_L' should be numeric")}
  if (!(is.numeric(x_hor))){
    stop("'x_hor' should be numeric")}
  if (random_intercept & random_slope == FALSE) {
    stop("At least one of random_intercept/random_slope should be TRUE")
  }

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

  if (!missing(cv_name)) {
    if (!(cv_name %in% names(data))) {
      stop(cv_name, " is not a column name in data")
    }
    if (any(is.na(data[[cv_name]]))){
      stop("The column ",cv_name, " contains NA values")
    }
  }

  if (!(length(fixed_effects_time) %in% c(length(fixed_effects),1))){
    stop("Length of fixed_effects_time should be equal to length of fixed_effects or 1")}

  if (length(fixed_effects_time)==1){fixed_effects_time<-rep(fixed_effects_time,times=length(fixed_effects))}

  if (length(random_effects_time)!=length(random_effects)) {
    stop("Length of random_effects_time should be equal to length of random_effects")
  }

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

  data[[patient_id]] <- as.character(data[[patient_id]])

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
        patient_id =patient_id,
        x_L =
          x_L
      )
    })
  data_LOCF <- Reduce(merge, LOCF_values_by_variable)
  data_LOCF<-data_LOCF[match(unique(data[[patient_id]]),data_LOCF[[patient_id]]),]
  if (!missing(cv_name)){
    data_LOCF <-
      dplyr::left_join(data_LOCF, unique(data[c(patient_id, cv_name)]), by =
                         patient_id)
  }

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
  if(!missing(cv_name)){
    data_fixed_effects<-
      do.call("rbind", replicate(
        n = length(random_effects),
        data_LME[, c(patient_id, fixed_effects, cv_name)],
        simplify = FALSE
      ))}else{
        data_fixed_effects<-
          do.call("rbind", replicate(
            n = length(random_effects),
            data_LME[, c(patient_id, fixed_effects)],
            simplify = FALSE
          ))
      }
  data_LME <-
    data.frame(data_fixed_effects, response_type, response, response_time)

  if(standardise_time==TRUE){
    mean_response_time<-mean(data_LME$response_time,na.rm=TRUE)
    sd_response_time<-stats::sd(data_LME$response_time,na.rm=TRUE)
    data_LME$response_time<-(data_LME$response_time-mean_response_time)/sd_response_time
    x_L<-(x_L-mean_response_time)/sd_response_time
    x_hor<-(x_hor-mean_response_time)/sd_response_time
  }

  data_LME_model_val <- data_LME[data_LME$response_time <= x_L,]
  data_LME_model_dev <- data_LME[!is.na(data_LME$response),]
  #####

  formula_random <- c()
  if (random_intercept == TRUE) {
    formula_random <- "response_type"
  }
  if (random_slope == TRUE) {
    formula_random <- c(formula_random, paste0("response_type:response_time"))
  }
  formula_random<-as.formula(paste0(
    "~ -1 + ",
    paste0(formula_random, collapse = "+"),
    " |",
    patient_id
  ))
  formula_fixed<-as.formula(paste0(c(
    paste0(c(paste0("response~ -1 + response_type"), c("response_time", fixed_effects)), collapse = "+"),
    paste0(paste0(paste0(
      c("response_time", fixed_effects), ":response_type")
    ), collapse = "+")
  ), collapse = "+"))

  if (!missing(cv_name)) {
    cv_numbers <- unique(data_LME_model_dev[[cv_name]])

    model_LME <- lapply(cv_numbers, function(cv_number) {

      data_dev_cv <- data_LME_model_dev[data_LME_model_dev[[cv_name]]!= cv_number,]
      model_LME_cv<-nlme::lme(formula_fixed,
                              random=formula_random,
                              data = data_dev_cv,
                              weights = nlme::varIdent(form =  ~ 1 | "response_type"),
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
          data_LME_model_val[data_LME_model_val[[cv_name]] == cv_number, ]
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
          names(data_LME_model_val_cv)[grep(name, names(data_LME_model_val_cv))] <-
            name
        }
        data_LME_model_val_cv[cv_name]<-cv_number
        data_LME_model_val_cv
      })
    data_LME <- do.call("rbind", data_LME)
  }
  else{
    model_LME<-nlme::lme(formula_fixed,
                         random=formula_random,
                         data = data_LME_model_dev,
                         weights = nlme::varIdent(form =  ~ 1 | "response_type"),
                         control = lme_control)
    model_LME$call$fixed<-formula_fixed

    response_type <- Reduce(c,lapply(random_effects,function(i){rep(i,dim(data_LOCF)[1])}))
    response <- as.numeric(rep(NA,length(response_type)))
    response_time<-as.numeric(rep(x_L,length(response_type)))
    data_fixed_effects<-do.call("rbind",replicate(n=length(random_effects),data_LOCF[,c(patient_id,fixed_effects)],simplify=FALSE))
    data_LOCF_model_val<-data.frame(data_fixed_effects,response_type,response,response_time,predict=1)
    data_LME_model_val$predict<-0
    data_LME_model_val<-dplyr::bind_rows(data_LOCF_model_val,data_LME_model_val)

    response_predictions <-
      which(data_LME_model_val$predict == 1)
    data_LME_model_val <-
      mixoutsamp(model = model_LME,
                 newdata = data_LME_model_val)$preddata[response_predictions, ][, c(patient_id,
                                                                                    fixed_effects,
                                                                                    "response_type",
                                                                                    "fitted")]
    data_LME_model_val <-
      stats::reshape(
        data_LME_model_val,
        timevar = "response_type",
        idvar = c(patient_id, fixed_effects),
        direction = "wide"
      )
    for (name in random_effects) {
      names(data_LME_model_val)[grep(name, names(data_LME_model_val))] <-
        name
    }
    data_LME <- data_LME_model_val
  }
  data_LME<-data_LME[match(unique(data[[patient_id]]),data_LME[[patient_id]]),]
  data_LME<-data_LME[,order(match(names(data_LME),names(data)))]
  rownames(data_LME)<-NULL

  return(list(
      data_longitudinal = data_LME,
      model_longitudinal = model_LME,
      call = call
    ))
}

fit_survival_model <- function(data,
                               patient_id,
                               cv_name,
                               covariates,
                               event_time,
                               event_status,
                               cr_method = c("standard_cox", "cause_specific", "fine_gray"),
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
  if (!missing(cv_name)) {
    if (!(cv_name %in% names(data))) {
      stop(cv_name, " is not a column name in data")
    }
    if (any(is.na(data[[cv_name]]))){
      stop("The column ",cv_name, " contains NA values")
    }
  }


  #####

  cr_method <- match.arg(cr_method)

  data[[patient_id]] <- as.character(data[[patient_id]])

  if (!missing(cv_name)) {
    cv_numbers <- unique(data[[cv_name]])
    model <- as.list(cv_numbers)
    names(model) <- cv_numbers

    data_cv <- lapply(cv_numbers, function(cv_number) {

      data_test <- data[data[[cv_name]] == cv_number, ]
      data_train <- data[data[[cv_name]] != cv_number, ]

      model_survival <- c()

      if (cr_method == "standard_cox") {
        formula_survival <-
          as.formula(paste0("Surv(", event_time, ", ", event_status, "==1) ~",
                            paste0(covariates, collapse = "+")))
        model_survival <- coxph(formula_survival, data = data_train)
        data_test$event_prediction <-pec::predictSurvProb.coxph(model_survival, times = x_hor, newdata = data_test)
      }

      if (cr_method == "cause_specific") {
        if(any(table(data_train[[event_status]]))==0){
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

      if (cr_method == "fine_gray") {
        if(any(table(data_train[[event_status]]))==0){
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
      Reduce(dplyr::bind_rows, lapply(data_cv, `[[`, 1))
    model_survival <- lapply(data_cv, `[[`, 2)
    names(model_survival)<-cv_numbers

  }

  else{
    data_test <- data
    data_train <- data

    model_survival <- c()

    if (cr_method == "standard_cox") {
      formula_survival <-
        as.formula(paste0("Surv(", event_time, ", ", event_status, "==1) ~",
                          paste0(covariates, collapse = "+")))
      model_survival <- coxph(formula_survival, data = data_train)
      data_test$event_prediction <-pec::predictSurvProb.coxph(model_survival, times = x_hor, newdata = data_test)
    }

    if (cr_method == "cause_specific") {
      if(any(table(data_train[[event_status]]))==0){
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

    if (cr_method == "fine_gray") {
      if(any(table(data_train[[event_status]]))==0){
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

    data_survival<-data_test
  }
  data_survival<-data_survival[,order(match(names(data_survival),names(data)))]
  data_survival<-data_survival[order(match(names(data_survival),names(data))),]
  rownames(data_survival)<-NULL
  list(data_survival = data_survival, model_survival = model_survival)
}

#' Fit a landmarking model using a linear mixed effects (LME) model for the longitudinal submodel
#'
#' This function performs the two-stage landmarking analysis. In the first stage, the longitudinal submodel is fitted using the LME model and in the
#' second stage the survival submodel is fitted.
#'
#' @param data Data frame containing longitudinal (repeat measurements) data and time to event data, there may be more than one row per individual
#' @template x_L
#' @template x_hor
#' @param standardise_time Boolean indicating whether to standarise the `covariates_time` (see Details section for more information)
#' @template event_status
#' @template event_time
#' @template cv_name
#' @template patient_id
#' @template fixed_effects
#' @template random_effects
#' @template fixed_effects_time
#' @template random_effects_time
#' @param random_intercept Boolean indicating whether to include a random intercept in the LME model
#' @param random_slope Boolean indicating whether to include a random slope in the LME model
#' @param lme_control Object create using nlme::lmeControl(), which contains the `control` argument in the `lme`
#' function
#' @template cr_method
#' @return List containing `data`, `model_longitudinal`, and `model_survival`.  The
#' data frame `data` contains one row for each individual and includes the `covariates`
#' that are calculated using the LME model (more information in Details section).
#'
#' `model_longitudinal` contains the output from
#' the `lme` function from package `nlme`. For a model using cross-validation,
#' `model_longitudinal` will contain a list of outputs with each
#' element in the list corresponds to a different cross-validation fold.
#'
#' `model_survival` contains the outputs from
#' the `coxph` function from package `survival`. This output contains the
#' estimated parameters for the survival submodel. For a model using cross-validation,
#' `model_survival` will contain a list of outputs with each
#' element in the list corresponds to a different cross-validation fold.
#'
#' @details
#'
#' There are two parts to fitting the landmark model: the longitudinal (LME) submodel and the survival submodel.
#'
#'
#' Firstly, the longitudinal submodel (LME) is described in more detail.
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
#' Extending this model to the case where there are multiple random effects, denoted \eqn{k}
#'
#' \deqn{Y_{ik} = X_{ik} \beta_k + Z_{ik} U_{ik} + \epsilon_{ik}}
#'
#' Using this model we can allow a covariance structure within the random effects term \eqn{U_{ik}}, i.e. it follows the
#' multivariate normal (MVN) distribution \eqn{MVN(0,\Sigma_u)}. This means information from one random effects variable informs about the
#' value of the other random effects variables, leading to more accurate predictions and allowing there to be missing data in the
#' random effects variables.
#'
#' This function uses the LME model with a covariance structure for the random effects that has been described
#' for the longitudinal submodel. This is fitted using the function \code{lme} from the package \code{nlme}.
#' The fixed effects are calculated as the LOCF for the variables \code{fixed_effects} at the landmark age \code{x_L} and the random effects
#' are those stated in \code{random_effects} and times \code{random_effects_time}. This model is used to predict the
#' random effects at the landmark time \code{x_L}.
#'#'
#' It is important to distinguish between the validation set and the development set for fitting the LME. The development set includes
#' all the repeat measurements (including those after the landmark age \code{x_L}). Conversely, the validation set only includes
#' the repeat measurements recorded up until and including the landmark age \code{x_L}.
#'
#'
#' There is an important consideration about fitting the linear mixed effects model. If the variable \code{random_effects_time}
#' is far from 0 this means that the random effects coefficients can be close to 0. This causes computational issues
#' as the variance of the random effects is constrained to
#' be greater than zero. Therefore the parameter \code{standard_time=TRUE} is it standardises the
#' \code{random_effects_time} by subtracting the mean of this variable and dividing by the standard deviation.
#'
#' The predictions of the random effects, in addition to the LOCF values for the fixed effects
#' are the covariates for the survival submodel. These covariates can be viewed in the data frame `data` that are returned by this function.
#'
#' For the survival submodel, there are three choices of model: the standard Cox model, the cause-specific model and the Fine Gray model.
#' The latter two models estimate the probability of the event of interest in the presence of competing events.
#'
#' This function uses the function: \code{coxph} from the pacakge \code{survival} to fit the standard Cox model;
#' the function \code{CSC} from package \code{riskRegression} to fit the cause-specific model; and the function \code{FGR} from package
#' \code{riskRegression} to fit the Fine Gray model.
#'
#'
#'
#' @author Isobel Barrott \email{isobel.barrott@@gmail.com}
#' @examples \dontrun{2^3}
#' @importFrom stats as.formula
#' @importFrom prodlim Hist
#' @export


fit_LME_landmark_model<-function(data,
                                 x_L,
                                 x_hor,
                                 fixed_effects,
                                 random_effects,
                                 fixed_effects_time,
                                 random_effects_time,
                                 random_intercept = TRUE,
                                 random_slope = TRUE,
                                 standardise_time=FALSE,
                                 cv_name,
                                 patient_id,
                                 lme_control = nlme::lmeControl(),
                                 event_time,
                                 event_status,
                                 cr_method){
  call <- match.call()

  #Checks
  #####
  if (!(is.data.frame(data))) {
    stop("data should be a dataframe")
  }
  for (col in c(fixed_effects,random_effects,
                event_time,
                event_status,
                patient_id,
                fixed_effects_time,
                random_effects_time)) {
    if (!(col %in% names(data))) {
      stop(col, " is not a column name in data")
    }
  }
  if (!missing(cv_name)) {
    if (!(cv_name %in% names(data))) {
      stop(cv_name, " is not a column name in data")
    }
    if (any(is.na(data[[cv_name]]))){
      stop("The column ",cv_name, " contains NA values")
    }
  }
  #####
  if (!(is.numeric(x_L))){
    stop("'x_L' should be numeric")}
  if (!(is.numeric(x_hor))){
    stop("'x_hor' should be numeric")}
  if (random_intercept & random_slope == FALSE) {
    stop("At least one of random_intercept/random_slope should be TRUE")
  }


  if (!(length(fixed_effects_time) %in% c(length(fixed_effects),1))){
    stop("Length of fixed_effects_time should be equal to length of fixed_effects or 1")}

  if (length(fixed_effects_time)==1){fixed_effects_time<-rep(fixed_effects_time,times=length(fixed_effects))}

  if (length(random_effects_time)!=length(random_effects)) {
    stop("Length of random_effects_time should be equal to length of random_effects")
  }

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

  data_events<-dplyr::distinct(data[,c(patient_id,event_time,event_status)])

  print("Fitting longitudinal submodel")
  data_model_longitudinal<-fit_LME_longitudinal_model(data=data,
                                                      x_L=x_L,
                                                      x_hor=x_hor,
                                                      fixed_effects=fixed_effects,
                                                      random_effects=random_effects,
                                                      fixed_effects_time=fixed_effects_time,
                                                      random_effects_time=random_effects_time,
                                                      random_intercept = TRUE,
                                                      random_slope = TRUE,
                                                      standardise_time=standardise_time,
                                                      cv_name=cv_name,
                                                      patient_id=patient_id,
                                                      lme_control = lme_control)
  print("Complete")

  data_longitudinal_events<-dplyr::left_join(data_model_longitudinal$data_longitudinal,data_events,by=patient_id)
  data_longitudinal_events<-data_longitudinal_events[,order(match(names(data_longitudinal_events),names(data))),]

  print("Fitting survival submodel")

  data_model_survival<-fit_survival_model(data=data_longitudinal_events,
                                          patient_id=patient_id,
                                          cv_name=cv_name,
                                          covariates=c(fixed_effects,random_effects),
                                          event_time=event_time,
                                          event_status=event_status,
                                          cr_method = cr_method,
                                          x_hor=x_hor)
  print("Complete")

  list(data=data_model_survival$data_survival,
       model_longitudinal=data_model_longitudinal$model_longitudinal,
       model_survival=data_model_survival$model_survival
  )
}

#' Fit a landmarking model using a last observation carried forward (LOCF) model for the longitudinal submodel
#'
#' This function performs the two-stage landmarking analysis. In the first stage, the longitudinal submodel is fitted using the LOCF model and in the
#' second stage the survival submodel is fitted.
#'
#' @param data Data frame containing longitudinal (repeat measurements) data and time to event data, there may be more than one row per individual
#' @template x_L
#' @template x_hor
#' @template event_status
#' @template event_time
#' @template cv_name
#' @template covariates
#' @template covariates_time
#' @template patient_id
#' @template cr_method
#' @return List containing `data`, `model_longitudinal`, and `model_survival`.  The
#' data frame `data` contains one row for each individual and includes the `covariates`
#' that are calculated using the LOCF model (more information in Details section).
#'
#' `model_longitudinal` indicates that the longitudinal submodel is LOCF.
#'
#' `model_survival` contains the outputs from
#' the `coxph` function from package `survival`. This output contains the
#' estimated parameters for the survival submodel. For a model using cross-validation,
#' `model_survival` will contain a list of outputs with each
#' element in the list corresponding to a different cross-validation fold.
#'
#' @details
#' #' There are two parts to fitting the landmark model: the longitudinal (LME) submodel and the survival submodel.
#'
#'
#' For the survival submodel, there are three choices of model: the standard Cox model, the cause-specific model and the Fine Gray model.
#' The latter two models estimate the probability of the event of interest in the presence of competing events.
#'
#' This function uses the function: \code{coxph} from the pacakge \code{survival} to fit the standard Cox model;
#' the function \code{CSC} from package \code{riskRegression} to fit the cause-specific model; and the function \code{FGR} from package
#' \code{riskRegression} to fit the Fine Gray model.
#'
#'
#' @author Isobel Barrott \email{isobel.barrott@@gmail.com}
#' @examples \dontrun{data(data_landmark_cv)
#' data_model_landmark_LME<-fit_LME_landmark_model(data=data_landmark_cv,x_L=60,x_hor=65,
#'   fixed_effects=c("ethnicity","smoking","diabetes","deprivation","atrial_fibrillation"),
#'   fixed_effects_time=rep("response_time_sbp_stnd",length(fixed_effects_col)),
#'   random_effects=random_effects_col,random_effects_time=random_effects_time_col,
#'   cv_name="cross_validation_number",patient_id=patient_id_col,
#'   lme_control = nlme::lmeControl(maxIter=100,msMaxIter=100,
#'   event_time="event_time",event_status="event_status",cr_method = "cause_specific")}
#' @importFrom stats as.formula
#' @importFrom survival Surv
#' @importFrom survival coxph
#' @importFrom prodlim Hist
#' @export

fit_LOCF_landmark_model<-function(data,
                                  x_L,
                                  x_hor,
                                  covariates,
                                  covariates_time,
                                  cv_name,
                                  patient_id,
                                  event_time,
                                  event_status,
                                  cr_method = c("standard_cox", "cause_specific", "fine_gray")){
  call <- match.call()
  if (!(is.data.frame(data))) {
    stop("data should be a dataframe")
  }
  if (!(is.numeric(x_L))){
    stop("'x_L' should be numeric")}
  if (!(is.numeric(x_hor))){
    stop("'x_hor' should be numeric")}

  for (col in c(
    covariates,
    covariates_time,
    event_time,
    event_status,
    patient_id
  )) {
    if (!(col %in% names(data))) {
      stop(col, " is not a column name in data")
    }
  }

  if (!missing(cv_name)) {
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
    stop("data contains individuals that do not have a LOCF for all covariates.
         Use function return_ids_with_LOCF to remove these individuals from the dataset data.")
  }

  call <- match.call()
  data_events<-dplyr::distinct(data[,c(patient_id,event_time,event_status)])

  print("Fitting longitudinal submodel")

  data_model_longitudinal<-fit_LOCF_longitudinal_model(data=data,
                                                       x_L=x_L,
                                                       x_hor=x_hor,
                                                       covariates=covariates,
                                                       covariates_time=covariates_time,
                                                       cv_name=cv_name,
                                                       patient_id=patient_id)

  print("Complete")
  data_longitudinal_events<-dplyr::left_join(data_model_longitudinal$data_longitudinal,data_events,by=patient_id)
  data_longitudinal_events<-data_longitudinal_events[,order(match(names(data_longitudinal_events),names(data))),]

  print("Fitting survival submodel")

  data_model_survival<-fit_survival_model(data=data_longitudinal_events,
                                          patient_id=patient_id,
                                          cv_name=cv_name,
                                          covariates=covariates,
                                          event_time=event_time,
                                          event_status=event_status,
                                          cr_method = cr_method,
                                          x_hor=x_hor)
  print("Complete")

  list(data=data_model_survival$data_survival,
       model_longitudinal=data_model_longitudinal$model_longitudinal,
       model_survival=data_model_survival$model_survival
  )
}

