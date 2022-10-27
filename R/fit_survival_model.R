#' Fit a survival sub-model as part of a landmarking analysis
#'
#' This function is a helper function for `fit_LOCF_landmark_model` and `fit_LME_landmark_model`.
#'
#' @param data Data frame containing covariates and time-to-event data, one row for each individual.
#' @template x_hor
#' @template covariates
#' @param cv_name Character string specifying the column name in `data` that indicates cross-validation fold. If no cross-validation is needed, set this parameter to `NA`.
#' @param event_time Character string specifying the column name in `data` which contains the event time
#' @param event_status	Character string specifying the column name in `data` which contains the event status (where 0=censoring, 1=event of interest, if there are competing events these are labelled 2 or above). Events at time x_hor should be labelled censored.
#' @param individual_id Character string specifying the column name in `data` which contains the individual identifiers
#' @template survival_submodel
#' @return List containing `data_survival` and `model_survival`
#'
#' `data_survival` contains the predicted risk of event by the horizon time `x_hor`.
#'
#' `model_survival` contains the outputs from the function used to fit the survival submodel, including the estimated parameters of the model.
#' For a model using cross-validation, `model_survival` contains a list of outputs with each
#' element in the list corresponding to a different cross-validation fold.
#'
#' @details
#'
#' This function fits the survival model from the landmark model framework. The individuals are censored at the time horizon `x_hor` and the survival model is fitted with
#' covariates specified in parameter `covariates`.
#'
#' For the survival model, there are three choices of model:
#' * the standard Cox model, this is a wrapper function for \code{coxph} from the package \code{survival}
#' * the cause-specific model, this is a wrapper function for \code{CSC} from package \code{riskRegression}
#' * the Fine Gray model, this is a wrapper function for \code{FGR} from package \code{riskRegression}
#'
#' The latter two models estimate the probability of the event of interest in the presence of competing events.
#'
#' @author Isobel Barrott \email{isobel.barrott@@gmail.com}
#' @export

fit_survival_model <- function(data,
                               individual_id,
                               cv_name=NA,
                               covariates,
                               event_time,
                               event_status,
                               survival_submodel = c("standard_cox", "cause_specific", "fine_gray"),
                               x_hor) {
  #Checks
  #####
  if (!(inherits(data,"data.frame"))) {
    stop("data should be a dataframe")
  }
  if (!(inherits(x_hor,"numeric"))) {
    stop("x_hor should be numeric")
  }
  for (col in c(covariates,
                event_time,
                event_status,
                individual_id)) {
    if (!(col %in% names(data))) {
      stop(col, " is not a column name in data")
    }
    if(any(is.na(data[[col]]))){
      stop(col, " contains NA values")
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

  if(!is.null(levels(data[[event_status]]))){
    data[[event_status]]<-as.numeric(levels(data[[event_status]]))[data[[event_status]]]
  }

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

  data[[individual_id]] <- as.factor(data[[individual_id]])

  if(is.na(cv_name)){data[["cross_validation_number"]]<-1;cv_name<-"cross_validation_number"}

    cv_numbers <- unique(data[[cv_name]])
    model <- as.list(cv_numbers)
    names(model) <- cv_numbers
    Surv<-survival::Surv

    #Censor at the time horizon
    data[[event_status]][data[[event_time]] > x_hor] <-
      0
    data[[event_time]][data[[event_time]] > x_hor] <-
      x_hor

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
        data_test$event_prediction <- as.numeric(riskRegression::predictRisk(model_survival, times = x_hor, newdata = data_test))
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

  data_survival<-data_survival[order(match(data_survival[[individual_id]],data[[individual_id]])),]
  rownames(data_survival)<-NULL
  list(data_survival = data_survival, model_survival = model_survival)
}
