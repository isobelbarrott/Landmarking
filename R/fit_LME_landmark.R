#' Fit a landmarking model using a linear mixed effects (LME) model for the longitudinal submodel
#'
#' This function is a helper function for `fit_LME_landmark`.
#'
#' @param data_long Data frame containing repeat measurement data and time-to-event data in long format.
#' @template x_L
#' @param standardise_time Boolean indicating whether to standardise the time variable by subtracting the mean
#' and dividing by the standard deviation (see Details section for more information)
#' @param cv_name Character string specifying the column name in `data_long` that indicates cross-validation fold
#' @template individual_id
#' @template fixed_effects
#' @template random_effects
#' @template fixed_effects_time
#' @template random_effects_time
#' @param random_slope_longitudinal Boolean indicating whether to include a random slope in the LME model
#' @param random_slope_survival Boolean indicating whether to include the random slope estimate from the LME model
#' as a covariate in the survival submodel.
#' @param include_data_after_x_L Boolean indicating whether to include all longitudinal data, including data after the landmark age x_L,
#' in the model development dataset. See Details for more information.
#' @param lme_control Object created using `nlme::lmeControl()`, which will be passed to the `control` argument of the `lme` function
#' @return List containing elements:
#' `data_longitudinal`, `model_longitudinal`, `model_LME`, and `model_LME_standardise_time`.
#'
#' `data_longitudinal` has one row for each individual in the risk set at `x_L` and
#' contains the value of the covariates at the landmark time `x_L` of the `fixed_effects` using the LOCF model and
#' `random_effects` using the LME model.
#'
#' `model_longitudinal` indicates that the LME approach is used.
#'
#' `model_LME` contains the output from
#' the `lme` function from package `nlme`. For a model using cross-validation,
#' `model_LME` contains a list of outputs with each
#' element in the list corresponds to a different cross-validation fold.
#'
#' `model_LME_standardise_time` contains a list of two objects `mean_response_time` and `sd_response_time` if the parameter `standardise_time=TRUE` is used. This
#' is the mean and standard deviation used to normalise times when fitting the LME model.
#'
#' @details For an individual \eqn{i}, the LME model can be written as
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
#' By using an LME model to fit repeat measures data, rather than a linear model, we can allow measurements from the same individuals to be
#' more similar than measurements from different individuals. This is done through the random intercept and/or
#' random slope.
#'
#' Extending this model to the case where there are multiple random effects, denoted \eqn{k}, we have
#'
#' \deqn{Y_{ik} = X_{ik} \beta_k + Z_{ik} U_{ik} + \epsilon_{ik}}
#'
#' Using this model we can allow a certain covariance structure within the random effects term \eqn{U_{ik}}, for example a sample from the
#' multivariate normal (MVN) distribution \eqn{MVN(0,\Sigma_u)}. This covariance structure means the value of one random effects variable informs about the
#' value of the other random effects variables, leading to more accurate predictions and allowing there to be missing data in the
#' random effects variables.
#'
#' The function \code{fit_LME_landmark} uses a unstructured covariance for the random effects when fitting the LME model (i.e. no constraints are imposed on the values).
#' To fit the LME model the function \code{lme} from the package \code{nlme} is used.
#'
#' It is important to distinguish between the validation set and the development set for fitting the LME model in this function.
#' The development dataset either includes all the repeat measurements (including those after the landmark age \code{x_L}), or only the repeat measurements
#' recorded up to and including the landmark age \code{x_L}. This is controlled using the parameter `include_data_after_x_L`.
#' The validation set only includes the repeat measurements recorded up until and including the landmark age \code{x_L},
#' i.e. it does not include future data in its predictions.
#'
#' The random intercept is always included in the LME model.
#' Additionally, the random slope can be included in the LME model using the parameter `random_slope_longitudinal=TRUE`.
#'
#' Using the fitted model, the values of the best linear unbiased predictions (BLUPs)
#' at the landmark age \code{x_L} are calculated. These are considered the predictions of the values of the \code{random_effects}
#' the landmark age \code{x_L}. The values of the fixed effects used within this prediction are calculated as the LOCF values for the \code{fixed_effects}
#' at the landmark age \code{x_L}. These predictions are used as covariates in the survival model along with the LOCF values of the fixed effects.
#' Additionally, the estimated value of the random slope can
#' be included as predictors in the survival model using the parameter `random_slope_survival=TRUE`.
#'
#' There is an important consideration about fitting the linear mixed effects model. As the variable \code{random_effects_time}
#' gets further from 0, the random effects coefficients get closer to 0. This causes computational issues
#' as the elements in the covariance matrix of the random effects, \eqn{\Sigma_u}, are constrained to
#' be greater than 0. Using parameter \code{standard_time=TRUE} can prevent this issue by standardising the
#' time variables to ensure that the \code{random_effects_time} values are not too close to 0.
#'
#' The LOCF values for the fixed effects and the prediction of the random effects at the landmark age
#' are used as the covariates for the survival submodel, in addition to the estimated random slopes
#' if option `random_effects_as_covariate` is selected.
#'
#' @export
fit_LME_longitudinal <- function(data_long,
                                 x_L,
                                 fixed_effects,
                                 random_effects,
                                 fixed_effects_time,
                                 random_effects_time,
                                 standardise_time = FALSE,
                                 random_slope_longitudinal = TRUE,
                                 random_slope_survival = TRUE,
                                 include_data_after_x_L=TRUE,
                                 cv_name = NA,
                                 individual_id,
                                 lme_control = nlme::lmeControl()) {
  call <- match.call()
  if (!(is.data.frame(data_long))) {
    stop("data_long should be a dataframe")
  }
  if (!(is.numeric(x_L))) {
    stop("'x_L' should be numeric")
  }

  for (col in c(
    fixed_effects,
    random_effects,
    fixed_effects_time,
    random_effects_time,
    individual_id
  )) {
    if (!(col %in% names(data_long))) {
      stop(col, " is not a column name in data_long")
    }
  }

  if (!is.na(cv_name)) {
    if (!(cv_name %in% names(data_long))) {
      stop(cv_name, " is not a column name in data_long")
    }
    if (any(is.na(data_long[[cv_name]]))) {
      stop("The column ", cv_name, " contains NA values")
    }
  }

  if (is.na(cv_name)) {
    data_long[["cross_validation_number"]] <-
      1
    cv_name <- "cross_validation_number"
  }

  if (!(length(fixed_effects_time) %in% c(length(fixed_effects), 1))) {
    stop("Length of fixed_effects_time should be equal to length of fixed_effects or 1")
  }

  if (length(fixed_effects_time) == 1) {
    fixed_effects_time <-
      rep(fixed_effects_time, times = length(fixed_effects))
  }

  if (!(length(random_effects_time) %in% c(length(random_effects), 1))) {
    stop("Length of random_effects_time should be equal to length of random_effects or 1")
  }

  if (length(random_effects_time) == 1) {
    random_effects_time <-
      rep(random_effects_time, times = length(random_effects))
  }

  if (dim(
    return_ids_with_LOCF(
      data_long = data_long,
      individual_id = individual_id,
      x_L = x_L,
      covariates = fixed_effects,
      covariates_time = fixed_effects_time
    )
  )[1] != dim(data_long)[1]) {
    stop(
      "data_long contains individuals that do not have a LOCF for all fixed_effects. Use function return_ids_with_LOCF to remove these individuals from the dataset data_long."
    )
  }

  data_long[[individual_id]] <-
    as.factor(data_long[[individual_id]])

  data_LOCF <- data_long
  data_LME <- data_long

  #Pick out LOCF for each variable
  LOCF_values_by_variable <-
    lapply(1:length(c(fixed_effects, random_effects)), function(x) {
      return_LOCF_by_variable(
        data_long = data_LOCF,
        i = x,
        covariates = c(fixed_effects, random_effects),
        covariates_time = c(fixed_effects_time, random_effects_time),
        individual_id = individual_id,
        x_L =
          x_L
      )
    })
  data_LOCF <- Reduce(merge, LOCF_values_by_variable)

  data_LOCF <-
    dplyr::left_join(data_LOCF, unique(data_long[c(individual_id, cv_name)]), by =
                       individual_id)

  #Create validation and development dataset
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
  data_fixed_effects <-
    do.call("rbind", replicate(
      n = length(random_effects),
      data_LME[, c(individual_id, fixed_effects, cv_name)],
      simplify = FALSE
    ))
  data_LME <-
    data.frame(data_fixed_effects, response_type, response, response_time)

  #Standardise response time
  if (standardise_time == TRUE) {
    mean_response_time <- mean(data_LME$response_time, na.rm = TRUE)
    sd_response_time <- stats::sd(data_LME$response_time, na.rm = TRUE)
  } else{
    mean_response_time <- 0
    sd_response_time <- 1
  }
  data_LME$response_time <-
    (data_LME$response_time - mean_response_time) / sd_response_time
  x_L <- (x_L - mean_response_time) / sd_response_time
  standardise_time <-
    list(mean_response_time = mean_response_time,
         sd_response_time = sd_response_time)

  data_LME_model_val <- data_LME[data_LME$response_time <= x_L, ] #validation set
  data_LME_model_dev <- data_LME[!is.na(data_LME$response), ] #development set

  if (include_data_after_x_L==FALSE){data_LME_model_dev <- data_LME_model_dev[data_LME_model_dev$response_time <= x_L, ]}

  if (length(random_effects) == 1){
    formula_weights <- NULL
    if (random_slope_longitudinal == FALSE) {
      formula_random <- as.formula(paste0(" ~ 1 |", individual_id))
    }else{
      formula_random <-
        as.formula(paste0(" ~ 1 + response_time |", individual_id))
    }
    formula_fixed <-
      as.formula(paste0(c(
        paste0("response ~ 1 "), c("response_time", fixed_effects)
      ), collapse = "+"))
  }
  if (length(random_effects) > 1) {
    formula_weights <- nlme::varIdent(form = ~ 1 | "response_type")
    if (random_slope_longitudinal == FALSE) {
      formula_random <-
        as.formula(paste0(" ~ -1 + response_type | ", individual_id))
    }
    else{
      formula_random <-
        as.formula(paste0(
          " ~ -1 + response_type + response_type:response_time | ",
          individual_id
        ))
    }
    formula_fixed <- as.formula(paste0(c(
      paste0(c(
        paste0("response~-1+ response_type"),
        c("response_time", fixed_effects)
      ), collapse = "+"),
      paste0(paste0(paste0(
        c("response_time", fixed_effects), ":response_type"
      )), collapse = "+")
    ), collapse = "+"))
  }

  cv_numbers <- unique(data_LME_model_dev[[cv_name]])
  model_LME <- lapply(cv_numbers, function(cv_number) {
    if (length(cv_numbers) > 1) {
      data_dev_cv <-
        data_LME_model_dev[data_LME_model_dev[[cv_name]] != cv_number, ]
    }
    if (length(cv_numbers) == 1) {
      data_dev_cv <- data_LME_model_dev
    }
    model_LME_cv <- nlme::lme(
      fixed = formula_fixed,
      random = formula_random,
      data = data_dev_cv,
      weights = formula_weights,
      control = lme_control
    )
    model_LME_cv$call$fixed <- formula_fixed
    model_LME_cv
  })
  names(model_LME) <- cv_numbers

  response_type <-
    Reduce(c, lapply(random_effects, function(i) {
      rep(i, dim(data_LOCF)[1])
    }))
  response <- as.numeric(rep(NA, length(response_type)))
  response_time <- as.numeric(rep(x_L, length(response_type)))
  data_fixed_effects <-
    do.call("rbind", replicate(
      n = length(random_effects),
      data_LOCF[, c(individual_id, fixed_effects, cv_name)],
      simplify = FALSE
    ))
  data_LOCF_model_val <-
    data.frame(data_fixed_effects,
               response_type,
               response,
               response_time,
               predict = 1)
  data_LME_model_val$predict <- 0
  data_LME_model_val <-
    dplyr::bind_rows(data_LOCF_model_val, data_LME_model_val)

  data_LME <-
    lapply(cv_numbers, function(cv_number) {
      data_LME_model_val_cv <-
        data_LME_model_val[which(data_LME_model_val[[cv_name]] == cv_number),]
      data_LME_model_val_cv <- droplevels(data_LME_model_val_cv)
      model_LME_cv <- model_LME[[as.character(cv_number)]]
      response_predictions <-
        which(data_LME_model_val_cv$predict == 1)

      mixoutsamp_LME_model_val_cv <-
        mixoutsamp(model = model_LME_cv,
                   newdata = data_LME_model_val_cv)

      data_LME_model_val_cv <-
        mixoutsamp_LME_model_val_cv$preddata[response_predictions,][, c(individual_id,
                                                                        fixed_effects,
                                                                        cv_name,
                                                                        "response_type",
                                                                        "fitted")]
      data_LME_model_val_cv <-
        stats::reshape(
          data_LME_model_val_cv,
          timevar = "response_type",
          idvar = c(individual_id, fixed_effects, cv_name),
          direction = "wide"
        )

      data_LME_model_val_cv<-dplyr::rename_with(data_LME_model_val_cv,~random_effects[which(paste0("fitted.",random_effects)==.x)],.cols=paste0("fitted.",random_effects))

      if (random_slope_survival == TRUE) {
        if (length(random_effects)==1){
          slopes_df<-mixoutsamp_LME_model_val_cv$random[,c("id","reffresponse_time")]
          slopes_df["reffresponse_time"]<-slopes_df["reffresponse_time"]+model_LME_cv$coefficients$fixed["response_time"]
        }
        if (length(random_effects)>1){
          slopes_df<-mixoutsamp_LME_model_val_cv$random[,c("id",paste0("reffresponse_type",random_effects,":response_time"))]
          for (i in 1:length(random_effects)){
            slopes_df[,paste0("reffresponse_type",random_effects[i],":response_time")]<-
               slopes_df[,paste0("reffresponse_type",random_effects[i],":response_time")]+model_LME_cv$coefficients$fixed[paste0("response_type",random_effects[i])]
          }
        }
      names(slopes_df)<-c("id",paste0(random_effects,"_slope"))
      data_LME_model_val_cv <- dplyr::left_join(data_LME_model_val_cv,
                                                slopes_df,
                                                by = individual_id)
      }

      data_LME_model_val_cv[[as.character(cv_name)]] <- cv_number
      data_LME_model_val_cv
    })
  data_LME <- do.call("rbind", data_LME)

  data_LME <-
    data_LME[match(unique(data_long[[individual_id]][data_long[[individual_id]] %in% data_LME[[individual_id]]]), data_LME[[individual_id]]), ]
  data_LME <-
    data_LME[, order(match(names(data_LME), names(data_long)))]
  rownames(data_LME) <- NULL

  if (length(cv_numbers) == 1) {
    model_LME <- model_LME[[1]]
  }
  if (length(unique(data_long[[cv_name]])) == 1) {
    data_LME[[cv_name]] <- NULL
  }

  return(
    list(
      data_longitudinal = data_LME,
      model_longitudinal = "LME",
      call = call,
      model_LME = model_LME,
      model_LME_standardise_time = standardise_time
    )
  )
}


#' Fit a landmarking model using a linear mixed effects (LME) model for the longitudinal submodel
#'
#' This function performs the two-stage landmarking analysis.
#'
#' @param data_long Data frame or list of data frames each corresponding to a landmark age `x_L` (each element of the list must be named the value of `x_L` it corresponds to).
#' Each data frame contains repeat measurements data and time-to-event data in long format.
#' @template x_L
#' @template x_hor
#' @param standardise_time Boolean indicating whether to standardise the time variables (`fixed_effects_time` and `random_effects_time`) by subtracting the mean
#' and dividing by the standard deviation (see Details section for more information)
#' @template event_status
#' @template event_time
#' @template k
#' @template cross_validation_df
#' @template individual_id
#' @template fixed_effects
#' @template random_effects
#' @template fixed_effects_time
#' @template random_effects_time
#' @param random_slope_longitudinal Boolean indicating whether to include a random slope in the LME model
#' @param random_slope_survival Boolean indicating whether to include the random slope estimate from the LME model
#' as a covariate in the survival submodel.
#' @param lme_control Object created using `nlme::lmeControl()`, which will be passed to the `control` argument of the `lme`
#' function
#' @template b
#' @template survival_submodel
#' @return List containing containing information about the landmark model at each of the landmark times.
#' Each element of this list is named the corresponding landmark time, and is itself a list containing elements:
#' `data`, `model_longitudinal`, `model_LME`, `model_LME_standardise_time`, `model_survival`, and `prediction_error`.
#'
#' `data` has one row for each individual in the risk set at `x_L` and
#' contains the value of the `fixed_effects` using the LOCF approach and predicted values of the
#' `random_effects` using the LME model at the landmark time `x_L`. It also includes the predicted
#' probability that the event of interest has occurred by time \code{x_hor}, labelled as \code{"event_prediction"}.
#' There is one row for each individual.
#'
#' `model_longitudinal` indicates that the longitudinal approach is LME.
#'
#' `model_LME` contains the output from
#' the `lme` function from package `nlme`. For a model using cross-validation,
#' `model_LME` contains a list of outputs with each
#' element in the list corresponds to a different cross-validation fold.
#' `prediction_error` contains a list indicating the c-index and Brier score at time `x_hor` and their standard errors if parameter `b` is used.
#' For more information on how the prediction error is calculated
#' please see `?get_model_assessment` which is the function used to do this within `fit_LME_landmark`.
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
#' please see `?fit_survival_model` which is the function used to do this within `fit_LME_landmark`.
#'
#' `prediction_error` contains a list indicating the c-index and Brier score at time `x_hor` and their standard errors if parameter `b` is used.
#' @details Firstly, this function selects the individuals in the risk set at the landmark time \code{x_L}.
#' Specifically, the individuals in the risk set are those that have entered the study before the landmark age
#' (there is at least one observation for each of the \code{fixed_effects} and\code{random_effects} on or before \code{x_L}) and
#' exited the study on after the landmark age (\code{event_time}
#' is greater than \code{x_L}).
#'
#' Secondly, if the option to use cross validation
#' is selected (using either parameter `k` or `cross_validation_df`), then an extra column `cross_validation_number` is added with the
#' cross-validation folds. If parameter `k` is used, then the function `add_cv_number`
#' randomly assigns these folds. For more details on this function see `?add_cv_number`.
#' If the parameter `cross_validation_df` is used, then the folds specified in this data frame are added.
#' If cross-validation is not selected then the landmark model is
#' fit to the entire group of individuals in the risk set (this is both the training and test dataset).
#'
#' Thirdly, the landmark model is then fit to each of the training data. There are two parts to fitting the landmark model: using the longitudinal data and using the survival data.
#' Using the longitudinal data is the first stage and is performed using `fit_LME_longitudinal`. See `?fit_LME_longitudinal` more for information about this function.
#' Using the survival data is the second stage and is performed using `fit_survival_model`. See `?fit_survival_model` more for information about this function.
#'
#' Fourthly, the performance of the model is then assessed on the set of predictions
#' from the entire set of individuals in the risk set by calculating Brier score and C-index.
#' This is performed using `get_model_assessment`. See `?get_model_assessment` more for information about this function.
#'
#' @author Isobel Barrott \email{isobel.barrott@@gmail.com}
#' @examples \donttest{
#' library(Landmarking)
#' data(data_repeat_outcomes)
#' data_model_landmark_LME <-
#'   fit_LME_landmark(
#'     data_long = data_repeat_outcomes,
#'     x_L = c(60, 61),
#'     x_hor = c(65, 66),
#'     k = 10,
#'     fixed_effects = c("ethnicity", "smoking", "diabetes"),
#'     fixed_effects_time = "response_time_sbp_stnd",
#'     random_effects = c("sbp_stnd", "tchdl_stnd"),
#'     random_effects_time = c("response_time_sbp_stnd", "response_time_tchdl_stnd"),
#'     individual_id = "id",
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


fit_LME_landmark <- function(data_long,
                             x_L,
                             x_hor,
                             fixed_effects,
                             random_effects,
                             fixed_effects_time,
                             random_effects_time,
                             random_slope_longitudinal=TRUE,
                             random_slope_survival=TRUE,
                             include_data_after_x_L=TRUE,
                             individual_id,
                             k,
                             cross_validation_df,
                             standardise_time = FALSE,
                             lme_control = nlme::lmeControl(),
                             event_time,
                             event_status,
                             survival_submodel = c("standard_cox", "cause_specific", "fine_gray"),
                             b) {
  call <- match.call()

  survival_submodel <- match.arg(survival_submodel)


  #Checks
  if (missing(k)) {
    k_add <- FALSE
  }
  else{
    k_add <- TRUE
    if (!(is.numeric(k))) {
      stop("k should be numeric")
    }
  }

  if (missing(cross_validation_df)) {
    cross_validation_df_add <- FALSE
  }else{
    cross_validation_df_add <- TRUE
    if (class(cross_validation_df) == "list") {
      if (!all(x_L %in% names(cross_validation_df))) {
        stop(
          "The names of elements in cross_validation_df list should be the landmark times in x_L"
        )
      }
      if (any(Reduce("c", lapply(cross_validation_df, function(x) {
        any(duplicated(dplyr::distinct(x[, c(individual_id, "cross_validation_number")])[, individual_id]))
      })))) {
        stop("Cross validation folds should be the same for the same individual")
      }
    }
    else if (class(cross_validation_df) == "data.frame") {
      if (any(duplicated(dplyr::distinct(cross_validation_df[, c(individual_id, "cross_validation_number")])[, individual_id]))) {
        stop("Cross validation folds should be the same for the same individual")
      }
    }
    else{
      stop("cross_validation_df should be either a data frame or a list")
    }
  }
  if (k_add == TRUE &&
      cross_validation_df_add == TRUE) {
    stop("Either use parameter k or cross_validation_df but not both")
  }
  if (k_add == FALSE &&
      cross_validation_df_add == FALSE) {
    cv_name <- NA
  } else{
    cv_name <- "cross_validation_number"
  }

  if (!(length(fixed_effects_time) %in% c(length(fixed_effects), 1))) {
    stop("Length of fixed_effects_time should be equal to length of fixed_effects or 1")
  }
  if (length(fixed_effects_time) == 1) {
    fixed_effects_time <-
      rep(fixed_effects_time, times = length(fixed_effects))
  }
  if (!(length(random_effects_time) %in% c(length(random_effects), 1))) {
    stop("Length of random_effects_time should be equal to length of random_effects or 1")
  }
  if (length(random_effects_time) == 1) {
    random_effects_time <-
      rep(random_effects_time, times = length(random_effects))
  }

  if (length(x_L) != length(x_hor)) {
    stop("Length of x_L should be the same as length of x_hor")
  }

  if (!(is.data.frame(data_long) ||
        is.list(data_long))) {
    stop("data_long should be a list or data.frame")
  }
  if (is.data.frame(data_long)) {
    data_long <- lapply(x_L, function(x_l) {
      data_long
    })
    names(data_long) <- x_L
  }
  if (is.list(data_long)) {
    if (!setequal(names(data_long), x_L)) {
      stop("Names of elements in data_long should be landmark ages x_L")
    }
  }

  if (missing(b)) {
    b <- NA
  }


  #Fit each landmark model
  data_long_x_L <- lapply(1:length(x_L), function(i) {
    x_l <- x_L[i]
    x_h <- x_hor[i]

    data_long_x_l <- data_long[[as.character(x_l)]]

    if (!is.null(levels(data_long_x_l[[event_status]]))) {
      data_long_x_l[[event_status]] <-
        as.numeric(levels(data_long_x_l[[event_status]]))[data_long_x_l[[event_status]]]
    }

    if (survival_submodel %in% c("cause_specific", "fine_gray")) {
      if (!(setequal(data_long_x_l[[event_status]], 0:2))) {
        stop(
          "event_status column should contain only values 0, 1, and 2 for cause_specific or fine_gray survival submodel,
          or values 0 and 1 for standard_cox survival submodel"
        )
      }
    }
    if (survival_submodel %in% c("standard_cox")) {
      if (!(setequal(data_long_x_l[[event_status]], 0:1))) {
        stop(
          "event_status column should contain only values 0, 1, and 2 for cause_specific or fine_gray survival submodel,
          or values 0 and 1 for standard_cox survival submodel"
        )
      }
    }

    if (!all(is.numeric(x_l))) {
      stop("'x_L' should be numeric")
    }
    if (!all(is.numeric(x_h))) {
      stop("'x_hor' should be numeric")
    }
    for (col in c(
      fixed_effects,
      fixed_effects_time,
      random_effects,
      random_effects_time,
      individual_id,
      event_time,
      event_status
    )) {
      if (!(col %in% names(data_long_x_l))) {
        stop(col, " is not a column name in data_long")
      }
    }
    data_long_x_l[[individual_id]] <-
      as.factor(data_long_x_l[[individual_id]])

    #Pull out individuals in the risk set
    data_long_x_l_risk_set <-
      return_ids_with_LOCF(
        data_long = data_long_x_l,
        individual_id = individual_id,
        x_L = x_l,
        covariates = fixed_effects,
        covariates_time = fixed_effects_time
      )
    data_long_x_l_risk_set <-
      data_long_x_l_risk_set[data_long_x_l_risk_set[[event_time]] > x_l,]
    n <-
      length(unique(data_long_x_l[[individual_id]])) - length(unique(data_long_x_l_risk_set[[individual_id]]))
    if (n >= 1) {
      warning(
        n,
        " individuals have been removed from the model building as they are not in the risk set at landmark age ",
        x_l
      )
    }
    data_long_x_l <- data_long_x_l_risk_set
    #Censor at the time horizon
    data_long_x_l[[event_status]][data_long_x_l[[event_time]] > x_h] <-
      0
    data_long_x_l[[event_time]][data_long_x_l[[event_time]] > x_h] <-
      x_h

    if (cross_validation_df_add == TRUE) {
      data_long_x_l <-
        dplyr::left_join(data_long_x_l, cross_validation_df[[as.character(x_l)]][, c(individual_id, "cross_validation_number")], by =
                           individual_id)
    }
    return(data_long_x_l)
  })
  names(data_long_x_L) <- x_L
  #Add cross-validation folds
  if (k_add == TRUE) {
    data_long_x_L_cv <-
      add_cv_number(
        data_long = Reduce("rbind", data_long_x_L),
        individual_id = individual_id,
        k = k
      )
    data_long_x_L <-
      lapply(data_long_x_L, function(x) {
        dplyr::left_join(x, dplyr::distinct(data_long_x_L_cv[, c(individual_id, "cross_validation_number")]), by =
                           individual_id)
      }
    )
  }

  out <- lapply(1:length(x_L), function(i) {
    x_l <- x_L[i]
    x_h <- x_hor[i]

    data_long <- data_long_x_L[[as.character(x_l)]]
    print(paste0("Fitting longitudinal submodel, landmark age ", x_l))
    data_model_longitudinal <-
      fit_LME_longitudinal(
        data_long = data_long,
        x_L = x_l,
        fixed_effects =
          fixed_effects,
        random_effects =
          random_effects,
        fixed_effects_time =
          fixed_effects_time,
        random_effects_time =
          random_effects_time,
        random_slope_longitudinal = random_slope_longitudinal,
        random_slope_survival =
          random_slope_survival,
        include_data_after_x_L=include_data_after_x_L,
        standardise_time =
          standardise_time,
        cv_name = cv_name,
        individual_id =
          individual_id,
        lme_control = lme_control
      )
    print(paste0("Complete, landmark age ", x_l))

    data_events <-
      dplyr::distinct(data_long[, c(individual_id, event_status, event_time)])
    data_longitudinal <-
      dplyr::left_join(data_model_longitudinal$data_longitudinal,
                       data_events,
                       by = individual_id)

    print(paste0("Fitting survival submodel, landmark age ", x_l))

    if (random_slope_survival) {
      random_effects <- c(random_effects, paste0(random_effects, "_slope"))

    }
    data_model_survival <- fit_survival_model(
      data = data_longitudinal,
      individual_id = individual_id,
      cv_name = cv_name,
      covariates = c(fixed_effects, random_effects),
      event_time = event_time,
      event_status = event_status,
      survival_submodel = survival_submodel,
      x_hor = x_h
    )
    print(paste0("Complete, landmark age ", x_l))

    prediction_error <-
      get_model_assessment(
        data = data_model_survival$data_survival,
        individual_id = individual_id,
        event_prediction = "event_prediction",
        event_status = event_status,
        event_time = event_time,
        x_hor =  max(setdiff(data_model_survival$data_survival[[event_time]],x_hor)),
        b = b
      )

    list(
      data = data_model_survival$data_survival,
      model_longitudinal = data_model_longitudinal$model_longitudinal,
      model_LME = data_model_longitudinal$model_LME,
      model_LME_standardise_time = data_model_longitudinal$model_LME_standardise_time,
      model_survival = data_model_survival$model_survival,
      prediction_error = prediction_error,
      call = call
    )
  })
  names(out) <- x_L
  class(out) <- "landmark"
  out
}
