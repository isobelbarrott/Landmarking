#' Create landmark dataset
#'
#' Selects the rows of individuals at risk at age `x_L`,
#' i.e individuals that started follow-up before
#' (or at) landmark time `x_L`, and ended follow-up after
#' (but not at) landmark time `x_L`. These entries comprise
#' the dataset for landmark time `x_L` in the landmark model.
#' Also censors individuals with event times over `x_hor`.
#'
#'
#' @param data Data frame with repeated measurements in long format; each row corresponds to an assessment entry
#' @template x_L
#' @template x_hor
#' @param assessment_time Character string specifying the column name in
#' `data` which contains the assessment time for each row
#' @template patient_id
#' @template event_status
#' @template event_time
#' @return Data frame `data` updated to contain only entries of
#' individuals at risk at age `x_L`
#' @author Isobel Barrott \email{isobel.barrott@@gmail.com}
#' @details An individual enters the risk set at their first assessment and exits the risk set
#' after their first event (this may be censoring, event of interest or a competing event).
#' @examples
#' data(data_repeat_outcomes)
#' data_repeat_outcomes <- create_landmark_dataset(data=data_repeat_outcomes,
#'   x_L=60,
#'   assessment_time="response_time_sbp_stnd",
#'   patient_id="id",
#'   event_time="event_time",
#'   event_status="event_status")
#'
#' @export

create_landmark_dataset <-
  function(data,
           x_L,
           x_hor,
           patient_id,
           assessment_time,
           event_status,
           event_time){

    ######

    for (col in c(patient_id,
                  assessment_time,
                  event_time,
                  event_status)){
      if (!(col %in% names(data))){
        stop(col, " is not a column name in data")
      }
    }

    for (col in c(patient_id,
                  assessment_time,
                  event_time,
                  event_status)){
      if (any(is.na(data[[col]]))){
        stop(col, " cannot contain NA values")
      }
    }

    for (col in c(assessment_time,
                  event_time)){
      if (!(is.numeric(data[[col]]))){
        stop(col, " should be a numeric column in data")
      }
    }


    if(length(assessment_time)!=1){
      stop("assessment_time should be length 1. If measurements are made at different assessment_times,
           then use the parameters start_time and end_time to indicate when each patient enters and leaves the risk set")
    }

    data_start_time<-stats::aggregate(stats::as.formula(paste0(assessment_time,"~",patient_id)),data,function(x){min(x)})
    data_end_time<-stats::aggregate(stats::as.formula(paste0(event_time,"~",patient_id)),data,function(x){min(x)})
    data[data[[patient_id]] %in% data_start_time[data_start_time[[assessment_time]]<=x_L & data_end_time[[event_time]]>x_L,patient_id],]

    data[[event_status]][data[[event_time]]>x_hor]<-0
    data[[event_time]][data[[event_time]]>x_hor]<-x_hor

    data
  }

