#' Simulated CVD risk assessment data
#' This is simulated data formed using an LME model.
#'
#' @format A dataset containing repeat measurements data from a CVD risk assessment from 10000 patients.
#'
#'
#' A dataset with 30018 rows and 12 columns:
#' \describe{
#'   \item{id}{Patient ID}
#'   \item{smoking}{Smoking status, 0 indicates the patient has 'never smoked', 1 indicated the patient has quit smoking, and 2 indicated the individual is a current smoker.}
#'   \item{diabetes}{Diabetes status, 0 indicates the patient does not have diabetes, 1 indicates the patient does have diabetes}
#'   \item{ethnicity}{Ethnicity, one of five ethnicities}
#'   \item{deprivation}{Deprivation score out of five}
#'   \item{index}{An index indicating assessment number for a patient}
#'   \item{response_date_sbp_stnd}{The date that systolic blood pressure was recorded, this is the same as the date that the fixed measured were recorded}
#'   \item{sbp_stnd}{Standardised systolic blood pressure}
#'   \item{response_date_tchdl_stnd}{The date that systolic total cholesterol to HDL ratio was recorded}
#'   \item{tchdl_stnd}{Standardised total cholesterol to HDL ratio}
#' }
#'
#'
#'
"data_repeat"
#' Simulated CVD event data
#'
#' A dataset containing CVD event data for 4100 patients at risk at landmark age 60. This is simulated data formed using an cause-specific survival model.
#'
#' A dataset with 4100 rows and 3 columns:
#' \describe{
#'   \item{id}{Patient ID}
#'   \item{event_status}{Event status, 0 indicates censoring, 1 indicates CVD event, 2 indicates death from other causes}
#'   \item{event_time}{Event time}
#' }
#'
"data_outcomes"
#' Simulated CVD event data
#'
#' A dataset containing CVD event data for 4100 patients at risk at landmark age 60. This is simulated data formed using an cause-specific survival model.
#'
#' A dataset with 4100 rows and 3 columns:
#' \describe{
#'   \item{id}{Patient ID}
#'   \item{event_status}{Event status, 0 indicates censoring, 1 indicates CVD event, 2 indicates death from other causes}
#'   \item{event_time}{Event time}
#' }
#'
"data_landmark"

#' Simulated CVD event data
#'
#' A dataset containing CVD event data for 4100 patients at risk at landmark age 60. This is simulated data formed using an cause-specific survival model.
#'
#' A dataset with 4100 rows and 3 columns:
#' \describe{
#'   \item{id}{Patient ID}
#'   \item{event_status}{Event status, 0 indicates censoring, 1 indicates CVD event, 2 indicates death from other causes}
#'   \item{event_time}{Event time}
#' }
#'
"data_landmark_cv"

#' Simulated CVD event data
#'
#' A dataset containing CVD event data for 4100 patients at risk at landmark age 60. This is simulated data formed using an cause-specific survival model.
#'
#' A dataset with 4100 rows and 3 columns:
#' \describe{
#'   \item{id}{Patient ID}
#'   \item{event_status}{Event status, 0 indicates censoring, 1 indicates CVD event, 2 indicates death from other causes}
#'   \item{event_time}{Event time}
#' }
#'
"data_model_landmark_LOCF"

#' Simulated CVD event data
#'
#' A dataset containing CVD event data for 4100 patients at risk at landmark age 60. This is simulated data formed using an cause-specific survival model.
#'
#' A dataset with 4100 rows and 3 columns:
#' \describe{
#'   \item{id}{Patient ID}
#'   \item{event_status}{Event status, 0 indicates censoring, 1 indicates CVD event, 2 indicates death from other causes}
#'   \item{event_time}{Event time}
#' }
#'
"data_repeat_outcomes"

#' Simulated CVD event data
#'
#' A dataset containing CVD event data for 4100 patients at risk at landmark age 60. This is simulated data formed using an cause-specific survival model.
#'
#' A dataset with 4100 rows and 3 columns:
#' \describe{
#'   \item{id}{Patient ID}
#'   \item{event_status}{Event status, 0 indicates censoring, 1 indicates CVD event, 2 indicates death from other causes}
#'   \item{event_time}{Event time}
#' }
#'
"data_model_longitudinal_LOCF"
