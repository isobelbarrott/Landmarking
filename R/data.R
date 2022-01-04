#' Simulated repeat measurement data
#' @description A simulated dataset formed using an LME model, whose parameters were based on CVD risk assessments
#' recorded at primary care practices in New Zealand. The dataset contains repeat measurements data from 3000 patients.
#'
#' @format A dataset with 9048 rows and 11 columns:
#' \describe{
#'   \item{id}{Patient ID}
#'   \item{smoking}{Smoking status, 0 indicates the patient has never smoked, 1 indicates the patient has quit smoking, and 2 indicates the patient is a current smoker}
#'   \item{diabetes}{Diabetes status, 0 indicates the patient is not diagnosed with diabetes, and 1 indicates the patient is diagnosed with diabetes}
#'   \item{ethnicity}{Ethnicity, one of five ethnicities}
#'   \item{deprivation}{Deprivation score, quintiles}
#'   \item{index}{An index indicating assessment number for a patient}
#'   \item{sbp_stnd}{Standardised systolic blood pressure}
#'   \item{tchdl_stnd}{Standardised total cholesterol to HDL ratio}
#'   \item{end_of_study_age}{Age that the patient left the study, either the age at event (CVD or death) or age at end of study (1st Jan 2010)}
#'   \item{response_time_tchdl_stnd}{Age that total cholesterol to HDL ratio was recorded}
#'   \item{response_time_sbp_stnd}{Age that systolic blood pressure and all other measurements were recorded}
#' }
#'
#'
#'
"data_repeat"

#' Simulated time-to-event data
#' @description This dataset contains time-to-event data simulated from the landmark model.
#' Specifically, the LME model was fitted to dataset \code{data_repeat} and values of \code{sbp_stnd}
#' and \code{tchdl_stnd} were estimated at landmark age 60. These values (along with the other measurements)
#' were used to simulate time-to-event data from a
#' cause specific model with parameters based on CVD events of patients
#' at primary care practices in New Zealand.
#'
#' @format A dataset with 3000 rows and 3 columns:
#'
#' \describe{
#'   \item{id}{Patient ID}
#'   \item{event_status}{Event status, 0 indicates censoring, 1 indicates CVD event, and 2 indicates death from other causes}
#'   \item{event_time}{Event time}
#' }
#'
"data_outcomes"

#' This dataset contains repeat measurements data and time-to-event data
#' from 3000 patients. See vignette "How to use the R package 'Landmarking'" to see how this dataset was formed.
#'
#' @format A dataset with 9048 rows and 14 columns:
#' \describe{
#'   \item{id}{Patient ID}
#'   \item{smoking}{Smoking status, 0 indicates the patient has never smoked, 1 indicates the patient has quit smoking, and 2 indicates the patient is a current smoker}
#'   \item{diabetes}{Diabetes status, 0 indicates the patient is not diagnosed with diabetes, and 1 indicates the patient is diagnosed with diabetes}
#'   \item{ethnicity}{Ethnicity, one of five ethnicities}
#'   \item{deprivation}{Deprivation score, quintiles}
#'   \item{index}{An index indicating assessment number for a patient}
#'   \item{sbp_stnd}{Standardised systolic blood pressure}
#'   \item{tchdl_stnd}{Standardised total cholesterol to HDL ratio}
#'   \item{end_of_study_age}{Age that individual left the study, either the age at event (CVD or death) or age at end of study (1st Jan 2010)}
#'   \item{response_time_tchdl_stnd}{Age that total cholesterol to HDL ratio was recorded}
#'   \item{response_time_sbp_stnd}{Age that systolic blood pressure was recorded, this is the same as the date that the fixed measures were recorded}
#'   \item{start_time}{Age the individual entered follow-up}
#'   \item{event_status}{Event status, 0 indicates censoring, 1 indicates CVD event, and 2 indicates death from other causes}
#'   \item{event_time}{Event time}
#' }
#'
"data_repeat_outcomes"
