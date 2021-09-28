#' Simulated repeat measurement data
#' @description A simulated dataset formed using an LME model, whose parameters were based on CVD risk assessments
#' recorded at primary care practices in New Zealand. The dataset contains repeat measurements data from 10000 patients.
#'
#' @format A dataset with 9048 rows and 11 columns:
#' \describe{
#'   \item{id}{Patient ID}
#'   \item{smoking}{Smoking status, 0 indicates the patient has 'never smoked', 1 indicates the patient has quit smoking, and 2 indicated the individual is a current smoker.}
#'   \item{diabetes}{Diabetes status, 0 indicates the patient does not have diabetes, 1 indicates the patient does have diabetes}
#'   \item{ethnicity}{Ethnicity, one of five ethnicities}
#'   \item{deprivation}{Deprivation score out of five}
#'   \item{index}{An index indicating assessment number for a patient}
#'   \item{sbp_stnd}{Standardised systolic blood pressure}
#'   \item{tchdl_stnd}{Standardised total cholesterol to HDL ratio}
#'   \item{end_of_study_age}{Age that individual left the study, either the age at event (CVD or death) or age at end of study (1st Jan 2010)}
#'   \item{response_time_tchdl_stnd}{Age that total cholesterol to HDL ratio was recorded}
#'   \item{response_time_sbp_stnd}{Age that systolic blood pressure was recorded, this is the same as the date that the fixed measures were recorded}
#' }
#'
#'
#'
"data_repeat"

#' Simulated time-to-event data
#' @description A simulated dataset formed using a cause specific model, whose parameters were based on CVD events of patients
#' at primary care practices in New Zealand. The dataset contains time-to-event data from 3000 patients.
#'
#' @format A dataset with 3000 rows and 3 columns:
#'
#' \describe{
#'   \item{id}{Patient ID}
#'   \item{event_status}{Event status, 0 indicates censoring, 1 indicates CVD event, 2 indicates death from other causes}
#'   \item{event_time}{Event time}
#' }
#'
"data_outcomes"

#' Simulated repeat measurement and time-to-event data
#' @description A simulated dataset formed using an LME model, whose parameters were based on CVD risk assessments
#' recorded at primary care practices in New Zealand. The dataset contains repeat measurements data and time-to-event data
#' from 10000 patients. See vignette "How to use the R package 'Landmarking'" to see how this dataset was formed.
#'
#' @format A dataset with 9048 rows and 14 columns:
#' \describe{
#'   \item{id}{Patient ID}
#'   \item{smoking}{Smoking status, 0 indicates the patient has 'never smoked', 1 indicates the patient has quit smoking, and 2 indicated the individual is a current smoker.}
#'   \item{diabetes}{Diabetes status, 0 indicates the patient does not have diabetes, 1 indicates the patient does have diabetes}
#'   \item{ethnicity}{Ethnicity, one of five ethnicities}
#'   \item{deprivation}{Deprivation score out of five}
#'   \item{index}{An index indicating assessment number for a patient}
#'   \item{sbp_stnd}{Standardised systolic blood pressure}
#'   \item{tchdl_stnd}{Standardised total cholesterol to HDL ratio}
#'   \item{end_of_study_age}{Age that individual left the study, either the age at event (CVD or death) or age at end of study (1st Jan 2010)}
#'   \item{response_time_tchdl_stnd}{Age that total cholesterol to HDL ratio was recorded}
#'   \item{response_time_sbp_stnd}{Age that systolic blood pressure was recorded, this is the same as the date that the fixed measures were recorded}
#' }
#'
"data_repeat_outcomes"
