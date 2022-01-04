#' @param event_status Character string specifying the column name in \code{data_long}
#' which contains the event status (where 0=censoring, 1=event of interest, if there are competing events these are labelled
#' 2 or above). Events at time `x_hor` should be labelled censored.
