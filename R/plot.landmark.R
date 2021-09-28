#' Create a calibration plot
#'
#' Creates a calibration plot for the landmark model fitted by `fit_LME_landmark_model` or `fit_LOCF_landmark_model`.
#' This function plots the observed frequencies of the event of interest against the predicted probabilities of the event of interest.
#'
#' @param x Object inheriting the class `landmark`, this should be the output from either `fit_LME_landmark_model` or `fit_LOCF_landmark_model`. It should contain a list
#' of landmark models corresponding to different landmark times `x_L`.
#' @param x_L Numeric specifying the landmark time. This indicates which landmark model in `x` to use.
#' @param n Numeric specifying the number of bins to use.
#' @param \dots Arguments passed on to the `aes()` function
#' @return Calibration plot showing the value of predicted probabilities against observed frequencies, with a `y=x` line.
#' @details This function bins the predicted probabilities of the event of interest into `n` bins. The event of interest is the event with
#' `event_status=1` when fitting the landmark model. For each of the `n` sets of individuals, the Aalen-Johansen estimator is fit to that set
#' and used to calculate the risk of an event at the horizon time. The predictions (from the landmark model) and the observed frequencies
#' (from the Aalen-Johansen estimator) are plotted against each other. For a perfect prediction model, the points will be plotted along the y=x line.
#' @examples
#' library(Landmarking)
#' data(data_repeat_outcomes)
#' data_repeat_outcomes <-
#'   return_ids_with_LOCF(
#'     data_long = data_repeat_outcomes,
#'     patient_id = "id",
#'     covariates =
#'       c("ethnicity", "smoking", "diabetes", "sbp_stnd", "tchdl_stnd"),
#'     covariates_time =
#'       c(rep("response_time_sbp_stnd", 4), "response_time_tchdl_stnd"),
#'     x_L = c(60,61)
#'   )
#' data_model_landmark_LOCF <-
#'   fit_LOCF_landmark_model(
#'     data_long = data_repeat_outcomes,
#'     x_L = c(60, 61),
#'     x_hor = c(65, 66),
#'     covariates =
#'       c("ethnicity", "smoking", "diabetes", "sbp_stnd", "tchdl_stnd"),
#'     covariates_time =
#'       c(rep("response_time_sbp_stnd", 4), "response_time_tchdl_stnd"),
#'     k = 10,
#'     start_study_time = "start_time",
#'     end_study_time = "event_time",
#'     patient_id = "id",
#'     event_time = "event_time",
#'     event_status = "event_status",
#'     survival_submodel = "cause_specific"
#'   )
#'  plot(x=data_model_landmark_LOCF,x_L=60,n=5)
#'  plot(x=data_model_landmark_LOCF,x_L=61,n=5)
#' @export

plot.landmark<-function(x,x_L,n,...){
  if(class(x)!="landmark"){
    stop("x must have class 'landmark'")}
  if(!is.numeric(x_L)){
    stop("x_L should be numeric")
  }
  if(!is.numeric(n)){stop("n must have class numeric")}
  if(!(as.character(x_L) %in% names(x))){stop("x_L must be the name of an element in 'x'")}
  x<-x[[as.character(x_L)]]
  event_status<-x$call$event_status

  strata<-survival::strata
  actual<-c()
  predicted<-c()
  data_survival<-x$data
  event_time<-x$call$event_time
  data_survival[["event_time"]]<-data_survival[[event_time]]
  events<-unique(x[["data"]][[event_status]])

  data_survival$quantile<-.bincode(data_survival$event_prediction, breaks=stats::quantile(data_survival$event_prediction,seq(0,1,by=1/n)),include.lowest=TRUE)
  if(length(table(data_survival$quantile))!=n){stop("Not enough data for number of quantiles, select lower value of n")}
  trans<-mstate::trans.comprisk(max(events),1:max(events))
  data_survival<-do.call("rbind",lapply(1:max(events),function(i){
    data_survival_i<-data.frame(data_survival,trans=i,status=as.numeric(+(data_survival[[event_status]]==i)))
    if(all(data_survival[["status"]]==i)==0){stop(paste0("Not enough events (for event ", i,")"))}
    data_survival_i
  }))
  for (i in 1:n){
    data_i<-data_survival[data_survival$quantile==i,]
    coxph_i<-survival::coxph(survival::Surv(event_time,status)~strata(trans),data=data_i,method="breslow")
    newdata<-data.frame(trans=1:max(events))
    msfit_i<-mstate::msfit(object=coxph_i,newdata=newdata,trans=trans)
    probtrans_i<-mstate::probtrans(msfit_i,predt=0,method="aalen")[[1]]

    actual<-c(actual, probtrans_i[dim(probtrans_i)[1],"pstate2"])
    predicted<-c(predicted,mean(data_i[,"event_prediction"]))
  }
  calibration_plot_df<-data.frame(actual,predicted)
  ggplot2::ggplot(calibration_plot_df,ggplot2::aes(x=predicted,y=actual,...))+
    ggplot2::geom_point()+ggplot2::geom_line()+ggplot2::labs(x="Predicted probability",y="Observed frequency")+
    ggplot2::geom_abline(intercept=0, slope=1,linetype="dashed")+
    ggplot2::scale_x_continuous(limits=c(min(calibration_plot_df),max(calibration_plot_df)))+
    ggplot2::scale_y_continuous(limits=c(min(calibration_plot_df),max(calibration_plot_df)))
}

