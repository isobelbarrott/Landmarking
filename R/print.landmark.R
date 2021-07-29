#' @export
print.landmark<-function(x,...){print(lapply(lapply(x,"[[",i=1),utils::head),...)}
