
<!-- README.md is generated from README.Rmd. Please edit that file -->

This package is designed to help the user perform dynamic prediction
using a landmark model. This model applies to repeat measurements data
and time-to-event data. There is the option to use the last observation
carried forward (LOCF) or linear mixed effects (LME) model to fit the
repeat measurements data. This package also allows the user to account
for competing risks by using either the Fine Gray or Cause-specific
model to fit the time-to-event data. Cross-validation with k-folds can
be applied easily using this package.

For more detailed examples and for an explanation of landmark models,
see <https://isobelbarrott.github.io/Landmarking>.

## Installation

You can install the development version of this package from
[GitHub](https://github.com/) with:

``` r
install.packages("devtools")
devtools::install_github("isobelbarrott/Landmarking")
```

## A simple example

Below is a simple example of how to use this package with the pbc2
dataset from package JM to predict the risk of death for a new patient.

``` r
library(Landmarking)
data(pbc2,package="JM")
levels(pbc2$status)<-c("0","2","1")
pbc2<-return_ids_with_LOCF(data=pbc2,
                           patient_id="id",
                           x_L=c(5),
                           covariates=c("drug","serBilir","serChol"),
                           covariates_time="year")
pbc2$years<-pbc2$years+pbc2$age
data_model_landmark_LOCF<-fit_LOCF_landmark_model(data=pbc2,
                                                  x_L=c(40),
                                                  x_hor=c(45),
                                                  covariates=c("drug","serBilir","serChol"),
                                                  covariates_time="year",
                                                  start_study_time="age",
                                                  end_study_time="years",
                                                  patient_id="id",
                                                  event_time="years",
                                                  event_status="status",
                                                  survival_submodel = "cause_specific",
                                                  b=50)
newdata<-rbind(data.frame(id=c(1,1,1),year=c(30,32,35),drug=c("placebo","placebo","placebo"),serBilir=c(2.4,2.7,2.6),serChol=c(220,234,234)))
predict(object=data_model_landmark_LOCF,x_L=40,x_hor=45,newdata=newdata)
```
