
<!-- README.md is generated from README.Rmd. Please edit that file -->

This package is designed to help the user perform dynamic prediction
using a landmark model. Dynamic prediction means that risk predictions
of an event are updated as more longitudinal data is collected. There is
the option to use the last observation carried forward (LOCF) or linear
mixed effects (LME) model in the first stage of the two-stage landmark
model. This package allows the user to account for competing risks by
using either the Fine Gray model or cause-specific model (or the Cox
model) in the second stage of the two-stage landmark model. k-fold
cross-validation can be performed using this package.

For more detailed examples and for an explanation of landmark models,
see <https://isobelbarrott.github.io/Landmarking/>.

## Installation

You can install the development version of this package from
[GitHub](https://github.com/) with:

``` r
install.packages("devtools")
devtools::install_github("isobelbarrott/Landmarking")
```

## A simple example

Below is a simple toy example of how to use this package with the pbc2
dataset from package JM to predict the risk of death for a new patient.

This uses the LOCF method for the first stage of the landmark model and
cause-specific model for the second stage.

``` r
#Load the library and dataset
library(Landmarking)
data(pbc2,package="JM")
#Change levels to make the death the event of interest (event_status=1), transplant the competing risks (event_status=2), and leave censoring (event_status=0)
levels(pbc2$status)<-c("0","2","1")
#Calculate the age of the patient at each assessment (as opposed to time since first assessment)
pbc2$years<-pbc2$years+pbc2$age
#Fit the landmark model
data_model_landmark_LOCF<-fit_LOCF_landmark(data=pbc2,
                                                  x_L=40,
                                                  x_hor=45,
                                                  covariates=c("serBilir","serChol"),
                                                  covariates_time="year",
                                                  individual_id="id",
                                                  event_time="years",
                                                  event_status="status",
                                                  survival_submodel = "cause_specific",
                                                  b=50)
#Define new dataset
newdata<-rbind(data.frame(id=c(313,313,313),year=c(30,32,35),serBilir=c(2.4,2.7,2.6),serChol=c(220,233,234)))
#Return event prediction and LOCF values
predict(object=data_model_landmark_LOCF,x_L=40,x_hor=45,newdata=newdata)
```

<!-- badges: start -->

[![R-CMD-check](https://github.com/isobelbarrott/Landmarking/workflows/R-CMD-check/badge.svg)](https://github.com/isobelbarrott/Landmarking/actions)
<!-- badges: end -->
