
<!-- README.md is generated from README.Rmd. Please edit that file -->

## Installation

You can install the development version of this package from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("isobelbarrott/Landmarking")
```

Below is an example of how to use this package with the pbc2 dataset.
For more detailed examples and for an explanation of Landmark models,
see <https://isobelbarrott.github.io/Landmarking>.

``` r
library(Landmarking)
data(pbc2)
new_status<-data.frame(status=c(0,1,2),status_cr=c(0,2,1))
pbc2<-dplyr::left_join(pbc2,new_status,by="status")
pbc2<-return_ids_with_LOCF(data=pbc2,
                           patient_id="id",
                           x_L=c(40),
                           covariates=c("drug","serBilir","serChol"),
                           covariates_time="year")
data_model_landmark_LOCF<-fit_LOCF_landmark_model(data=pbc2,
                                                  x_L=c(40),
                                                  x_hor=c(45),
                                                  covariates=c("drug","serBilir","serChol"),
                                                  covariates_time="year",
                                                  start_study_time="age",
                                                  end_study_time="years",
                                                  patient_id="id",
                                                  event_time="years",
                                                  event_status="status_cr",
                                                  survival_submodel = "cause_specific",
                                                  b=50)
#> [1] "Fitting longitudinal submodel, landmark age 40"
#> [1] "Complete, landmark age 40"
#> [1] "Fitting survival submodel, landmark age 40"
#> [1] "Complete, landmark age 40"

newdata<-rbind(data.frame(id=c(1,1,1),year=c(30,32,35),drug=c("placebo","placebo","placebo"),serBilir=c(2.4,2.7,2.6),serChol=c(220,234,234)))
predict(object=data_model_landmark_LOCF,x_L=40,x_hor=45,newdata=newdata)
#>    id    drug serBilir serChol event_prediction
#> 1:  1 placebo      2.6     234        0.1834694
```
