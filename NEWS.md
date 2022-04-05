# Landmarking 1.0.0

* Added a `NEWS.md` file to track changes to the package.

# Landmarking 1.0.1

* Fix bug that returned prediction of survival, rather than prediction of event, for the coxph survival model.

# Landmarking 1.0.2

* Fix bug for non-competing risks scenario in get_model_assessment()

* Change the slopes in the survival model to be fixed slopes plus random slopes in function fit_LME_longitudinal()

* Add feature to use all the data or only data up to and incluing the landmark age in the validation dataset of function fit_LME_longitudinal()