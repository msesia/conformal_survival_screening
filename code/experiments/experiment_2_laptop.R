## Load required libraries
suppressMessages(library(tidyverse))
library(survival)
library(confsurv)

## Source utility functions for data generation and analysis
source("../utils/utils_misc.R")

## Source utility functions for loading and processing real data for experiments described in paper
source("utils_data.R")

######################
## Input parameters ##
######################

run_analysis <- function(data.train, data.cal, data.test, surv_model_type, cens_model_type,
                         screening_times, screening_probs, time_points=NULL, move_censored_obs=FALSE) {
  
    ###################################################
    ## Instantiate the survival and censoring models ##
    ###################################################

    surv_model <- init_surv_model(surv_model_type)
    surv_model_large <- init_surv_model(surv_model_type)

    ## List of covariates to use for censoring model    
    num_feat_censor <- num_features ## Use all features
    use.covariates <- paste("X", 1:min(num_features, num_feat_censor), sep="")

    ## Instantiate censoring model based on the specified type
    cens_base_model <- init_censoring_model(cens_model_type, use_covariates=use.covariates)
    ## Create an instance of the CensoringModel class with the model
    cens_model <- CensoringModel$new(model = cens_base_model)

    #######################################
    ## Define function to analyze the data #
    #######################################

    analyze_data <- function(data.train, data.cal, data.test, surv_model, cens_model) {
        ## Fit the Kaplan-Meier survival model
        surv_object <- Surv(time = data.cal$time, event = data.cal$status)
        km_fit <- survival::survfit(surv_object ~ 1)

        ## Fit the survival model on the training data
        surv_model$fit(Surv(time, status) ~ ., data = data.train)

        ## Fit the censoring model on all training data
        idx.train.cens <- 1:nrow(data.train)
        cens_model$fit(data = data.train[idx.train.cens,])
        
        ## Apply CSB method
        csb <- conformal_survival_band(data.test, data.cal, surv_model, cens_model$model, time_points=time_points)

        ## Fit the survival model on all training and calibration data
        data.supervised <- rbind(data.train, data.cal)
        surv_model_large$fit(Surv(time, status) ~ ., data = data.supervised)

        time_points <- csb$time_points
        model_pred <- matrix(surv_model_large$predict(data.test, time_points)$predictions, nrow(data.test), ncol = length(time_points))
        colnames(model_pred) <- time_points


        ## Define parameter grid for screening analysis
        param_grid <- expand.grid(
            screening_time = screening_times,
            screening_prob = screening_probs,
            screening_crit = c("low risk", "high risk"),
            stringsAsFactors = FALSE
        )
       
        ## Iterate through all combinations
        results_all <- pmap_dfr(param_grid, function(screening_time, screening_prob, screening_crit) {

            ## Selection with black-box model
            time.points <- as.numeric(colnames(model_pred))
            sel_model <- select_patients_point(time.points, model_pred, screening_time, screening_prob, screening_crit)$selected

            ## Selections with CSB
            sel_csb <- select_patients_band(time.points, csb$lower, csb$upper,
                                            screening_time, screening_prob, screening_crit)$selected

            ## Selections with KM
            km_vec <- summary(km_fit, times = time.points, extend = TRUE)$surv
            n_patients <- nrow(data.test)
            km_matrix <- matrix(rep(km_vec, each = n_patients), nrow = n_patients, byrow = FALSE)
            colnames(km_matrix) <- time.points
            sel_km <- select_patients_point(time.points, km_matrix, screening_time, screening_prob, screening_crit)$selected

            ## Combine all selections
            selections <- list("model" = sel_model,
                               "KM" = sel_km,
                               "CSB" = sel_csb
                               )

            ## Evaluate and format results
            evaluated <- map2_dfr(selections, names(selections), function(selected, method_name) {
                res.raw <- evaluate_selections_without_oracle(data.test, selected, screening_time, screening_prob, screening_crit)
                res.raw %>%
                    transmute(
                        Method = method_name,
                        Time = Screening.time,
                        Criterion = screening_crit,
                        Probability = screening_prob,
                        Screened = Screened,
                        Survival_lower = Proportion.survived.lower,
                        Survival_upper = Proportion.survived.upper
                    )
            })
            return(evaluated)
        })

        return(results_all)
    }


    ## This calibration method only uses observations with events.
    ## Therefore, we can move all censored observations from the calibration to the training set
    if(move_censored_obs) {
        idx.cal.censored <- which(data.cal$status==0)
        idx.cal.event <- which(data.cal$status==1)
        if(length(idx.cal.censored)>0) {
            data.train <- rbind(data.train, data.cal[idx.cal.censored,])
        }
        if(length(idx.cal.censored)>0) {
            data.cal <- data.cal[idx.cal.event,]
        } else {
            warning("Warning! No events in calibration data set")
            data.cal <- head(data.cal,0)
        }
    }

    ## Run analysis
    results <- analyze_data(data.train, data.cal, data.test, surv_model, cens_model) %>%
        mutate(n_train = nrow(data.train), n_cal = nrow(data.cal), n_test = nrow(data.test))
}
