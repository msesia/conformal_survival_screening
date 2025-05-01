## Load required libraries
suppressMessages(library(tidyverse))
library(survival)
library(confsurv)

## Source utility functions for data generation and analysis
source("../utils/utils_synthetic_data.R")
source("../utils/utils_misc.R")

######################
## Input parameters ##
######################

run_all_experiments <- function(setting, surv_model_type, cens_model_type,
                                num_samples_train, num_samples_cal, num_samples_test,
                                time_points, screening_prob,
                                batch, batch_size=10, prop_shift=0, setup=1) {
    
    ######################
    ## Fixed parameters ##
    ######################
    
    ## Do not use censored observations for calibration
    move_censored_obs <- FALSE

    ## Number of independent data points used to estimate interesting time points
    num_samples_ref <- 10000

    ## Number of time points (for adaptive time points)
    num_time_points <- 100

    ## Parameters controlling training of censoring model
    num_samples_train_cens <- num_samples_train
    num_feat_censor <- 10

    ####################
    ## Prepare output ##
    ####################
    
    ## Store important parameters including model types
    header <- tibble(setup = setup,
                     setting = setting,
                     surv_model_type = surv_model_type,
                     cens_model_type = cens_model_type,
                     n_train = num_samples_train,
                     n_train_cens = num_samples_train_cens,
                     n_cal = num_samples_cal,
                     n_test = num_samples_test,
                     num_feat_censor = num_feat_censor,
                     prop_shift = prop_shift,
                     batch = batch)
    
    ## Generate a unique and interpretable file name based on the input parameters
    output_file <- paste0("results/setup_", setup, "/",
                          "setting", setting,
                          "_surv_", surv_model_type,
                          "_cens_", cens_model_type,
                          "_train", num_samples_train,
                          "_trainc", num_samples_train_cens,
                          "_cal", num_samples_cal,
                          "_nfc", num_feat_censor,
                          "_shift", prop_shift,
                          "_batch", batch, ".txt")
    
    ## Print the output file name to verify
    cat("Output file name:", output_file, "\n")
   
    
    ##############################
    ## Define data distribution ##
    ##############################
    
    ## Initialize the data generator
    datagen <- init_data_generator(setting)
    num_features <- datagen$num_features
    covariate_generator <- datagen$covariate_generator
    survival_generator <- datagen$survival_generator
    censoring_generator <- datagen$censoring_generator
    survival_generator_shift <- datagen$survival_generator_shift
    generator <- SurvivalDataGenerator$new(covariate_generator, survival_generator, censoring_generator)
    generator_shift <- SurvivalDataGenerator$new(covariate_generator, survival_generator_shift, censoring_generator)

    ## Generate independent reference data set to find meaningful time grid for this survival distribution
    if(is.null(time_points)) {
        set.seed(1)
        data.ref <- generator$sample(num_samples_ref)
        max_event_time <- max(data.ref$event_time)
        time_points <- get_pretty_quantiles(c(0,data.ref$time), n=num_time_points)
    }
            
    ###################################################
    ## Instantiate the survival and censoring models ##
    ###################################################

    surv_model <- init_surv_model(surv_model_type)
    surv_model_large <- init_surv_model(surv_model_type)

    ## List of covariates to use for censoring model
    use.covariates <- paste("X", 1:min(num_features, num_feat_censor), sep="")

    ## Instantiate censoring model based on the specified type
    cens_base_model <- init_censoring_model(cens_model_type, use_covariates=use.covariates)
    ## Create an instance of the CensoringModel class with the model
    cens_model <- CensoringModel$new(model = cens_base_model)

    
    #########################################
    ## Define function to analyze the data ##
    #########################################

    analyze_data <- function(data.train, data.cal, data.test, surv_model, cens_model, data.test.oracle) {
        ## Calculate oracle survival curves (using true population model)
        oracle_pred <- survival_generator$predict(data.test, time_points)$predictions
        colnames(oracle_pred) <- time_points
        
        ## Fit the Kaplan-Meier survival model
        surv_object <- Surv(time = data.cal$time, event = data.cal$status)
        km_fit <- survival::survfit(surv_object ~ 1)

        ## Fit the survival model on the training data
        surv_model$fit(Surv(time, status) ~ ., data = data.train)

        ## Fit the survival model on all training and calibration data
        data.supervised <- rbind(data.train, data.cal)
        surv_model_large$fit(Surv(time, status) ~ ., data = data.supervised)
        model_pred <- matrix(surv_model_large$predict(data.test, time_points)$predictions, nrow(data.test), ncol = length(time_points))
        colnames(model_pred) <- time_points    
        
        ## Fit the censoring model on a subset of the training data
        if(num_samples_train_cens < num_samples_train) {
            idx.train.cens <- sort(sample(1:nrow(data.train), num_samples_train_cens))
        } else {
            idx.train.cens <- 1:nrow(data.train)
        }
        cens_model$fit(data = data.train[idx.train.cens,])
        
        ## Apply CSB method
        csb <- conformal_survival_band(data.test, data.cal, surv_model, cens_model$model, time_points=time_points)

        ## Define parameter grid for screening analysis
        param_grid <- expand.grid(
            screening_time = time_points,
            screening_prob = screening_prob,
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

            ## Selections with oracle
            time.points.oracle <- as.numeric(colnames(oracle_pred))
            sel_oracle <- select_patients_point(time.points.oracle, oracle_pred, screening_time, screening_prob, screening_crit)$selected

            ## Selections with KM
            km_vec <- summary(km_fit, times = time.points.oracle, extend = TRUE)$surv
            n_patients <- nrow(data.test)
            km_matrix <- matrix(rep(km_vec, each = n_patients), nrow = n_patients, byrow = FALSE)
            colnames(km_matrix) <- time.points.oracle
            sel_km <- select_patients_point(time.points.oracle, km_matrix, screening_time, screening_prob, screening_crit)$selected
            
            ## Combine all selections
            selections <- list("model" = sel_model,
                               "KM" = sel_km,
                               "CSB" = sel_csb,
                               "oracle" = sel_oracle)
            ## Evaluate and format results
            evaluated <- map2_dfr(selections, names(selections), function(selected, method_name) {
                res.raw <- evaluate_selections(data.test.oracle, selected, screening_time, screening_prob, screening_crit)
                res.raw %>%
                    transmute(
                        Method = method_name,
                        Time = Screening.time,
                        Criterion = screening_crit,
                        Probability = screening_prob,
                        Screened = Screened,
                        Survival = Proportion.survived
                    )
            })

            ## Evaluate classification performance
            evaluated_class <- evaluate_classification(selections, oracle=selections$oracle, n_test=num_samples_test)
            evaluated <- full_join(evaluated, evaluated_class, by="Method")
            
            return(evaluated)
        })

        return(results_all)
    }
   
    
    #######################################
    ## Define function to run experiment ##
    #######################################
    
    run_experiment <- function(random.state) {
        set.seed(random.state)

        ## Generate training, calibration, and test data
        ## (including true event and censoring times)
        data.train.oracle <- sample_mixture(generator, generator_shift, num_samples_train, prop_shift = prop_shift)
        data.cal.oracle <- generator$sample(num_samples_cal)
        data.test.oracle <- generator$sample(num_samples_test)

        ## This calibration method only uses observations with events.
        ## Therefore, we can move all censored observations from the calibration to the training set
        if(move_censored_obs) {
            idx.cal.censored <- which(data.cal.oracle$status==0)
            idx.cal.event <- which(data.cal.oracle$status==1)
            if(length(idx.cal.censored)>0) {
                data.train.oracle <- rbind(data.train.oracle, data.cal.oracle[idx.cal.censored,])
            }
            if(length(idx.cal.censored)>0) {
                data.cal.oracle <- data.cal.oracle[idx.cal.event,]
            } else {
                warning("Warning! No events in calibration data set")
                data.cal.oracle <- head(data.cal.oracle,0)
            }
        }
        
        ## Remove true event and censoring times from the data (right-censoring)
        data.train <- data.train.oracle |> select(-event_time, -censoring_time)
        data.cal <- data.cal.oracle |> select(-event_time, -censoring_time)
        data.test <- data.test.oracle |> select(-event_time, -censoring_time)
        
        ## Run analysis
        results <- analyze_data(data.train, data.cal, data.test, surv_model, cens_model, data.test.oracle)

        return(results)
    }

    ## Function to run multiple experiments and gather results
    ## Args:
    ##   batch_size: Number of repetitions for each experimental setting
    ## Returns:
    ##   A tibble containing the combined results of all experiments
    run_multiple_experiments <- function(batch_size) {
        results_df <- data.frame()  ## Initialize an empty data frame to store cumulative results

                                        ## Print a progress bar header
        cat("Running experiments\n")
        pb <- txtProgressBar(min = 0, max = batch_size, style = 3)  ## Initialize progress bar

                                        ## Loop over each repetition
        for (i in 1:batch_size) {
            random.state <- batch*1000 + i
            res <- run_experiment(random.state)  ## Run experiment and get the result

            ## Combine the results with experiment metadata
            result_df <- tibble(Seed = random.state) |> cbind(header) |> cbind(res)

                                        ## Add the result to the cumulative data frame
            results_df <- rbind(results_df, result_df) %>% as_tibble()

                                        ## Write the cumulative results to the CSV file
            write.csv(results_df, output_file, row.names = FALSE)

            setTxtProgressBar(pb, i)  ## Update progress bar
        }

        close(pb)  ## Close the progress bar

        return(results_df)  ## Return the cumulative results data frame
    }
    
    #####################
    ## Run experiments ##
    #####################
    
    ## Run the experiments with specified parameters
    results <- run_multiple_experiments(batch_size)
    
    return(results)
}
