# Load required libraries
suppressMessages(library(tidyverse))
suppressMessages(library(survival))

library(confsurv)

## Source utility functions for data generation and analysis
source("../utils/utils_misc.R")

## Source utility functions for loading and processing real data for experiments described in paper
source("utils_data.R")

######################
## Input parameters ##
######################

## Flag to determine if input should be parsed from command line
parse_input <- TRUE

if(parse_input) {
    ## Reading command line arguments
    args <- commandArgs(trailingOnly = TRUE)

    ## Checking if the correct number of arguments is provided
    if (length(args) < 5) {
        stop("Insufficient arguments provided. Expected 4 arguments.")
    }

    ## Assigning command line arguments to variables
    dataset <- args[1]
    surv_model_type <- args[2]
    cens_model_type <- args[3]
    train_prop_sub <- as.double(args[4])
    batch <- as.integer(args[5])

} else {
    dataset <- "PBC"
    surv_model_type <- "grf"
    cens_model_type <- "grf"
    train_prop_sub <- 1
    batch <- 1
}

######################
## Fixed parameters ##
######################

## Relative calibration sample size
cal_prop = 0.2

## Relative test sample size
test_prop = 0.2

## Relative training sample size
train_prop = 1 - (cal_prop + test_prop)

## Number of repetitions
batch_size <- 10

## Do not use censored observations for calibration
move_censored_obs <- FALSE

## Number of time points (for adaptive time points)
num_time_points <- 100

## Fixed time points and screening probabilities
time_points <- NULL
time_points_screen_q <- c(0.1, 0.5, 0.75, 0.9)
screening_prob <- c(0.25, 0.5, 0.75, 0.8, 0.9)

save_plots <- FALSE

####################
## Prepare output ##
####################

## Store important parameters including model types
header <- tibble(dataset = dataset,
                 surv_model_type = surv_model_type,
                 cens_model_type = cens_model_type,
                 train_prop_sub = train_prop_sub,
                 batch = batch)

## Generate a unique and interpretable file name based on the input parameters
output_file <- paste0("results/data/", dataset,
                      "_surv_", surv_model_type,
                      "_cens_", cens_model_type,
                      "_train_", train_prop_sub,
                      "_batch", batch, ".txt")

## Print the output file name to verify
cat("Output file name:", output_file, "\n")


##################
# Load the data ##
##################

data <- load_csv(dataset)

## Compute screening time point
time.points_values <- data$time

## Data features
num_features <- ncol(data) - 2

# Rename the columns
col.names <- c("time", "status", paste("X", 1:num_features, sep = ""))
colnames(data) <- col.names

## Use all features
num_feat_censor <- num_features

## Find meaningful time grid for this survival distribution
if(is.null(time_points)) {
    max_event_time <- max(data$time)
    time_points <- get_pretty_quantiles(c(0,data$time), n=num_time_points)
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


#######################################
# Define function to analyze the data #
#######################################

analyze_data <- function(data.train, data.cal, data.test, surv_model, cens_model) {      
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
    
    ## Fit the censoring model on all training data
    idx.train.cens <- 1:nrow(data.train)
    cens_model$fit(data = data.train[idx.train.cens,])
      
    ## Apply CSB method
    csb <- conformal_survival_band(data.test, data.cal, surv_model, cens_model$model, time_points=time_points, boost_scores=FALSE)
    csb.plus <- conformal_survival_band(data.test, data.cal, surv_model, cens_model$model, time_points=time_points, boost_scores=TRUE)
    
    ## Define parameter grid for screening analysis
    param_grid <- expand.grid(
        screening_time = time_points_screen_q,
        screening_prob = screening_prob,
        screening_crit = c("low risk", "high risk"),
        stringsAsFactors = FALSE
    )

    ## Iterate through all combinations
    results_all <- pmap_dfr(param_grid, function(screening_time_q, screening_prob, screening_crit) {
        screening_time <- round(as.numeric(quantile(time.points_values, screening_time_q)))

        ## Selection with black-box model
        time.points <- as.numeric(colnames(model_pred))
        sel_model <- select_patients_point(time.points, model_pred, screening_time, screening_prob, screening_crit)$selected
        
        ## Selections with CSB
        sel_csb <- select_patients_band(time.points, csb$lower, csb$upper,
                                         screening_time, screening_prob, screening_crit)$selected
        sel_csb_plus <- select_patients_band(time.points, csb.plus0$lower, csb.plus$upper,
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
                           "CSB" = sel_csb,
                           "CSB+" = sel_csb_plus
                           )
        if(save_plots) {
            patient.list <- c(1:4)
            title.str <- sprintf("Screening for: %s.\n", interpret_screening_rule(screening_time, screening_prob, screening_crit))
            pp <- plot_survival_panel(
                pred_list = list(
                    Model = model_pred[patient.list, ],
                    KM = km_matrix[patient.list, ]
                ),
                band_list = list(`CSB+` = list(
                                     lower = csb_plus$lower[patient.list, ],
                                     upper = csb_plus$upper[patient.list, ]
                                 )),
                screening_time = screening_time,
                screening_prob = screening_prob,
                screening_crit = screening_crit,
                patient_names = paste0("Patient ", seq_along(patient.list)),
                x_label="Time",
                save_fig = FALSE
            )
            pp <- pp + ggtitle(title.str) + theme(plot.title = element_text(margin = margin(b = -5)))
            fname <- sprintf("figures/data%s_T%.2f_P%.2f_%s.pdf",
                             dataset, screening_time_q, screening_prob,
                             gsub(" ", "_", screening_crit))
            ggsave(fname, plot = pp, width = 12, height = 4.25, device = cairo_pdf)
        }

        ## Evaluate and format results
        evaluated <- map2_dfr(selections, names(selections), function(selected, method_name) {           
            res.raw <- evaluate_selections_without_oracle(data.test, selected, screening_time, screening_prob, screening_crit)
            res.raw %>%
                transmute(
                    Method = method_name,
                    Time_q = screening_time_q,
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

#######################################
## Define function to run experiment ##
#######################################

run_experiment <- function(random.state) {
    set.seed(random.state)

    ## Split the dataset
    split_result <- split_data(data, train_prop = train_prop, cal_prop = cal_prop, test_prop = test_prop, seed = random.state)

    ## Access subsets
    data.train <- split_result$train
    data.cal <- split_result$cal
    data.test <- split_result$test

    ## Subsample the training set
    data.train = data.train[sample(1:nrow(data.train), ceiling(train_prop_sub*nrow(data.train))),]

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
        results_df <- rbind(results_df, result_df)

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
