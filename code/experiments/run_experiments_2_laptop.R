suppressMessages(library(tidyverse))

source("experiment_2_laptop.R")

## Available data sets:
## | Dataset name    | Obs.  | Vars | Source                 | Citation                            |
## |-----------------|-------|------|------------------------|-------------------------------------|
## | **COLON**       | 1,858 | 11   | `survival` package     | Moertel et al. (1990)               |
## | **GBSG**        | 2,232 | 6    | GitHub (`DeepSurv`)    | Katzman et al. (2018)               |
## | **HEART**       | 172   | 4    | `survival` package     | Crowley and Hu (1977)               |
## | **METABRIC**    | 1,981 | 41   | cBioPortal             | Curtis et al. (2012)                |
## | **PBC**         | 418   | 17   | `survival` package     | Therneau et al. (2000)              |
## | **RETINOPATHY** | 394   | 5    | `survival` package     | Blair et al. (1980)                 |
## | **VALCT**       | 137   | 6    | `survival` package     | Kalbfleisch and Prentice (2002)     |

##############################################
## Dataset loading and splitting parameters ##
##############################################

dataset <- "COLON"

## Relative calibration sample size
cal_prop = 0.2

## Relative test sample size
test_prop = 0.2

## Relative training sample size
train_prop = 1 - (cal_prop + test_prop)

#########################
## Analysis parameters ##
#########################

## Types of model
surv_model_type <- "grf"
cens_model_type <- "grf"

## Random seed for reproducibility
random.state <- 2025

## Number of time points (for adaptive time points, used to construct CSB)
num_time_points <- 100

## Fixed time points and screening probabilities
time_points_screen_q <- c(0.1, 0.5, 0.75, 0.9)   # Quantiles of observed time distribution at which to do screening
screening_probs <- c(0.25, 0.5, 0.75, 0.8, 0.9)   # Probability levels for screening

## Do not use censored observations for calibration
move_censored_obs <- FALSE

## Make sure directory for saving results exists
out.dir <- sprintf("results/data")
if(!dir.exists(file.path(out.dir))) {
    dir.create(file.path(out.dir))
} else {
    cat(sprintf("Directory %s exists.\n", out.dir))
}


##################
# Load the data ##
##################

data <- load_csv(dataset)

########################
## Data preprocessing ##
########################

## Data features
num_features <- ncol(data) - 2

# Rename the columns
col.names <- c("time", "status", paste("X", 1:num_features, sep = ""))
colnames(data) <- col.names

## Compute screening time point
time.points_values <- data$time
screening_times <- round(as.numeric(quantile(time.points_values, time_points_screen_q)))

## Split the dataset
split_result <- split_data(data, train_prop = train_prop, cal_prop = cal_prop, test_prop = test_prop, seed = random.state)

## Access subsets
data.train <- split_result$train
data.cal <- split_result$cal
data.test <- split_result$test

##############
## Analysis ##
##############

results <- run_analysis(data.train, data.cal, data.test, surv_model_type, cens_model_type,
                        screening_times, screening_probs,
                        move_censored_obs=move_censored_obs)
