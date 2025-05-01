suppressMessages(library(tidyverse))

source("experiment_1_laptop.R")

setup <- 1 ## The only thing this does is determine the sub-directory where the results are saved

setting <- 1
surv_model_type <- "grf"
cens_model_type <- "grf"
num_samples_cal <- 500
prop_shift <- 0         ## proportion of distribution shift in training data (default: 0, other values only used for synthetic data setting 4)
batch <- 1
batch_size <- 1
num_samples_test <- 1000


## Fixed time points and screening probabilities
if(setting==1) {
    time_points <- c(6,12)
} else if(setting==2) {
    time_points <- c(2,3)
} else if(setting==3) {
    time_points <- c(3,10)
} else if(setting==4) {
    time_points <- c(3)
} else {
    time_points <- NULL
}
screening_prob <- c(0.25, 0.5, 0.8, 0.9)

## Make sure directory for saving results exists
out.dir <-sprintf("results/setup_%d", setup)   
if(!dir.exists(file.path(out.dir))) {
    dir.create(file.path(out.dir))
} else {
    cat(sprintf("Directory %s exists.\n", out.dir))
}

## Run all the experiments for different values of training sample size
for(num_samples_train in c(100,1000)) {
    run_all_experiments(setting, surv_model_type, cens_model_type,
                        num_samples_train, num_samples_cal, num_samples_test,
                        time_points, screening_prob,
                        batch, batch_size=batch_size, prop_shift=prop_shift, setup=setup)
}

## Load the results from sub-directory corresponding to given 'setup' number
load_data <- function(setup) {
    idir <- sprintf("results/setup_%d", setup)
    ifile.list <- list.files(out.dir)
    results <- do.call("rbind", lapply(ifile.list, function(ifile) {
        df <- read_delim(sprintf("%s/%s", idir, ifile), delim=",", col_types=cols(), guess_max=2)
    }))
    return(results)
}

results <- load_data(setup)
