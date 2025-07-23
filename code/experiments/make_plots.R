options(width = 300)

library(tidyverse)
library(kableExtra)

interpret_screening_rule <- function(time, prob, crit, latex=FALSE) {
    stopifnot(crit %in% c("low risk", "high risk"))
    stopifnot(time>=0)
    stopifnot(prob>=0 & prob<=1)
    if(crit=="high risk") {
        if(latex) {
            interpretation <- sprintf("high-risk patients with P($T> %.2f$)$ < %.2f$", time, prob)
        } else {
            interpretation <- sprintf("high-risk patients with P(T> %.2f) < %.2f", time, prob)
        }
    } else {
        if(latex) {
            interpretation <- sprintf("low-risk patients with P($T>%.2f$)$ > %.2f$", time, prob)
        } else {
            interpretation <- sprintf("low-risk patients with P(T> %.2f) > %.2f", time, prob)
        }    }
    return(interpretation)
}

save_fig <- TRUE

variable.values <- c("Screened", "Survival", "Precision", "Recall")
variable.labels <- c("Screened Proportion", "Survival Rate", "Precision", "Recall")

method.values <- c("model", "KM", "CSB+", "oracle")
method.labels <- c("Model", "KM", "CSB", "Oracle")

#method.values <- c("model", "KM", "CSB+", "oracle")
#method.labels <- c("Model", "KM", "CSB+", "Oracle")


load_results <- function(setup) {
    idir <- sprintf("results_hpc/setup_%s", setup)
    ifile.list <- list.files(idir)
    results <- do.call("rbind", lapply(ifile.list, function(ifile) {
        df <- read_delim(sprintf("%s/%s", idir, ifile), delim=",", col_types=cols(), guess_max=2)
    }))
    results <- results %>%
        mutate(
            across(c(Precision, Recall), ~replace(., is.na(.), 0))
        )
    return(results)
}

source("../utils/utils_plotting.R")
custom_colors <- get_custom_colors_model()
custom_shapes <- get_custom_shapes_model()

#############
## SETUP 1 ##
#############

make_figures_setup_1 <- function(results) {
    ## Define specific parameter combinations
    configs <- list(
        list(setting = 1, screening_time = 6, screening_prob = 0.8, screening_crit = "low risk"),
        list(setting = 1, screening_time = 12, screening_prob = 0.8, screening_crit = "high risk"),
        list(setting = 2, screening_time = 2, screening_prob = 0.8, screening_crit = "low risk"),
        list(setting = 2, screening_time = 3, screening_prob = 0.25, screening_crit = "high risk"),
        list(setting = 3, screening_time = 3, screening_prob = 0.9, screening_crit = "low risk"),
        list(setting = 3, screening_time = 10, screening_prob = 0.5, screening_crit = "high risk")
        ## Add more custom combinations as needed
    )

    df <- results %>%
        filter(setup == 1, Method %in% method.values, n_train_cens==10000) %>%
        mutate(Method = factor(Method, method.values, method.labels)) %>%
        group_by(Seed, setting, n_train, n_cal, Criterion, Time, Probability) %>%
        mutate(Oracle_Screened = Screened[Method == "Oracle"]) %>%        
        group_by(setting, n_train, n_cal, Method, Criterion, Time, Probability) %>%
        mutate(Screened = Screened / n_test) %>%
        mutate(Precision = ifelse(Screened==0, NA, Precision),
               Recall = ifelse(Oracle_Screened==0, NA, Recall)) %>%
        summarise(
            across(c(Screened, Survival, Precision, Recall),
                   list(mean = ~mean(., na.rm=TRUE), se = ~sd(., na.rm=TRUE) / sqrt(n())),
                   .names = "{.col}_{.fn}"),
            N = n(),
            .groups = 'drop'
        ) %>%
        mutate(
            across(where(is.numeric), ~ifelse(is.nan(.), NA, .))
        )

    ## Function to filter df based on one config
    filter_one_config <- function(cfg) {
        df %>%
            filter(
                setting == cfg$setting,
                Time == cfg$screening_time,
                Probability == cfg$screening_prob,
                Criterion == cfg$screening_crit
            )
    }
    
    ## Apply the filtering to each config
    filtered_dfs <- map(configs, filter_one_config)

    plot_summary <- function(idx) {

        df.plot <- filtered_dfs[[idx]] %>% pivot_longer(
                                                    cols = c(Screened_mean, Survival_mean, Precision_mean, Recall_mean,
                                                             Screened_se, Survival_se, Precision_se, Recall_se),
                                                    names_to = c("Variable", ".value"),
                                                    names_pattern = "(.*)_(mean|se)"
                                                ) %>%
            mutate(Variable = factor(Variable, variable.values, variable.labels))

        df.shade <- df.plot %>%
            filter(Variable == "Survival Rate") %>%
            distinct(Probability, Variable, Criterion) %>%
            mutate(
                ymin = ifelse(Criterion == "high risk", 0, Probability),
                ymax = ifelse(Criterion == "high risk", Probability, 1)
            ) 

        if(configs[[idx]]$screening_crit == "low risk") {
            title.str <- sprintf("(a) Screening for: %s.\n",
                                 interpret_screening_rule(configs[[idx]]$screening_time,
                                                          configs[[idx]]$screening_prob,
                                                          configs[[idx]]$screening_crit))
        } else {
            title.str <- sprintf("(b) Screening for: %s.\n",
                                 interpret_screening_rule(configs[[idx]]$screening_time,
                                                          configs[[idx]]$screening_prob,
                                                          configs[[idx]]$screening_crit))
        }

        pp <- df.plot %>%
            ggplot(aes(x=n_train, y=mean, color=Method, shape=Method)) +
            geom_rect(
                data = df.shade,
                aes(xmin = min(df.plot$n_train), xmax = max(df.plot$n_train), ymin = ymin, ymax = ymax),
                inherit.aes = FALSE,
                fill = "gray80", alpha = 0.5
            ) +
            geom_point() +
            geom_line() +
            geom_errorbar(aes(ymin = mean - 2*se, ymax = mean + 2*se), width = 0.1) +
            facet_wrap(.~Variable, nrow=1, scales="free") +
            scale_x_log10() +
            ylim(0,1) +
            geom_hline(
                data = df.plot %>% filter(Variable == "Survival Rate") %>% distinct(Probability, Variable),
                aes(yintercept = Probability),
                linetype = "dashed", color = "gray50"
            ) +
            scale_color_manual(name = "Method", values = custom_colors,
                               guide = guide_legend(order = 1)) +
            scale_shape_manual(name = "Method", values = custom_shapes,
                               guide = guide_legend(order = 1)) +   
            labs(
                x = "Training Sample Size",
                y = "Mean ± 2 SE",
                subtitle = title.str
            ) +
            theme_bw() +
            theme(legend.position = "right",
                  plot.subtitle = element_text(margin = margin(b = -5)))

        if (save_fig) {
            fname <- sprintf("figures/setup1_setting%d_T%d_P%.2f_%s.pdf",
                             configs[[idx]]$setting,
                             configs[[idx]]$screening_time,
                             configs[[idx]]$screening_prob,
                             gsub(" ", "_", configs[[idx]]$screening_crit))
            ggsave(fname, pp, width = 9, height = 2.25)
        }
    }

    for(idx in 1:length(configs)) {
        cat(sprintf("Making plot %d of %d...\n", idx, length(configs)))
        plot_summary(idx)
    }
    
}

results <- load_results(1)
make_figures_setup_1(results)


#############
## SETUP 2 ##
#############

### Tables for Setup 2
make_tables_setup_2 <- function(results) {
    ## Define specific parameter combinations
    configs <- list(
        list(setting = 4, screening_time = 3, screening_prob = 0.8, screening_crit = "low risk"),
        list(setting = 4, screening_time = 3, screening_prob = 0.5, screening_crit = "high risk")
    )

    df <- results %>%
        filter(setup == 2, Method %in% method.values) %>%
        mutate(Method = factor(Method, method.values, method.labels)) %>%
        group_by(Seed, setting, n_train, n_cal, prop_shift, Criterion, Time, Probability) %>%
        mutate(Oracle_Screened = Screened[Method == "Oracle"]) %>%        
        group_by(setting, n_train, n_cal, prop_shift, Method, Criterion, Time, Probability) %>%
        mutate(Screened = Screened / n_test) %>%
        mutate(Precision = ifelse(Screened==0, NA, Precision),
               Recall = ifelse(Oracle_Screened==0, NA, Recall)) %>%
        summarise(
            across(c(Screened, Survival, Precision, Recall),
                   list(mean = ~mean(., na.rm=TRUE), se = ~sd(., na.rm=TRUE) / sqrt(n())),
                   .names = "{.col}_{.fn}"),
            N = n(),
            .groups = 'drop'
        ) %>%
        mutate(
            across(where(is.numeric), ~ifelse(is.nan(.), NA, .))
        )

    ## Function to filter df based on one config
    filter_one_config <- function(cfg) {
        df %>%
            filter(
                setting == cfg$setting,
                Time == cfg$screening_time,
                Probability == cfg$screening_prob,
                Criterion == cfg$screening_crit
            )
    }

    ## Apply the filtering to each config
    filtered_dfs <- map(configs, filter_one_config)

    prop.shift.values <- c(1,0)
    model.quality.labels <- c("Low-Quality Model", "High-Quality Model")

    make_table <- function(idx) {
        df.tab <- filtered_dfs[[idx]] %>%
            mutate(
                color = case_when(
                    Criterion == "low risk" & Survival_mean < Probability ~ "red",
                    Criterion == "high risk" & Survival_mean > Probability ~ "red",
                    TRUE ~ "black"
                ),                
                Screened = ifelse(
                    is.na(Screened_mean) | is.na(Screened_se),
                    "NA",
                    sprintf("%.3f ± %.3f", Screened_mean, 2 * Screened_se)
                ),
                Survival = ifelse(Screened_mean==0, "NA", 
                           ifelse(is.na(Survival_mean),
                                  "NA",
                                  sprintf("\\textcolor{%s}{%.3f ± %.3f}", color, Survival_mean, 2 * Survival_se)
                                  )),
                Precision = ifelse(
                    is.na(Precision_mean),
                    "NA",
                    sprintf("%.3f ± %.3f", Precision_mean, 2 * Precision_se)
                ),
                Recall = ifelse(
                    is.na(Recall_mean),
                    "NA",
                    sprintf("%.3f ± %.3f", Recall_mean, 2 * Recall_se)
                )
            ) %>%
            mutate(`Model quality` = factor(prop_shift, prop.shift.values, model.quality.labels)) %>%
            select(`Model quality`, Method, Screened, Survival, Precision, Recall) %>%
            arrange(`Model quality`, Method)

        target.str <- sprintf("%s.",
                             interpret_screening_rule(configs[[idx]]$screening_time,
                                                      configs[[idx]]$screening_prob,
                                                      configs[[idx]]$screening_crit,
                                                      latex=TRUE))
        
        print(df.tab)
        caption.str <- sprintf("Performance metrics for %s", target.str)

        latex_table <- df.tab %>%
            select(-`Model quality`) %>%
            kable(format = "latex", booktabs = TRUE, longtable = FALSE,
                  escape = FALSE,
                  align = "lccccc",
                  caption = NULL) %>%  # <- no caption inside
            collapse_rows(columns = 1, valign = "top", latex_hline = "none") %>%    # no lines between collapsed rows
            group_rows(index = table(df.tab$`Model quality`), escape = FALSE)
        
        fname <- sprintf("tables/setup2_%s.txt",
                         gsub(" ", "_", configs[[idx]]$screening_crit))
        writeLines(latex_table, fname)
        cat(sprintf("Table saved on %s\n.", fname))
        
    }
    for(idx in 1:length(configs)) {
        cat(sprintf("Making table %d of %d...\n", idx, length(configs)))
        make_table(idx)
    }   
}
results <- load_results(2)
make_tables_setup_2(results)


#############
## SETUP 3 ##
#############

make_figures_setup_3 <- function(results) {
    ## Define specific parameter combinations
    configs <- list(
        list(setting = 1, screening_time = 6, screening_prob = 0.8, screening_crit = "low risk"),
        list(setting = 1, screening_time = 12, screening_prob = 0.8, screening_crit = "high risk"),
        list(setting = 2, screening_time = 2, screening_prob = 0.8, screening_crit = "low risk"),
        list(setting = 2, screening_time = 3, screening_prob = 0.25, screening_crit = "high risk"),
        list(setting = 3, screening_time = 3, screening_prob = 0.9, screening_crit = "low risk"),
        list(setting = 3, screening_time = 10, screening_prob = 0.5, screening_crit = "high risk")
        ## Add more custom combinations as needed
    )

    df <- results %>%
        filter(setup == 3, Method %in% method.values) %>%
        mutate(Method = factor(Method, method.values, method.labels)) %>%
        group_by(Seed, setting, n_train, n_cal, Criterion, Time, Probability) %>%
        mutate(Oracle_Screened = Screened[Method == "Oracle"]) %>%        
        group_by(setting, n_train, n_cal, Method, Criterion, Time, Probability) %>%
        mutate(Screened = Screened / n_test) %>%
        mutate(Precision = ifelse(Screened==0, NA, Precision),
               Recall = ifelse(Oracle_Screened==0, NA, Recall)) %>%
        summarise(
            across(c(Screened, Survival, Precision, Recall),
                   list(mean = ~mean(., na.rm=TRUE), se = ~sd(., na.rm=TRUE) / sqrt(n())),
                   .names = "{.col}_{.fn}"),
            N = n(),
            .groups = 'drop'
        ) %>%
        mutate(
            across(where(is.numeric), ~ifelse(is.nan(.), NA, .))
        )

    ## Function to filter df based on one config
    filter_one_config <- function(cfg) {
        df %>%
            filter(
                setting == cfg$setting,
                Time == cfg$screening_time,
                Probability == cfg$screening_prob,
                Criterion == cfg$screening_crit
            )
    }
    
    ## Apply the filtering to each config
    filtered_dfs <- map(configs, filter_one_config)

    plot_summary <- function(idx) {

        df.plot <- filtered_dfs[[idx]] %>% pivot_longer(
                                                    cols = c(Screened_mean, Survival_mean, Precision_mean, Recall_mean,
                                                             Screened_se, Survival_se, Precision_se, Recall_se),
                                                    names_to = c("Variable", ".value"),
                                                    names_pattern = "(.*)_(mean|se)"
                                                ) %>%
            mutate(Variable = factor(Variable, variable.values, variable.labels))

        df.shade <- df.plot %>%
            filter(Variable == "Survival Rate") %>%
            distinct(Probability, Variable, Criterion) %>%
            mutate(
                ymin = ifelse(Criterion == "high risk", 0, Probability),
                ymax = ifelse(Criterion == "high risk", Probability, 1)
            ) 

        if(configs[[idx]]$screening_crit == "low risk") {
            title.str <- sprintf("(a) Screening for: %s.\n",
                                 interpret_screening_rule(configs[[idx]]$screening_time,
                                                          configs[[idx]]$screening_prob,
                                                          configs[[idx]]$screening_crit))
        } else {
            title.str <- sprintf("(b) Screening for: %s.\n",
                                 interpret_screening_rule(configs[[idx]]$screening_time,
                                                          configs[[idx]]$screening_prob,
                                                          configs[[idx]]$screening_crit))
        }

        pp <- df.plot %>%
            ggplot(aes(x=n_cal, y=mean, color=Method, shape=Method)) +
            geom_rect(
                data = df.shade,
                aes(xmin = min(df.plot$n_cal), xmax = max(df.plot$n_cal), ymin = ymin, ymax = ymax),
                inherit.aes = FALSE,
                fill = "gray80", alpha = 0.5
            ) +
            geom_point() +
            geom_line() +
            geom_errorbar(aes(ymin = mean - 2*se, ymax = mean + 2*se), width = 0.1) +
            facet_wrap(.~Variable, nrow=1, scales="free") +
            scale_x_log10() +
            ylim(0,1) +
            geom_hline(
                data = df.plot %>% filter(Variable == "Survival Rate") %>% distinct(Probability, Variable),
                aes(yintercept = Probability),
                linetype = "dashed", color = "gray50"
            ) +
            scale_color_manual(name = "Method", values = custom_colors,
                               guide = guide_legend(order = 1)) +
            scale_shape_manual(name = "Method", values = custom_shapes,
                               guide = guide_legend(order = 1)) +   
            labs(
                x = "Calibration Sample Size",
                y = "Mean ± 2 SE",
                subtitle = title.str
            ) +
            theme_bw() +
            theme(legend.position = "right",
                  plot.subtitle = element_text(margin = margin(b = -5)))

        if (save_fig) {
            fname <- sprintf("figures/setup3_setting%d_T%d_P%.2f_%s_cal.pdf",
                             configs[[idx]]$setting,
                             configs[[idx]]$screening_time,
                             configs[[idx]]$screening_prob,
                             gsub(" ", "_", configs[[idx]]$screening_crit))
            ggsave(fname, pp, width = 9, height = 2.25)
        }
    }

    for(idx in 1:length(configs)) {
        cat(sprintf("Making plot %d of %d...\n", idx, length(configs)))
        plot_summary(idx)
    }
    
}

results <- load_results(3)
make_figures_setup_3(results)



#############
## SETUP 4 ##
#############

### Tables for Setup 4
make_tables_setup_4 <- function(results, plot.n_train) {
    ## Define specific parameter combinations
    configs <- list(
        list(setting = 1, screening_time = 6, screening_prob = 0.8, screening_crit = "low risk"),
        list(setting = 1, screening_time = 12, screening_prob = 0.8, screening_crit = "high risk"),
        list(setting = 2, screening_time = 2, screening_prob = 0.8, screening_crit = "low risk"),
        list(setting = 2, screening_time = 3, screening_prob = 0.25, screening_crit = "high risk"),
        list(setting = 3, screening_time = 3, screening_prob = 0.9, screening_crit = "low risk"),
        list(setting = 3, screening_time = 10, screening_prob = 0.5, screening_crit = "high risk")
#        list(setting = 4, screening_time = 3, screening_prob = 0.8, screening_crit = "low risk"),
#        list(setting = 4, screening_time = 3, screening_prob = 0.5, screening_crit = "high risk")
    )

    df <- results %>%
        filter(setup == 4, Method %in% method.values, n_train==plot.n_train) %>%
        mutate(Method = factor(Method, method.values, method.labels)) %>%
        group_by(Seed, setting, n_train, n_cal, surv_model_type, cens_model_type, Criterion, Time, Probability) %>% 
        mutate(Oracle_Screened = Screened[Method == "Oracle"]) %>%        
        group_by(setting, n_train, n_cal, surv_model_type, cens_model_type, Method, Criterion, Time, Probability) %>%
        mutate(Screened = Screened / n_test) %>% 
        mutate(Precision = ifelse(Screened==0, NA, Precision),
               Recall = ifelse(Oracle_Screened==0, NA, Recall)) %>%
        summarise(
            across(c(Screened, Survival, Precision, Recall),
                   list(mean = ~mean(., na.rm=TRUE), se = ~sd(., na.rm=TRUE) / sqrt(n())),
                   .names = "{.col}_{.fn}"),
            N = n(),
            .groups = 'drop'
        ) %>%
        mutate(
            across(where(is.numeric), ~ifelse(is.nan(.), NA, .))
        )

    ## Function to filter df based on one config
    filter_one_config <- function(cfg) {
        df %>%
            filter(
                setting == cfg$setting,
                Time == cfg$screening_time,
                Probability == cfg$screening_prob,
                Criterion == cfg$screening_crit
            )
    }

    ## Apply the filtering to each config
    filtered_dfs <- map(configs, filter_one_config)

    make_table <- function(idx) {
        df.tab <- filtered_dfs[[idx]] %>%
            mutate(
                color = case_when(
                    Criterion == "low risk" & Survival_mean < Probability ~ "red",
                    Criterion == "high risk" & Survival_mean > Probability ~ "red",
                    TRUE ~ "black"
                ),                
                Screened = ifelse(
                    is.na(Screened_mean) | is.na(Screened_se),
                    "NA",
                    sprintf("%.3f ± %.3f", Screened_mean, 2 * Screened_se)
                ),
                Survival = ifelse(Screened_mean==0, "NA", 
                           ifelse(is.na(Survival_mean),
                                  "NA",
                                  sprintf("\\textcolor{%s}{%.3f ± %.3f}", color, Survival_mean, 2 * Survival_se)
                                  )),
                Precision = ifelse(
                    is.na(Precision_mean),
                    "NA",
                    sprintf("%.3f ± %.3f", Precision_mean, 2 * Precision_se)
                ),
                Recall = ifelse(
                    is.na(Recall_mean),
                    "NA",
                    sprintf("%.3f ± %.3f", Recall_mean, 2 * Recall_se)
                )
            ) %>%
            mutate(`Survival Model` = factor(surv_model_type, c("grf", "survreg", "rf", "cox"),
                                             c("grf", "survreg", "rf", "Cox"))) %>%
            mutate(`Censoring Model` = factor(cens_model_type, c("grf", "cox"),
                                                   c("Censoring Model: grf", "Censoring Model: Cox"))) %>%
            select(`Censoring Model`, `Survival Model`, Method, Screened, Survival, Precision, Recall) %>%
            arrange(`Censoring Model`, `Survival Model`, Method)

        target.str <- sprintf("%s.",
                             interpret_screening_rule(configs[[idx]]$screening_time,
                                                      configs[[idx]]$screening_prob,
                                                      configs[[idx]]$screening_crit,
                                                      latex=TRUE))
        
        print(df.tab)
        caption.str <- sprintf("Performance metrics for %s", target.str)

        latex_table <- df.tab %>%
            select(-`Censoring Model`) %>%
            kable(format = "latex", booktabs = TRUE, longtable = FALSE,
                  escape = FALSE,
                  align = "lccccc",
                  caption = NULL) %>%  # <- no caption inside
            collapse_rows(columns = 1, valign = "top", latex_hline = "none") %>%
            group_rows(index = table(df.tab$`Censoring Model`), escape = FALSE)
        
        fname <- sprintf("tables/setup4_setting%d_ntrain%d_%s.txt",
                         configs[[idx]]$setting,
                         plot.n_train,                         
                         gsub(" ", "_", configs[[idx]]$screening_crit))
        writeLines(latex_table, fname)
        cat(sprintf("Table saved on %s\n.", fname))
        
    }
    for(idx in 1:length(configs)) {
        cat(sprintf("Making table %d of %d...\n", idx, length(configs)))
        make_table(idx)
    }   
}
results <- load_results(4)
make_tables_setup_4(results, plot.n_train=1000)


#############
## SETUP 5 ##
#############

make_figures_setup_5 <- function(results) {
    ## Define specific parameter combinations
    configs <- list(
        list(setting = 1, screening_time = 6, screening_prob = 0.8, screening_crit = "low risk"),
        list(setting = 1, screening_time = 12, screening_prob = 0.8, screening_crit = "high risk"),
        list(setting = 2, screening_time = 2, screening_prob = 0.8, screening_crit = "low risk"),
        list(setting = 2, screening_time = 3, screening_prob = 0.25, screening_crit = "high risk"),
        list(setting = 3, screening_time = 3, screening_prob = 0.9, screening_crit = "low risk"),
        list(setting = 3, screening_time = 10, screening_prob = 0.5, screening_crit = "high risk")
        ## Add more custom combinations as needed
    )

    df <- results %>%
        filter(setup == 5, Method %in% method.values) %>%
        mutate(Method = factor(Method, method.values, method.labels)) %>%
        group_by(Seed, setting, n_train, n_train_cens, n_cal, Criterion, Time, Probability) %>%
        mutate(Oracle_Screened = Screened[Method == "Oracle"]) %>%        
        group_by(setting, n_train, n_train_cens, n_cal, Method, Criterion, Time, Probability) %>%
        mutate(Screened = Screened / n_test) %>% 
        mutate(Precision = ifelse(Screened==0, NA, Precision),
               Recall = ifelse(Oracle_Screened==0, NA, Recall)) %>%
        summarise(
            across(c(Screened, Survival, Precision, Recall),
                   list(mean = ~mean(., na.rm=TRUE), se = ~sd(., na.rm=TRUE) / sqrt(n())),
                   .names = "{.col}_{.fn}"),
            N = n(),
            .groups = 'drop'
        ) %>%
        mutate(
            across(where(is.numeric), ~ifelse(is.nan(.), NA, .))
        )

    ## Function to filter df based on one config
    filter_one_config <- function(cfg) {
        df %>%
            filter(
                setting == cfg$setting,
                Time == cfg$screening_time,
                Probability == cfg$screening_prob,
                Criterion == cfg$screening_crit
            )
    }
    
    ## Apply the filtering to each config
    filtered_dfs <- map(configs, filter_one_config)

    plot_summary <- function(idx) {

        df.plot <- filtered_dfs[[idx]] %>% pivot_longer(
                                                    cols = c(Screened_mean, Survival_mean, Precision_mean, Recall_mean,
                                                             Screened_se, Survival_se, Precision_se, Recall_se),
                                                    names_to = c("Variable", ".value"),
                                                    names_pattern = "(.*)_(mean|se)"
                                                ) %>%
            mutate(Variable = factor(Variable, variable.values, variable.labels))

        df.shade <- df.plot %>%
            filter(Variable == "Survival Rate") %>%
            distinct(Probability, Variable, Criterion) %>%
            mutate(
                ymin = ifelse(Criterion == "high risk", 0, Probability),
                ymax = ifelse(Criterion == "high risk", Probability, 1)
            ) 

        title.str <- sprintf("Screening for: %s.\n",
                             interpret_screening_rule(configs[[idx]]$screening_time,
                                                      configs[[idx]]$screening_prob,
                                                      configs[[idx]]$screening_crit))

        pp <- df.plot %>%
            ggplot(aes(x=n_train_cens, y=mean, color=Method, shape=Method)) +
            geom_rect(
                data = df.shade,
                aes(xmin = min(df.plot$n_train_cens), xmax = max(df.plot$n_train_cens), ymin = ymin, ymax = ymax),
                inherit.aes = FALSE,
                fill = "gray80", alpha = 0.5
            ) +
            geom_point() +
            geom_line() +
            geom_errorbar(aes(ymin = mean - 2*se, ymax = mean + 2*se), width = 0.1) +
            facet_wrap(.~Variable, nrow=1, scales="free") +
            scale_x_log10() +
            ylim(0,1) +
            geom_hline(
                data = df.plot %>% filter(Variable == "Survival Rate") %>% distinct(Probability, Variable),
                aes(yintercept = Probability),
                linetype = "dashed", color = "gray50"
            ) +
            scale_color_manual(name = "Method", values = custom_colors,
                               guide = guide_legend(order = 1)) +
            scale_shape_manual(name = "Method", values = custom_shapes,
                               guide = guide_legend(order = 1)) +   
            labs(
                x = "Training Sample Size for Censoring Model",
                y = "Mean ± 2 SE",
                subtitle = title.str
            ) +
            theme_bw() +
            theme(legend.position = "right",
                  plot.subtitle = element_text(margin = margin(b = -5)))

        if (save_fig) {
            fname <- sprintf("figures/setup5_setting%d_T%d_P%.2f_%s.pdf",
                             configs[[idx]]$setting,
                             configs[[idx]]$screening_time,
                             configs[[idx]]$screening_prob,
                             gsub(" ", "_", configs[[idx]]$screening_crit))
            ggsave(fname, pp, width = 9, height = 2.25)
        }
    }

    for(idx in 1:length(configs)) {
        cat(sprintf("Making plot %d of %d...\n", idx, length(configs)))
        plot_summary(idx)
    }
    
}

results <- load_results(5)
make_figures_setup_5(results)

#############
## SETUP 6 ##
#############

make_figures_setup_6 <- function(results) {
    ## Define specific parameter combinations
    configs <- list(
        list(setting = 1, screening_time = 6, screening_prob = 0.8, screening_crit = "low risk"),
        list(setting = 1, screening_time = 12, screening_prob = 0.8, screening_crit = "high risk"),
        list(setting = 2, screening_time = 2, screening_prob = 0.8, screening_crit = "low risk"),
        list(setting = 2, screening_time = 3, screening_prob = 0.25, screening_crit = "high risk"),
        list(setting = 3, screening_time = 3, screening_prob = 0.9, screening_crit = "low risk"),
        list(setting = 3, screening_time = 10, screening_prob = 0.5, screening_crit = "high risk")
        ## Add more custom combinations as needed
    )

    df <- results %>%
        filter(setup == 6, Method %in% method.values) %>%
        mutate(Method = factor(Method, method.values, method.labels)) %>%
        group_by(Seed, setting, n_train, n_train_cens, num_feat_censor, n_cal, Criterion, Time, Probability) %>%
        mutate(Oracle_Screened = Screened[Method == "Oracle"]) %>%        
        group_by(setting, n_train, n_train_cens, num_feat_censor, n_cal, Method, Criterion, Time, Probability) %>%
        mutate(Screened = Screened / n_test) %>% 
        mutate(Precision = ifelse(Screened==0, NA, Precision),
               Recall = ifelse(Oracle_Screened==0, NA, Recall)) %>%
        summarise(
            across(c(Screened, Survival, Precision, Recall),
                   list(mean = ~mean(., na.rm=TRUE), se = ~sd(., na.rm=TRUE) / sqrt(n())),
                   .names = "{.col}_{.fn}"),
            N = n(),
            .groups = 'drop'
        ) %>%
        mutate(
            across(where(is.numeric), ~ifelse(is.nan(.), NA, .))
        )

    ## Function to filter df based on one config
    filter_one_config <- function(cfg) {
        df %>%
            filter(
                setting == cfg$setting,
                Time == cfg$screening_time,
                Probability == cfg$screening_prob,
                Criterion == cfg$screening_crit
            )
    }
    
    ## Apply the filtering to each config
    filtered_dfs <- map(configs, filter_one_config)

    plot_summary <- function(idx) {

        df.plot <- filtered_dfs[[idx]] %>% pivot_longer(
                                                    cols = c(Screened_mean, Survival_mean, Precision_mean, Recall_mean,
                                                             Screened_se, Survival_se, Precision_se, Recall_se),
                                                    names_to = c("Variable", ".value"),
                                                    names_pattern = "(.*)_(mean|se)"
                                                ) %>%
            mutate(Variable = factor(Variable, variable.values, variable.labels))

        df.shade <- df.plot %>%
            filter(Variable == "Survival Rate") %>%
            distinct(Probability, Variable, Criterion) %>%
            mutate(
                ymin = ifelse(Criterion == "high risk", 0, Probability),
                ymax = ifelse(Criterion == "high risk", Probability, 1)
            ) 

        title.str <- sprintf("Screening for: %s.\n",
                             interpret_screening_rule(configs[[idx]]$screening_time,
                                                      configs[[idx]]$screening_prob,
                                                      configs[[idx]]$screening_crit))

        pp <- df.plot %>%
            ggplot(aes(x=num_feat_censor, y=mean, color=Method, shape=Method)) +
            geom_rect(
                data = df.shade,
                aes(xmin = min(df.plot$num_feat_censor), xmax = max(df.plot$num_feat_censor), ymin = ymin, ymax = ymax),
                inherit.aes = FALSE,
                fill = "gray80", alpha = 0.5
            ) +
            geom_point() +
            geom_line() +
            geom_errorbar(aes(ymin = mean - 2*se, ymax = mean + 2*se), width = 0.1) +
            facet_grid(n_train ~ Variable, scales = "free", 
                       labeller = labeller(n_train = function(n) paste0("Train sample size = ", n))) +
            scale_x_log10() +
            ylim(0,1) +
            geom_hline(
                data = df.plot %>% filter(Variable == "Survival Rate") %>% distinct(Probability, Variable),
                aes(yintercept = Probability),
                linetype = "dashed", color = "gray50"
            ) +
            scale_color_manual(name = "Method", values = custom_colors,
                               guide = guide_legend(order = 1)) +
            scale_shape_manual(name = "Method", values = custom_shapes,
                               guide = guide_legend(order = 1)) +   
            labs(
                x = "Number of Features for Censoring Model",
                y = "Mean ± 2 SE",
                subtitle = title.str
            ) +
            theme_bw() +
            theme(legend.position = "right",
                  plot.subtitle = element_text(margin = margin(b = -5)))

        if (save_fig) {
            fname <- sprintf("figures/setup6_setting%d_T%d_P%.2f_%s.pdf",
                             configs[[idx]]$setting,
                             configs[[idx]]$screening_time,
                             configs[[idx]]$screening_prob,
                             gsub(" ", "_", configs[[idx]]$screening_crit))
            ggsave(fname, pp, width = 9, height = 4)
        }
    }

    for(idx in 1:length(configs)) {
        cat(sprintf("Making plot %d of %d...\n", idx, length(configs)))
        plot_summary(idx)
    }
    
}

results <- load_results(6)
make_figures_setup_6(results)


##########
## DATA ##
##########

load_results_data <- function() {
    idir <- sprintf("results_hpc/data")
    ifile.list <- list.files(idir)
    results <- do.call("rbind", lapply(ifile.list, function(ifile) {
        df <- read_delim(sprintf("%s/%s", idir, ifile), delim=",", col_types=cols(), guess_max=2)
    }))
    return(results)
}


### Tables for real data
make_tables_data <- function(results, surv_model_type.plot) {
    ## Define specific parameter combinations
    configs <- list(
        list(time_q = 0.1, screening_prob = 0.8, screening_crit = "low risk"),
        list(time_q = 0.1, screening_prob = 0.8, screening_crit = "high risk"),
        list(time_q = 0.9, screening_prob = 0.25, screening_crit = "low risk"),
        list(time_q = 0.9, screening_prob = 0.25, screening_crit = "high risk")
    )

    df <- results %>%
        filter(Method %in% method.values, surv_model_type==surv_model_type.plot, train_prop_sub==1) %>%
        mutate(Method = factor(Method, method.values, method.labels)) %>%
        group_by(dataset, n_train, train_prop_sub, n_cal, n_test, Method, Criterion, Time_q, Time, Probability) %>%
        mutate(Screened = Screened / n_test) %>%
        summarise(
            across(c(Screened, Survival_lower, Survival_upper),
                   list(mean = ~mean(., na.rm=TRUE), se = ~sd(., na.rm=TRUE) / sqrt(n())),
                   .names = "{.col}_{.fn}"),
            N = n(),
            .groups = 'drop'
        ) %>%
        mutate(
            across(where(is.numeric), ~ifelse(is.nan(.), NA, .))
        )

    ## Function to filter df based on one config
    filter_one_config <- function(cfg) {
        df %>%
            filter(
                Time_q == cfg$time_q,
                Probability == cfg$screening_prob,
                Criterion == cfg$screening_crit
            )
    }

    ## Apply the filtering to each config
    filtered_dfs <- map(configs, filter_one_config)

    ##prop.shift.values <- c(1,0)
    ##model.quality.labels <- c("Low-Quality Model", "High-Quality Model")

    make_table <- function(idx) {
        df.tab <- filtered_dfs[[idx]] %>%
            mutate(
                color = case_when(
                    Criterion == "low risk" & Survival_upper_mean < Probability ~ "red",
                    Criterion == "low risk" & Survival_lower_mean - 2 * Survival_lower_se < Probability ~ "orange",
                    Criterion == "high risk" & Survival_lower_mean > Probability ~ "red",
                    Criterion == "high risk" & Survival_upper_mean + 2 * Survival_upper_se > Probability ~ "orange",
                    TRUE ~ "darkgreen"
                ),                
                Screened = ifelse(
                    is.na(Screened_mean) | is.na(Screened_se),
                    "NA",
                    sprintf("%.3f ± %.3f", Screened_mean, 2 * Screened_se)
                ),
                `Survival (lower bound)` = ifelse(Screened_mean==0, "NA", 
                                           ifelse(is.na(Survival_lower_mean),
                                                  "NA",
                                                  sprintf("\\textcolor{%s}{%.3f ± %.3f}",
                                                          color, Survival_lower_mean, 2 * Survival_lower_se)
                                                  )),
                `Survival (upper bound)` = ifelse(Screened_mean==0, "NA", 
                                           ifelse(is.na(Survival_upper_mean),
                                                  "NA",
                                                  sprintf("\\textcolor{%s}{%.3f ± %.3f}",
                                                          color, Survival_upper_mean, 2 * Survival_upper_se)
                                                  )),
                Dataset = dataset
            ) %>%
            ##mutate(`Model quality` = factor(prop_shift, prop.shift.values, model.quality.labels)) %>%
            select(Dataset, Method, Screened, `Survival (lower bound)`, `Survival (upper bound)`) %>%
            arrange(Dataset, Method)

        target.str <- sprintf("%s.",
                             interpret_screening_rule(configs[[idx]]$screening_time,
                                                      configs[[idx]]$screening_prob,
                                                      configs[[idx]]$screening_crit,
                                                      latex=TRUE))
        
        print(df.tab)
        caption.str <- sprintf("Performance metrics for %s", target.str)

        latex_table <- df.tab %>%
            select(-`Dataset`) %>%
            kable(format = "latex", booktabs = TRUE, longtable = FALSE,
                  escape = FALSE,
                  align = "lccccc",
                  caption = NULL) %>%  # <- no caption inside
            collapse_rows(columns = 1, valign = "top", latex_hline = "none") %>%    # no lines between collapsed rows
            group_rows(index = table(df.tab$Dataset), escape = FALSE)
        
        fname <- sprintf("tables/data_%s_t%s_p%s_%s.txt",
                         surv_model_type.plot,
                         configs[[idx]]$time_q,
                         configs[[idx]]$screening_prob,
                         gsub(" ", "_", configs[[idx]]$screening_crit))
        writeLines(latex_table, fname)
        cat(sprintf("Table saved on %s\n.", fname))
        
    }
    for(idx in 1:length(configs)) {
        cat(sprintf("Making table %d of %d...\n", idx, length(configs)))
        make_table(idx)
    }   
}

results <- load_results_data()
make_tables_data(results, surv_model_type.plot="grf")
make_tables_data(results, surv_model_type.plot="cox")
make_tables_data(results, surv_model_type.plot="survreg")

## Meta table ##

make_metatable_data <- function(results) {
    ## Define specific parameter combinations
    configs <- list(
        list(time_q = 0.1, screening_prob = 0.8, screening_crit = "low risk"),
        list(time_q = 0.1, screening_prob = 0.8, screening_crit = "high risk"),
        list(time_q = 0.9, screening_prob = 0.25, screening_crit = "low risk"),
        list(time_q = 0.9, screening_prob = 0.25, screening_crit = "high risk")
    )

    df <- results %>%
        filter(Method %in% method.values, train_prop_sub==1) %>%
        mutate(Method = factor(Method, method.values, method.labels)) %>%
        group_by(dataset, n_train, train_prop_sub, n_cal, n_test, surv_model_type, Method, Criterion, Time_q, Time, Probability) %>%
        mutate(Screened = Screened / n_test) %>%
        summarise(
            across(c(Screened, Survival_lower, Survival_upper),
                   list(mean = ~mean(., na.rm=TRUE), se = ~sd(., na.rm=TRUE) / sqrt(n())),
                   .names = "{.col}_{.fn}"),
            N = n(),
            .groups = 'drop'
        ) %>%
        mutate(
            across(where(is.numeric), ~ifelse(is.nan(.), NA, .))
        )
           
    ## Function to filter df based on one config
    filter_one_config <- function(cfg) {
        df %>%
            filter(
                Time_q == cfg$time_q,
                Probability == cfg$screening_prob,
                Criterion == cfg$screening_crit
            )
    }

    ## Apply the filtering to each config
    filtered_dfs <- map(configs, filter_one_config)

    make_summary <- function(idx) {
        df.tab <- filtered_dfs[[idx]] %>%
            mutate(                
                Validity = case_when(
                    Criterion == "low risk" & Survival_upper_mean < Probability ~ "invalid",
                    Criterion == "low risk" & Survival_lower_mean - 2 * Survival_lower_se < Probability ~ "dubious",
                    Criterion == "high risk" & Survival_lower_mean > Probability ~ "invalid",
                    Criterion == "high risk" & Survival_upper_mean + 2 * Survival_upper_se > Probability ~ "dubious",
                    TRUE ~ "valid"
                ),                
                Screened = Screened_mean,
                Dataset = dataset,
                Task = idx,
                `Survival Model` = factor(surv_model_type, c("grf", "survreg", "rf", "cox"),
                                                 c("grf", "survreg", "rf", "Cox"))
            ) %>%
            ##mutate(`Model quality` = factor(prop_shift, prop.shift.values, model.quality.labels)) %>%
            select(Dataset, `Survival Model`, Method, Task, Screened, Validity) %>%
            arrange(Dataset, `Survival Model`, Method)
        return(df.tab)
    }
    df.summary <- do.call("rbind", lapply(1:length(configs), function(idx) make_summary(idx)))

    df.tab <- df.summary %>%
        mutate(placeholder = 1) %>%  # dummy variable
        pivot_wider(
            names_from = Validity,
            values_from = placeholder,
            values_fill = 0
        ) %>%
        group_by(Method, `Survival Model`) %>%
        select(Method, `Survival Model`, everything()) %>%
        summarise(`Screened Proportion`=mean(Screened), Valid=mean(valid), Dubious=mean(dubious), Invalid=mean(invalid))
    
    latex_table <- df.tab %>%
        mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
      kable(format = "latex", booktabs = TRUE, longtable = FALSE,
              escape = FALSE,
              align = "lccccc",
              caption = NULL) %>%  # <- no caption inside
        collapse_rows(columns = 1, valign = "top", latex_hline = "none") %>%
    add_header_above(c(" " = 3, "Verification Outcome" = 3)) # <- add group title
    
    fname <- sprintf("tables/data_summary.txt")
    writeLines(latex_table, fname)
    cat(sprintf("Table saved on %s\n.", fname))   
    
}


make_metatable_data(results)
