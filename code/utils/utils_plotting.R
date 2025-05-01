library(tidyverse)

get_custom_colors_model <- function() {
    custom_colors <- c(
        "Model" = "black",
        "Oracle" = "forestgreen",
        "KM" = "blue",
        "CSB" = "red"
    )
}
get_custom_shapes_model <- function() {
    custom_shapes <- c(
        "Model" = 16,
        "Oracle" = 8,
        "KM" = 18,
        "CSB" = 17
    )
}

plot_survival_curves <- function(predictions, individuals = NULL) {
    # Predict survival curves
    survival_curves <- as.data.frame(predictions$predictions)
    time_points <- predictions$time.points

    df <- survival_curves
    colnames(df) <- time_points
    df$individual <- seq(1, nrow(df))
    df <- df %>%
        pivot_longer(cols = -individual,
                     names_to = "time",
                     values_to = "survival_probability") %>%
        mutate(time = as.numeric(time))

    # Filter for specific individuals if provided
    if (!is.null(individuals)) {
        df <- df %>% filter(individual %in% individuals)
    }

    # Plot using ggplot
    ggplot(df, aes(x = time, y = survival_probability, color = factor(individual))) +
        geom_line(linewidth = 1) +
        labs(title = "Survival Curves",
             x = "Time",
             y = "Survival Probability",
             color = "Individual") +
        theme_minimal()
}

plot_csb <- function(csb, patient.list, nrow=2, oracle_pred=NULL,
                             screening_time=NULL, screening_prob=NULL, screening_crit=NULL) {
    df <- make_calibration_tidy(csb, external_preds=list("oracle"=oracle_pred))

    if(!is.null(screening_time) & !is.null(screening_prob) & !is.null(screening_crit)) {
        stopifnot(screening_crit %in% c("low risk", "high risk"))

        ## Selection with black-box model
        time.points <- as.numeric(colnames(csb$model_pred))
        selections.model <- select_patients_point(time.points, csb$model_pred, screening_time, screening_prob, screening_crit)$selected

        ## Selections with conformal method
        selections.calib <- select_patients_band(time.points, csb$lower, csb$upper, screening_time, screening_prob, screening_crit)$selected

        if(!is.null(oracle_pred)) {
            time.points.oracle <- as.numeric(colnames(oracle_pred))
            selections.oracle <- select_patients_point(time.points.oracle, oracle_pred, screening_time, screening_prob, screening_crit)$selected
        } else {
            selections.oracle <- c()
        }
        axis_add <- c(0.25, 0.05)
        title.str <- sprintf("Screening for: %s.\n", interpret_screening_rule(screening_time, screening_prob, screening_crit))
    } else {
        selections.model <- c()
        selections.calib <- c()
        selections.oracle <- c()
        axis_add <- c(0,0)
        title.str <- NULL
    }

    if(!is.null(screening_crit)) {
        if(screening_crit=="low risk") {
            selected.model.str <- sprintf("* (low risk, S > %.2f, black-box)", screening_prob)
            selected.calib.str <- sprintf("* (low risk, S > %.2f, conformal)", screening_prob)
        } else {
            selected.model.str <- sprintf("* (high risk, S < %.2f, black-box)", screening_prob)
            selected.calib.str <- sprintf("* (high risk, S < %.2f, conformal)", screening_prob)
        }
    }

    min.time <- min(df$time)
    max.time <- max(df$time)

    pp <- df |>
        filter(patient_id %in% patient.list, source %in% c("model", "oracle")) |>
        mutate(Survival = factor(source, c("model", "oracle"), c("Model", "Oracle")),
               Selected.model = ifelse(patient_id %in% selections.model, selected.model.str, ""),
               Selected.calib = ifelse(patient_id %in% selections.calib, selected.calib.str, "")) |>
        ggplot(aes(x = time, y = value, color = Survival, linetype = Survival, alpha = Survival)) +
        geom_ribbon(data = df |> filter(patient_id %in% patient.list, source %in% c("lower", "upper")) |>
                          pivot_wider(names_from = source, values_from = value),
                    aes(x = time, ymin = lower, ymax = upper, fill = "Survival Band"),
                    inherit.aes = FALSE, alpha = 0.5) +  # Shaded area between lower and upper curves
        geom_line() +
        geom_vline(xintercept=screening_time, linetype=3) +
        geom_hline(yintercept=screening_prob, linetype=3) +
        geom_text(x = 0.75*max.time, y = -0.05, aes(label = Selected.model), check_overlap = TRUE, color="darkorange", show.legend=FALSE) +
        geom_text(x = 0.75*max.time, y = -0.15, aes(label = Selected.calib), check_overlap = TRUE, color="red2", show.legend=FALSE) +
        facet_wrap(patient_id ~ ., nrow = nrow, labeller = label_both) +
        ylab("Survival Probability") +
        xlab("Time") +
        scale_y_continuous(limits=c(0,1), breaks=c(0,0.25,0.5,0.75,1), expand = expansion(add = axis_add)) +
        scale_color_manual(values = c("black", "darkgreen")) +
        scale_linetype_manual(values = c(1, 2)) +
        scale_alpha_manual(values = c(1, 1)) +
        scale_fill_manual(values = c("gray"), name = "Calibration") +  # Color for shading
        theme_bw() +
        ggtitle(title.str) +
        theme(
            text = element_text(size = 16),
            axis.title = element_text(size = 18),
            axis.text = element_text(size = 14),
            plot.title = element_text(size = 20, face = "bold")
        )
    pp
}


plot_survival_panel <- function(pred_list,
                                band_list = NULL,
                                screening_time = NULL,
                                screening_prob = NULL,
                                screening_crit = NULL,
                                patient_names = NULL,
                                custom_colors = NULL,
                                custom_linetypes = NULL,
                                custom_shapes = NULL,
                                model_order = NULL,
                                facet_label = NULL,
                                x_label = "Time",
                                save_fig = FALSE) {

    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(RColorBrewer)
    library(ggnewscale)

    if(is.null(custom_colors)) {
        custom_colors <- get_custom_colors_model()
    }
    if(is.null(custom_shapes)) {
        custom_shapes <- get_custom_shapes_model()
    }
    if(is.null(custom_linetypes)) {
        custom_linetypes <- c(
            "Model" = 1,
            "Oracle" = 2,
            "KM" = 4
        )
    }
    if(is.null(model_order)) {
        model_order <- c("Oracle", "CSB", "Model", "KM")
    }
    
    model_names <- names(pred_list)
    n_models <- length(model_names)
    n_patients <- nrow(pred_list[[1]])
    time_points <- as.numeric(colnames(pred_list[[1]]))
    axis_add <- c(0.05, 0.05)

    ## Patient IDs
    patient_ids <- if (!is.null(patient_names)) patient_names else 1:n_patients
    patient_ids <- factor(patient_ids, levels = patient_ids)

    if (!is.null(facet_label)) {
        facet_group <- facet_label
    } else {
        facet_group <- NULL
    }
    
    ## Long-format survival predictions
    df_list <- lapply(model_names, function(name) {
        pred_list[[name]] |>
        as.data.frame() |>
        mutate(patient_id = patient_ids) |>
        pivot_longer(-patient_id, names_to = "time", values_to = "value") |>
        mutate(model_line = name, time = as.numeric(time))
    })
    df <- bind_rows(df_list)
    df$model_line <- factor(df$model_line, levels = model_order)
    
    ## Ribbon (uncertainty band)
    ribbon_df <- NULL
    if (!is.null(band_list)) {
        stopifnot(length(band_list) == 1)
        band_name <- names(band_list)[1]
        band_obj <- band_list[[1]]
        lower <- band_obj$lower
        upper <- band_obj$upper

        ribbon_df <- data.frame(
            patient_id = rep(patient_ids, each = length(time_points)),
            time = rep(time_points, times = n_patients),
            lower = as.vector(t(lower)),
            upper = as.vector(t(upper)),
            band_label = band_name
        )
        ribbon_df$patient_id <- factor(ribbon_df$patient_id, levels = patient_ids)
    }
    if (!is.null(ribbon_df)) ribbon_df$facet_group <- facet_group
    
    ## Flagged patients
    flagged_df <- NULL
    if (!is.null(screening_time) & !is.null(screening_prob) & !is.null(screening_crit)) {
        all_flag_names <- model_order  # assume this includes band name too
        offset_seq <- seq(0, 0.1 * (length(all_flag_names) - 1), length.out = length(all_flag_names))
        offset_map <- setNames(offset_seq, all_flag_names)

        flagged_df_list <- list()

        ## Point-based flagging
        for (name in model_names) {
            mat <- pred_list[[name]]
            sel <- select_patients_point(time_points, mat, screening_time, screening_prob, screening_crit)$selected
            if (length(sel) > 0) {
                flagged_df_list[[paste0(name, "_point")]] <- data.frame(
                    patient_id = factor(patient_ids[sel], levels = patient_ids),
                    x = offset_map[[name]] * max(time_points),
                    y = 1.06,
                    model_flag = name
                )
            }
        }

        ## Band-based flagging (if provided)
        if (!is.null(band_list)) {
            band_name <- names(band_list)[1]
            band_obj <- band_list[[1]]
            lower_mat <- band_obj$lower
            upper_mat <- band_obj$upper

            sel_band <- select_patients_band(
                time_points, lower_mat, upper_mat,
                screening_time, screening_prob, screening_crit
            )$selected

            if (length(sel_band) > 0) {
                flagged_df_list[["band"]] <- data.frame(
                    patient_id = factor(patient_ids[sel_band], levels = patient_ids),
                    x = offset_map[[band_name]] * max(time_points),
                    y = 1.06,
                    model_flag = band_name
                )
            }
        }

        flagged_df <- bind_rows(flagged_df_list)
        if ("model_flag" %in% colnames(flagged_df) && nrow(flagged_df) > 0) {
            flagged_df$model_flag <- factor(flagged_df$model_flag, levels = model_order)
        } else {
            flagged_df <- NULL  # nothing to plot
        }        
    }
    if (!is.null(flagged_df)) flagged_df$facet_group <- facet_group
        
    ## Color & shape settings
    band_name <- names(band_list)[1]
    full_flag_names <- c(model_names, band_name)
    
    ## Legend title for flags
    screen_str <- if (!is.null(screening_crit)) {
                      if (screening_crit == "low risk") {
                          sprintf("Flagged as 'low risk'\n(i.e., P[T>%.2f] > %.2f)", screening_time, screening_prob)
                      } else {
                          sprintf("Flagged as 'high risk'\n(i.e., P[T>%.2f] < %.2f)", screening_time, screening_prob)
                      }
                  } else {
                      "Flagged"
                  }

    ## Build plot
    p <- ggplot() +
        ## Survival model lines
        geom_line(data = df,
                  aes(x = time, y = value, color = model_line, linetype = model_line, group = interaction(patient_id, model_line)),
                  linewidth = 1) +        
        scale_color_manual(name = "Survival Model", values = custom_colors,
                           guide = guide_legend(order = 1, override.aes = list(shape = NA))) +
        scale_linetype_manual(name = "Survival Model", values = custom_linetypes,
                              guide = guide_legend(order = 1))

    ## Ribbon band
    if (!is.null(ribbon_df)) {
        p <- p +
            geom_ribbon(data = ribbon_df,
                        aes(x = time, ymin = lower, ymax = upper, fill = band_label),
                        inherit.aes = FALSE, alpha = 0.5) +
            scale_fill_manual(name = "Uncertainty",
                              values = c("gray"),
                              guide = guide_legend(order = 2, override.aes = list(shape = NA, linetype = 0)))
    }

    ## Separate color scale for flagged points
    if (!is.null(flagged_df) && nrow(flagged_df) > 0) {
        p <- p + new_scale_color()
        
        p <- p +
            geom_point(data = flagged_df,
                       aes(x = x, y = y, color = model_flag, shape = model_flag),
                       size = 3, inherit.aes = FALSE) +
            scale_color_manual(name = screen_str, values = custom_colors,
                               guide = guide_legend(order = 3, override.aes = list(linetype = 0))) +
            scale_shape_manual(name = screen_str, values = custom_shapes,
                               guide = guide_legend(order = 3))
    }

    ## Screening threshold lines
    if (!is.null(screening_time)) {
        p <- p + geom_vline(xintercept = screening_time, linetype = 3)
    }
    if (!is.null(screening_prob)) {
        p <- p + geom_hline(yintercept = screening_prob, linetype = 3)
    }

    ## Final layout
    if(is.null(facet_group)) {
        p <- p + facet_grid(. ~ patient_id)
    } else {
        p <- p + facet_grid(facet_group ~ patient_id)
    }                
    
    p <- p +
        ylab("Survival Probability") +
        xlab(x_label) +
        scale_x_continuous(limits = c(0, NA), breaks = scales::pretty_breaks(n = 4)) +
        scale_y_continuous(limits = c(0, 1.07),
                           expand = expansion(add = axis_add),
                           breaks = c(0, 0.25, 0.5, 0.75, 1)) +
        theme_bw() +
        theme(
            text = element_text(size = 14),
            axis.title = element_text(size = 16),
            axis.text = element_text(size = 12),
            strip.text = element_text(size = 14),
            legend.title = element_text(size = 14),
            legend.text = element_text(size = 12),
            legend.box = "vertical",
            legend.key.width = unit(1, "cm")
        )

    ## Save if requested
    if (save_fig) {
        fname <- sprintf("figures/survival_panel_T%d_P%.2f_%s.pdf",
                         screening_time, screening_prob,
                         gsub(" ", "_", screening_crit))
        ggsave(fname, p, width = 12, height = 4)
    }

    return(p)
}

plot_survival_two_panels <- function(csb_list, facet_labels, oracle_pred,
                                     patient.list, screening_time, screening_prob, screening_crit,
                                     save_fig = FALSE, fname=NULL) {

    res.1 <- csb_list[[1]]

    pp.1 <- plot_survival_panel(
        pred_list = list(
            Model = res.1$model_pred[patient.list, ],
            Oracle = oracle_pred[patient.list, ],
            KM = res.1$km_pred[patient.list, ]
        ),
        band_list = list(`CSB` = list(
                             lower = res.1$lb_dr[patient.list, ],
                             upper = res.1$ub_dr[patient.list, ]
                         )),
        screening_time = screening_time,
        screening_prob = screening_prob,
        screening_crit = screening_crit,
        patient_names = paste0("Patient ", seq_along(patient.list)),
        facet_label = facet_labels[1],
        save_fig = FALSE
    )

    res.2 <- csb_list[[2]]

    pp.2 <- plot_survival_panel(
        pred_list = list(
            Model = res.2$model_pred[patient.list, ],
            Oracle = oracle_pred[patient.list, ],
            KM = res.2$km_pred[patient.list, ]
        ),
        band_list = list(`CSB` = list(
                             lower = res.2$lb_dr[patient.list, ],
                             upper = res.2$ub_dr[patient.list, ]
                         )),
        screening_time = screening_time,
        screening_prob = screening_prob,
        screening_crit = screening_crit,
        patient_names = paste0("Patient ", seq_along(patient.list)),
        facet_label = facet_labels[2],
        save_fig = FALSE
    )

    library(patchwork)

    ## Remove legends from both plots
    pp.1_noleg <- pp.1 + theme(legend.position = "none")
    pp.2_noleg <- pp.2 + theme(legend.position = "none")

    ## Extract legend from pp.2
    library(cowplot)
    legend_shared <- get_legend(pp.2 + theme(legend.position = "right"))

    ## Stack plots and attach legend
    final_plot <- plot_grid(
        plot_grid(pp.1_noleg, pp.2_noleg, ncol = 1, align = "v"),
        legend_shared,
        ncol = 2,
        rel_widths = c(1, 0.3)
    )

    ## Save if requested
    if (save_fig) {
        if(is.null(fname)) {
            fname <- sprintf("figures/survival_panel_setting%d_T%d_P%.2f_%s.pdf",
                             setting, screening_time, screening_prob,
                             gsub(" ", "_", screening_crit))
        }
        ggsave(fname, plot = final_plot, width = 12, height = 5.25, device = cairo_pdf)
    }
    return(final_plot)
}

make_calibration_tidy <- function(csb, external_preds = NULL) {
    df_model <- as_tibble(csb$model_pred) %>%
      mutate(source = "model", patient_id = row_number()) %>%
      pivot_longer(-c(source, patient_id), names_to = "time", values_to = "value")

    df_lower <- as_tibble(csb$lower) %>%
      mutate(source = "lower", patient_id = row_number()) %>%
      pivot_longer(-c(source, patient_id), names_to = "time", values_to = "value")

    df_upper <- as_tibble(csb$upper) %>%
      mutate(source = "upper", patient_id = row_number()) %>%
      pivot_longer(-c(source, patient_id), names_to = "time", values_to = "value")

    # External prediction curves
    df_external <- list()
    if (!is.null(external_preds)) {
        for (name in names(external_preds)) {
            df_external[[name]] <- as_tibble(external_preds[[name]]) %>%
              mutate(source = name, patient_id = row_number()) %>%
              pivot_longer(-c(source, patient_id), names_to = "time", values_to = "value")
        }
    }

    combined_df <- bind_rows(
        df_lower,
        df_upper,
        df_model,
        bind_rows(df_external)  # include external preds here
    ) %>%
      mutate(time = as.numeric(time))

    return(combined_df)
}
