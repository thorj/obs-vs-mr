library(tidyverse)
library(patchwork)

box::use(data.table[fread, fwrite])
source("R/utils.R")

## Constants
other_lty <- 0

## ============================================================
## Load data
## ============================================================

primary_mr_summary <- readRDS(data_paths$primary_mr_summary)
primary_mr_agreement <- primary_mr_summary$primary_mr_agreement

## ============================================================
## Plot: Agreement between observationl and primary MR estimates 
## ============================================================

plot_df <- 
    primary_mr_agreement |>
    filter(include == 1)

observed_agreement <-
    plot_df |>
    ggplot(aes(x = order)) +
    geom_col(aes(y = rf, fill = "Forward MR"), alpha = 0.6, color = "black") +
    geom_col(aes(y = rr, fill = "Reverse MR"), 
             lty = other_lty, alpha = 0.6, color = "black") +
    geom_hline(yintercept = 0, lty = 1) +
    geom_hline(yintercept = -0.5, lty = 2, alpha = 0.7) +
    geom_hline(yintercept = 0.5, lty = 2, alpha = 0.7) +
    geom_label(aes(y = rf, label = rftext), size = 3, hjust = 1.1) +
    geom_label(aes(y = rr, label = rrtext), size = 3, hjust = -0.1) +
    labs(y = "Observed agreement ratio", x = "", fill = "") +
    scale_x_discrete(label = primary_mr_agreement$phenotype) +
    scale_y_continuous(labels = function(x) scales::percent(abs(x)),
                       breaks = -10:10/10) +
    scale_fill_manual(values = c("Forward MR" = "#882255",
                                 "Reverse MR" = "#332288"), drop = F) +
    coord_flip(ylim = c(-1.05, 1.15)) +
    theme_bw(base_size = 12) +
    theme(axis.title = element_text(face = "bold"),
          legend.position = "bottom",
          legend.box = "vertical",
          legend.margin = margin())

## ============================================================
## Export
## ============================================================

export_image(plot = observed_agreement, 
             fig_name = "figure4_agreement_ratios", 
             width = 12, height = 7, dpi = 300)
