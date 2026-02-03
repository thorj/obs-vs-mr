library(tidyverse)
library(patchwork)
library(ggrepel)

box::use(data.table[fread, fwrite])
source("R/utils.R")

### Load data
h1 <- fread("data/primary_mr_overview.csv")
f <- fread("data/phenotype_overview.csv")

hmf_df <-
    h1 |>
    mutate(hmf_n = ntile(hmf, n = 10)) |>
    mutate(ac = factor(any_coding, levels = c(0,1), labels = c("No", "Yes"))) |>
    mutate(hmf_n2 = factor(hmf_n, levels = 1:10, labels = paste0(10*1:10, "%"))) |>
    group_by(hmf_n2, ac) |>
    summarize(m = mean(agree),
              n = n())

hmf_agree_cod <-
    hmf_df |>
    ggplot(aes(x = hmf_n2, y = m, color = ac, group = ac)) +
    scale_color_manual(values=c("No" = "#648FFF", "Yes" = "#FFB000")) +
    scale_y_continuous(breaks = 0:10/10, labels = scales::percent, limit = c(0, 1)) +
    geom_point(aes(shape = ac), size = 3) +
    geom_line() +
    ggrepel::geom_label_repel(aes(label = n), size = 3) +
    labs(x = "Harmonic mean F-statistic decile", 
         y = "Observed agreement rate", 
         color = "Contains coding variant",
         shape = "Contains coding variant") +
    theme_bw(base_size = 12) +
    theme(axis.title = element_text(face = "bold"), 
          legend.position = "bottom")

group_dumb <- 
    h1 |>
    group_by(groupf, group, any_coding) |>
    summarize(m = mean(agree),
              n = n()) |>
    ungroup() |>
    mutate(group = factor(group), 
           any_coding = factor(any_coding, levels = c(0, 1), 
                               labels = c("No", "Yes")))

gd_plt <- 
    group_dumb |>
    ggplot(aes(x = m, y = group)) +
    geom_vline(xintercept = 0.5, lty = 2, alpha = 0.5) +
    geom_line() +
    geom_point(aes(color = any_coding, shape = any_coding), size = 3) +
    geom_label_repel(aes(color = any_coding, label = n), size = 3) +
    scale_color_manual(values = c("#648FFF", "#FFB000")) +
    scale_x_continuous(breaks = 0:10/10, labels = scales::percent) +
    labs(x = "Agreement ratio", y = "", 
         color = "Contains coding variant", 
         shape = "Contains coding variant") +
    coord_cartesian(xlim = c(0, 1)) +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom",
          axis.title = element_text(face = "bold"))


pheno_dumb <- 
    f |> 
    filter(N > 0) |>
    select(phenotype, group, event) |>
    inner_join(h1 |>
                   group_by(event, any_coding) |>
                   summarize(n = n(),
                             m = mean(agree))) |>
    ungroup() |>
    mutate(phenotype = factor(phenotype), 
           any_coding = factor(any_coding, levels = c(0, 1), 
                               labels = c("No", "Yes"))) 

pd_plt <- 
    pheno_dumb |>
    ggplot(aes(x = m, y = phenotype)) +
    geom_vline(xintercept = 0.5, lty = 2, alpha = 0.5) +
    geom_line() +
    geom_point(aes(color = any_coding, shape = any_coding), size = 3) +
    geom_label_repel(aes(color = any_coding, label = n), size = 3) +
    scale_color_manual(values = c("#648FFF", "#FFB000")) +
    scale_x_continuous(breaks = 0:10/10, labels = scales::percent) +
    labs(x = "Agreement ratio", 
         y = "", 
         color = "Contains coding variant", 
         shape = "Contains coding variant") +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom",
          axis.title = element_text(face = "bold"))

finalp <-
    (gd_plt |  pd_plt | hmf_agree_cod) + 
    plot_annotation(tag_levels = "A", tag_prefix = "Fig. ")

export_image(plot = finalp, 
             fig_name = "figure5_coding_vs_covars",
             width = 20, height = 6, dpi = 300)

## Export tables
#fwrite(x = group_dumb |> select(-groupf), file = "tables/tables_new/figure6A_table.csv")
#fwrite(x = pheno_dumb, file = "tables/tables_new/figure6B_table.csv")
#fwrite(x = hmf_df, file = "tables/tables_new/figure6C_table.csv")
