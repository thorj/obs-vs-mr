library(tidyverse)
library(patchwork)
library(ggrepel)
box::use(data.table[fread, fwrite])

### Load data
data_for_bayes <- readRDS(file = "data/data_for_analysis.rds")
f <- data_for_bayes$f
f2 <- data_for_bayes$f2
group <- data_for_bayes$map_phenotype_to_group
protein_map <- data_for_bayes$map_somamer_to_egs
h1 <- data_for_bayes$h1
vp <- data_for_bayes$variant_profile
harmon <- fread("../../papers/02_obs_vs_mr/data/forward_mr_harmonized.csv")
obs <- fread("../../papers/02_obs_vs_mr/data/observational_results_full.csv")

f <- 
    f |>
    mutate(group = 
               case_when(
                   event == "did" ~ "Blood pressure",
                   event == "s_glu" ~ "Diabetes",
                   TRUE ~ group
               )) |>
    rename(phenotype = phenotype2)

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
    #geom_col(position = position_dodge(width = 1)) +
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
    labs(x = "Agreement ratio", y = "", color = "Contains coding variant", shape = "Contains coding variant") +
    coord_cartesian(xlim = c(0, 1)) +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom",
          axis.title = element_text(face = "bold"))


pheno_dumb <- 
    f2 |> 
    select(phenotype2, group, event) |>
    inner_join(h1 |>
                   group_by(event, any_coding) |>
                   summarize(n = n(),
                             m = mean(agree))) |>
    ungroup() |>
    mutate(phenotype2 = factor(phenotype2), 
           any_coding = factor(any_coding, levels = c(0, 1), 
                               labels = c("No", "Yes"))) 

pd_plt <- 
    pheno_dumb |>
    ggplot(aes(x = m, y = phenotype2)) +
    geom_vline(xintercept = 0.5, lty = 2, alpha = 0.5) +
    geom_line() +
    geom_point(aes(color = any_coding, shape = any_coding), size = 3) +
    geom_label_repel(aes(color = any_coding, label = n), size = 3) +
    scale_color_manual(values = c("#648FFF", "#FFB000")) +
    scale_x_continuous(breaks = 0:10/10, labels = scales::percent) +
    labs(x = "Agreement ratio", y = "", color = "Contains coding variant", shape = "Contains coding variant") +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom",
          axis.title = element_text(face = "bold"))

finalp <-
    (gd_plt |  pd_plt | hmf_agree_cod) + plot_annotation(tag_levels = "A", tag_prefix = "Fig. ")

#ggsave(filename = "img/highres/figure6_coding_vs_covars.pdf", plot = finalp, width = 20, height = 6)
ggsave(filename = "img/lowres/figure6_coding_vs_covars.png", plot = finalp, width = 20, height = 6, dpi = 300)

## Export tables
fwrite(x = group_dumb |> select(-groupf), file = "tables/tables_new/figure6A_table.csv")
fwrite(x = pheno_dumb, file = "tables/tables_new/figure6B_table.csv")
fwrite(x = hmf_df, file = "tables/tables_new/figure6C_table.csv")
