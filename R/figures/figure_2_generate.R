library(tidyverse)
library(patchwork)
box::use(data.table[fread, fwrite])
source("R/utils.R")

### Load data
f <- fread("data/phenotype_overview.csv")
obs <- fread("data/observational_results_full.csv")

### CONST
other_lty <- 0
####

sig_df <- fread("data/reoccurring_data/bonferroni_sig_obs.csv")

sig_pan <- 
    sig_df |>
    ggplot(aes(x = order2, y = n)) +
    geom_col(aes(fill = trait, linetype = trait),
             color = "black", alpha = 0.6) + 
    scale_y_continuous(breaks = scales::extended_breaks(n = 10), 
                       limits = c(0, 4500)) +
    scale_x_discrete(labels = sig_df$phenotype) +
    scale_fill_manual(values = c("Disease" = "#1E88E5", 
                                 "Risk factor" = "orange"), 
                      drop = T) +
    scale_linetype_manual(values = c("Disease" = 1, 
                                     "Risk factor" = other_lty)) +
    geom_label(aes(label = n), size = 3, hjust = -0.1) +
    coord_flip() +
    labs(x = "", 
         y = "Number of phenotype-level significant SOMAmers", 
         fill = "", linetype = "") +
    theme_bw(base_size = 12) +
    theme(axis.title = element_text(face = "bold"),
          legend.position = "bottom")

### cis info
cis_info <- readRDS("data/reoccurring_data/cis_protein_data.rds")
prop_w_cis <- cis_info$overall_cis_prop
cis_df <- cis_info$bonferroni_with_cis
has_cis <- cis_info$has_cis

background_df <- data.frame(
    xmin = as.numeric(levels(has_cis$order)) - 0.5,
    xmax = as.numeric(levels(has_cis$order)) + 0.5,
    ymin = -Inf,
    ymax = Inf,
    fill = rep(c("grey90", "white"), length.out = length(levels(has_cis$order)))
)

cis_pan <- 
    cis_df |>
    mutate(ln = paste0(ciss, " / ", n)) |>
    ggplot(aes(x = order2, y = p)) +
    geom_hline(yintercept = prop_w_cis, lty = 2, alpha = 0.7) +
    geom_col(aes(fill = trait, linetype = trait), color = "black", alpha = 0.6) +
    scale_y_continuous(labels = scales::percent, breaks = 0:20/20, limits = c(0, 0.75)) +
    scale_x_discrete(labels = cis_df$phenotype) +
    scale_fill_manual(values = c("Disease" = "#1E88E5", "Risk factor" = "orange"), drop = T) +
    scale_linetype_manual(values = c("Disease" = 1, "Risk factor" = other_lty)) +
    geom_label(aes(label = ln), size = 3, hjust = -0.1) +
    coord_flip() +
    labs(x = "", y = "Proportion of SOMAmers with cis-instruments", fill = "", linetype = "") +
    theme_bw(base_size = 12) +
    theme(axis.title = element_text(face = "bold"),
          legend.position = "bottom")

### Bonferroni tail enrichment
tail_res <-
    lapply(obs$event |> unique(),
           FUN = function(x) {
               df <- obs |> filter(event == {{ x }}) |> mutate(intail = pval < 0.05/7288)
               m1 <- 
                   glm(intail ~ has_cis, data = df, family = binomial())
               res <- 
                   broom::tidy(m1) |> 
                   filter(term == "has_cis") |> 
                   mutate(term = {{ x }})
               return(res)
           }) |> 
    data.table::rbindlist()


tail_res_table <- 
    f |> 
    select(trait, group, phenotype, term = event, order) |>
    inner_join(tail_res) |>
    select(-order) |>
    mutate(conf.low = signif(estimate - 1.96 * std.error, digits = 3),
           conf.high = signif(estimate + 1.96 * std.error, digits = 3),
           p.value = signif(p.value, digits = 3),
           estimate = signif(estimate, digits = 3),
           std.error = signif(std.error, digits = 3),
           statistic = signif(statistic, digits = 3),
           significant = as.integer(p.value < 0.05)) |>
    relocate(conf.low, .before = "statistic") |>
    relocate(conf.high, .before = "statistic") 

#fwrite(x = tail_res_table, file = "tables/tables_new/results_cis_enrichment_in_tail.csv")

cis_tail <-
    tail_res |>
    rename(event = term) |>
    inner_join(has_cis |> select(order2, event, phenotype, trait), by = join_by("event")) |>
    arrange(order2) |>
    mutate(sig = if_else(p.value < 0.05, 1, 0),
           sig = factor(sig, levels = c(0, 1), labels = c("P ≥ 0.05", "P < 0.05")),
           conf.low = estimate - 1.96 * std.error,
           conf.high = estimate + 1.96 * std.error) |>
    ggplot(aes(x = order2, y = estimate)) +
    geom_rect(data = background_df, 
              aes(xmin = xmin, xmax = xmax, 
                  ymin = ymin, ymax = ymax), 
              alpha = 0.5, 
              inherit.aes = F,
              fill = background_df$fill,
              show.legend = F) +
    geom_crossbar(aes(ymin = conf.low, 
                      ymax = conf.high, 
                      fill = trait, 
                      linetype = trait), 
                  middle.linetype = "blank",
                  middle.color = "black", 
                  middle.linewidth = 0.5,
                  width = 0.6, 
                  alpha = 0.6) +
    geom_point(aes(shape = sig), size = 2) +
    geom_hline(yintercept = 0, lty = 2) +
    scale_y_continuous(labels = scales::number, breaks = scales::pretty_breaks(n = 7)) +
    scale_x_discrete(labels = has_cis$phenotype) +
    scale_fill_manual(values = c("Disease" = "#1E88E5", "Risk factor" = "orange"), drop = T) +
    scale_shape_manual(values = c("P ≥ 0.05" = 1, "P < 0.05" = 16)) +
    scale_linetype_manual(values = c("Disease" = 1, "Risk factor" = other_lty)) +
    coord_flip() +
    guides(
        fill = guide_legend(nrow = 2, order = 1,
                            title = "",
                            theme = theme(legend.key.spacing.y = unit(0, "cm"),
                                          legend.margin = margin())),
        shape = guide_legend(nrow = 2, title = "", order = 2,
                             theme = theme(legend.key.spacing.y = unit(0, "cm"),
                                           legend.margin = margin())),
        linetype = guide_legend(nrow = 2, title = "", order = 1,
                                theme = theme(legend.key.spacing.y = unit(0, "cm"),
                                              legend.margin = margin()))
    ) +
    labs(x = "", 
         y = "Log-odds ratio for cis-instrument presence",
         linetype = "") +
    theme_bw(base_size = 12) +
    theme(axis.title = element_text(face = "bold"),
          legend.position = "bottom",
          legend.box = "horizontal")


layout <- "
AABB
CCC#
"
subpan_1 <- (sig_pan | cis_pan | free(cis_tail)) +
    plot_layout(design = layout) +
    plot_annotation(tag_levels = "A", 
                    tag_prefix = "Fig. ") &
    theme(plot.tag.position = c(0, 1),
          plot.tag = element_text(size = 12, hjust = 0, vjust = 0))

export_image(plot = subpan_1, 
             fig_name = "figure2_obs_and_cis", 
             width = 15, height = 12, dpi = 300)