library(tidyverse)
library(patchwork)
source("R/utils.R")
box::use(data.table[fread, fwrite])

### Load data
f <- fread("data/phenotype_overview.csv")
obs <- fread("data/observational_results_full.csv")

### CONST
other_lty <- 0
####

sig_df <-
    obs |>
    filter(p.bon < 0.05) |>
    group_by(event) |>
    summarize(n = n()) |>
    mutate(p = n/7288) |>
    inner_join(f) |>
    arrange(p) |>
    ungroup() |>
    mutate(order2 = row_number(),
           order2 = factor(order2))

sig_pan <- 
    sig_df |>
    ggplot(aes(x = order2, y = n)) +
    geom_col(aes(fill = trait, linetype = trait), color = "black", alpha = 0.6) + 
    scale_y_continuous(breaks = scales::extended_breaks(n = 10), limits = c(0, 4500)) +
    scale_x_discrete(labels = sig_df$phenotype) +
    #scale_fill_manual(values = c("Disease" = "#882255", "Risk factor" = "#332288")) +
    scale_fill_manual(values = c("Disease" = "#1E88E5", "Risk factor" = "orange"), drop = T) +
    scale_linetype_manual(values = c("Disease" = 1, "Risk factor" = other_lty)) + ## <-linetype
    geom_label(aes(label = n), size = 3, hjust = -0.1) +
    coord_flip() +
    labs(x = "", y = "Number of phenotype-level significant SOMAmers", fill = "", linetype = "") +
    theme_bw(base_size = 12) +
    theme(axis.title = element_text(face = "bold"),
          legend.position = "bottom")

sig_pan
### cis
prop_w_cis <- 
    obs |> 
    filter(has_cis == 1) |> 
    pull(somamer) |>
    unique() |> 
    length() / 7288
cis_df <-
    obs |>
    filter(p.bon < 0.05) |>
    group_by(event) |>
    summarize(n = n(), ciss = sum(has_cis)) |>
    mutate(p = ciss/n) |>
    inner_join(f) |>
    ungroup() |>
    inner_join(sig_df |> select(event, order2)) |>
    arrange(order2)


cis_pan <- 
    cis_df |>
    mutate(ln = paste0(ciss, " / ", n)) |>
    ggplot(aes(x = order2, y = p)) +
    geom_hline(yintercept = prop_w_cis, lty = 2, alpha = 0.7) +
    geom_col(aes(fill = trait, linetype = trait), color = "black", alpha = 0.6) +
    scale_y_continuous(labels = scales::percent, breaks = 0:20/20, limits = c(0, 0.75)) +
    scale_x_discrete(labels = cis_df$phenotype) +
    scale_fill_manual(values = c("Disease" = "#1E88E5", "Risk factor" = "orange"), drop = T) +
    #scale_fill_manual(values = c("Disease" = "#882255", "Risk factor" = "#332288")) +
    scale_linetype_manual(values = c("Disease" = 1, "Risk factor" = other_lty)) +
    geom_label(aes(label = ln), size = 3, hjust = -0.1) +
    coord_flip() +
    labs(x = "", y = "Proportion of SOMAmers with cis-instruments", fill = "", linetype = "") +
    theme_bw(base_size = 12) +
    theme(axis.title = element_text(face = "bold"),
          legend.position = "bottom")

### total fdr significant
forwardsig <-
    obs |>
    filter(p.bon < 0.05, has_cis == 1) |>
    group_by(event) |>
    mutate(fdr = p.adjust(pval.mr, "fdr", n())) |>
    summarize(n = n(), a = sum(fdr < 0.05), r = a/n) |>
    mutate(rftext = paste0(a, " / ", n))

reversesig <-
    obs |>
    filter(p.bon < 0.05) |>
    group_by(event) |>
    mutate(fdr.rev = p.adjust(pval.rev.mr, "fdr", n())) |>
    summarize(n = n(), a = sum(fdr.rev < 0.05), r = a/n) |>
    mutate(rrtext = paste0(a, " / ", n))

agreementsig <- 
    full_join(forwardsig |> select(event, rftext, rf = r), reversesig |> select(event, rrtext, rr = r)) |>
    mutate(rf = -rf) |>
    #select(-c("rftext", "rrtext")) |>
    arrange(rf) |>
    mutate(order = factor(row_number())) |>
    mutate(rf = replace_na(rf, 0),
           rr = replace_na(rr, 0)) |>
    #filter(!(rr == 0 & rf == 0)) |>
    full_join(f |> select(-order, -N, -Y))


fdr_sig <-
    agreementsig |>
    # mutate(rrtext = ifelse(rr == 0, NA, rrtext),
    #        rftext = ifelse(rf == 0, NA, rftext)) |>
    ggplot(aes(x = order)) +
    geom_col(aes(y = rf, fill = "Forward MR"), alpha = 0.6, color = "black") +
    geom_col(aes(y = rr, fill = "Reverse MR"), lty = other_lty, alpha = 0.6, color = "black") +
    geom_hline(yintercept = 0, lty = 1) +
    geom_label(aes(y = rf, label = rftext), size = 3, hjust = 1.1) +
    geom_label(aes(y = rr, label = rrtext), size = 3, hjust = -0.1) +
    labs(y = "Proportion of MR significant SOMAmers", x = "", fill = "") +
    scale_x_discrete(label = agreementsig$phenotype) +
    scale_y_continuous(labels = function(x) scales::percent(abs(x)),
                       breaks = -10:10/10) +
    scale_fill_manual(values = c("Forward MR" = "#882255", "Reverse MR" = "#332288"), drop = F) +
    coord_flip(ylim = c(-0.3, 0.55)) +
    theme_bw(base_size = 12) +
    theme(axis.title = element_text(face = "bold"),
          legend.position = "bottom",
          legend.box = "verticle",
          legend.margin = margin())

## same plot but do it for secondary analysis
allcis_fdrsig <-
    obs |>
    filter(has_cis == 1) |>
    group_by(event) |>
    mutate(fdrnew = p.adjust(pval.mr, "fdr", n())) |>
    summarize(n = n(), a = sum(fdrnew < 0.05), r = a/n) |>
    mutate(rftext = paste0(a, " / ", n))

allrev_fdrsig <-
    obs |>
    group_by(event) |>
    mutate(fdrnew = p.adjust(pval.rev.mr, "fdr", n())) |>
    summarize(n = n(), a = sum(fdrnew < 0.05), r = a/n) |>
    mutate(rrtext = paste0(a, " / ", n)) 

agreementsecsig <- 
    full_join(allcis_fdrsig |> select(event, rftext, rf = r), allrev_fdrsig |> select(event, rrtext, rr = r)) |>
    mutate(rf = -rf) |>
    inner_join(agreementsig |> select(event, order)) |>
    mutate(rf = replace_na(rf, 0),
           rr = replace_na(rr, 0)) |>
    full_join(f |> select(-order, -N, -Y)) |>
    arrange(order)

fdr_sec_sig <-
    agreementsecsig |>
    # mutate(rrtext = ifelse(rr == 0, NA, rrtext),
    #        rftext = ifelse(rf == 0, NA, rftext)) |>
    ggplot(aes(x = order)) +
    geom_col(aes(y = rf, fill = "Forward MR"), alpha = 0.6, color = "black") +
    geom_col(aes(y = rr, fill = "Reverse MR"), lty = other_lty, alpha = 0.6, color = "black") +
    geom_hline(yintercept = 0, lty = 1) +
    geom_label(aes(y = rf, label = rftext), size = 3, hjust = 1.1) +
    geom_label(aes(y = rr, label = rrtext), size = 3, hjust = -0.1) +
    labs(y = "Proportion of MR significant SOMAmers", x = "", fill = "") +
    scale_x_discrete(label = agreementsecsig$phenotype) +
    #scale_x_discrete(label = NULL) +
    scale_y_continuous(labels = function(x) scales::percent(abs(x)),
                       breaks = -10:10/20) +
    scale_fill_manual(values = c("Forward MR" = "#882255", "Reverse MR" = "#332288"), drop = F) +
    coord_flip(ylim = c(-0.17, 0.25)) +
    theme_bw(base_size = 12) +
    theme(axis.title = element_text(face = "bold"),
          legend.position = "bottom",
          legend.box = "verticle",
          legend.margin = margin())



### adding 
has_cis <-
    obs |>
    filter(p.bon < 0.05) |>
    group_by(event) |>
    summarize(n = n(), 
              nc = sum(has_cis == 1),
              m = mean(has_cis),
              s = sd(has_cis)) |>
    ungroup() |>
    arrange(desc(m)) |>
    mutate(order = row_number(),
           order = factor(order),
           ne = floor(0.28 * n)) |>
    mutate(label = paste0(nc, "/", n)) |>
    inner_join(cis_df |> select(event, phenotype, order2, trait)) |>
    arrange(order2)

background_df <- data.frame(
    xmin = as.numeric(levels(has_cis$order2)) - 0.5,
    xmax = as.numeric(levels(has_cis$order2)) + 0.5,
    ymin = -Inf,
    ymax = Inf,
    fill = rep(c("grey90", "white"), length.out = length(levels(has_cis$order)))
)

#### Bonferroni tail enrichment
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
    #scale_fill_manual(values = c("Disease" = "#882255", "Risk factor" = "#332288"), drop = T) +
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

export_image(plot = subpan_1, fig_name = "figure2_obs_and_cis", width = 15, height = 12, dpi = 300)