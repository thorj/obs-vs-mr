library(tidyverse)
library(patchwork)
box::use(data.table[fread, fwrite])
source("R/utils.R")

## Constants
other_lty <- 0
base_theme <-
    theme_bw(base_size = 12) +
    theme(axis.title = element_text(face = "bold"),
          legend.position = "bottom",
          legend.box = "vertical",
          legend.margin = margin())

mr_fill <- c("Forward MR" = "#882255", "Reverse MR" = "#332288")

## ============================================================
## Load data
## ============================================================

f <- fread(data_paths$phenotype_overview)
obs <- fread(data_paths$observational_full)
sig_df <- fread(data_paths$observational_significant)

primary_mr_summary <- readRDS(data_paths$primary_mr_summary)
primary_mr_significant <- primary_mr_summary$primary_mr_significant

secondary_mr_summary <- readRDS(data_paths$secondary_mr_summary)
secondary_mr_significant <- secondary_mr_summary$secondary_mr_significant

cis_info <- readRDS(data_paths$cis_info)
prop_w_cis <- cis_info$overall_cis_prop
cis_df <- cis_info$bonferroni_with_cis
has_cis <- cis_info$has_cis

## ============================================================
## Helpers
## ============================================================

make_mr_barplot <- function(df, ylims, y_breaks) {
    ggplot(df, aes(x = order)) +
        geom_col(aes(y = rf, fill = "Forward MR"), alpha = 0.6, color = "black") +
        geom_col(aes(y = rr, fill = "Reverse MR"),
                 linetype = other_lty, alpha = 0.6, color = "black") +
        geom_hline(yintercept = 0, lty = 1) +
        geom_label(aes(y = rf, label = rftext), size = 3, hjust = 1.1) +
        geom_label(aes(y = rr, label = rrtext), size = 3, hjust = -0.1) +
        scale_x_discrete(labels = df$phenotype) +
        scale_y_continuous(
            labels = function(x) scales::percent(abs(x)),
            breaks = y_breaks
        ) +
        scale_fill_manual(values = mr_fill, drop = FALSE, name = "") +
        coord_flip(ylim = ylims) +
        labs(x = "", y = "Proportion of MR significant SOMAmers") +
        base_theme
}

crossbar_guides <- function() {
    guides(
        fill = guide_legend(
            nrow = 2, order = 1, title = "",
            theme = theme(legend.key.spacing.y = unit(0, "cm"),
                          legend.margin = margin())
        ),
        shape = guide_legend(
            nrow = 2, order = 2, title = "",
            theme = theme(legend.key.spacing.y = unit(0, "cm"),
                          legend.margin = margin())
        ),
        linetype = guide_legend(
            nrow = 2, order = 1, title = "",
            theme = theme(legend.key.spacing.y = unit(0, "cm"),
                          legend.margin = margin())
        )
    )
}

## ============================================================
## TOP PANELS: MR significance (primary + secondary)
## ============================================================

primary_mr_sig <-
    make_mr_barplot(
        df       = primary_mr_significant,
        ylims    = c(-0.3, 0.55),
        y_breaks = (-10:10) / 10
    )

secondary_mr_sig <-
    make_mr_barplot(
        df       = secondary_mr_significant,
        ylims    = c(-0.17, 0.25),
        y_breaks = (-10:10) / 20
    )

## ============================================================
## LOWER PANELS: enrichment of MR significance in observational tail
## ============================================================

tester <- 
    obs |>
    filter(has_cis == 1) |>
    group_by(event)|>
    mutate(mrfdr = p.adjust(pval.mr, "fdr", n()),
           sec_analysis = as.integer(mrfdr < 0.05),
           obs_sig.mr2 = if_else(is.na(obs_sig.mr), 0, obs_sig.mr),
           mr_significant = if_else(obs_sig.mr2 + sec_analysis > 0, 1, 0)) |>
    ungroup()

mrsig_res <-
    lapply(tester$event |> unique(),
           FUN = function(x) {
               df <- tester |> 
                   filter(event == {{ x }}) |> 
                   mutate(intail = as.integer(pval < 0.05/7288))
               m1 <- 
                   glm(mr_significant ~ intail, data = df, family = binomial())
               res <- 
                   broom::tidy(m1) |>
                   filter(term == "intail") |> 
                   mutate(term = {{ x }})
               return(res)
           }) |> 
    data.table::rbindlist()

fisher_test_on_bad <-
    tester |>
    filter(event %in% mrsig_res$term[mrsig_res$std.error > 5]) |>
    mutate(intail = as.integer(p.bon < 0.5)) |>
    group_by(event) |>
    group_nest() |>
    mutate(tabx = map(data, ~table(.x$mr_significant, .x$intail)),
           fisher = map(tabx, fisher.test),
           fisher_p = map_dbl(fisher, "p.value"),
           fisher_est = log(map_dbl(fisher, "estimate"))) |>
    select(event, fisher_p, fisher_est)

fisher_anno <-
    mrsig_res |>
    filter(std.error > 5) |>
    mutate(estimate2 = -Inf,
           event = term) |>
    select(term, estimate2, event)  |>
    inner_join(fisher_test_on_bad |> select(-fisher_est)) |>
    inner_join(primary_mr_significant |> 
                   select(order, event, phenotype, trait), 
               by = join_by("event")) |>
    mutate(pp = round(fisher_p, 2),
           label = paste0("FET P-value: ", pp)) |>
    select(event, order, trait, label, estimate2)

mrsig_df <-
    mrsig_res |>
    mutate(estimate2 = if_else(std.error > 5, NA, estimate),
           conf.low = estimate2 - 1.96 * std.error,
           conf.high = estimate2 + 1.96 * std.error) |>
    rename(event = term) |>
    inner_join(primary_mr_significant |>
                   select(order, event, phenotype, trait), 
               by = join_by("event")) |>
    arrange(order) |>
    mutate(sig = if_else(p.value < 0.05, 1, 0),
           sig = factor(sig, 
                        levels = c(0, 1), 
                        labels = c("P ≥ 0.05", "P < 0.05")))


background_df <- data.frame(
    xmin = as.numeric(levels(has_cis$order2)) - 0.5,
    xmax = as.numeric(levels(has_cis$order2)) + 0.5,
    ymin = -Inf,
    ymax = Inf,
    fill = rep(c("grey90", "white"), length.out = length(levels(has_cis$order)))
)

mrsig_plot <- 
    mrsig_df |>
    ggplot(aes(x = order, y = estimate2)) +
    geom_rect(data = background_df, 
              aes(xmin = xmin, xmax = xmax, 
                  ymin = ymin, ymax = ymax,), 
              fill = background_df$fill,
              alpha = 0.5, 
              inherit.aes = F, 
              show.legend = F) +
    geom_crossbar(aes(ymin = conf.low, 
                      ymax = conf.high, 
                      fill = trait,
                      linetype = trait), 
                  middle.linetype = "blank",
                  width = 0.6, 
                  alpha = 0.6) +
    geom_point(aes(shape = sig), size = 2) +
    geom_hline(yintercept = 0, lty = 2) +
    geom_label(data = fisher_anno, aes(label = label),
               hjust = -0.1, size = 2.5, alpha = 0.3, linewidth = 0.3) +
    scale_y_continuous(labels = scales::number, breaks = scales::pretty_breaks(n = 7)) +
    scale_x_discrete(labels = mrsig_df$phenotype) +
    scale_fill_manual(values = c("Disease" = "#1E88E5", "Risk factor" = "orange"), drop = T) +
    scale_shape_manual(values = c("P ≥ 0.05" = 1, "P < 0.05" = 16)) +
    scale_linetype_manual(values = c("Disease" = 1, "Risk factor" = other_lty)) + 
    coord_flip() +
    crossbar_guides() +
    labs(x = "", 
         y = "Log-odds ratio", 
         linetype = "") +
    base_theme

### REVERSE

rev_mr_obs <-
    obs |>
    group_by(event) |>
    filter(p.bon < 0.05) |>
    mutate(revfdr = p.adjust(pval.rev.mr, "fdr", n())) |>
    filter(revfdr < 0.05) |>
    mutate(prim_rev_sig = 1) |>
    ungroup() |>
    select(event, somamer, prim_rev_sig)

tester2 <-
    obs |>
    left_join(rev_mr_obs, by = join_by("event", "somamer")) |>
    mutate(prim_rev_sig2 = if_else(is.na(prim_rev_sig), 0, prim_rev_sig)) |>
    group_by(event) |>
    mutate(mrfdr = p.adjust(pval.rev.mr, "fdr", n()),
           sec_analysis = as.integer(mrfdr < 0.05),
           mr_significant = if_else(prim_rev_sig2  + sec_analysis > 0, 1, 0)) |>
    ungroup()


revmrsig_res <-
    lapply(tester2$event |> unique(),
           FUN = function(x) {
               df <- tester2 |> 
                   filter(event == {{ x }}) |> 
                   mutate(intail = as.integer(pval < 0.05/7288))
               wmsg <- NULL
               m1 <- withCallingHandlers(
                   glm(mr_significant ~ intail, data = df, family = binomial()),
                   warning = function(w) {
                       wmsg <<- conditionMessage(w)
                       invokeRestart("muffleWarning")
                   }
               )
               if (!is.null(wmsg)) {
                   message("Warning for outcome ", x, ": ", wmsg)
               }
               res <- 
                   broom::tidy(m1) |>
                   filter(term == "intail") |> 
                   mutate(term = {{ x }},
                          wmsg = ifelse(is.null(wmsg), 0, 1))
               return(res)
           }) |> 
    data.table::rbindlist()

revmrsig_res 

mrsig_plot2_df <- 
    revmrsig_res |>
    mutate(include = if_else(wmsg == 1 | abs(estimate) > 10, 0, 1),
           estimate2 = if_else(include == 0, NA, estimate), 
           conf.low = estimate2 - 1.96 * std.error,
           conf.high = estimate2 + 1.96 * std.error) |>
    rename(event = term) |>
    inner_join(primary_mr_significant |> select(order, event, phenotype, trait), by = join_by("event")) |>
    arrange(order) |>
    mutate(sig = if_else(p.value < 0.05, 1, 0),
           sig = factor(sig, levels = c(0, 1), labels = c("P ≥ 0.05", "P < 0.05")))

## Eight phenotypes are annoying need to fix explicitly 

tester2 |>
    filter(event %in% mrsig_plot2_df$event[mrsig_plot2_df$include == 0]) |>
    mutate(intail = as.integer(p.bon < 0.5)) |>
    group_by(event) |>
    summarize(a = sum(mr_significant))

tmptmp <- 
    tester2 |>
    filter(event %in% mrsig_plot2_df$event[mrsig_plot2_df$include == 0]) |>
    mutate(intail = as.integer(p.bon < 0.5)) |>
    group_by(event) |>
    group_nest() |>
    mutate(tabx = map(data, ~table(.x$mr_significant, .x$intail)),
           totalsigs = map_dbl(data, ~sum(.x$mr_significant)))

fisher_secanal <-
    tmptmp |>
    filter(totalsigs > 0) |>
    mutate(fisher = map(tabx, fisher.test),
           fisher_p = map_dbl(fisher, "p.value"),
           fisher_estimate = map_dbl(fisher, "estimate")) |>
    select(event, fisher_p, fisher_estimate)

mrsig_plot2_df |>
    filter(include == 0) |>
    mutate(note = 
               case_when(
                   event %in% fisher_secanal$event ~ "Tested with Fisher's exact test",
                   TRUE ~ "No significant hits"
               ))


fisher_secanno <-
    mrsig_plot2_df |>
    filter(include == 0) |>
    mutate(estimate2 = -Inf) |>
    select(estimate2, event)  |>
    left_join(fisher_secanal) |>
    inner_join(primary_mr_significant |> select(order, event, phenotype, trait), by = join_by("event")) |>
    mutate(pp = round(fisher_p, 2),
           label = if_else(is.na(fisher_p),
                           "No significant hits",
                           paste0("FET P-value: ", pp))) |>
    select(event, order, trait, label, estimate2)

mrsig_2 <-
    mrsig_plot2_df |>
    ggplot(aes(x = order, y = estimate2)) +
    geom_rect(data = background_df,
              aes(xmin = xmin, xmax = xmax,
                  ymin = ymin, ymax = ymax,),
              fill = background_df$fill,
              alpha = 0.5,
              inherit.aes = F,
              show.legend = F) +
    geom_crossbar(aes(ymin = conf.low, 
                      ymax = conf.high, 
                      fill = trait,
                      linetype = trait), 
                  middle.linetype = "blank",
                  width = 0.6, 
                  alpha = 0.6) +
    geom_point(aes(shape = sig), size = 2) +
    geom_hline(yintercept = 0, lty = 2) +
    geom_label(data = fisher_secanno, aes(label = label),
               hjust = -0.1, size = 2.5, alpha = 0.3, linewidth = 0.3) +
    scale_y_continuous(labels = scales::number, breaks = scales::pretty_breaks(n = 7)) +
    scale_x_discrete(labels = mrsig_plot2_df$phenotype) +
    scale_fill_manual(values = c("Disease" = "#1E88E5", "Risk factor" = "orange"), drop = T) +
    scale_shape_manual(values = c("P ≥ 0.05" = 1, "P < 0.05" = 16)) +
    scale_linetype_manual(values = c("Disease" = 1, "Risk factor" = other_lty)) + 
    coord_flip() +
    guides(
        fill = guide_legend(nrow = 2, order = 1,
                            title = "",
                            theme = theme(legend.key.spacing.y = unit(0, "cm"),
                                          legend.margin = margin())),
        shape = guide_legend(nrow = 2, 
                             title = "", 
                             order = 2,
                             theme = theme(legend.key.spacing.y = unit(0, "cm"),
                                           legend.margin = margin())),
        linetype = guide_legend(nrow = 2, 
                                title = "", 
                                order = 1,
                                theme = theme(
                                    legend.key.spacing.y = unit(0,"cm"),
                                    legend.margin = margin()))
        ) +
    labs(x = "", 
         y = "Log-odds ratio", 
         linetype = "") +
    theme_bw(base_size = 12) +
    theme(axis.title = element_text(face = "bold"),
          legend.position = "bottom",
          legend.box = "horizontal")

subpan_3 <- 
    (fdr_sig | fdr_sec_sig) / (mrsig_plot | mrsig_2) + 
    plot_annotation(tag_levels = "A", tag_prefix = "Fig. ") &
    theme(plot.tag.position = c(0, 1),
          plot.tag = element_text(size = 12, hjust = 0, vjust = 0))

export_image(plot = subpan_3,
             fig_name = "figure3_mr_res_all", 
             width = 18, height = 12, dpi = 300)