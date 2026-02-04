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
          legend.box = "horizontal",
          legend.margin = margin())

mr_fill <- c("Forward MR" = "#882255", "Reverse MR" = "#332288")
crossbar_fill <- c("Disease" = "#1E88E5", "Risk factor" = "orange")

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
has_cis <- cis_info$has_cis

background_df <- data.frame(
    xmin = as.numeric(levels(has_cis$order2)) - 0.5,
    xmax = as.numeric(levels(has_cis$order2)) + 0.5,
    ymin = -Inf,
    ymax = Inf,
    fill = rep(c("grey90", "white"), length.out = length(levels(has_cis$order)))
)

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

make_crossbar_plot <- function(plot_df, background_df) {
    plot_df |>
        ggplot(aes(x = order, y = estimate2)) +
        geom_rect(data = background_df, 
                  aes(xmin = xmin, xmax = xmax, 
                      ymin = ymin, ymax = ymax,), 
                  fill = background_df$fill,
                  alpha = 0.5, 
                  inherit.aes = F, 
                  show.legend = F
        ) +
        geom_crossbar(aes(ymin = conf.low, 
                          ymax = conf.high, 
                          fill = trait,
                          linetype = trait), 
                      middle.linetype = "blank",
                      width = 0.6, 
                      alpha = 0.6) +
        geom_point(aes(shape = sig), size = 2) +
        geom_hline(yintercept = 0, lty = 2) +
        geom_label(aes(label = label),
                   hjust = -0.1, 
                   size = 2.5, 
                   alpha = 0.3, 
                   linewidth = 0.3) +
        scale_y_continuous(labels = scales::number, 
                           breaks = scales::pretty_breaks(n = 7)) +
        scale_x_discrete(labels = plot_df$phenotype) +
        scale_fill_manual(values = crossbar_fill, drop = T) +
        scale_shape_manual(values = c("P ≥ 0.05" = 1, "P < 0.05" = 16)) +
        scale_linetype_manual(values = c("Disease" = 1, "Risk factor" = other_lty)) + 
        coord_flip() +
        crossbar_guides() +
        labs(x = "", 
             y = "Log-odds ratio", 
             linetype = "") +
        base_theme 
}

perform_enrichment_analysis <- function(df, order_df) {
    nested_df <- construct_enrichment_df(df)
    events <- nested_df$event
    lapply(X = events, enrichment_analysis, df = nested_df) |>
        data.table::rbindlist() |>
        enrichment_annotation(order_df = order_df)
}

construct_enrichment_df <- function(df) {
    df |> 
        group_by(event) |>
        group_nest() |>
        mutate(tabx = map(data, ~table(.x$mr_significant, .x$intail)),
               totalsigs = map_dbl(data, ~sum(.x$mr_significant)),
               any_bad_sep = map_dbl(tabx, check_table_for_bad_sep),
               include = if_else(totalsigs == 0 | any_bad_sep > 0, 0, 1))
}

## Bad events have either zero significant proteins OR zero cells in frequency table

check_table_for_bad_sep <- function(tbl) {
    sapply(X = tbl, FUN = function(x) { x == 0}) |> sum()
}

enrichment_analysis <- function(df, event) {
    event_df <- df |> filter(event == {{ event }})
    if (event_df$totalsigs[1] == 0) {
        enrichment_res <- no_significant_fit(event_df)
    }
    else if (event_df$include[1] == 0) {
        enrichment_res <- fisher_enrichment_fit(event_df)
    } else {
        event_df <- event_df |> select(data) |> unnest(cols = c(data))
        enrichment_res <- glm_enrichment_fit(df = event_df, event = {{ event }})
    }
    return(enrichment_res)
}

no_significant_fit <- function(df) {
    tibble(
        event = df$event[1],
        estimate = NA,
        std.error = NA,
        statistic = NA,
        p.value = NA,
        method = "no_fit"
    )
}

glm_enrichment_fit <- function(df, independent_var = "intail", response_var = "mr_significant", event) {
    model_formula <- reformulate(termlabels = independent_var, response = response_var)
    model <- glm(model_formula, data = df, family = binomial())
    broom::tidy(model) |>
        filter(term == {{ independent_var }}) |> 
        mutate(term = {{ event }},
               method = "glm") |>
        rename(event = term)
}

fisher_enrichment_fit <- function(df) {
    df |>
        mutate(
            fisher = map(tabx, fisher.test),
            p.value = map_dbl(fisher, "p.value"),
            estimate = log(map_dbl(fisher, "estimate")),
            std.error = NA,
            statistic = NA,
            method = "fisher"
        ) |>
        select (
            event,
            estimate,
            std.error,
            statistic,
            p.value,
            method
        )
}

enrichment_annotation <- function(df, order_df) {
    df |>
        mutate(estimate2 = 
                   case_when(
                       method == "fisher" ~ -Inf, 
                       method == "no_fit" ~ -Inf,
                       method == "glm" ~ estimate,
                       .default = NA),
               conf.low = estimate2 - 1.96 * std.error,
               conf.high = estimate2 + 1.96 * std.error,
               label = 
                   case_when(
                       method == "fisher" ~ paste0("FET P-value: ", round(p.value, 2)),
                       method == "no_fit" ~ "No significant hits",
                       method == "glm" ~ "",
                       .default = "ERROR"),
               sig = if_else(p.value < 0.05, 1, 0),
               sig = factor(sig, 
                            levels = c(0, 1), 
                            labels = c("P ≥ 0.05", "P < 0.05"))
        ) |>
        inner_join(order_df |>
                       select(order, event, phenotype, trait), 
                   by = join_by("event")) |>
        arrange(order) 
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

forward_df <- 
    obs |>
    filter(has_cis == 1) |>
    group_by(event)|>
    mutate(mrfdr = p.adjust(pval.mr, "fdr", n()),
           sec_analysis = as.integer(mrfdr < 0.05),
           obs_sig.mr2 = if_else(is.na(obs_sig.mr), 0, obs_sig.mr),
           mr_significant = if_else(obs_sig.mr2 + sec_analysis > 0, 1, 0),
           intail = as.integer(p.bon < 0.05))

forward_enrich <-
    perform_enrichment_analysis(df = forward_df, 
                                order_df = primary_mr_significant)


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

reverse_df <-
    obs |>
    left_join(rev_mr_obs, by = join_by("event", "somamer")) |>
    mutate(prim_rev_sig2 = if_else(is.na(prim_rev_sig), 0, prim_rev_sig)) |>
    group_by(event) |>
    mutate(mrfdr = p.adjust(pval.rev.mr, "fdr", n()),
           sec_analysis = as.integer(mrfdr < 0.05),
           mr_significant = if_else(prim_rev_sig2  + sec_analysis > 0, 1, 0),
           intail = as.integer(p.bon < 0.05))

reverse_enrich <-
    perform_enrichment_analysis(df = reverse_df, 
                                order_df = primary_mr_significant)

## ============================================================
## LOWER PANELS: enrichment of MR significance in observational tail
## ============================================================

mrsig_plot <- 
    make_crossbar_plot(plot_df = forward_enrich, 
                       background_df = background_df)
mrsig_2 <-
    make_crossbar_plot(plot_df = reverse_enrich, 
                       background_df = background_df)

## ============================================================
## Combine + export
## ============================================================

subpan_3 <- 
    (fdr_sig | fdr_sec_sig) / (mrsig_plot | mrsig_2) + 
    plot_annotation(tag_levels = "A", tag_prefix = "Fig. ") &
    theme(plot.tag.position = c(0, 1),
          plot.tag = element_text(size = 12, hjust = 0, vjust = 0))

export_image(plot = subpan_3,
             fig_name = "figure3_mr_res_all", 
             width = 18, height = 12, dpi = 300)