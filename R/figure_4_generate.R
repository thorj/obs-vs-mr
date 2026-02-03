library(tidyverse)
library(rstan)
library(patchwork)
library(Cairo)
box::use(data.table[fread, fwrite])
source("R/utils.R")

### CONST
other_lty <- 0

### Load data
f <- fread("data/phenotype_overview.csv")
obs <- fread("data/observational_results_full.csv")

forward <-
    obs |>
    filter(obs_sig.mr == 1) |>
    mutate(sig = as.integer(sign(beta) == sign(beta.mr))) |>
    group_by(event) |>
    summarize(n = n(), a = sum(sig), r = a/n) |>
    mutate(rftext = paste0(a, " / ", n))

fmissevent <- base::setdiff(obs$event, forward$event)

forward_missing <-
    obs |>
    filter(event %in% fmissevent) |>
    filter(has_cis == 1, p.bon < 0.05) |>
    group_by(event) |>
    mutate(mrfdr = p.adjust(pval.mr, "fdr", n())) |>
    summarize(n = n(), a = sum(mrfdr < 0.05)) |>
    mutate(r = 0,
           rftext = "No significant forward hits")

forward <- rbind(forward, forward_missing)


reverse <-
    obs |>
    filter(p.bon < 0.05) |>
    group_by(event) |>
    mutate(fdr.rev = p.adjust(pval.rev.mr, "fdr", n())) |>
    filter(fdr.rev < 0.05) |>
    mutate(sig = as.integer(sign(beta) == sign(beta.rev.mr))) |>
    summarize(n = n(), a = sum(sig), r = a/n) |>
    mutate(rrtext = paste0(a, " / ", n))

rmissevent <- base::setdiff(obs$event, reverse$event)

reverse_missing <-
    obs |>
    filter(event %in% rmissevent) |>
    filter(p.bon < 0.05) |>
    group_by(event) |>
    mutate(mrfdr = p.adjust(pval.rev.mr, "fdr", n())) |>
    summarize(n = n(), a = sum(mrfdr < 0.05)) |>
    mutate(r = 0,
           rrtext = "No significant reverse hits")

reverse <- rbind(reverse, reverse_missing)


agreement <- 
    full_join(forward |> select(event, rftext, rf = r), reverse |> select(event, rrtext, rr = r)) |>
    ungroup() |>
    mutate(rf = -rf) |>
    arrange(rf) |> 
    mutate(order = factor(row_number())) |>
    mutate(rf = replace_na(rf, 0),
           rr = replace_na(rr, 0)) |>
    full_join(f |> select(-order)) |>
    mutate(include = if_else(N == 0 & rr == 0, 0, 1))


observed_agreement <-
    agreement |>
    filter(include == 1) |>
    ggplot(aes(x = order)) +
    geom_col(aes(y = rf, fill = "Forward MR"), alpha = 0.6, color = "black") +
    geom_col(aes(y = rr, fill = "Reverse MR"), lty = other_lty, alpha = 0.6, color = "black") +
    geom_hline(yintercept = 0, lty = 1) +
    geom_hline(yintercept = -0.5, lty = 2, alpha = 0.7) +
    geom_hline(yintercept = 0.5, lty = 2, alpha = 0.7) +
    geom_label(aes(y = rf, label = rftext), size = 3, hjust = 1.1) +
    geom_label(aes(y = rr, label = rrtext), size = 3, hjust = -0.1) +
    labs(y = "Observed agreement ratio", x = "", fill = "") +
    scale_x_discrete(label = agreement$phenotype) +
    scale_y_continuous(labels = function(x) scales::percent(abs(x)),
                       breaks = -10:10/10) +
    scale_fill_manual(values = c("Forward MR" = "#882255", "Reverse MR" = "#332288"), drop = F) +
    coord_flip(ylim = c(-1.05, 1.15)) +
    theme_bw(base_size = 12) +
    theme(axis.title = element_text(face = "bold"),
          legend.position = "bottom",
          legend.box = "verticle",
          legend.margin = margin())

export_image(plot = observed_agreement, fig_name = "figure4_agreement_ratios", width = 12, height = 7, dpi = 300)

observed_agreement
