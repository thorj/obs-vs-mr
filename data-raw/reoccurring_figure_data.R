library(tidyverse)
box::use(data.table[fread, fwrite])

### Load data
f <- fread("data/phenotype_overview.csv")
obs <- fread("data/observational_results_full.csv")

### CONSTANTS
n_somamers <- unique(obs$somamer) |> length()

## Observational significance
sig_df <-
    obs |>
    filter(p.bon < 0.05) |>
    group_by(event) |>
    summarize(n = n()) |>
    mutate(p = n/n_somamers) |>
    inner_join(f) |>
    arrange(p) |>
    ungroup() |>
    mutate(order2 = row_number(),
           order2 = factor(order2))

fwrite(x = sig_df, file = "data/reoccurring_data/bonferroni_sig_obs.csv")

## Figure 2
## CIS 
## How many SOMAmers have a cis-instrument?
prop_w_cis <- 
    obs |> 
    filter(has_cis == 1) |> 
    pull(somamer) |>
    unique() |> 
    length() / n_somamers

cis_df <-
    obs |>
    filter(p.bon < 0.05) |>
    group_by(event) |>
    summarize(n = n(), 
              ciss = sum(has_cis)) |>
    mutate(p = ciss/n) |>
    inner_join(f, by = join_by(event)) |>
    ungroup() |>
    inner_join(sig_df |> select(event, order2), by = join_by(event)) |>
    arrange(order2)

## This data can probably be removed but we won't do that just eyt
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
           ne = floor(prop_w_cis * n)) |>
    mutate(label = paste0(nc, "/", n)) |>
    inner_join(cis_df |> select(event, phenotype, order2, trait),
               by = join_by(event)) |>
    arrange(order2)

saveRDS(object = list(overall_cis_prop = prop_w_cis,
                      bonferroni_with_cis = cis_df,
                      has_cis = has_cis), 
        file = "data/reoccurring_data/cis_protein_data.rds")

### total fdr significant <000 you are here (figure 2 and 3)
## How many bonferroni significant proteins are significant in forward MR 
primary_forward_mr_significant <-
    obs |>
    filter(p.bon < 0.05, has_cis == 1) |>
    group_by(event) |>
    mutate(fdr = p.adjust(pval.mr, "fdr", n())) |>
    summarize(n = n(), a = sum(fdr < 0.05), r = a/n) |>
    mutate(rftext = paste0(a, " / ", n))

## Agreement ratios for primary forward MR 
primary_mr_forward_agreement <- 
    obs |>
    filter(obs_sig.mr == 1) |>
    mutate(sig = as.integer(sign(beta) == sign(beta.mr))) |>
    group_by(event) |>
    summarize(n = n(), a = sum(sig), r = a/n) |>
    mutate(rftext = paste0(a, " / ", n)) |>
    rbind(
        primary_forward_mr_significant |>
            filter(a == 0) |>
            mutate(rftext = "No significant forward hits")
        )

## How many bonferroni significant proteins are significant in reverse MR
primary_reverse_mr_significant <-
    obs |>
    filter(p.bon < 0.05) |>
    group_by(event) |>
    mutate(fdr.rev = p.adjust(pval.rev.mr, "fdr", n())) |>
    summarize(n = n(), a = sum(fdr.rev < 0.05), r = a/n) |>
    mutate(rrtext = paste0(a, " / ", n))

## How many bonferroni & MR significant proteins are in agreement
primary_reverse_mr_agreement <-
    obs |>
    filter(p.bon < 0.05) |>
    group_by(event) |>
    mutate(fdr.rev = p.adjust(pval.rev.mr, "fdr", n())) |>
    filter(fdr.rev < 0.05) |>
    mutate(sig = as.integer(sign(beta) == sign(beta.rev.mr))) |>
    summarize(n = n(), a = sum(sig), r = a/n) |>
    mutate(rrtext = paste0(a, " / ", n)) |>
    rbind(
        primary_reverse_mr_significant |>
            filter(a == 0) |>
            mutate(rrtext = "No significant reverse hits")
    )

## Total bonferroni and MR significant proteins in primary analysis
primary_mr_significant <- 
    full_join(primary_forward_mr_significant |> select(event, rftext, rf = r), 
              primary_reverse_mr_significant |> select(event, rrtext, rr = r),
              by = join_by(event)) |>
    mutate(rf = -rf) |>
    arrange(rf) |>
    mutate(order = factor(row_number())) |>
    mutate(rf = replace_na(rf, 0),
           rr = replace_na(rr, 0)) |>
    full_join(f |> select(-order, -N, -Y),
              by = join_by(event))

## Total agreements 
primary_mr_agreement <- 
    full_join(primary_mr_forward_agreement |> select(event, rftext, rf = r),
              primary_reverse_mr_agreement |> select(event, rrtext, rr = r),
              by = join_by(event)) |>
    ungroup() |>
    mutate(rf = -rf) |>
    arrange(rf) |> 
    mutate(order = factor(row_number())) |>
    mutate(rf = replace_na(rf, 0),
           rr = replace_na(rr, 0)) |>
    full_join(f |> select(-order),
              by = join_by(event)) |>
    mutate(include = if_else(N == 0 & rr == 0, 0, 1))

saveRDS(object = 
            list(primary_mr_forward_significant = primary_forward_mr_significant ,
                 primary_mr_reverse_significant = primary_reverse_mr_significant,
                 primary_mr_significant = primary_mr_significant,
                 primary_mr_forward_agreement = primary_mr_forward_agreement,
                 primary_mr_reverse_agreement = primary_reverse_mr_agreement,
                 primary_mr_agreement = primary_mr_agreement),
        file = "data/reoccurring_data/primary_mr_summary.rds")

## Secondary MR significance
secondary_forward_mr_base <-
    obs |>
    filter(has_cis == 1) |>
    group_by(event) |>
    mutate(fdrnew = p.adjust(pval.mr, "fdr", n()),
           sig = as.integer(sign(beta) == sign(beta.mr)))

secondary_forward_mr_significant <-
    secondary_forward_mr_base |>
    summarize(n = n(), a = sum(fdrnew < 0.05), r = a/n) |>
    mutate(rftext = paste0(a, " / ", n)) 

secondary_forward_mr_agreement <- 
    secondary_forward_mr_base |>
    filter(fdrnew < 0.05) |>
    summarize(n = n(), a = sum(sig), r = a/n) |>
    mutate(rftext = paste0(a, " / ", n)) |>
    rbind(
        secondary_forward_mr_significant |>
            filter(a == 0) |>
            mutate(rftext = "No significant forward hits")
    )
### reverse
secondary_reverse_mr_base <-
    obs |>
    group_by(event) |>
    mutate(fdrnew = p.adjust(pval.rev.mr, "fdr", n()),
           sig = as.integer(sign(beta) == sign(beta.rev.mr)))

secondary_reverse_mr_significant <-
    secondary_reverse_mr_base |>
    summarize(n = n(), a = sum(fdrnew < 0.05), r = a/n) |>
    mutate(rrtext = paste0(a, " / ", n)) 

secondary_reverse_mr_agreement <-
    secondary_reverse_mr_base |>
    filter(fdrnew < 0.05) |>
    summarize(n = n(), a = sum(sig), r = a/n) |>
    mutate(rrtext = paste0(a, " / ", n)) |>
    rbind(
        secondary_reverse_mr_significant |>
            filter(a == 0) |>
            mutate(rrtext = "No significant forward hits")
    )

secondary_mr_significant <- 
    full_join(secondary_forward_mr_significant |> 
                  select(event, rftext, rf = r), 
              secondary_reverse_mr_significant |> 
                  select(event, rrtext, rr = r),
              by = join_by(event)) |>
    mutate(rf = -rf) |>
    inner_join(primary_mr_significant |> select(event, order), by = join_by(event)) |>
    mutate(rf = replace_na(rf, 0),
           rr = replace_na(rr, 0)) |>
    full_join(f |> select(-order, -N, -Y), by = join_by(event)) |>
    arrange(order)

secondary_mr_agreement <- 
    full_join(secondary_forward_mr_agreement |> select(event, rftext, rf = r),
              secondary_reverse_mr_agreement |> select(event, rrtext, rr = r),
              by = join_by(event)) |>
    ungroup() |>
    mutate(rf = -rf) |>
    arrange(rf) |> 
    mutate(order = factor(row_number())) |>
    mutate(rf = replace_na(rf, 0),
           rr = replace_na(rr, 0)) |>
    full_join(f |> select(-order,),
              by = join_by(event)) |>
    mutate(include = if_else(N == 0 & rr == 0, 0, 1))

fread("data/00b_no_obs_res.csv") |> filter(event == "a2dm2")

saveRDS(object = 
            list(secondary_mr_forward_significant = secondary_forward_mr_significant ,
                 secondary_mr_reverse_significant = secondary_reverse_mr_significant,
                 secondary_mr_significant = secondary_mr_significant), 
        file = "data/reoccurring_data/secondary_mr_significant.rds")
