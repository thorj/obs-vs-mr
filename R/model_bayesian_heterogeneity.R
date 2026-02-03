library(tidyverse)
library(rstan)
options(mc.cores = parallel::detectCores())
box::use(data.table[fread, fwrite])

### Load data
data_for_bayes <- readRDS(file = "data/data_for_analysis.rds")
f <- data_for_bayes$f
f2 <- data_for_bayes$f2
group <- data_for_bayes$map_phenotype_to_group
protein_map <- data_for_bayes$map_somamer_to_egs

obs <- fread("../../papers/02_obs_vs_mr/data/observational_results_full.csv")

## Onto coding and instruments
harmon <- fread("../../papers/02_obs_vs_mr/data/forward_mr_harmonized.csv")
harmon_snps <- harmon$snp_id |> unique() |> sort()

dbsnp <- readRDS("data-raw/view_snp_info.rds")
dbsnp <- 
    dbsnp |> 
    select(snp_id, 
           coding_variant) |> 
    filter(snp_id %in% harmon_snps) |>
    group_by(snp_id) |> 
    filter(row_number() == 1) |> 
    ungroup()

## hetero
get_dosages <- 
    harmon |>
    inner_join(h1 |> select(gwas, somamer)) |>
    pull(snp_id) |>
    unique() |> sort()

saveRDS(object = get_dosages, file = "data-raw/dosages_for_forward_mr_snp_ids.rds")

matched_ids <- data.table::fread("data-raw/old_id_to_new_id.csv") |>
    rename(new_db_snpid = snp_id,
           snp_id = old_snp_id,
           dosage_id = id)

wald_delta <- function(wr, beta.out, beta.exp, se.out, se.exp) {
    exposure_variance <- (se.exp/beta.exp)^2
    outcome_variance <- (se.out/beta.out)^2
    return(wr^2 * (exposure_variance + outcome_variance))
}

forward_harmon <-     
    harmon |>
    inner_join(h1 |> select(gwas, somamer)) |>
    group_by(gwas, somamer, snp_id) |>
    filter(row_number() == 1) |>
    ungroup() |>
    inner_join(matched_ids) |>
    inner_join(data_for_bayes$variant_profile  |> select(somamer, gwas, snp_id, coding_variant)) |>
    mutate(beta_hat = beta.outcome/beta.exposure)

forward_harmon$var_beta_hat <- 
    wald_delta(wr = forward_harmon$beta_hat, 
               beta.out = forward_harmon$beta.outcome, 
               beta.exp = forward_harmon$beta.exposure,
               se.out = forward_harmon$se.outcome,
               se.exp = forward_harmon$se.exposure)

forward_harmon$se_beta_hat <- sqrt(forward_harmon$var_beta_hat)

forward_harmon_multi <- 
    forward_harmon |> 
    filter(!is.nan(se_beta_hat)) |>
    group_by(gwas, somamer) |> filter(n() > 1) |> ungroup()

dosages <- readRDS("data-raw/dosages_from_database.rds")
genotype_matrix <- dosages$dosages
soma_pheno_pairs <- forward_harmon_multi |> distinct(gwas, somamer)
model <- stan_model("stan/heterogeneity/snp_heterogeneity.stan")

hetergen_test <- function(i, pair_data, snp_data, dosage_matrix, model) {
    print(i)
    ## Get somamer x phenotype pair
    pair <- pair_data[i, ]
    somamer <- pair$somamer
    gwas <- pair$gwas
    d <- snp_data |> filter(gwas == {{ gwas }}, somamer == {{ somamer }}) |> arrange(snp_id)
    
    ## Get dosages, match order
    dosage_ids <- d$dosage_id
    dosage_idx <- which(colnames(dosage_matrix) %in% dosage_ids)
    dosage_submatrix <- dosage_matrix[, dosage_idx]
    match_idx <- match(dosage_ids, colnames(dosage_submatrix))
    dosage_submatrix <- dosage_submatrix[, match_idx]
    
    ## Prepare data
    N <- nrow(d)
    beta_hat <- d$beta_hat
    se_beta_hat <- d$se_beta_hat
    R <- cor(dosage_submatrix)
    Sigma <- diag(se_beta_hat) %*% R %*% diag(se_beta_hat)
    
    stan_data <-
        list(N = N,
             beta_hat = beta_hat,
             Sigma = Sigma
        )
    fit <- sampling(model, 
                    data = stan_data, 
                    chains = 4, 
                    iter = 10000, 
                    seed = 99, 
                    control = list(adapt_delta = 0.99))
    posterior_fit <- extract(fit, pars = c("theta", "tau"))
    
    save_path <- sprintf("data/heterogeneity/full_posteriors/%s_posterior.rds", i)
    saveRDS(object = posterior_fit, file = save_path)
}

res <- 
    lapply(1:nrow(soma_pheno_pairs), 
           hetergen_test, 
           pair_data = soma_pheno_pairs, 
           snp_data = forward_harmon_multi, 
           dosage_matrix = genotype_matrix, 
           model = model)