library(tidyverse)
library(rstan)
options(mc.cores = parallel::detectCores())
box::use(data.table[fread])

stan_data <- readRDS("data/stan_agreement_model_data.rds")

run_model <- function(model, L, seed) {
    model_path <- sprintf("stan/%s.stan", model)
    fit <- 
        stan(
            file = model_path,
            data = stan_data,
            chains = 4,
            iter = L,
            seed = seed,
            control = list(adapt_delta = 0.99)
        )
    save_path <- sprintf("data/stan-models/%s.rds", model)
    saveRDS(fit, save_path)
}

L <- 10000

run_model(model = "model_base_fstat_total_coding_proteincorr", L = L, seed = 99)
