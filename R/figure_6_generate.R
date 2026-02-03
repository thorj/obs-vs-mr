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
h1 <- fread("data/primary_mr_overview.csv")
group <- readRDS("data/stan_map_pheno_to_group.rds")
protein_map <-readRDS("data/stan_map_soma_to_egs.rds")
stan_data <- readRDS("data/stan_agreement_model_data.rds")
bayes_res <- readRDS("data/model_base_fstat_total_coding_proteincorr.rds")

posterior_p <- 
    rstan::extract(bayes_res, 
                   pars = names(bayes_res)[grepl("^alpha$|alpha_p|gamma_g|beta_s|theta_f|theta_tc", names(bayes_res))])

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


## Helper function to get variable from posterior
retrieve_post_col <- function(posterior, para, idx) {
    cname <- sprintf("%s[%s]", para, idx)
    return(posterior[[cname]])
}

## Compute q for each SOMAmer per phenotype, then average over phenotype to 
## get phenotypic effect
compute_tilde_q_p_full <- function(posterior, group_vec, L, P, input_data) {
    set.seed(1)
    tilde_q_p <- matrix(NA_real_, nrow = L, ncol = P)
    
    # Extract relevant posterior columns
    S <- sample(1:length(posterior$alpha), size = L, replace = FALSE)
    alpha_draws    <- posterior$alpha[S]
    theta_f <- posterior$theta_f[S]
    theta_tc <- posterior$theta_tc[S]
    
    for (p in 1:P) {
        # Get input data for phenotype p
        df  <- input_data |> filter(phenotype_id == {{ p }})
        somamers <- df$somamer_id
        hmf_vec <- df$hmf_s
        tc_vec <- df$log_tc
        K <- length(somamers)
        
        ## Get posterior for phenotype p
        alpha_p <- retrieve_post_col(posterior, "alpha_p", p)[S]
        
        ## Get posterior for group g
        g <- group_vec[p]
        gamma_p <- retrieve_post_col(posterior, "gamma_g", g)[S]
        
        ## Initialize eta matrix
        linear_pred <- matrix(NA_real_, nrow = L, ncol = K)
        
        for (j in 1:K) {
            s_idx <- somamers[j]
            beta_s <- retrieve_post_col(posterior, "beta_s", s_idx)[S]
            linear_pred_j <- 
                alpha_draws + alpha_p + gamma_p + beta_s +
                theta_f * hmf_vec[j] + 
                theta_tc * tc_vec[j]
            linear_pred[, j] <- plogis(linear_pred_j)
        }
        tilde_q_p[, p] <- rowMeans(linear_pred)
    }
    
    colnames(tilde_q_p) <- paste0("phenotype_", 1:P)
    return(tilde_q_p)
}

compute_tilde_q_p_full(posterior_p, group$groupf, L = 10, P = stan_data$P, input_data = h1)
q_t_p <- compute_tilde_q_p_full(posterior_p, group$groupf, L = 10000, P = stan_data$P, input_data = h1)

q_t_p_df <- 
    h1 |> group_by(phenotype_id, event, groupf) |> count() |> select(-n) |>
    inner_join(f |> filter(N > 0) |> mutate(r = Y/N), by = join_by("event")) |>  ### careful... any f |> filter(N > 0) might be wrong after update
    inner_join(data.frame(phenotype_id = 1:21, 
                          m = apply(q_t_p, MARGIN = 2, median),
                          q025 = apply(q_t_p, MARGIN = 2, quantile, probs = 0.025),
                          q975 = apply(q_t_p, MARGIN = 2, quantile, probs = 0.975))) |>
    arrange(m) |>
    ungroup() |>
    mutate(order = n() - row_number()) |>
    inner_join(
        fread("data/00b_no_obs_res.csv") |>
            mutate(r_no_obs = total_agreef/total_mr_for) |>
            select(event, r_no_obs) 
    )


global_mean <- rowMeans(q_t_p) |> median()
global_q95 <- rowMeans(q_t_p) |> quantile(c(0.025, 0.975))

qt_area <- 
    data.frame(order2 = c(0.5, 22))

point_size <- 2
q_t_p_df <-
    q_t_p_df |> inner_join(agreement |> select(event, order2 = order)) |>
    arrange(order2)

hibayes_plt <- 
    q_t_p_df |>
    ggplot(aes(x = order2)) +
    geom_ribbon(data = qt_area, aes(x = order2, ymin = global_q95[1], ymax = global_q95[2]), fill = "gray30", alpha = 0.2) +
    geom_hline(yintercept = global_mean, lty = 2, alpha = 0.5) +
    geom_point(aes(y = m, color = "Predicted agreement", shape = "Predicted agreement"), size = point_size) +
    geom_errorbar(aes(ymin = q025, ymax = q975), width = 0.1) +
    geom_point(aes(y = r, color = "Observed agreement (Primary)", shape = "Observed agreement (Primary)"), size = point_size) +
    geom_point(aes(y = r_no_obs, color = "Observed agreement (Secondary)", shape = "Observed agreement (Secondary)"), size = point_size) +
    scale_y_continuous(breaks = 0:10/10, labels = scales::percent) +
    scale_x_discrete(labels = q_t_p_df$phenotype) +
    scale_color_manual(
        name = NULL,
        values = c("Predicted agreement" = "black", 
                   "Observed agreement (Primary)" = "#1E88E5", 
                   "Observed agreement (Secondary)" = "orange")
    ) +
    scale_shape_manual(
        name = NULL,
        values = c("Predicted agreement" = 19, 
                   "Observed agreement (Primary)" = 15, 
                   "Observed agreement (Secondary)" = 17)
    ) +
    labs(x = "", y = "Agreement ratio") +
    theme_bw(base_size = 12) +
    theme(axis.title = element_text(face = "bold"), legend.position = "bottom") +
    coord_flip()

#hibayes_plt

#### Posterior density
## Start with posterior distribution of tc
tc_med <- posterior_p$theta_tc |> median()

ci <- posterior_p$theta_tc |> quantile(probs = c(0.025, 0.975))

dens <- density(posterior_p$theta_tc)
dens_df <- tibble(x = dens$x, y = dens$y)

tc_post <- 
    ggplot(dens_df, aes(x = x, y = y)) +
    # Shade left tail
    geom_area(data = filter(dens_df, x <= ci[1]),
              aes(x = x, y = y), fill = "gray80", alpha = 0.6) +
    # Shade 95% CI region
    geom_area(data = filter(dens_df, x >= ci[1], x <= ci[2]),
              aes(x = x, y = y), fill = "#648FFF", alpha = 0.6) +
    # Shade right tail
    geom_area(data = filter(dens_df, x >= ci[2]),
              aes(x = x, y = y), fill = "gray80", alpha = 0.6) +
    # Density line
    geom_line(color = "black") +
    # Reference line at 0
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    geom_vline(xintercept = tc_med, linetype = "dotted", color = "red", linewidth = 1.2) +
    labs(
        x = "Log-odds ratio for coding variant burden",
        y = "Posterior density"
    ) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 9)) +
    theme_bw(base_size = 12) +
    theme(axis.title = element_text(face = "bold"))


agreement_plots <-
    free(tc_post, side = "l") + hibayes_plt +
    plot_annotation(tag_levels = "A", tag_prefix = "Fig. ") &
    theme(plot.tag.position = c(0, 1),
          plot.tag = element_text(size = 12, hjust = 0, vjust = 0))

export_image(plot = agreement_plots, fig_name = "figure6_agreement_ratios_posterior", width = 18, height = 6, dpi = 300)

#fwrite(x = q_t_p_df, file = "tables/raw-tables/predicted_ratios_phenotypes.csv")

