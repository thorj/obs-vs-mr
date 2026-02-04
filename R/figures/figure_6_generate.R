library(tidyverse)
library(rstan)
library(patchwork)
library(Cairo)
box::use(data.table[fread, fwrite])
source("R/utils.R")

## Constants
other_lty <- 0
point_size <- 2

## ============================================================
## Load data
## ============================================================
f <- fread(data_paths$phenotype_overview)
obs <- fread(data_paths$observational_full)
h1 <- fread(data_paths$primary_mr_overview)
group <- readRDS(data_paths$group_map)
protein_map <-readRDS(data_paths$protein_map)
stan_data <- readRDS(data_paths$stan_data)
bayes_res <- readRDS(data_paths$bayes_results)
primary_mr_summary <- readRDS(data_paths$primary_mr_summary)
primary_mr_agreement <- primary_mr_summary$primary_mr_agreement
secondary_mr_summary <- readRDS(data_paths$secondary_mr_summary)
secondary_mr_agreement <- secondary_mr_summary$secondary_mr_agreement

## ============================================================
## Posterior extraction
## ============================================================
posterior_p <- 
    rstan::extract(bayes_res, 
                   pars = names(bayes_res)[grepl("^alpha$|alpha_p|gamma_g|beta_s|theta_f|theta_tc", names(bayes_res))])

## ============================================================
## Helpers
## ============================================================

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
## ============================================================
## Compute predicted agreement per phenotype
## ============================================================
q_t_p <- compute_tilde_q_p_full(posterior_p, 
                                group$groupf, 
                                L = 10000, 
                                P = stan_data$P, 
                                input_data = h1)

q_summary <- 
    data.frame(phenotype_id = 1:ncol(q_t_p), 
                        m = apply(q_t_p, MARGIN = 2, median),
                        q025 = apply(q_t_p, MARGIN = 2, quantile, probs = 0.025),
                        q975 = apply(q_t_p, MARGIN = 2, quantile, probs = 0.975))

global_mean <- rowMeans(q_t_p) |> median()
global_q95 <- rowMeans(q_t_p) |> quantile(c(0.025, 0.975))

qt_area <- data.frame(order2 = c(0.5, 22))

## ============================================================
## Assemble plotting data frame
## ============================================================

q_t_p_df <- 
    h1 |> 
    group_by(phenotype_id, event, groupf) |> 
    count() |> 
    select(-n) |>
    inner_join(f |> 
                   filter(N > 0) |> 
                   mutate(r = Y/N), 
               by = join_by("event")) |> 
    inner_join(q_summary,
               by = join_by(phenotype_id)) |>
    arrange(m) |>
    ungroup() |>
    mutate(order = n() - row_number()) |>
    inner_join(
        secondary_mr_agreement |>
            mutate(r_no_obs = abs(rf)) |>
            select(event, r_no_obs),
        by = join_by(event)
    ) |>
    inner_join(primary_mr_agreement |> 
                   select(event, order2 = order),
               by = join_by(event)
    ) |>
    arrange(order2)

## ============================================================
## Plot: agreement ratios
## ============================================================

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

## ============================================================
## Plot: posterior density for theta_tc
## ============================================================

tc_med <- posterior_p$theta_tc |> median()
ci <- posterior_p$theta_tc |> quantile(probs = c(0.025, 0.975))

dens <- density(posterior_p$theta_tc)
dens_df <- tibble(x = dens$x, y = dens$y)

tc_post <- 
    ggplot(dens_df, aes(x = x, y = y)) +
    geom_area(data = filter(dens_df, x <= ci[1]),
              aes(x = x, y = y), fill = "gray80", alpha = 0.6) +
    geom_area(data = filter(dens_df, x >= ci[1], x <= ci[2]),
              aes(x = x, y = y), fill = "#648FFF", alpha = 0.6) +
    geom_area(data = filter(dens_df, x >= ci[2]),
              aes(x = x, y = y), fill = "gray80", alpha = 0.6) +
    geom_line(color = "black") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    geom_vline(xintercept = tc_med, linetype = "dotted", color = "red", linewidth = 1.2) +
    labs(
        x = "Log-odds ratio for coding variant burden",
        y = "Posterior density"
    ) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 9)) +
    theme_bw(base_size = 12) +
    theme(axis.title = element_text(face = "bold"))

## ============================================================
## Combine + export
## ============================================================

agreement_plots <-
    free(tc_post, side = "l") + hibayes_plt +
    plot_annotation(tag_levels = "A", tag_prefix = "Fig. ") &
    theme(plot.tag.position = c(0, 1),
          plot.tag = element_text(size = 12, hjust = 0, vjust = 0))

export_image(plot = agreement_plots, 
             fig_name = "figure6_agreement_ratios_posterior", 
             width = 18, height = 6, dpi = 300)

#fwrite(x = q_t_p_df, file = "tables/raw-tables/predicted_ratios_phenotypes.csv")

