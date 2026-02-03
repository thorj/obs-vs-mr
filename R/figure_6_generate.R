stan_data <- data_for_bayes$stan_data


bayes_res <- readRDS("data/stan-models/model_base_fstat_total_coding_proteincorr.rds")

posterior_p <- 
    rstan::extract(bayes_res, 
                   pars = names(bayes_res)[grepl("^alpha$|alpha_p|gamma_g|beta_s|theta_f|theta_tc", names(bayes_res))])


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
    inner_join(f2, by = join_by("event")) |>
    inner_join(data.frame(phenotype_id = 1:21, 
                          m = apply(q_t_p, MARGIN = 2, median),
                          q025 = apply(q_t_p, MARGIN = 2, quantile, probs = 0.025),
                          q975 = apply(q_t_p, MARGIN = 2, quantile, probs = 0.975))) |>
    arrange(m) |>
    ungroup() |>
    mutate(order = n() - row_number()) |>
    inner_join(
        fread("tables/00b_no_obs_res.csv") |>
            mutate(r_no_obs = total_agreef/total_mr_for) |>
            select(event, r_no_obs) 
    )


global_mean <- rowMeans(q_t_p) |> median()
global_q95 <- rowMeans(q_t_p) |> quantile(c(0.025, 0.975))

pred_table <- 
    q_t_p_df |>
    select(trait, group, phenotype = phenotype2, event, N, obs_ratio = r, 
           no_obs_ratio = r_no_obs, bayes_r = m, bayes_low = q025, bayes_high = q975) |>
    mutate(global_median = global_mean,
           global_low = global_q95[1], 
           global_high = global_q95[2]) |>
    inner_join(f2 |> select(event, order)) |>
    arrange(order) |>
    select(-order) |>    
    mutate(across(where(is.numeric), ~round(.x, 2))) 


pred_table2 <- 
    pred_table |>
    mutate(within_ci_obs = as.numeric(ifelse(obs_ratio >= bayes_low & obs_ratio <= bayes_high, 1, 0)),
           within_ci_nobs = as.numeric(ifelse(no_obs_ratio >= bayes_low & no_obs_ratio <= bayes_high, 1, 0)),
           within_glob_obs = as.numeric(ifelse(obs_ratio >= global_low & obs_ratio <= global_high, 1, 0)),
           within_glob_nobs = as.numeric(ifelse(no_obs_ratio >= global_low & no_obs_ratio <= global_high, 1, 0)))

pred_table2$within_glob_obs |> sum(na.rm = T)

pred_table2 |> filter(within_glob_nobs != 1)

write_csv(x = pred_table2, file = "tables/tables_new/predicted_ratios.csv")

## no_obs_filter
q_t_p_df |> select(phenotype2, m, q025, q975, r, r_no_obs) |> print(n=50)

q_t_p_df |> select(phenotype2, m, q025, q975, r, r_no_obs) |>
    filter(r >= global_q95[1] & r <= global_q95[2])

q_t_p_df |> select(phenotype2, m, q025, q975, r, r_no_obs) |>
    filter(r_no_obs >= global_q95[1] & r_no_obs <= global_q95[2])

q_t_p_df |> select(phenotype2, m, q025, q975, r, r_no_obs) |>
    filter(r_no_obs >= q025 & r_no_obs <= q975)

q_t_p_df |> select(phenotype2, m, q025, q975, r, r_no_obs) |>
    filter(r < q025 | r > q975)

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
    # geom_hline(yintercept = global_q95[1], lty = 2, alpha = 0.5) +
    # geom_hline(yintercept = global_q95[2], lty = 2, alpha = 0.5) +
    # geom_crossbar(aes(y = m,
    #                   ymin = q025, 
    #                   ymax = q975), 
    #               middle.linetype = "blank",
    #               width = 0.6, 
    #               alpha = 0.6) +
    geom_point(aes(y = m, color = "Predicted agreement", shape = "Predicted agreement"), size = point_size) +
    geom_errorbar(aes(ymin = q025, ymax = q975), width = 0.1) +
    geom_point(aes(y = r, color = "Observed agreement (Primary)", shape = "Observed agreement (Primary)"), size = point_size) +
    geom_point(aes(y = r_no_obs, color = "Observed agreement (Secondary)", shape = "Observed agreement (Secondary)"), size = point_size) +
    scale_y_continuous(breaks = 0:10/10, labels = scales::percent) +
    scale_x_discrete(labels = q_t_p_df$phenotype2) +
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
    coord_flip() #+
# guides(
#     color = guide_legend(nrow = 2, order = 1,
#                         title = "",
#                         theme = theme(legend.key.spacing.y = unit(0, "cm"),
#                                       legend.margin = margin())),
#     shape = guide_legend(nrow = 2, title = "", order = 1,
#                          theme = theme(legend.key.spacing.y = unit(0, "cm"),
#                                        legend.margin = margin()))
# )

hibayes_plt

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
    #((observed_agreement | tc_post) / hibayes_plt) +
    free(tc_post, side = "l") + hibayes_plt +
    plot_annotation(tag_levels = "A", tag_prefix = "Fig. ") &
    theme(plot.tag.position = c(0, 1),
          plot.tag = element_text(size = 12, hjust = 0, vjust = 0))
ggsave(filename = "img/highres/figure4_agreement_ratios_posterior.pdf", 
       plot = agreement_plots, 
       width = 18, height = 6, device = cairo_pdf)

ggsave(filename = "img/lowres/figure4_agreement_ratios_posterior.png", 
       plot = agreement_plots, 
       width = 18, height = 6, dpi = 300)

ggsave(filename = "img/highres/figure4_agreement_ratios.pdf", 
       plot = observed_agreement, 
       width = 12, height = 7, device = cairo_pdf)

ggsave(filename = "img/lowres/figure4_agreement_ratios.png", 
       plot = observed_agreement, 
       width = 12, height = 7, dpi = 300)

#ggsave(filename = "img/highres/figure4_predicted_ratio_agree.pdf", plot = hibayes_plt, width = 12, height = 8, device = cairo_pdf)

fwrite(x = q_t_p_df, file = "tables/raw-tables/predicted_ratios_phenotypes.csv")


### coding variant
posterior_p$theta_tc |> quantile(probs = c(0.025, 0.5, 0.975))
