data_paths <- list(
    phenotype_overview = "data/phenotype_overview.csv",
    observational_full = "data/observational_results_full.csv",
    observational_significant = "data/reoccurring_data/bonferroni_sig_obs.csv",
    primary_mr_overview = "data/primary_mr_overview.csv",
    group_map = "data/stan_map_pheno_to_group.rds",
    protein_map = "data/stan_map_soma_to_egs.rds",
    stan_data = "data/stan_agreement_model_data.rds",
    bayes_results = "data/model_base_fstat_total_coding_proteincorr.rds",
    primary_mr_summary = "data/reoccurring_data/primary_mr_summary.rds",
    secondary_mr_summary = "data/reoccurring_data/secondary_mr_summary.rds",
    cis_info = "data/reoccurring_data/cis_protein_data.rds"
)

export_image <- function(plot, fig_name, width, height, dpi = 300) {
    path <- 
        file.path("figures", 
                  c("low_res", "high_res"), 
                  paste0(fig_name, ".", c("png", "pdf")))
    ggplot2::ggsave(
        filename = path[1],
        plot = plot,
        width = width,
        height = height,
        dpi = dpi
    )
    ggplot2::ggsave(
        filename = path[2],
        plot = plot,
        width = width,
        height = height,
        device = cairo_pdf
    )
    invisible(path)
}

create_folder_structure <- function() {
    ## Figures
    figures_path <- file.path("figures", c("low_res", "high_res"))
    for (p in figures_path) {
        dir.create(p, recursive = TRUE, showWarnings = FALSE)
    }
    invisible(figures_path)
}


