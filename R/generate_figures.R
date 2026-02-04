source("R/utils.R")
generate_figures <- function() {
    create_folder_structure()
    source("R/figures/figure_2_generate.R")
    source("R/figures/figure_3_generate.R")
    source("R/figures/figure_4_generate.R")
    source("R/figures/figure_5_generate.R")
    source("R/figures/figure_6_generate.R")
}