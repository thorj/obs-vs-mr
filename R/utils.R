
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
