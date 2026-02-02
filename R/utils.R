
export_image <- function(plot, fig_name, extension, width, height, dpi = 300) {
    extension <- match.arg(extension, c("png", "pdf"))
    res_dir <- if (extension == "png") "low_res" else "high_res"
    path <- file.path("images", res_dir, paste0(fig_name, ".", extension))
    if (extension == "png") {
        ggplot2::ggsave(
            filename = path,
            plot = plot,
            width = width,
            height = height,
            dpi = dpi
        )
    } else {
        ggplot2::ggsave(
            filename = path,
            plot = plot,
            width = width,
            height = height
        )
    }
    invisible(path)
}