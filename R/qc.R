


#' Elbow Plot
#'
#' @param seurat object
#' @param group name
#' @return data frame
#' @export
#'
ElbowPlot2 <- function(obj, reduc='pca')
{ 
    p <- Seurat::ElbowPlot(obj, ndims = 50, reduction = reduc)
    p <- p & geom_hline(yintercept=1, linetype="dashed", color = "red")
    ggplot2::ggsave(paste0('ElbowPlot.', reduc, '.pdf'), w=4, h=2.5)
}






