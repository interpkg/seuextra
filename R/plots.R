#' @include plot_basic.R
#'
NULL


#' Elbow Plot
#'
#' @param seurat object
#' @param group name
#'
#' @import ggplot2
#'
#' @export
#'
#' @examples
#' library(Seurat)
#' ElbowPlot2(seu_obj = pbmc_small, reduc='pca')
#'
ElbowPlot2 <- function(seu_obj, reduc='pca')
{ 
    p <- Seurat::ElbowPlot(seu_obj, ndims = 50, reduction = reduc)
    p <- p & geom_hline(yintercept=1, linetype="dashed", color = "blue")
    p <- p & geom_hline(yintercept=2, linetype="dashed", color = "orange")
    p <- p & geom_hline(yintercept=3, linetype="dashed", color = "red")
    ggplot2::ggsave(paste0('ElbowPlot.', reduc, '.pdf'), w=4, h=2.5)
}




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#  Feature Plot
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Seurat FeaturePlot
#'
#' @param seurat object
#' @param feature name
#'
#' @return data frame
#'
#' @import stringr ggplot2
#'
#' @export
#'
Seurat_FeaturePlot <- function(
    obj=NULL, 
    slot="data", 
    reduction='umap', 
    features=NULL, 
    colorset='lightgrey,red',
    xlab='UMAP 1', 
    ylab='UMAP 2'
){
    features <- stringr::str_split(features, ',')[[1]]
    colorset <- stringr::str_split(colorset, ',')[[1]]

    p <- Seurat::FeaturePlot(object = obj, slot=slot, features = features, reduction = reduction, order=T, pt.size=0.01, cols=colorset)
    p <- p & CustThemeOption1() & labs(x =xlab, y = ylab) 
            
    return(p)
}








