#' @include plot_parameter.R
#'
NULL


#' Elbow Plot
#'
#' @param obj object
#' @param reduc name
#'
#' @import ggplot2
#'
#' @export
#'
ElbowPlot2 <- function(obj, reduc='pca')
{ 
    p <- Seurat::ElbowPlot(obj, ndims = 50, reduction = reduc)
    p <- p & geom_hline(yintercept=1, linetype="dashed", color = "blue")
    p <- p & geom_hline(yintercept=2, linetype="dashed", color = "orange")
    p <- p & geom_hline(yintercept=3, linetype="dashed", color = "red")
    
    p
}




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#  Feature Plot
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Seurat FeaturePlot
#'
#' @param obj object
#' @param features gene name
#'
#' @return data frame
#'
#' @import stringr ggplot2
#'
#' @export
#'
Seurat_FeaturePlot <- function(
    obj=NULL, 
    item='',
    slot="data", 
    reduction='umap', 
    features=NULL, 
    colorset='lightgrey,red',
    title = '',
    xlab='UMAP 1', 
    ylab='UMAP 2',
    min_cutoff = NA,
    max_cutoff = NA
){
    features <- stringr::str_split(features, ',')[[1]]
    colorset <- stringr::str_split(colorset, ',')[[1]]

    if (item == motif){
        colorset <- c("lightgrey", "#FAAD61", "#BA2724")
        min_cutoff = 'q10'
        max_cutoff = 'q90'
    }

    p <- Seurat::FeaturePlot(object = obj, slot=slot, features = features, 
            reduction = reduction, order=T, pt.size=0.01, cols=colorset,
            min.cutoff=min_cutoff, max.cutoff=max_cutoff)

    if (item == 'motif'){
        p <- p & CustThemeOption2() & labs(x =xlab, y = ylab) & ggtitle(title)
    } else {
        p <- p & CustThemeOption1() & labs(x =xlab, y = ylab)
    }
    
            
    return(p)
}







