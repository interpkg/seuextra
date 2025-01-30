
#' ggplot2 theme setting 1
#'
#' @return plot
#'
#' @import ggplot2
#'
#' @export
#'
CustThemeOption1 <- function()
{ 
    p <- theme(plot.title = element_text(size = 8),
            text=element_text(size=6), 
            axis.text=element_text(size=5),
            axis.line = element_line(colour = 'black', size = 0.3),
            axis.ticks = element_line(linewidth = 0.3),
            axis.ticks.length=unit(1, "mm"),
            legend.title = element_text(size=6),
            legend.text = element_text(size=5),
            legend.key.width = unit(3, 'mm'),
            legend.key.height = unit(3, 'mm')
        )
        
    p
}



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
    slot="data", 
    reduction='umap', 
    features=NULL, 
    colorset='#E0E0E0,red',
    title = '',
    xlab='UMAP 1', 
    ylab='UMAP 2',
    min_cutoff = NA,
    max_cutoff = NA
){
    features <- stringr::str_split(features, ',')[[1]]
    colorset <- stringr::str_split(colorset, ',')[[1]]

    p <- Seurat::FeaturePlot(object = obj, slot=slot, features = features, 
            reduction = reduction, order=T, pt.size=0.01, cols=colorset,
            min.cutoff=min_cutoff, max.cutoff=max_cutoff)

    p <- p & CustThemeOption1() & labs(title=title, x =xlab, y = ylab)
     
    return(p)
}







