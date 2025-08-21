


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#  Feature Plot
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
            text=element_text(size=6, face="bold"), 
            axis.text=element_text(size=5, face="bold"),
            axis.line = element_line(colour = 'black', size = 0.3),
            axis.ticks = element_line(linewidth = 0.3),
            axis.ticks.length=unit(1, "mm"),
            legend.title = element_text(size=7, face="bold"),
            legend.text = element_text(size=6, face="bold"),
            legend.key.width = unit(3, 'mm'),
            legend.key.height = unit(3, 'mm')
        )
        
    p
}





#' Seurat FeaturePlot
#'
#' @param obj object
#' @param features gene name
#'
#' @return plot
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
    colors=c('#E0E0E0', 'red'),
    title = '',
    xlab='UMAP 1', 
    ylab='UMAP 2',
    order=TRUE,
    min_cutoff = NA,
    max_cutoff = NA
){
    #features <- stringr::str_split(features, ',')[[1]]
    #colors <- stringr::str_split(colors, ',')[[1]]

    p <- Seurat::FeaturePlot(object = obj, slot=slot, features = features, 
            reduction = reduction, order=order, pt.size=0.01, cols=colors,
            min.cutoff=min_cutoff, max.cutoff=max_cutoff)

    p <- p & CustThemeOption1() & labs(title=title, x =xlab, y = ylab)
     
    return(p)
}




#' scCustom FeaturePlot
#'
#' #gene exp: colors_use = viridis_inferno_dark_high
#' #gene act: colors_use = viridis_light_high
#' #motif: colors_use = viridis_plasma_dark_high
#'
#' @param obj object
#' @param features gene name
#'
#' @return plot
#'
#' @import scCustomize ggplot2
#'
#' @export
#'
scCustom_FeaturePlot <- function(
    obj=NULL, slot="data", 
    gene=NULL, features=NULL, 
    reduction='wnn.umap', 
    title='',
    xlab='UMAP 1', ylab='UMAP 2',
    order=TRUE,
    colors=viridis_plasma_dark_high
){
    p <- scCustomize::FeaturePlot_scCustom(seurat_object = obj, 
            features = features, reduction=reduction, order=order, 
            pt.size=0.01, colors_use = colors, na_color = "lightgray")

    p <- p & CustThemeOption1() & labs(x =xlab, y = ylab)

    return(p)
}







