

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#  UMAP Plot
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Seurat DimPlot
#'
#' @param obj object
#' @param reduction
#'
#' @return plot
#'
#' @import ggplot2
#'
#' @export
#'
Seurat_DimPlot <- function(obj=NULL,
    reduction='umap', group_by='cell_type2',
    title='', x_lab='UMAP 1', y_lab='UMAP 2',
    colors='', label=FALSE, repel=TRUE
){
    p <- Seurat::DimPlot(obj, reduction=reduction, group.by=group_by, label.size = 3, label=label, pt.size=0.01, repel=repel) 

    if (length(colors) > 1){
        p <- p + scale_color_manual(values=colors)
    }
    
    # theme_void() + 
    p <- p + ggtitle(title) + labs(x =x_lab, y =y_lab) +
            theme_bw() +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                    axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank()) +
            theme(plot.title = element_text(hjust=0.5, size=11),
                axis.title=element_text(size=9),
                legend.text=element_text(size=7),
                legend.key.size = unit(0.6, 'mm')
            )
            
    p <- p + guides(colour = guide_legend(override.aes=list(size=2)))

    return(p)
}



#' Seurat DimPlot with arrow
#'
#' @param obj object
#' @param reduction
#'
#' @return plot
#'
#' @import ggplot2 patchwork grid
#'
#' @export
#'
Seurat_DimPlot2 <- function(obj=NULL,
    reduction='umap', group_by='cell_type2',
    title='', x_lab='UMAP 1', y_lab='UMAP 2',
    colors='', label=FALSE, legend=FALSE, repel=TRUE,
    xa=1.2, xb=.3, ya=1.1, yb=.25
){

    p <- Seurat::DimPlot(obj, reduction=reduction, group.by=group_by, label.size=3, label=label, pt.size=0.01, repel=repel)
    
    p <- p + patchwork::plot_layout(guides = "collect") +
            theme(plot.title = element_text(size = 10),
                text=element_text(size=8), axis.text=element_text(size=6)) +
            theme(axis.text=element_blank(), axis.ticks=element_blank()) +
            labs(title=title, x = x_lab, y = y_lab) 

    if (!legend){
        p <- p + theme(legend.position = 'none')
    }
    
    if (length(colors) > 1){
        p <- p + scale_color_manual(values=colors)
    }

    # customized umap
    #print(colnames(obj@reductions$umap@cell.embeddings))
    xmin <- min(obj@reductions[[reduction]]@cell.embeddings[,1]) # UMAP-1
    xmax <- max(obj@reductions[[reduction]]@cell.embeddings[,1])

    ymin <- min(obj@reductions[[reduction]]@cell.embeddings[,2]) # UMAP-2
    ymax <- max(obj@reductions[[reduction]]@cell.embeddings[,2])

    # (optional) arrow = arrow(length = unit(2, "mm"), type = "closed")
    p <- p + theme(panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(),
                axis.line = element_blank()) +
            # x
            annotation_custom(grob = grid::linesGrob(), xmin = xmin*xa, xmax = xmin + abs(xmin)*xb, ymin = ymin*ya, ymax = ymin*ya) +
            # y
            annotation_custom(grob = grid::linesGrob(), xmin = xmin*xa, xmax = xmin*xa, ymin = ymin*ya, ymax = ymin + abs(ymin)*yb) +
            coord_cartesian(xlim=c(xmin, xmax), ylim = c(ymin, ymax), clip = "off") +
            theme(axis.title = element_text(hjust = 0))

    return(p)
}




#' scCustom DimPlot
#'
#' @param obj object
#' @param reduction
#'
#' @return plot
#'
#' @import scCustomize 
#'
#' @export
#'
scCustom_DimPlot <- function(obj=NULL,
    reduction='umap', group_by='cell_type2',
    title='', x_lab='UMAP 1', y_lab='UMAP 2',
    colors=NULL, label=FALSE, repel=TRUE, figure_plot=TRUE
){
    
    p <- scCustomize::DimPlot_scCustom(seurat_object = obj, group.by=group_by, reduction = reduction, 
                colors_use=colors, label=label, repel=repel, figure_plot = figure_plot)
    p <- p + labs(title=title, x=x_lab, y=y_lab) +
        theme(plot.title = element_text(size = 10), text=element_text(size=8), axis.text=element_text(size=6)) +
        theme(legend.text=element_text(size=6), legend.key.size = unit(0.6, 'mm'))

    return(p)
}







