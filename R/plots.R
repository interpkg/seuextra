
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
    colors='#E0E0E0,red',
    title = '',
    xlab='UMAP 1', 
    ylab='UMAP 2',
    min_cutoff = NA,
    max_cutoff = NA
){
    features <- stringr::str_split(features, ',')[[1]]
    colors <- stringr::str_split(colors, ',')[[1]]

    p <- Seurat::FeaturePlot(object = obj, slot=slot, features = features, 
            reduction = reduction, order=T, pt.size=0.01, cols=colors,
            min.cutoff=min_cutoff, max.cutoff=max_cutoff)

    p <- p & CustThemeOption1() & labs(title=title, x =xlab, y = ylab)
     
    return(p)
}





#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% scCustomize


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
    x_label='UMAP 1', y_label='UMAP 2',
    colors=viridis_plasma_dark_high
){

    x_label <- 'UMAP 1'
    y_label <- 'UMAP 2'

    if (reduction=='wnn.umap'){
        x_label <- 'UMAP 1 (WNN)'
        y_label <- 'UMAP 2 (WNN)'
    }

    p <- scCustomize::FeaturePlot_scCustom(seurat_object = obj, 
            features = features, reduction=reduction, order=T, 
            pt.size=0.01, colors_use = colors, na_color = "lightgray")

    p <- p + theme(plot.title = element_text(size = 10),
                text=element_text(size=8), 
                axis.text=element_text(size=6)) +
            theme(legend.key.width = unit(3, 'mm')) +
            theme(axis.text.x=element_blank(),
                axis.ticks.x=element_blank(),
                axis.text.y=element_blank(),
                axis.ticks.y=element_blank()) +
            ggtitle(title) +
            labs(x = x_label, y = y_label) +
            theme(axis.line = element_line(colour = 'black', size = 0.3),
                axis.ticks = element_line(linewidth = 0.3),
                axis.ticks.length=unit(1, "mm"))

    return(p)
}



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

    if (length(color_set) > 1){
        p <- p + scale_color_manual(values=color_set)
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
#' @import ggplot2
#'
#' @export
#'
Seurat_DimPlot2 <- function(obj=NULL,
    reduction='umap', group_by='cell_type2',
    title='', x_lab='UMAP 1', y_lab='UMAP 2',
    colors='', label=FALSE, repel=TRUE
){

    p <- Seurat::DimPlot(obj, reduction=reduction, group.by=group_by, label.size=3, label=label, pt.size=0.01, repel=repel)

    if (length(colors) > 2){
        p <- p + scale_color_manual(values=colors)
    }
    
    p <- p + plot_layout(guides = "collect") & 
            theme(plot.title = element_text(size = 10),
                text=element_text(size=8), 
                axis.text=element_text(size=6)) &
            theme(legend.position = 'none') &
            theme(axis.text=element_blank(),
                axis.ticks=element_blank()) &
            labs(title=title, x = x_lab, y = y_lab) 


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
            annotation_custom(grob = linesGrob(), xmin = xmin*1.5, xmax = 0.8, ymin = ymin*1.3, ymax = ymin*1.3) +
            # y
            annotation_custom(grob = linesGrob(), xmin = xmin*1.5, xmax = xmin*1.5, ymin = ymin*1.3, ymax = 0) +
            coord_cartesian(xlim=c(xmin, xmax), ylim = c(ymin, ymax), clip = "off") +
            theme(axis.title = element_text(hjust = 0))

    return(p)
}














