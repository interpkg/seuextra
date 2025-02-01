
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
#  Dot Plot
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' scCustomize ClusteredDotPlot 
#'
#' @param obj object
#' @param group
#'
#' @return plot
#'
#' @import ggplot2
#'
#' @export
#'
scCustomize_ClusteredDotPlot <- function(obj=NULL, group='', d_markers='', n=7){
    group_size <- length(unique(obj@meta.data[[group]]))

    top_markers <- scCustomize::Extract_Top_Markers(marker_dataframe = d_markers, num_genes = n, 
            group_by = group, rank_by = "avg_log2FC",
            named_vector = FALSE, make_unique = TRUE)
    
    p <- scCustomize::Clustered_DotPlot(seurat_object = obj, features = top_markers, group.by = group, k = group_size, plot_km_elbow=FALSE)

    return(p)
}





#' Seurat DotPlot
#'
#' @param obj object
#' @param group
#'
#' @return plot
#'
#' @import ggplot2
#'
#' @export
#'
Seurat_DotPlotMarkers <- function(obj=NULL, assay_use='SCT', d_markers=NULL, n=5, group='')
{
    topX.markers <- data.frame(d_markers %>% dplyr::group_by(cluster) %>% dplyr::slice(1:n))

    meta <- obj@meta.data

    if (group != ''){
        meta$cluster_plus <- paste0('C', meta$seurat_clusters, '.', meta[[group]])
    } else {
        meta$cluster_plus <- paste0('C', meta$seurat_clusters)
    }

    sorted_name <- unique(meta[,c('seurat_clusters', 'cluster_plus')])
    sorted_name <- as.vector(sorted_name$cluster_plus[order(as.numeric(sorted_name$seurat_clusters))])

    obj@meta.data <- meta

    color.scheme <- rev(brewer.pal(9,"RdBu"))
    p <- Seurat::DotPlot(obj, assay = assay_use, group.by = 'cluster_plus', features = unique(topX.markers$gene)) +
            theme_bw(base_line_size=0.1) +
            scale_size_area(max_size = 3) +
            scale_color_gradientn(colors=color.scheme, limits = c(-2.5, 2.5)) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
            labs(title='', x='', y='') +
            scale_y_discrete(limits=sorted_name) +
            theme(
                text=element_text(size=7), 
                axis.text=element_text(colour="black", size=7), 
                axis.title.x=element_blank()) +
            labs(y="") +
            theme(
                    legend.position="bottom", 
                    legend.key.width = unit(4, 'mm'),
                    legend.key.height = unit(2.5, 'mm'),
                    legend.title = element_text(size=7),
                    legend.text=element_text(size=7)
            )

    p

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
    p <- p + labs(title=title, x=x_lab, y=y_lab)

    return(p)
}






