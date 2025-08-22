
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
scCustomize_ClusteredDotPlot <- function(obj=NULL, group='', d_markers='', top_n=7, k_size=7){

    top_markers <- scCustomize::Extract_Top_Markers(marker_dataframe = d_markers, num_genes = top_n, 
            group_by = "cluster", rank_by = "avg_log2FC",
            named_vector = FALSE, make_unique = TRUE)
    
    p <- scCustomize::Clustered_DotPlot(seurat_object = obj, features = top_markers, group.by = group, 
            row_label_size = 6, column_label_size = 8, legend_label_size = 6, legend_title_size = 8,
            k = k_size, plot_km_elbow=FALSE)

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
    sorted_name <- unique(topX.markers$cluster)

    color.scheme <- rev(brewer.pal(9,"RdBu"))
    p <- Seurat::DotPlot(obj, assay = assay_use, group.by = group, features = unique(topX.markers$gene)) +
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
            theme(
                    legend.position="bottom", 
                    legend.key.width = unit(4, 'mm'),
                    legend.key.height = unit(2.5, 'mm'),
                    legend.title = element_text(size=7),
                    legend.text=element_text(size=7)
            )

    p

}






#' Seurat DotPlot - default motif
#'
#' @param obj object
#' @param assay_use assay
#'
#' @return plot
#'
#' @import ggplot2
#'
#' @export
#'
Seurat_DotPlot <- function(
    obj=NULL, 
    assay_use='chromvar', 
    features=NULL, 
    col_min = 0,
    cex = 6,
    colors="viridis::mako"
) {

    p <- DotPlot(obj,
            features = features,
            dot.scale = cex,
            col.min = col_min, 
            scale = F, 
            assay = assay_use) +  
        RotatedAxis()   +
        paletteer::scale_colour_paletteer_c(colors, direction = -1) +
        labs(x='', y='')

    p <- p + theme(plot.title = element_text(size = 8),
            text=element_text(size=6, face="bold"), 
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








#' Seurat DotPlot - default motif
#'
#' @param obj object
#' @param assay_use assay
#'
#' @return plot
#'
#' @import ggplot2
#'
#' @export
#'
Seurat_DotPlot_test <- function(
    obj=NULL, 
    assay_use='chromvar', 
    features=NULL, 
    col_min = 0,
    cex = 6,
    title="",
    colors="magma"
) {

    p <- DotPlot(obj,
            features = features,
            dot.scale = cex,
            col.min = col_min, 
            scale = F, 
            assay = assay_use)

    p <- p + geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.4) +
        guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"), title='Percentage')) +
        scale_colour_viridis(option="magma") +
        labs(title=title, x='', y='') +
        coord_flip() +
        #theme(text=element_text(size=7, face='bold'), axis.text=element_text(size=6, face='bold')) +
        theme(plot.title = element_text(size=8, hjust=0.5), text=element_text(size=7), axis.text=element_text(size=6)) +
        theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1)) +
        theme(legend.key.size = unit(4, 'mm')) +
        theme(axis.line=element_line(size=.3), axis.ticks = element_line(size = .3)) 

    p

}















