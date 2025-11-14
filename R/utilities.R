

#' Convert Gene To Mouse
#'
#' @param gene list
#'
#' @import gprofiler2
#'
#' @export
#'
ConvertGeneToMouse <- function(gene){
    mouse_gene <- gprofiler2::gorth(gene, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name
    return(mouse_gene)
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





#' Export UMAP
#'
#' @param obj seurat object
#' @param umap str
#'
#' @return data frame
#'
#' @export
#'
ExportUMAP <- function(obj, umap)
{
    obj$UMAP_1 <- obj@reductions[[umap]]@cell.embeddings[,1]
    obj$UMAP_2 <- obj@reductions[[umap]]@cell.embeddings[,2]
    d_umap <- as.data.frame(obj@meta.data[,c('UMAP_1', 'UMAP_2')])

    return(d_umap)
}





#' Seurat Add Score
#'
#' @param obj seurat object
#' @param name score title 
#' @param features gene vectors 
#'
#' @return seurat object
#'
#' @import Seurat
#'
#' @export
#'
SeuratAddScore <- function(
    obj=NULL, 
    assay='SCT',
    slot='data',
    name='', 
    features=''
){
    all_gene <- rownames(obj[[assay]]@data)
    gene_intersect <- intersect(features, all_gene)

    obj <- AddModuleScore(
        object = obj,
        assay = assay,
        slot = slot,
        features = list(gene_intersect),
        name = name
    )

    return(obj)
}





#' EstimatedSignal
#'
#' @param obj seurat object
#' @param features genes
#' @param umap umap
#'
#' @return data frame
#'
#' @export
#'
EstimatedSignal <- function(obj=NULL, assay='SCT', name='', features='', min_cutoff=0.1, max_cutoff=0.2, umap='umap')
{
    obj <- SeuratAddScore(obj=obj, assay=assay, name=name, features=features)

    df <- obj@meta.data[,c('orig.ident', 'seurat_clusters', paste0(name, '1'))]
    df$signal <- 'Middle'
    df$signal[ df[[paste0(name, '1')]] > max_cutoff ] <- 'High'
    df$signal[ df[[paste0(name, '1')]] < min_cutoff ] <- 'Low'

    # umap
    df$UMAP_1 <- obj@reductions[[umap]]@cell.embeddings[,1]
    df$UMAP_2 <- obj@reductions[[umap]]@cell.embeddings[,2]

    return(df)
}




#' Find markers for all group
#'
#' @param obj object
#' @param assay assay
#' @return data frame
#' @export
#'
RunFindAllMarkers <- function(
    obj=NULL, 
    assay='SCT', 
    idents='seurat_clusters', 
    logfc=.25, 
    min_pct=0.1, 
    min_diff_pct=0, 
    recor_umi=TRUE
){ 
    Idents(obj) <- idents

    diff_markers <- Seurat::FindAllMarkers(
        object = obj,
        only.pos = TRUE,
        logfc.threshold = logfc,
        min.pct = min_pct,
        min.diff.pct = min_diff_pct,
        recorrect_umi = recor_umi
    )

    diff_markers$diff_pct <- diff_markers$pct.1 - diff_markers$pct.2

    # filter genes not in SCT
    all_genes <- rownames(obj@assays[[assay]]@data)
    # wrong gene name from duplicated gene markers
    wrong_gene_name <- setdiff(rownames(diff_markers), all_genes)

    diff_markers$gene <- rownames(diff_markers)
    diff_markers$gene[diff_markers$gene %in% wrong_gene_name] <- str_sub(diff_markers$gene[diff_markers$gene %in% wrong_gene_name], end = -2)

    return(diff_markers)
}





#' Find markers for two group
#'
#' @param obj object
#' @param assay assay
#' @return data frame
#' @export
#'
RunFindMarkers <- function(
    obj=NULL, 
    assay='SCT', 
    idents='seurat_clusters',
    logfc=0.25, 
    only_pos=TRUE, 
    ct_target=NULL, 
    ct_bg=NULL, 
    recor_umi=TRUE
){
    diff_markers <- FindMarkers(
        object = obj,
        assay = assay,
        group.by = idents,
        ident.1 = ct_target,
        ident.2 = ct_bg,
        only.pos = only_pos,
        logfc.threshold = logfc,
        recorrect_umi = recor_umi
    )
    #logfc.threshold = 0.25,

    diff_markers$diff_pct = diff_markers$pct.1 - diff_markers$pct.2

    return(diff_markers)
}






#' Find Conserved CellType Markers
#'
#' @param obj object
#' @param assay assay
#' @return data frame
#' @export
#'
RunFindConservedCellTypeMarkers <- function(obj, assay='SCT', group='cell_type2')
{ 
    Idents(obj) <- group
    meta <- obj@meta.data
    celltype_set <- unique(meta[[group]])

    d_markers <- NULL
    i = 1
    for (x in celltype_set){
        temp_markers <- FindConservedMarkers(obj, ident.1 = x, grouping.var = "orig.ident", assay=assay, verbose = FALSE)
        temp_markers[[group]] <- x
        
        if (i == 1){
            d_markers <- temp_markers
            i = 2
        } else {
            d_markers <- rbind(d_markers, temp_markers)
        }
    }
    
    return(d_markers)
}






