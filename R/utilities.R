

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




#' Add UMAP
#'
#' @param obj seurat object
#' @param umap str
#'
#' @return data frame
#'
#' @export
#'
AddUMAP <- function(obj=NULL, df=NULL, umap=NULL)
{
    obj$UMAP_1 <- obj@reductions[[umap]]@cell.embeddings[,1]
    obj$UMAP_2 <- obj@reductions[[umap]]@cell.embeddings[,2]
    d_umap <- as.data.frame(obj@meta.data[,c('UMAP_1', 'UMAP_2')])

    d_merged <- merge(df, d_umap, by='row.names', all=T)

    return(d_merged)
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
SeuratAddScore <- function(obj=NULL, assay_use='SCT', name='', features=''){
    all_gene <- rownames(obj[[assay_use]]@data)
    gene_intersect <- intersect(features, all_gene)

    obj <- AddModuleScore(
        object = obj,
        assay = assay_use,
        slot = "data",
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
EstimatedSignal <- function(obj=NULL, assay_use='SCT', name='', features='', min_cutoff=0.1, max_cutoff=0.2, umap='umap')
{
    obj <- SeuratAddScore(obj=obj, assay_use=assay_use, name=name, features=features)

    df <- obj@meta.data[,c('orig.ident', 'seurat_clusters', paste0(name, '1'))]
    df$signal <- 'Middle'
    df$signal[ df[[paste0(name, '1')]] > max_cutoff ] <- 'High'
    df$signal[ df[[paste0(name, '1')]] < min_cutoff ] <- 'Low'

    df$group <- 2
    df$group[df$signal == 'High'] <- 1
    df$group[df$signal == 'Low'] <- 3

    # umap
    df$UMAP_1 <- obj@reductions[[umap]]@cell.embeddings[,1]
    df$UMAP_2 <- obj@reductions[[umap]]@cell.embeddings[,2]

    return(df)
}




#' Gene expression markers for all group
#'
#' @param obj object
#' @param assay_use assay
#' @return data frame
#' @export
#'
GenxFindAllMarkers <- function(obj, assay_use='SCT', idents='seurat_clusters', logfc=.25, min_pct=0.05, min_diff_pct=0, recorrect_umi=TRUE)
{ 
    Idents(obj) <- idents

    diff_exp <- Seurat::FindAllMarkers(
        object = obj,
        only.pos = TRUE,
        logfc.threshold = logfc,
        min.pct = min_pct,
        min.diff.pct = min_diff_pct,
        recorrect_umi = recorrect_umi
    )

    diff_exp$diff_pct <- diff_exp$pct.1 - diff_exp$pct.2

    # filter genes not in SCT
    all_genes <- rownames(obj@assays[[assay_use]]@data)
    # wrong gene name from duplicated gene markers
    wrong_gene_name <- setdiff(rownames(diff_exp), all_genes)

    diff_exp$gene <- rownames(diff_exp)
    diff_exp$gene[diff_exp$gene %in% wrong_gene_name] <- str_sub(diff_exp$gene[diff_exp$gene %in% wrong_gene_name], end = -2)

    return(diff_exp)
}





#' Genx Conserved CellType Markers
#'
#' @param obj object
#' @param assay_use assay
#' @return data frame
#' @export
#'
GenxConservedCellTypeMarkers <- function(obj, assay_use='SCT', group='cell_type2')
{ 
    Idents(obj) <- group
    meta <- obj@meta.data
    celltype_set <- unique(meta[[group]])

    d_markers <- NULL
    i = 1
    for (x in celltype_set){
        temp_markers <- FindConservedMarkers(obj, ident.1 = x, grouping.var = "orig.ident", assay=assay_use, verbose = FALSE)
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






