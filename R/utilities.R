

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






