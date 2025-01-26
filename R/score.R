

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





#' Seurat_EstimateCycling
#'
#' @param obj seurat object
#' @param ref reference
#' @param s.genes gene vectors 
#' @param g2m.genes gene vectors 
#'
#' @return data frame
#'
#' @export
#'
SeuratCellCycleScoring <- function(obj=NULL, ref='mouse', s.genes=NULL, g2m.genes=NULL, umap='umap')
{
    # CellCycleScoring
    obj <- Seurat::CellCycleScoring(
        object = obj,
        s.features = s.genes,
        g2m.features = g2m.genes
    )

    # CellCycleScoring2 for all cell cycle genes
    cc.genes <- c(s.genes, g2m.genes)
    obj <- SeuratAddScore(obj=obj, name='cc_score', features=cc.genes, assay_use='SCT')

    df <- obj@meta.data[,c('cc_score1', 'S.Score', 'G2M.Score', 'Phase', 'seurat_clusters')]
    
    # define cycling
    df$signal <- '0.1-0.2'
    df$signal[df$cc_score1 > 0.2 & (df$S.Score > 0.2 | df$G2M.Score > 0.2) ] <- '> 0.2'
    df$signal[df$cc_score1 < 0.1 & (df$S.Score < 0.1 | df$G2M.Score < 0.1) ] <- '< 0.1'

    # group
    df$group <- 'Intermediate'
    df$group[df$signal == '> 0.2'] <- 'Cycling'
    df$group[df$signal == '< 0.1'] <- 'Non-cycling'

    # umap
    df$UMAP_1 <- obj@reductions[[umap]]@cell.embeddings[,1]
    df$UMAP_2 <- obj@reductions[[umap]]@cell.embeddings[,2]

    return(df)
}









