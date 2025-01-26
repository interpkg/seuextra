

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




