# 2024-12-27

library(dplyr)


#' Export UMAP
#'
#' @param reference anno
#' @return anotation
#' @export
#'
ExportUMAP <- function(object, umap)
{
    object$UMAP_1 <- object@reductions[[umap]]@cell.embeddings[,1]
    object$UMAP_2 <- object@reductions[[umap]]@cell.embeddings[,2]
    d_umap <- as.data.frame(object@meta.data[,c('UMAP_1', 'UMAP_2')])

    return(d_umap)
}




