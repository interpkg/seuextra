

#' Export UMAP
#'
#' @param reference anno
#' @return anotation
#' @export
#'
#' @examples
#' d_umap <- ExportUMAP(seu, 'wnn.umap')
#'
ExportUMAP <- function(object, umap)
{
    object$UMAP_1 <- object@reductions[[umap]]@cell.embeddings[,1]
    object$UMAP_2 <- object@reductions[[umap]]@cell.embeddings[,2]
    d_umap <- as.data.frame(object@meta.data[,c('UMAP_1', 'UMAP_2')])

    return(d_umap)
}




