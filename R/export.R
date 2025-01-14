

#' Export UMAP
#'
#' @param reference anno
#' @return anotation
#' @export
#'
#' @examples
#' d_umap <- ExportUMAP(seu, 'wnn.umap')
#'
ExportUMAP <- function(seu_obj, umap)
{
    seu_obj$UMAP_1 <- seu_obj@reductions[[umap]]@cell.embeddings[,1]
    seu_obj$UMAP_2 <- seu_obj@reductions[[umap]]@cell.embeddings[,2]
    d_umap <- as.data.frame(seu_obj@meta.data[,c('UMAP_1', 'UMAP_2')])

    return(d_umap)
}




