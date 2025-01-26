

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





