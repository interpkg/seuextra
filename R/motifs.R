# 2024-12-27

library(dplyr)
options(scipen = 999)



#' Annotations
#'
#' @param object anno
#' @return anotation
#' @export
#'
CalMeanMotifSig <- function(seu, motifs='all', group='')
{   
    d_motif <- NULL

    if (motifs == 'all'){
        d_motif <- data.frame(t(seu@assays$chromvar@data))
    } else {
        tar_motifs <- stringr::str_split(motifs, ',')[[1]]
        d_motif <- data.frame(t(seu@assays$chromvar@data[tar_motifs,]))
    }
    #                               MA1615.1    MA1548.1    MA0163.1
    # H_02_2138_AAACAGCCAAGTGAAC.1  0.9646257  0.03137187  0.23462118

    metadata <- seu@meta.data
    d_motif[group] <- metadata[[group]][match(rownames(metadata), rownames(d_motif))]

    d_motif <- d_motif %>% 
            group_by(get(group)) %>%
            summarise(across(everything(), mean))

    return(d_motif)
}





