


#' Call mean
#'
#' @param data frame
#' @param group name
#' @return data frame
#' @import dplyr
#' @export
#'
CallMeanByGroup <- function(df, group='')
{ 
    # remove group name
    col_name <- colnames(df)
    col_name <- col_name[ !col_name == group]

    # call means
    df <- df %>% 
            group_by(.data[[group]]) %>%
            summarise(across(col_name, mean))

    df <- as.data.frame(df)

    return(df)
}




#' Calculate mean motif signal
#'
#' @param object anno
#' @param motif ids
#' @param group name
#' @param split by sample name yes/no
#' @return data frame
#' @export
#'
CalMeanMotifSig <- function(obj, motifs='all', group='', split=FALSE)
{   
    d_motif <- NULL

    if (motifs == 'all'){
        d_motif <- data.frame(t(obj@assays$chromvar@data))
    } else {
        tar_motifs <- stringr::str_split(motifs, ',')[[1]]
        d_motif <- data.frame(t(obj@assays$chromvar@data[tar_motifs,]))
    }
    #                               MA1615.1    MA1548.1    MA0163.1
    # H_02_2138_AAACAGCCAAGTGAAC.1  0.9646257  0.03137187  0.23462118

    metadata <- obj@meta.data
    d_motif[group] <- metadata[[group]][match(rownames(metadata), rownames(d_motif))]

    # calculare
    d_avg_score <- NULL

    if (split){
        print(table(metadata$orig.ident))

        sampleset <- unique(metadata$orig.ident)
        i <- 1
        for (s in sampleset){
            cells <- rownames(metadata[metadata$orig.ident==s,])
            dsub <- d_motif[cells,]
            print(paste0(s, ':', nrow(dsub)))
            dsub_score <- CallMeanByGroup(dsub, group)
            dsub_score$orig.ident <- s

            if (i == 1){
                d_avg_score <- dsub_score
                i = 2
            } else {
                d_avg_score <- rbind(d_avg_score, dsub_score)
            }
        }   
    } else {
        d_avg_score <- CallMeanByGroup(d_motif, group)
    }
    
    return(d_avg_score)
}




