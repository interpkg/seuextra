




#' Extract motif signal
#'
#' @param obj data
#' @param motifs ids
#' 
#' @return data frame
#' @export
#'
ExtractMotifSig <- function(obj, motifs='all')
{
    d_motif <- NULL

    if (motifs == 'all'){
        d_motif <- data.frame(obj@assays$chromvar@data)
    } else {
        tar_motifs <- stringr::str_split(motifs, ',')[[1]]
        d_motif <- data.frame(obj@assays$chromvar@data[tar_motifs,])
    }

    #           H_02_2138_AAACAGCCAAGTGAAC.1   H_02_2138_AAACAGCCAAGTGATT.1 
    # MA1615.1   0.9646257  0.23462118
    # MA1548.1   0.0100001  0.03137187

    return(d_motif)
}





#' Extract motif signal
#'
#' @param obj data
#' @param motifs ids
#' @param group name
#' 
#' @return data frame
#' @export
#'
ExtractMotifSigTranspose <- function(obj, motifs='all', group='')
{
    d_motif <- NULL

    if (motifs == 'all'){
        d_motif <- data.frame(t(obj@assays$chromvar@data))
    } else {
        tar_motifs <- stringr::str_split(motifs, ',')[[1]]
        d_motif <- data.frame(t(obj@assays$chromvar@data[tar_motifs,]))
    }

    metadata <- obj@meta.data
    d_motif$orig.ident <- metadata$orig.ident[match(rownames(metadata), rownames(d_motif))]
    d_motif[group] <- metadata[[group]][match(rownames(metadata), rownames(d_motif))]

    #                               MA1615.1    MA1548.1    MA0163.1   orig.ident  <group>
    # H_02_2138_AAACAGCCAAGTGAAC.1  0.9646257  0.03137187  0.23462118   x             x

    return(d_motif)
}






#' Call mean
#'
#' @param df data
#' @param group name
#' @return data frame
#' @import dplyr
#' @export
#'
CallMeanByGroup <- function(df, group='')
{ 
    # remove group name
    col_name <- colnames(df)
    col_name <- col_name[ !col_name %in% c('orig.ident', group)]

    # call means
    df <- df %>% 
            group_by(.data[[group]]) %>%
            summarise(across(all_of(col_name), mean))

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
    d_motif <- ExtractMotifSigTranspose(obj, motifs, group)

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
            dsub_score$orig.ident <- s
            dsub_score <- CallMeanByGroup(dsub, group)
            
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





#' Call count/ratio
#'
#' @param df data
#' @param group name
#' @param cutoff signal cutoff
#' @return data frame
#' @import dplyr
#' @export
#'
CallSignalCountRatioByGroup <- function(df, columns, group='', subgroup='', cutoff=0)
{   
    if (group != '' & subgroup != ''){
        df <- df %>% 
            group_by(.data[[group]], .data[[subgroup]]) %>%
            summarise(across(all_of(columns), ~ sum(.x > cutoff )))
    } else {
        df <- df %>% 
            group_by(.data[[group]]) %>%
            summarise(across(all_of(columns), ~ sum(.x > cutoff )))
    }
    
    df <- as.data.frame(df)

    return(df)
}




#' Calculate cell count/ratio for motif signal
#'
#' @param object anno
#' @param motif ids
#' @param group name
#' @param subgroup name
#' @param split by sample name yes/no
#' @return data frame
#' @export
#'
CalCellRatioForMotifSig <- function(obj, motifs='all', group='', subgroup='', split=FALSE)
{   
    d_motif <- ExtractMotifSigTranspose(obj, motifs, group)

    if (motifs != 'all'){
        motifs <- stringr::str_split(motifs, ',')[[1]]
    }

    # calculare
    d_final <- NULL

    if (split){
        print(table(metadata$orig.ident))

        sampleset <- unique(metadata$orig.ident)
        i <- 1
        for (s in sampleset){
            cells <- rownames(metadata[metadata$orig.ident==s,])
            dsub <- d_motif[cells,]
            print(paste0(s, ':', nrow(dsub)))
            dsub_score$orig.ident <- s

            dsub_score <- CallSignalCountRatioByGroup(df=dsub, columns=motifs, group=group, subgroup=subgroup, cutoff=0)

            if (i == 1){
                d_final <- dsub_score
                i = 2
            } else {
                d_final <- rbind(d_final, dsub_score)
            }
        }   
    } else {
        d_final <- CallSignalCountRatioByGroup(df=d_motif, columns=motifs, group=group, subgroup=subgroup, cutoff=0)
    }
    
    return(d_final)
}







