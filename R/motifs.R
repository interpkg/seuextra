


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
CalMeanByGroup <- function(df, group='')
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
CalMotifMeanSig <- function(obj, motifs='all', group='')
{   
    d_motif <- ExtractMotifSigTranspose(obj, motifs, group)

    # calculare
    d_avg_score <- CalMeanByGroup(d_motif, group)
    
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
CalSignalCountByGroup <- function(df, columns, group='', subgroup='', cutoff=0)
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




#' Test Calculate cell count/ratio for motif signal
#'
#' @param object anno
#' @param motif ids
#' @param group name
#' @param subgroup name
#' @return data frame
#' @export
#'
CalCellRatioForMotifSig_test <- function(obj, motifs='all', group='', sample_group='', cutoff=0)
{   
    d_motif <- ExtractMotifSigTranspose(obj, motifs, group)

    if (motifs != 'all'){
        motifs <- stringr::str_split(motifs, ',')[[1]]
    }

    # calculare count
    d_final <- CalSignalCountByGroup(df=d_motif, columns=motifs, group=group, subgroup=sample_group, cutoff=cutoff)
    #cell_type2  orig.ident  MA1615.1    MA1548.1
    #CycProg-Like    M_EPN_IUE0  114 112
    #CycProg-Like    M_EPN_IUE1  193 198

    # calculare ratio
    d_ratio <- d_final[, c(group, sample_group)]
    for (x in motifs){
        d_temp <- d_final %>%
                group_by(.data[[sample_group]]) %>%
                mutate(percentage = round(.data[[x]] / sum(.data[[x]]) * 100, 2))
        d_ratio[[x]] <- d_temp$percentage
    }
    
    return(d_ratio)
}






#' Calculate cell count/ratio for motif signal
#'
#' @param object anno
#' @param motif ids
#' @param group name
#' @param subgroup name
#' @return data frame
#' @export
#'
CalCellRatioForMotifSig <- function(obj, motifs='all', group='', sample_group='', cutoff=0)
{   
    d_motif <- ExtractMotifSigTranspose(obj, motifs, group)

    if (motifs != 'all'){
        motifs <- stringr::str_split(motifs, ',')[[1]]
    }

    # calculare count
    d_count <- CalSignalCountByGroup(df=d_motif, columns=motifs, group=group, subgroup=sample_group, cutoff=cutoff)
    #cell_type2  orig.ident  MA1615.1    MA1548.1
    #CycProg-Like    M_EPN_IUE0  114 112
    #CycProg-Like    M_EPN_IUE1  193 198

    # wide to long
    d_long <- d_count %>% tidyr::pivot_longer(cols=motifs, names_to='motifs', values_to='motifSig_count')
    #cell_type2  orig.ident  motifs    motifSig_count

    # total cell count of cell typ
    d_totalCellTypeCount <- data.frame(table(obj@meta.data[,c(sample_group, group)]))
    colnames(d_totalCellTypeCount) <- c(sample_group, group, 'ct_cells')

    d_totalCells <- data.frame(table(obj@meta.data[[sample_group]]))
    colnames(d_totalCells) <- c(sample_group, 'sample_cells')

    d_ratio <- merge(d_totalCellTypeCount, d_totalCells, by=sample_group, all.x=T)
    d_ratio <- merge(d_totalCellTypeCount, d_long, by=c(sample_group, group), all.x=T)

    # calculare ratio
    d_ratio$motifSig_ratio <- round(d_ratio$motifSig_count / d_ratio$sample_cells * 100, 2)
    d_ratio$ct_motifSig_ratio <- round(d_ratio$motifSig_count / d_ratio$ct_cells * 100, 2)
    d_ratio$ct_ratio <- round(d_ratio$ct_cells / d_ratio$sample_cells * 100, 2)
    #cell_type2  orig.ident  ct_cells sample_cells motifs motifSig_count motifSig_ratio  ct_motifSig_ratio ct_ratio
    
    return(d_ratio)
}







