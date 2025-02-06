


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





#' Extract multiple motif set signal
#'
#' @param obj data
#' @param motifs ids
#' @param group add group name
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

        if (length(tar_motifs) == 1){
            d_motif <- as.data.frame(seu@assays$chromvar@data[tar_motifs,])
            colnames(d_motif) <- tar_motifs
        } else {
            d_motif <- data.frame(t(obj@assays$chromvar@data[tar_motifs,]))
        }
    }

    metadata <- obj@meta.data
    d_motif$orig.ident <- metadata$orig.ident[match(rownames(metadata), rownames(d_motif))]
    #                               MA1615.1    MA1548.1    MA0163.1   orig.ident
    # H_02_2138_AAACAGCCAAGTGAAC.1  0.9646257  0.03137187  0.23462118   x        

    if (group != ''){
        d_motif[group] <- metadata[[group]][match(rownames(metadata), rownames(d_motif))]
        #                               MA1615.1    MA1548.1    MA0163.1   orig.ident  <group>
        # H_02_2138_AAACAGCCAAGTGAAC.1  0.9646257  0.03137187  0.23462118   x             x
    }
    

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
CalMeanByGroup <- function(df, group='', with_sample=FALSE)
{ 
    # remove group name
    col_name <- colnames(df)
    col_name <- col_name[ !col_name %in% c('orig.ident', group)]

    # call means
    if (with_sample){
        df <- df %>% group_by(orig.ident, .data[[group]]) %>% summarise(across(all_of(col_name), mean))
    } else {
        df <- df %>% group_by(.data[[group]]) %>% summarise(across(all_of(col_name), mean))
    }
    

    df <- as.data.frame(df)

    return(df)
}




#' Call median
#'
#' @param df data
#' @param group name
#' @return data frame
#' @import dplyr
#' @export
#'
CalMedianByGroup <- function(df, group='', with_sample=FALSE)
{ 
    # remove group name
    col_name <- colnames(df)
    col_name <- col_name[ !col_name %in% c('orig.ident', group)]

    # call median
    if (with_sample){
        df <- df %>% group_by(orig.ident, .data[[group]]) %>% summarise(across(all_of(col_name), median))
    } else {
        df <- df %>% group_by(.data[[group]]) %>% summarise(across(all_of(col_name), median))
    }

    df <- as.data.frame(df)

    return(df)
}




#' Calculate mean motif signal
#'
#' @param obj anno
#' @param motif ids
#' @param group name
#' @return data frame
#' @export
#'
CalMotifMeanOrMedianSig <- function(obj, motifs='all', group='', method='mean', with_sample=FALSE, add_cellcount=FALSE)
{   
    d_motif <- ExtractMotifSigTranspose(obj, motifs, group)
    # <motif> ... orig.ident  <group>

    # calculare
    d_score <- NULL
    # mean
    if (method == 'mean'){
        d_score <- CalMeanByGroup(d_motif, group, with_sample)
    }
    if (method == 'median'){
        d_score <- CalMedianByGroup(d_motif, group, with_sample)
    }
    
    if (add_cellcount){
        if (with_sample){
            dcount <- as.data.frame(table(d_motif[,group]))
            colnames(dcount) <- c(group, 'cell_count')
            d_score <- merge(d_score, dcount, by=group)
            # <motif> ...  <group>  cell_count
        } else {
            dcount <- as.data.frame(table(d_motif[,c('orig.ident', group)]))
            colnames(dcount) <- c('orig.ident', group, 'cell_count')
            d_score <- merge(d_score, dcount, by=c('orig.ident', group))
            # <motif> ... orig.ident  <group>  cell_count
        }
    }
    
    return(d_score)
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




#' Calculate cell count/ratio for motif signal
#'
#' @param obj anno
#' @param motif ids
#' @param group name
#' @param sample_group name
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
    d_ratio <- merge(d_ratio, d_long, by=c(sample_group, group), all.x=T)

    # calculare ratio
    d_ratio$ct_ratio <- round(d_ratio$ct_cells / d_ratio$sample_cells * 100, 2)
    d_ratio$motifSig_ratio <- round(d_ratio$motifSig_count / d_ratio$sample_cells * 100, 2)
    d_ratio$motifSig_ratio_in_ct <- round(d_ratio$motifSig_count / d_ratio$ct_cells * 100, 2)
    #cell_type2  orig.ident  ct_cells sample_cells motifs motifSig_count motifSig_ratio  ct_motifSig_ratio ct_ratio
    
    return(d_ratio)
}




#' Motif SetSignalGroup
#'
#' @param obj object
#' @param motif id
#' @param group name 
#' @param min_cutoff min signal cutoff
#' @param max_cutoff max signal cutoff
#' @return data frame
#' @export
#'
Motif_SetSignalGroup <- function(obj, motif=NULL, group='cell_type2', min_cutoff=0.5, max_cutoff=1)
{ 
    d_motif <- ExtractMotifSigTranspose(obj=obj, motifs=motif, group=group)
    d_motif$signal <- 'Middle'
    d_motif$signal[ d_motif[[motif]] > max_cutoff ] <- 'High'
    d_motif$signal[ d_motif[[motif]] < min_cutoff ] <- 'Low'
    d_motif$signal[ d_motif[[motif]] <= 0 ] <- 'No'

    return(d_motif)
}





#' Motif find markers
#'
#' @param obj object
#' @param query name
#' @param avg_fc cutoff avg_diff
#' @param pct_1 cutoff pct.1
#' @return data frame
#' @export
#'
Motif_FindMarkers <- function(obj, query='', avg_fc=.585, pct_1=0.2, diff_pct=0.1)
{ 
    DefaultAssay(obj) <- 'chromvar'

    all_items <- levels(obj)
    bkg <- all_items[!all_items %in% query]

    d_markers <- Seurat::FindMarkers(
      object = obj,
      ident.1 = query,
      ident.2 = bkg,
      only.pos = TRUE,
      mean.fxn = rowMeans,
      fc.name = "avg_diff"
    )

    d_markers$diff_pct <- d_markers$pct.1 - d_markers$pct.2
    d_filtered <- dplyr::filter(d_markers, avg_diff > avg_fc & p_val_adj<0.005 & pct.1 > pct_1 & diff_pct > diff_pct)

    return(d_filtered)
}




#' Motif find all markers
#'
#' @param obj object
#' @param avg_fc cutoff avg_diff
#' @param pct_1 cutoff pct.1
#' @return data frame
#' @export
#'
Motif_FindAllMarkers <- function(obj, avg_fc=.585, pct_1=0.2, diff_pct=0.1)
{ 
    DefaultAssay(obj) <- 'chromvar'

    d_markers <- Seurat::FindAllMarkers(
      object = obj,
      only.pos = TRUE,
      mean.fxn = rowMeans,
      fc.name = "avg_diff"
    )

    d_markers$diff_pct <- d_markers$pct.1 - d_markers$pct.2
    d_filtered <- dplyr::filter(d_markers, avg_diff > avg_fc & p_val_adj<0.005 & pct.1 > pct_1 & diff_pct > diff_pct)

    return(d_filtered)
}





