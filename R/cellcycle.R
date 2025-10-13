#' @include utilities.R
#'
NULL




#' Seurat_EstimateCycling
#'
#' @param obj seurat object
#' @param ref reference
#' @param s.genes gene vectors 
#' @param g2m.genes gene vectors 
#'
#' @return data frame
#'
#' @export
#'
SeuratCellCycleScoring <- function(obj=NULL, ref='mouse', s.genes=NULL, g2m.genes=NULL, umap='umap')
{
    # CellCycleScoring
    obj <- Seurat::CellCycleScoring(
        object = obj,
        s.features = s.genes,
        g2m.features = g2m.genes
    )

    # CellCycleScoring2 for all cell cycle genes
    cc.genes <- c(s.genes, g2m.genes)
    obj <- SeuratAddScore(obj=obj, name='cc_score', features=cc.genes, assay='SCT')

    df <- obj@meta.data[,c('cc_score1', 'S.Score', 'G2M.Score', 'Phase', 'seurat_clusters')]
    
    # define cycling
    df$signal <- '0.1-0.2'
    df$signal[df$cc_score1 > 0.2 & (df$S.Score > 0.2 | df$G2M.Score > 0.2) ] <- '> 0.2'
    df$signal[df$cc_score1 < 0.1 & (df$S.Score < 0.1 | df$G2M.Score < 0.1) ] <- '< 0.1'

    # group
    df$group <- 'Intermediate'
    df$group[df$signal == '> 0.2'] <- 'Cycling'
    df$group[df$signal == '< 0.1'] <- 'Non-cycling'

    # umap
    df$UMAP_1 <- obj@reductions[[umap]]@cell.embeddings[,1]
    df$UMAP_2 <- obj@reductions[[umap]]@cell.embeddings[,2]

    return(df)
}





#' CC_UMAPPlot_Group
#'
#' @param df dataframe
#'
#' @export
#'
CC_UMAPPlot_Group <- function(df=NULL, label_name='Cell-Cycle\nSignal'){
    df <- df[order(df$group, decreasing=TRUE), ]

    p <- ggplot(df,aes(x=UMAP_1, y=UMAP_2, color=group)) + 
            geom_point(size=0.1) +
            theme_classic(base_line_size=0.2) +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank()) +
            labs(x = "UMAP 1", y="UMAP 2") +
            scale_color_manual(name = label_name, breaks = c("Cycling", "Intermediate", "Non-cycling"), values=c("#CB4335", "#F5B7B1", "#D7DBDD")) +
            guides(color = guide_legend(override.aes = list(size = 2))) +
            theme(legend.title = element_text(size=8)) +
            theme(text=element_text(size=8))

    p
}



#' CC_UMAPPlot_Signal
#'
#' @param df dataframe
#'
#' @export
#'
CC_UMAPPlot_Signal <- function(df=NULL, label_name='Cell-Cycle\nSignal'){
    df <- df[order(df$group, decreasing=TRUE), ]
    
    p <- ggplot(df,aes(x=UMAP_1, y=UMAP_2, color=signal)) + 
            geom_point(size=0.1) +
            theme_classic(base_line_size=0.2) +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank()) +
            labs(x = "UMAP 1", y="UMAP 2") +
            scale_color_manual(name = label_name, breaks = c("> 0.2", "0.1-0.2", "< 0.1"), values=c("#CB4335", "#2E86C1", "#D7DBDD")) +
            guides(color = guide_legend(override.aes = list(size = 2))) +
            theme(legend.title = element_text(size=8)) +
            theme(text=element_text(size=8))

    p
}





#' CC_SignalUMAP_perSample
#'
#' @param df dataframe
#'
#' @export
#'
CC_SignalUMAP_perSample <- function(df, sorted_sample='', label_name='Cell-Cycle\nSignal')
{
    df <- df[order(df$signal, decreasing=FALSE), ]
    

    if (nchar(sorted_sample) > 2){ 
        df$Sample <- factor(df$Sample, levels=sorted_sample) 
    }
    
    p <- ggplot(df,aes(x=UMAP_1, y=UMAP_2, color=group)) + 
            geom_point(size=0.2) +
            theme_classic(base_line_size=0.2) +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank()) +
            labs(x = "UMAP 1", y="UMAP 2") +
            scale_color_manual(name = label_name, breaks = c("> 0.2", "0.1-0.2", "< 0.1"), values=c("#CB4335", "#2E86C1", "#D7DBDD")) +
            guides(color = guide_legend(override.aes = list(size = 2))) +
            theme(legend.title = element_text(size=9))

    # per sample
    p <- p + facet_wrap(~Sample, ncol = 5) + 
        theme(strip.background = element_blank(), strip.text = element_text(size=10))
    
    p
}




#' CC_Barplot_nonCyclining
#'
#' @param df dataframe
#'
#' @export
#'
CC_Barplot_nonCyclining <- function(df){
    # 1.order cell type
    report <- as.data.frame(table(df$clusters))
    colnames(report) <- c('clusters', 'n')
    order_clusters <- (dplyr::arrange(report, desc(n)))$clusters

    # 2.
    table <- df %>% count(clusters, group, sort = TRUE)
    colnames(table) <- c('clusters', 'group', 'n')
    print(head(table, n=3))
    
    #colors <- c("> 0.2"="#CB4335", "0.1-0.2"="#2E86C1", "< 0.1"="#D7DBDD")
    colors <- c("Cycling"="#CB4335", "Intermediate"="#F5B7B1", "Non-cycling"="#D7DBDD")
    p <- ggplot(table) +
            aes(x = clusters, y=n, fill=group) +
            geom_bar(stat = "identity") + 
            coord_flip() +
            scale_x_discrete(limits = rev(order_clusters)) +
            scale_fill_manual(breaks=c('Cycling', 'Intermediate', 'Non-cycling'), values=colors) +
            theme_classic() + 
            labs(title='') + xlab('') + ylab('The number of cells') +
            theme(
                axis.title = element_text(size = 8),
                text = element_text(size=7),
                axis.text.y = element_text(size = 7),
                axis.title.x = element_text(size = 7),

                legend.title = element_blank(),
                legend.text = element_text(),
                legend.key.size=unit(3,"mm"),
                    
                axis.line.y = element_blank(),
                axis.line.x = element_line(color="grey"),
                axis.ticks = element_line(color="grey"),

                panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(), 
                panel.background = element_blank(), 
                panel.border = element_blank()
            )

    return(p)
}






