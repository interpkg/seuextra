



#' ConvertPeakToBed
#'
#' @param peaks region
#' 
#' @return bed form
#' @export
#'
ConvertPeakToBed <- function(peaks){
      # format to GRanges
      bed <- as.data.frame(stringr::str_split_fixed(peaks, "-", 3))
      colnames(bed) <- c('chr', 'start', 'end')
      # chr start end

      bed
}



#' ConvertPeakToBed4
#'
#' @param peaks region
#' 
#' @return bed form
#' @export
#'
ConvertPeakToBed4 <- function(peaks){
      # format to GRanges
      bed <- as.data.frame(stringr::str_split_fixed(peaks, "-", 3))
      colnames(bed) <- c('chr', 'start', 'end')
      bed$name <- paste0('peak', 1:nrow(bed))
      # # chr start end name

      bed
}



#' ConvertBedToGRange
#'
#' @param bed form
#' 
#' @return gr_region form
#' @export
#'
ConvertBedToGRange <- function(bed){
      gr_region <- GRanges(seqnames=bed$chr, 
                  ranges=IRanges(start=as.numeric(bed$start), end=as.numeric(bed$end)))

      gr_region
}




#' BedtoolsIntersectQuery
#'
#' @param query query regions
#' @param db bed data
#' @param cutoff cutoff overlap
#' 
#' @return data frame
#' @export
#'
BedtoolsIntersectQuery <- function(query=NULL, db=NULL, cutoff=0.5){
      command <- paste0("bedtools intersect -a ", query, " -b ", db, " -f ", cutoff, " -wa -wb | cut -f1-3 | sed 's/\t/-/g' ")
      intersected_bed <- system(command, intern = TRUE)
      
      unique(intersected_bed) 
      # chr-start-end
}






#' CalIntersectPeaks
#'
#' @param peaks bed file
#' @param fbed bed file
#' @param cutoff overlap cutoff
#' @param temp_dir temporary dir
#' 
#' @return data frame
#' @export
#'
CalIntersectPeaks <- function(
      peaks=NULL, 
      fbed='', 
      cutoff=0.5, 
      temp_dir='temp'
){
      dir.create(temp_dir)

      # peaks
      # chr-start-end
      peaks_bed <- ConvertPeakToBed(peaks)
      write.table(peaks_bed, file=paste0(temp_dir, '/temp.query_peaks.bed'), sep='\t', quote=F, row.names=F, col.names=F)
      fquery_bed <- paste0(temp_dir, '/temp.query_peaks.bed')
      
      # query overlap peaks (low cutoff 0.3)
      overlap_peaks <- BedtoolsIntersectQuery(query=fquery_bed, db=fbed, cutoff=cutoff)
      # chr-start-end

      return(overlap_peaks)
}




#' CalIntersectMotifOpenPeaksPerCell
#'
#' @param obj multiome object
#' @param fquery bed file
#' @param fdb bed file
#' @param cutoff overlap cutoff
#' @param count_cutoff count cutoff
#' @param temp_dir temporary dir
#' 
#' @return data frame
#' @export
#'
CalIntersectMotifOpenPeaksPerCell <- function(obj=NULL, motif='', fdb='', cutoff=0.5, count_cutoff=2, temp_dir='temp'){
      dir.create(temp_dir)

      # motif related peaks
      # get motif
      peaks_in_motif <- Signac::GetMotifData(object=obj)[, motif]
      # filter NA
      peaks_in_motif <- names(Filter(any, peaks_in_motif))

      # intesect peaks
      #peaks_bed <- ConvertPeakToBed(peaks_in_motif)
      #write.table(peaks_bed, file=paste0(temp_dir, '/temp.query_peaks.bed'), sep='\t', quote=F, row.names=F, col.names=F)
      #fquery <- paste0(temp_dir, '/temp.query_peaks.bed')
      
      # query overlap peaks (low cutoff 0.3)
      #overlap_peaks <- BedtoolsIntersectQuery(query=fquery, db=fdb, cutoff=cutoff)

      overlap_peaks <- CalIntersectPeaks(peaks=peaks_in_motif, fbed=fdb, cutoff=cutoff)


      # peak count
      mtx_peak_count <- as.matrix(obj@assays$peaks@counts[overlap_peaks,])
      # region cell1 cell2 ...

      # count cutoff 3
      d_open_peak_in_cell <- data.frame(apply(mtx_peak_count, 2, function(col) sum(col > count_cutoff)))
      colnames(d_open_peak_in_cell) <- 'open_peaks'
      # <cell>  open_peaks

      return(d_open_peak_in_cell)
}





#' CalIntersectOpenPeaksPerCell
#'
#' @param obj multiome object
#' @param fquery bed file
#' @param fdb bed file
#' @param cutoff overlap cutoff
#' @param count_cutoff count cutoff
#' @param temp_dir temporary dir
#' 
#' @return data frame
#' @export
#'
CalIntersectOpenPeaksPerCell <- function(
      obj=NULL, 
      fdb='', 
      cutoff=0.5, 
      count_cutoff=2, 
      temp_dir='temp'
){
      dir.create(temp_dir)

      all_peaks <- rownames(obj@assays$peaks)

      # intesect peaks
      #peaks_bed <- ConvertPeakToBed(all_peaks)
      #write.table(peaks_bed, file=paste0(temp_dir, '/temp.all_peaks.bed'), sep='\t', quote=F, row.names=F, col.names=F)
      #fquery <- paste0(temp_dir, '/temp.all_peaks.bed')
      # query overlap peaks (low cutoff 0.3)
      #overlap_peaks <- BedtoolsIntersectQuery(query=fquery, db=fdb, cutoff=cutoff)

      overlap_peaks <- CalIntersectPeaks(peaks=all_peaks, fbed=fdb, cutoff=cutoff)

      # peak count
      mtx_peak_count <- as.matrix(obj@assays$peaks@counts[overlap_peaks,])
      # region cell1 cell2 ...

      # count cutoff 3
      d_open_peak_in_cell <- data.frame(apply(mtx_peak_count, 2, function(col) sum(col > count_cutoff)))
      colnames(d_open_peak_in_cell) <- 'open_peaks'
      # <cell>  open_peaks

      return(d_open_peak_in_cell)
}








#' FragmentDensityForGene
#'
#' @param obj seurat object
#' @param fragment fragment file
#' @param peaks peak name
#' @param gene gene name
#' @param colors color set 
#' 
#' @return data frame
#' @export
#'
FragmentDensityForGene <- function(
      obj=NULL, 
      fragment='atac_fragments.tsv.gz', 
      peaks='peaks',
      gene=NULL, 
      colors=c('RGC'='#BC243C', 'Neuron'='#0F4C81')
) {

    # set new path for fragment
    if (fragment != ''){
        print(fragment)
        fg <- CreateFragmentObject(path = fragment, cells = colnames(obj), validate.fragments = TRUE)
        Fragments(obj@assays[[peaks]]) <- NULL
        Fragments(obj@assays[[peaks]]) <- fg
    }

    p <- CoveragePlot(
            object = obj,
            region = gene,
            peaks = TRUE,
            links = FALSE,
            extend.upstream = 1000,
            extend.downstream = 1000
        ) &
        theme(text = element_text(size = 7, face="bold", color = "black"),
            axis.title.y = element_text(size = 6, face="bold")) &
        scale_fill_manual(values=colors)

    p <- p & plot_layout(ncol=1) & NoLegend() &
            theme(plot.title = element_text(hjust = 0.5, face="bold")) 

    p
}





