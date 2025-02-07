



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





#' CalIntersectOpenPeaksPerCell
#'
#' @param obj multiome object
#' @param fquery bed file
#' @param fdb bed file
#' @param cutoff cutoff overlap
#' @param temp_dir temporary dir
#' 
#' @return data frame
#' @export
#'
CalIntersectOpenPeaksPerCell <- function(obj=NULL, assay='peaks', motif='', fdb='', cutoff=0.5, temp_dir='temp'){
      dir.create(temp_dir)

      DefaultAssay(obj) <- assay

      # motif related peaks
      # get motif
      peaks_in_motif <- Signac::GetMotifData(object=obj)[, motif]
      # filter NA
      peaks_in_motif <- names(Filter(any, peaks_in_motif))

      # intesect peaks
      peaks_bed <- ConvertPeakToBed(peaks_in_motif)
      write.table(peaks_bed, file=paste0(temp_dir, '/temp.query_peaks.bed'), sep='\t', quote=F, row.names=F, col.names=F)
      fquery <- paste0(temp_dir, '/temp.query_peaks.bed')
      
      # query overlap peaks (low cutoff 0.3)
      overlap_peaks <- BedtoolsIntersectQuery(query=fquery, db=fdb, cutoff=cutoff)

      # peak count
      mtx_peak_count <- as.matrix(obj@assays$peaks@counts[overlap_peaks,])
      # region cell1 cell2 ...

      # count cutoff 3
      d_open_peak_in_cell <- data.frame(apply(mtx_peak_count, 2, function(col) sum(col > 2)))
      colnames(d_open_peak_in_cell) <- 'open_peaks'
      # <cell>  open_peaks

      return(d_open_peak_in_cell)
}






