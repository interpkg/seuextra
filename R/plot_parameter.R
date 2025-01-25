


#' ggplot2 theme setting
#'
#' @return plot
#'
#' @import ggplot2
#'
#' @export
#'
CustThemeOption1 <- function()
{ 
    p <- theme(plot.title = element_text(size = 8),
            text=element_text(size=6), 
            axis.text=element_text(size=5),
            axis.line = element_line(colour = 'black', size = 0.3),
            axis.ticks = element_line(linewidth = 0.3),
            axis.ticks.length=unit(1, "mm"),
            legend.key.width = unit(3, 'mm'),
            legend.key.height = unit(3, 'mm')
        )
        
    p
}










