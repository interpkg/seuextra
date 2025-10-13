

#' Color Box
#'
#' @return color set
#'
#' @export
#'
colbox <- function(x)
{   
    col <- switch(x,
            'motif' = c('#D3D3D3','#BFC9BD','#22A884FF','#2A788EFF','#414487FF','#440154FF'),
            c("gray")
        )

    return(col)
}




