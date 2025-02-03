

#' Mouse Brain Features
#'
#' @return vector
#'
#' @export
#'
MouseBrainFeatures <- function()
{
    vec <- c('Non-neuronal' = 'Gli3,Aqp4,Gja1,Gfap,Gpr50,Tmem212,Pdgfra,Cldn11,Slc47a1,Apod,Kcnj8,Acta2,Ly6c1,C1qa,Mrc1',
            'Midbrain-neuron' = 'Pax7,Otx2os1,Tfap2b,Pax5,Gata3',
            'Hindbrain-neuron' = 'Irx3,Pax7,Gata3,Vsx2',
            'GABA-neuron' = 'Grpr,Egln3,Mctp2,Rbp4,Arx,Sp9'
            )

    return(vec)
}







