

#' Mouse Brain Features
#'
#' @return vector
#'
#' @export
#'
MouseBrainFeatures <- function()
{
    vec <- c('Non-neuronal' = c('Gli3','Aqp4','Gja1','Gfap','Gpr50','Tmem212','Pdgfra','Cldn11','Slc47a1','Apod','Kcnj8','Acta2','Ly6c1','C1qa','Mrc1'),
            'Midbrain-neuron' = c('Pax7','Otx2os1','Tfap2b','Pax5','Gata3'),
            'Hindbrain-neuron' = c('Irx3','Pax7','Gata3','Vsx2'),
            'GABA-neuron' = c('Grpr','Egln3','Mctp2','Rbp4','Arx','Sp9')
            )

    return(vec)
}



#' Mouse Cell Cyle 
#'
#' @return vector
#'
#' @export
#'
MouseCCFeatures2019 <- function()
{
    vec <- c('s' = c('Mcm5','Pcna','Tyms','Fen1','Mcm7','Mcm4','Rrm1','Ung','Gins2','Mcm6','Cdca7','Dtl','Prim1','Uhrf1','Cenpu','Hells','Rfc2','Polr1b','Nasp','Rad51ap1','Gmnn','Wdr76','Slbp','Ccne2','Ubr7','Msh2','Rad51','Rrm2','Cdc45','Cdc6','Exo1','Tipin','Dscc1','Blm','Casp8ap2','Usp1','Clspn','Pola1','Chaf1b','Mrpl36','E2f8'),
            'g2m' = c('Hmgb2','Cdk1','Nusap1','Ube2c','Birc5','Tpx2','Top2a','Ndc80','Cks2','Nuf2','Cks1b','Mki67','Tmpo','Cenpf','Tacc3','Pimreg','Smc4','Ccnb2','Ckap2l','Ckap2','Aurkb','Bub1','Kif11','Anp32e','Tubb4b','Gtse1','Kif20b','Hjurp','Cdca3','Jpt1','Cdc20','Ttk','Cdc25c','Kif2c','Rangap1','Ncapd2','Dlgap5','Cdca2','Cdca8','Ect2','Kif23','Hmmr','Aurka','Psrc1','Anln','Lbr','Ckap5','Cenpe','Ctcf','Nek2','G2e3','Gas2l3','Cbx5','Cenpa')
            )

    return(vec)
}





