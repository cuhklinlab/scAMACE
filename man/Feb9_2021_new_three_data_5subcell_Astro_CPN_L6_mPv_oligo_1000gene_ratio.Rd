\name{Feb9_2021_new_three_data_5subcell_Astro_CPN_L6_mPv_oligo_1000gene_ratio}
\alias{Feb9_2021_new_three_data_5subcell_Astro_CPN_L6_mPv_oligo_1000gene_ratio}
\docType{data}
\title{Application 2: Mouse neocortex dataset}
\description{
Application 2: Mouse neocortex dataset, including scCAS data, scRNA-Seq data and sc-methylation data.
}
\usage{data("Feb9_2021_new_three_data_5subcell_Astro_CPN_L6_mPv_oligo_1000gene_ratio")}
\format{
  large list (11 elements, 171.1 Mb)
  \describe{
    \item{f10ratio_acc}{f1/f0 for scCAS data.}
    \item{f1}{the known pdf f1 of scCAS data when o=1.}
    \item{f0}{the known pdf f0 of scCAS data when o=0.}
    \item{atac_cell}{the true cell labels for scCAS data.}

    \item{f10ratio_exp}{g1/g0 for scRNA-Seq data.}
    \item{g1}{the known pdf g1 of scRNA-Seq data when v=1.}
    \item{g0}{the known pdf g0 of scRNA-Seq data when v=0.}
    \item{rna_cell}{the true cell labels for scRNA-Seq data.}

    \item{met_data}{sc-methylation data, cells by genes.}
    \item{f10ratio_met}{h1/h0 for sc-methylation data.}
    \item{met_cell}{the true cell labels for sc-methylation data.}



}


\keyword{datasets}
