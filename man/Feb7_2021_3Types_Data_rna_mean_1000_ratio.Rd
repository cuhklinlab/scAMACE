\name{Feb7_2021_3Types_Data_rna_mean_1000_ratio}
\alias{Feb7_2021_3Types_Data_rna_mean_1000_ratio}
\docType{data}
\title{Application 1: K562-GM12878 dataset}
\description{
Application 1: K562-GM12878 dataset, including scCAS data, scRNA-Seq data and sc-methylation data.
}
\usage{data("Feb7_2021_3Types_Data_rna_mean_1000_ratio")}
\format{
  large list (13 elements, 35.1 Mb)
  \describe{
    \item{atac_data}{scCAS data, cells by genes.}
    \item{f10ratio_acc}{f1/f0 for scCAS data.}
    \item{f1}{the known pdf f1 of scCAS data when o=1.}
    \item{f0}{the known pdf f0 of scCAS data when o=0.}
    \item{atac_cell}{the true cell labels for scCAS data, 1: HL60 cell line, 2: K562 cell line.}
    \item{rna_data}{scRNA-Seq data (FPKM/TPM), cells by genes.}
    \item{f10ratio_exp}{g1/g0 for scRNA-Seq data.}
    \item{g1}{the known pdf g1 of scRNA-Seq data when v=1.}
    \item{g0}{the known pdf g0 of scRNA-Seq data when v=0.}
    \item{rna_cell}{the true cell labels for scRNA-Seq data, 1: HL60 cell line, 2: K562 cell line.}

    \item{met_data}{sc-methylation data, cells by genes.}
    \item{f10ratio_met}{h1/h0 for sc-methylation data.}
    \item{met_cell}{the true cell labels for sc-methylation data, 1: HL60 cell line, 2: K562 cell line.}



}


\keyword{datasets}
