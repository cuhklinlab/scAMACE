\name{App2_Feb9_2021_new_three_data_5subcell_Astro_CPN_L6_mPv_oligo_complicated_rna_vtilde_all_link_real_3data_ini_square}
\alias{App2_Feb9_2021_new_three_data_5subcell_Astro_CPN_L6_mPv_oligo_complicated_rna_vtilde_all_link_real_3data_ini_square}
\docType{data}
\title{Application 2: Mouse neocortex dataset initialization data}
\description{
Application 2: Mouse neocortex dataset initialization dataset, including initialization data for scCAS data, scRNA-Seq data and sc-methylation data.
}
\usage{data("App2_Feb9_2021_new_three_data_5subcell_Astro_CPN_L6_mPv_oligo_complicated_rna_vtilde_all_link_real_3data_ini_square")}
\format{
  large list (15 elements, 121.1 Mb)
  \describe{
    \item{f10ratio_acc}{f1/f0 for scCAS data.}
    \item{f10ratio_exp}{g1/g0 for scRNA-Seq data.}
    \item{f10ratio_met}{h1/h0 for sc-methylation data.}

    \item{omega_acc}{the initialization for omega_acc in scCAS data.}
    \item{qi_acc}{the initialization of the probability that scCAS data has high gene score when gene g is active in cell i (i.e. pi_{i1}).}

    \item{omega_exp}{the initialization for omega_rna in scRNA-Seq data.}
    \item{pi_rna}{the initialization of the probability for scRNA-Seq data that gene g is expressed in cell l.}

    \item{omega_met}{the initialization for omega_met in sc-methylation data.}
    \item{pi_met}{the initialization of the probability for sc-methylation data that gene g is methylated in cell d.}

    \item{z_acc}{the cluster assignment matrix z_acc in scCAS data.}
    \item{z_exp}{the cluster assignment matrix z_rna.}
    \item{z_met}{the cluster assignment matrix z_met.}

    \item{u_acc_tilda}{gene activity score level o in scCAS data.}
    \item{v_exp_tilda}{gene expression status v in scRNA-Seq data.}
    \item{v_met}{gene methylated status m in sc-methylation data.}



}


\keyword{datasets}
