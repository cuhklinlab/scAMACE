\name{cal_post}
\alias{cal_post}
\title{cal_post}
\description{
Calculate the posterior probability for one iteration in the EM algorithm.
}
\usage{
cal_post(phi_atac_new,w_acc_new,qi_new,phi_rna_new,w_exp_new,pi_rna_new,phi_met_new,w_met_new,pi_met_new)
}
\arguments{
  \item{phi_atac}{the probability for scCAS data belonging to cluster k (i.e. z[i,k]=1).}

  \item{w_acc}{the probability for scCAS data that gene g is active in cluster k.}

  \item{qi}{the probability that scCAS data has high gene score when gene g is active in cell i (i.e. pi_{i1}).}


  \item{phi_rna}{the probability for scRNA-Seq data belonging to cluster k (i.e. z[l,k]=1).}

  \item{w_exp}{the probability for scRNA-Seq data that gene g is active in cluster k.}

  \item{pi_rna}{the probability for scRNA-Seq data that gene g is expressed in cell l.}


  \item{phi_met}{the probability for sc-methylation data belonging to cluster k (i.e. z[d,k]=1).}

  \item{w_met}{the probability for sc-methylation data that gene g is active in cluster k.}

  \item{pi_met}{the probability for sc-methylation data that gene g is methylated in cell d.}

}

\value{
\item{pst}{the posterior probability for one iteration.}
}

\author{
Jiaxuan WANGWU
}
