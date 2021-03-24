\name{cal_M_step}
\alias{cal_M_step}
\title{cal_M_step}
\description{
scAMACE expectation-maximization (EM) implementation, a model-based approach to the joint clustering of single-cell data on chromatin accessibility, gene expression and methylation.
}
\usage{
cal_M_step(phi_atac,w_acc,qi,phi_rna,w_exp,pi_rna,phi_met,w_met,pi_met)
}


\arguments{
  \item{phi_atac}{the initialization of the probability for scCAS data belonging to cluster k (i.e. z[i,k]=1).}

  \item{w_acc}{the initialization of the probability for scCAS data that gene g is active in cluster k.}

  \item{qi}{the initialization of the probability that scCAS data has high gene score when gene g is active in cell i (i.e. pi_{i1}).}


  \item{phi_rna}{the initialization of the probability for scRNA-Seq data belonging to cluster k (i.e. z[l,k]=1).}

  \item{w_exp}{the initialization of the probability for scRNA-Seq data that gene g is active in cluster k.}

  \item{pi_rna}{the initialization of the probability for scRNA-Seq data that gene g is expressed in cell l.}


  \item{phi_met}{the initialization of the probability for sc-methylation data belonging to cluster k (i.e. z[d,k]=1).}

  \item{w_met}{the initialization of the probability for sc-methylation data that gene g is active in cluster k.}

  \item{pi_met}{the initialization of the probability for sc-methylation data that gene g is methylated in cell d.}

}


\value{
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


\author{
Jiaxuan WANGWU
}
