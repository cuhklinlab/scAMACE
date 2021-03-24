\name{cal_E_acc}
\alias{cal_E_acc}
\title{cal_E_acc}
\description{
Perform E-step (i.e. calculate the expectations of missing data) for one iteration of scCAS data in the EM algorithm.
}
\usage{
cal_E_acc(rij, replacef0, phi_atac, w_acc, qi)
}

\arguments{
  \item{rij}{ratio f1/(f1+f0) for scCAS data.}

  \item{replacef0}{(1-rij) for scCAS data.}


  \item{phi_atac}{the initialization of the probability for scCAS data belonging to cluster k (i.e. z[i,k]=1).}

  \item{w_acc}{the initialization of the probability for scCAS data that gene g is active in cluster k.}

  \item{qi}{the initialization of the probability that scCAS data has high gene score when gene g is active in cell i (i.e. pi_{i1}).}
}



\value{
  \item{prob_z_atac}{Expectation of z.}

  \item{prob_z_u_atac}{Expectation of z*u.}

  \item{prob_z_u_ut_atac}{Expectation of z*u*o.}


}

\author{
Jiaxuan WANGWU
}
