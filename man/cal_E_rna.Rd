\name{cal_E_rna}
\alias{cal_E_rna}
\title{cal_E_rna}
\description{
Perform E-step (i.e. calculate the expectations of missing data) for one iteration of scRNA-Seq or sc-methylation data in the EM algorithm.
}
\usage{
cal_E_rna(rlg,replaceg0,phi_rna, w_exp, pi_rna)
}

\arguments{
  \item{rlg}{ratio g1/(g1+g0) for scRNA-Seq data or ratio h1/h0 for sc-methylation data.}

  \item{replaceg0}{(1-rlg) for scRNA-Seq data or a matrix of constant 1 for sc-methylation data.}

  \item{phi_rna}{the probability for scRNA-Seq or sc-methylation data belonging to cluster k (i.e. z[l,k]=1).}

  \item{w_exp}{the probability for scRNA-Seq or sc-methylation data that gene g is active in cluster k.}

  \item{pi_rna}{the probability for scRNA-Seq or sc-methylation data that gene g is expressed in cell l.}
}
}

\value{

  \item{prob_z_rna}{Expectation of z.}

  \item{prob_z_u_rna}{Expectation of z*u.}

  \item{prob_z_1_u_rna}{Expectation of z*(1-u).}

  \item{prob_z_u_v_rna}{Expectation of z*u*v for scRNA-Seq data or expectation of z*u*m for sc-methylation data.}

  \item{prob_z_1_u_v_rna}{Expectation of z*(1-u)*v for scRNA-Seq data or expectation of z*(1-u)*m for sc-methylation data.}

  \item{prob_z_v_rna}{Expectation of z*v for scRNA-Seq data or expectation of z*m for sc-methylation data.}



}

\author{
Jiaxuan WANGWU
}
