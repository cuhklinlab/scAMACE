\name{simData_3data}
\alias{simData_3data}
\title{simData_3data}
\description{
Generate simulation data $\boldsymbol x$, $\boldsymbol y$ and $\boldsymbol z$.
}
\usage{
simData_3data(n1, n2, n3, p, Xprobs1, Xprobs2, Xprobs3,
                           shape_rna, scale_rna, shape_met1, shape_met2,
                           qi,pi_rna,pi_met,cutoff=10^-6,
                           alpha1,beta1,phi_1,phi_2,omega,prop)
}
\arguments{
  \item{n1}{sample size of data type $\boldsymbol x$.}
  \item{n2}{sample size of data type $\boldsymbol y$.}
  \item{n3}{sample size of data type $\boldsymbol z$.}
  \item{p}{number of features.}

  \item{Xprobs1}{prior probability of k clusters in data type $\boldsymbol x$.}
  \item{Xprobs2}{prior probability of k clusters in data type $\boldsymbol y$.}
  \item{Xprobs3}{prior probability of k clusters in data type $\boldsymbol z$.}

  \item{shape_rna}{shape parameter in data type $\boldsymbol y$ gamma mixture distribution.}
  \item{scale_rna}{scale parameter in data type $\boldsymbol y$ gamma mixture distribution.}

  \item{shape_met1}{shape 1 parameter in data type $\boldsymbol z$ beta mixture distribution.}
  \item{shape_met2}{shape 2 parameter in data type $\boldsymbol z$ beta mixture distribution.}

  \item{qi}{the probability that scCAS data has high gene score when gene g is active in cell i (i.e. pi_{i1}).}}
  \item{pi_rna}{the probability for scRNA-Seq data that gene g is expressed in cell l.}
  \item{pi_met}{the probability for sc-methylation data that gene g is methylated in cell d.}

  \item{alpha1}{$\alpha$ in prior probability (beta distribution) of $\omega_{kg}$.}
  \item{beta1}{$\beta$ in prior probability (beta distribution) of $\omega_{kg}$.}

  \item{phi_1}{precision parameter in the logit link between scCAS data and scRNA-Seq data.}
  \item{phi_2}{precision parameter in the logit link between sc-methylation data and scRNA-Seq data.}

  \item{omega}{differential $\omega_{kg}$ in scRNA-Seq data.}
  \item{prop}{proportion of differential $\omega_{kg}$ in scRNA-Seq data.}


}

\value{
\item{}{Simulation data $\boldsymbol x$, $\boldsymbol y$ and $\boldsymbol z$ and parameters.}
}

\author{
Jiaxuan WANGWU
}
