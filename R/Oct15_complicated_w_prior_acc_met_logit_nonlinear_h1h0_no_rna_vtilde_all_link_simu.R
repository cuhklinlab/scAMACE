
library(magic)

get_mu_quadratic <- function(w){
  f <- eta+gamma*w+tau*w^2
  mu1 <- 1/(1+exp(-f))
  return(mu1)
}

get_mu_linear <- function(w){
  g <- delta+theta*w
  mu2 <- 1/(1+exp(-g))
  return(mu2)
}




# add constant omega to rna, acc & methy.
# allow missing cluster (k are different in different dataset) i.e. Xprobs1 has 0 and are different
simData_3data  <- function(n1, n2, n3, p, Xprobs1, Xprobs2, Xprobs3,
                           shape_rna, scale_rna, shape_met1, shape_met2,
                           qi,pi_rna,pi_met,cutoff=10^-6,
                           alpha1,beta1,phi_1,phi_2,omega,prop){
  # K are different across different dataset
  nComp1 <- length(Xprobs1) # number of clusters: k_atac
  nComp2 <- length(Xprobs2) # number of clusters: k_rna
  nComp3 <- length(Xprobs3) # number of clusters: k_met

  #########################################################################################################################
  ## Simulate some "fixed" omega in w_Exp
  div <- nComp2%/%2
  if(nComp2%%2==1){
    # nComp2 is odd number
    omega_seq1 <- sort(c(seq(0,omega,omega/div),0.5,(1-seq(0,omega,omega/div))),decreasing =T)
  }else{
    # nComp2 is even number
    omega_seq1 <- sort(c(seq(0,omega,omega/div),(1-seq(0,omega,omega/div))),decreasing =T)
   }


  omega_seq <- omega_seq1[2:(length(omega_seq1)-1)]
  w <- matrix(0,nrow=nComp2,ncol=p*prop*nComp2)
  for (j in 1:nComp2){
    w[,((j-1)*p*prop+1):(j*p*prop)] <- shift(omega_seq,nComp2+1-j)
    print(shift(omega_seq,nComp2+1-j))
  }



  ## simulate pi_pk_rna --> generate w_kj^exp
  w_cons <- rbeta(n=p, shape1=alpha1, shape2=beta1)
  pi_pk_rna <- matrix(rep(w_cons,nComp2), nrow=nComp2, ncol=p,byrow = T)
  pi_pk_rna[,1:(p*prop*nComp2)] <- w


  ## simulate pi_pk_atac --> generate w_kj^acc
  mu_atac <- get_mu_quadratic(w)
  mu_atac_cons <- get_mu_quadratic(w_cons)
  pi_pk_atac_diff <- matrix(rbeta(n=p*prop*nComp1*nComp1, shape1=as.numeric(mu_atac*phi_1),shape2=as.numeric(-mu_atac*phi_1+phi_1)), nrow=nComp1, ncol=p*prop*nComp1)
  w_acc_cons <- rbeta(n=p, shape1=as.numeric(mu_atac_cons*phi_1),shape2=as.numeric(-mu_atac_cons*phi_1+phi_1))
  pi_pk_atac <- matrix(rep(w_acc_cons,nComp1), nrow=nComp1, ncol=p,byrow = T)
  pi_pk_atac[,1:(p*prop*nComp1)] <- pi_pk_atac_diff


  ## simulate pi_pk_met --> generate w_km^met
  mu_met <- get_mu_linear(w)
  mu_met_cons <- get_mu_linear(w_cons)
  pi_pk_met_diff <- matrix(rbeta(n=p*prop*nComp3*nComp3, shape1=as.numeric(mu_met*phi_2),shape2=as.numeric(-mu_met*phi_2+phi_2)), nrow=nComp3, ncol=p*prop*nComp3)
  w_met_cons <- rbeta(n=p, shape1=as.numeric(mu_met_cons*phi_2),shape2=as.numeric(-mu_met_cons*phi_2+phi_2))
  pi_pk_met <- matrix(rep(w_met_cons,nComp3), nrow=nComp3, ncol=p,byrow = T)
  pi_pk_met[,1:(p*prop*nComp3)] <- pi_pk_met_diff

  ## implement the cutoff
  pi_pk_atac[which(pi_pk_atac<=cutoff)] <- cutoff
  pi_pk_rna[which(pi_pk_rna<=cutoff)] <- cutoff
  pi_pk_met[which(pi_pk_met<=cutoff)] <- cutoff
  pi_pk_atac[which(pi_pk_atac>=(1-cutoff))] <- 1-cutoff
  pi_pk_rna[which(pi_pk_rna>=(1-cutoff))] <- 1-cutoff
  pi_pk_met[which(pi_pk_met>=(1-cutoff))] <- 1-cutoff


  #################################################################################################################
  ## simulate x_atac --> generate z^acc
  x_atac <- t(rmultinom(n=n1, size=1, prob=Xprobs1))
  cluster_atac <- rep(0, n1)
  for (i in 1:nComp1){
    cluster_atac[which(x_atac[,i]==1)] <- i
  }

  ## simulate x_rna --> generate z^exp
  x_rna <- t(rmultinom(n=n2, size=1, prob=Xprobs2))
  cluster_rna <- rep(0, n2)
  for (i in 1:nComp2){
    cluster_rna[which(x_rna[,i]==1)] <- i
  }

  ## simulate x_met --> generate z^met
  x_met <- t(rmultinom(n=n3, size=1, prob=Xprobs3))
  cluster_met <- rep(0, n3)
  for (i in 1:nComp3){
    cluster_met[which(x_met[,i]==1)] <- i
  }

  ###############################################################################################
  ## simulate u_atac --> generate u^acc -->> use >runif to generate Bernoulli(w_kj^acc)
  pi_pk_atac_x <- x_atac%*%pi_pk_atac
  # pi_pk_atac_x2 <- x_atac[,which(Xprobs1!=0)]%*%pi_pk_atac[which(Xprobs1!=0),]
  # which(pi_pk_atac_x!=pi_pk_atac_x2) #0

  u_atac <- (pi_pk_atac_x>runif(length(pi_pk_atac_x))) + 0

  ## simulate u_rna --> generate u^exp
  pi_pk_rna_x <- x_rna%*%pi_pk_rna
  # pi_pk_rna_x2 <- x_rna[,which(Xprobs2!=0)]%*%pi_pk_rna[which(Xprobs2!=0),]
  # which(pi_pk_rna_x!=pi_pk_rna_x2) #0

  u_rna <- (pi_pk_rna_x>runif(length(pi_pk_rna_x))) + 0

  ## simulate u_met --> generate u^met
  pi_pk_met_x <- x_met%*%pi_pk_met
  # pi_pk_met_x2 <- x_met[,which(Xprobs3!=0)]%*%pi_pk_met[which(Xprobs3!=0),]
  # which(pi_pk_met_x!=pi_pk_met_x2) #0

  u_met <- (pi_pk_met_x>runif(length(pi_pk_met_x))) + 0

  ############################################################################################3
  ## simulate u_t_atac --> generate u tilde^acc
  pi_gA_1 <- matrix(rep(qi[1,], p), nrow=ncol(qi))
  pi_gA_0 <- matrix(rep(qi[2,], p), nrow=ncol(qi))
  u_t_atac <- u_atac
  u_t_atac[which(u_atac==1)] <- (runif(length(which(u_atac==1))) < pi_gA_1[which(u_atac==1)]) + 0
  u_t_atac[which(u_atac==0)] <- (runif(length(which(u_atac==0))) < pi_gA_0[which(u_atac==0)]) + 0


  ############################################################################################################################################################################
  ## simulate v_rna --> generate v^exp -->> include Gp (p*overlap_prop_rna) and G_p
  pi_gA_1 <- matrix(rep(pi_rna[1,], p), nrow=ncol(pi_rna))
  pi_gA_0 <- matrix(rep(pi_rna[2,], p), nrow=ncol(pi_rna))
  v_rna <- u_rna
  v_rna[which(u_rna==1)] <- (runif(length(which(u_rna==1))) < pi_gA_1[which(u_rna==1)]) + 0
  v_rna[which(u_rna==0)] <- (runif(length(which(u_rna==0))) < pi_gA_0[which(u_rna==0)]) + 0


  ############################################################################################################################################################################
  ## simulate v_met --> generate v^met -->> include Gm (p*overlap_prop_met) and G_m
  pi_gA_1 <- matrix(rep(pi_met[1,], p), nrow=ncol(pi_met))
  pi_gA_0 <- matrix(rep(pi_met[2,], p), nrow=ncol(pi_met))
  v_met <- u_met
  v_met[which(u_met==1)] <- (runif(length(which(u_met==1))) < pi_gA_1[which(u_met==1)]) + 0
  v_met[which(u_met==0)] <- (runif(length(which(u_met==0))) < pi_gA_0[which(u_met==0)]) + 0


  #######################################################################################################
  ## simulate atac data --> generate x_r: 0/1 ### Modified!!!
  Data <- matrix(nrow=n1, ncol=p)
  Data[which(u_t_atac==1)] <- 1
  Data[which(u_t_atac==0)] <- 0
  f1 <- f0 <- matrix(0,nrow=n1, ncol=p)
  f1[which(u_t_atac==1)] <- 1
  f0[which(u_t_atac==0)] <- 1
  Data_atac <- Data

  ## simulate rna data --> generate y_g: 2 component gamma ### Modified!!!
  Data <- matrix(nrow=n2, ncol=p)
  Data[which(v_rna==1)] <- rgamma(length(which(v_rna==1)), shape=shape_rna[1], scale=scale_rna[1])
  Data[which(v_rna==0)] <- rgamma(length(which(v_rna==0)), shape=shape_rna[2], scale=scale_rna[2])

  g1 <- dgamma(Data, shape=shape_rna[1], scale=scale_rna[1])
  g0 <- dgamma(Data, shape=shape_rna[2], scale=scale_rna[2])
  Data_rna <- Data

  ## simulate methy. data --> generate t_m: 2 component beta ### Modified!!!
  Data <- matrix(nrow=n3, ncol=p)
  Data[which(v_met==1)] <- rbeta(length(which(v_met==1)), shape1=shape_met1[1], shape2=shape_met2[1])
  Data[which(v_met==0)] <- rbeta(length(which(v_met==0)), shape1=shape_met1[2], shape2=shape_met2[2])
  h1 <- dbeta(Data, shape1=shape_met1[1], shape2=shape_met2[1])
  h0 <- dbeta(Data, shape1=shape_met1[2], shape2=shape_met2[2])
  Data_met <- Data

  ##
  return(list(f1=f1, f0=f0,
              g1=g1, g0=g0,
              h1=h1, h0=h0,
              cluster_atac=cluster_atac, cluster_rna=cluster_rna, cluster_met=cluster_met,
              x_atac=x_atac[,which(Xprobs1!=0)],
              x_rna=x_rna[,which(Xprobs2!=0)],
              x_met=x_met[,which(Xprobs3!=0)],
              x_atac_original=x_atac,
              x_rna_original=x_rna,
              x_met_original=x_met,
              u_atac=u_atac, u_rna=u_rna, u_met=u_met,
              v_rna=v_rna, v_met=v_met,
              u_t_atac=u_t_atac, #to generate v^exp and v tilde^exp in Gp
              Data_atac=Data_atac, Data_rna=Data_rna, Data_met=Data_met,
              w=w,pi_pk_atac=pi_pk_atac[which(Xprobs1!=0),],
              pi_pk_rna=pi_pk_rna[which(Xprobs2!=0),],
              pi_pk_met=pi_pk_met[which(Xprobs3!=0),],
              pi_pk_atac_original=pi_pk_atac,
              pi_pk_rna_original=pi_pk_rna,
              pi_pk_met_original=pi_pk_met,
              qi=qi, pi_rna=pi_rna, pi_met=pi_met,
              n1=n1,n2=n2,n3=n3, p=p,k1=nComp1,k2=nComp2,k3=nComp3,
              Xprobs1=Xprobs1,Xprobs2=Xprobs2,Xprobs3=Xprobs3,
              shape_rna=shape_rna, scale_rna=scale_rna,
              hape_met1=shape_met1, shape_met2=shape_met2,
              alpha1=alpha1,beta1=beta1,phi_1=phi_1,phi_2=phi_2,
              eta=eta,gamma=gamma,tau=tau,delta=delta,theta=theta,
              omega=omega,prop=prop
  ))
}

