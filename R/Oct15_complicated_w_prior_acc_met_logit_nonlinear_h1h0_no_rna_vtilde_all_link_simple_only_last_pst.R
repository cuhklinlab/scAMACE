### calculate ALL expectation for RNA/methy. in one function (since rna & methy. have same model)
cal_E_rna <- function(rlg,replaceg0,phi_rna, w_exp, pi_rna){
  n <- ncol(pi_rna)
  ### calculate E_z_lk   n*k
  product_1 <- product_2 <- array(0,dim=c(n,k0,p0))
  for(k in c(1:k0)){
    product_1[,k,] <- t(t(rlg*pi_rna[1,]+replaceg0*(1-pi_rna[1,]))*w_exp[k,]) #n*k*p
    product_2[,k,] <- t(t(rlg*pi_rna[2,]+replaceg0*(1-pi_rna[2,]))*(1-w_exp[k,]))
  }

  prob_z_rna <- t(t(apply(log(product_1+product_2),c(1,2),sum))+log(phi_rna)) #n*k + k*1
  prob_z_rna2 <- prob_z_rna

  for(k in c(1:k0)){
    if (k0>2){
      prob_z_rna[,k]<-1/(1+(rowSums(exp(prob_z_rna2[,-k]-prob_z_rna2[,k])))) # for k0>2
      #print(k0)
    }else if (k0==2){
      prob_z_rna[,k]<-1/(1+((exp(prob_z_rna2[,-k]-prob_z_rna2[,k])))) # for k0=2
    }else if (k0==1){
      prob_z_rna <- matrix(1,nrow=n,ncol=k0) # l0*k0
    }
  }


  ### calculate E_z_lk_u_lg n*k*p
  prob_z_u_rna <- array(0,dim=c(n,k0,p0))
  divi_part = product_1 / (product_1 + product_2) #n*k*p
  for(k in c(1:k0)){
    prob_z_u_rna[,k,]<-prob_z_rna[,k]*divi_part[,k,]
  }


  ### calculate E_z_lk_(1-u_lg) n*k*p
  prob_z_1_u_rna <- array(0,dim=c(n,k0,p0))
  for(k in c(1:k0)){
    prob_z_1_u_rna[,k,]<-prob_z_rna[,k] - prob_z_u_rna[,k,]
  }


  ### # calculate E_z_lk_u_lg_v_lg n*k*p
  prob_z_u_v_rna <- array(0,dim=c(n,k0,p0))
  div_1 = pi_rna[1,]*rlg #n*p
  div_2 = (1 - pi_rna[1,])*replaceg0
  tmp_E = div_1 / (div_1 + div_2)
  for(k in c(1:k0)){
    prob_z_u_v_rna[,k,]<-prob_z_u_rna[,k,]*tmp_E
  }


  ### calculate E_z_lk_(1-u_lg)_v_lg n*k*p
  prob_z_1_u_v_rna <- array(0,dim=c(n,k0,p0))
  div_1 = pi_rna[2,]*rlg #n*p
  div_2 = (1 - pi_rna[2,])*replaceg0
  tmp_E = div_1 / (div_1 + div_2)
  for(k in c(1:k0)){
    prob_z_1_u_v_rna[,k,]<-(prob_z_rna[,k]-prob_z_u_rna[,k,])*tmp_E
  }


  ### calculate E_z_lk_u_lg_v_lg n*k*p
  prob_z_v_rna<-prob_z_u_v_rna+prob_z_1_u_v_rna

  E_rna <- list(prob_z_rna=prob_z_rna,prob_z_u_rna=prob_z_u_rna,
                prob_z_1_u_rna=prob_z_1_u_rna,prob_z_u_v_rna=prob_z_u_v_rna,
                prob_z_1_u_v_rna=prob_z_1_u_v_rna,prob_z_v_rna=prob_z_v_rna)
  return(E_rna)

}


### calculate ALL expectation for acc in one function
cal_E_acc <- function(rij, replacef0, phi_atac, w_acc, qi){
  n <- length(qi)
  # calculate E_z_ik n*k, which is also P_old(z_ik = 1)
  product_1 <- product_2 <- array(0,dim=c(n,k0,p0))
  for(k in c(1:k0)){
    product_1[,k,] <- t(t(rij*qi+replacef0*(1-qi))*w_acc[k,]) #n*k*p
    product_2[,k,] <- t(t(replacef0)*(1-w_acc[k,]))
  }

  prob_z_atac <- t(t(apply(log(product_1+product_2),c(1,2),sum))+log(phi_atac)) #n*k + k*1
  prob_z_atac2 <- prob_z_atac

  for(k in c(1:k0)){
    if(k0>2){
      prob_z_atac[,k]<-1/(1+(rowSums(exp(prob_z_atac2[,-k]-prob_z_atac2[,k])))) # for k0>2
      #print(k0)
    }else if (k0==2){
      prob_z_atac[,k]<-1/(1+((exp(prob_z_atac2[,-k]-prob_z_atac2[,k])))) # for k0=2
    }else if (k0==1){
      prob_z_atac <- matrix(1,nrow=n,ncol=k0) # i0*k0
    }
  }


  ### calculate  E_z_ik_u_ir  dimension n*k*p
  prob_z_u_atac <- array(0,dim=c(n,k0,p0))
  divi_part = product_1 / (product_1 + product_2) #n*k*p
  for(k in c(1:k0)){
    prob_z_u_atac[,k,]<-prob_z_atac[,k]*divi_part[,k,]
  }


  ###  calculate E_z_u_u_t dimension n*k*p
  prob_z_u_ut_atac <- array(0,dim=c(n,k0,p0))
  div_1 = qi*rij #n*p
  div_2 = (1 - qi)*replacef0
  tmp_E = div_1 / (div_1 + div_2)
  for(k in c(1:k0)){
    prob_z_u_ut_atac[,k,]<-prob_z_u_atac[,k,]*tmp_E
  }


  E_acc <- list(prob_z_atac=prob_z_atac,prob_z_u_atac=prob_z_u_atac,
                prob_z_u_ut_atac=prob_z_u_ut_atac)
  return(E_acc)

}



## update common w (i.e. linked part) ### grid search: NEW!!! --> modified!!!
# common w: only exists in rna, prior of atac & met and itself
update_w_linked <- function(sum_prob_z_rna,sum_prob_z_u_rna,w_acc,w_met){
  # w seq
  temp <- seq(from=0.01, to=1, by=0.01)
  ntemp <- length(temp)
  pmin <- p0

  ## all linked part (1:k0*pmin)
  w_temp <- matrix(rep(temp,each=k0*pmin),nrow=k0*pmin,byrow = FALSE) # a (k0*po)*ntemp matrix i.e. a 160*100 matrix
  w_acc_linked <- as.vector(w_acc[,(1:pmin)])
  w_met_linked <- as.vector(w_met[,(1:pmin)])

  #w_exp_linked <- as.vector(w_exp[,(1:pmin)])
  sum_prob_z_u_rna_linked <- as.vector(sum_prob_z_u_rna[,(1:pmin)])
  diff_z_and_z_u_rna <- as.vector(sum_prob_z_rna-sum_prob_z_u_rna[,(1:pmin)])

  # rna part
  pt1 <- outer(sum_prob_z_u_rna_linked,log(temp)) # a (k0*po)*ntemp matrix i.e. a 160*100 matrix
  pt2 <- outer(diff_z_and_z_u_rna,log(1-temp))

  # met part
  mu2 <- 1/(1+exp(-(delta+theta*temp)))
  w_temp_met <- 1/(1+exp(-(delta+theta*w_temp)))
  pt3 <- outer(log(w_met_linked),(mu2*phi_2-1))
  pt4 <- outer(log(1-w_met_linked),(phi_2-1-phi_2*mu2))
  pt5 <- beta(w_temp_met*phi_2,(phi_2-phi_2*w_temp_met))

  # atac part
  mu <- 1/(1+exp(-(eta+gamma*temp+tau*temp^2)))
  w_temp_acc <- 1/(1+exp(-(eta+gamma*w_temp+tau*w_temp^2)))
  pt6 <- outer(log(w_acc_linked),(mu*phi_1-1))
  pt7 <- outer(log(1-w_acc_linked),(phi_1-1-phi_1*mu))
  pt8 <- beta(w_temp_acc*phi_1,(phi_1-phi_1*w_temp_acc))

  # prior part
  pt9 <- (alpha1-1)*log(w_temp)+(beta1-1)*log(1-w_temp)

  y <- -pt1-pt2-pt3-pt4+log(pt5)-pt6-pt7+log(pt8)-pt9
  w_index <- apply(y,1, function(x) which.min(x))
  w_index <- unlist(w_index) ### modified
  w_all <- temp[w_index]
  w_exp <- matrix(w_all, nrow=k0, byrow = FALSE)

  return(w_exp)
}


### M-step
cal_M_step <- function(phi_atac,w_acc,qi,
                       phi_rna,w_exp,pi_rna,
                       phi_met,w_met,pi_met){
  # prob_z_atac: n*k
  # prob_z_u_atac: n*k*p

  acc_E <- cal_E_acc(rij, replacef0,phi_atac, w_acc, qi)

  prob_z_atac <- acc_E$prob_z_atac
  prob_z_u_atac <- acc_E$prob_z_u_atac
  prob_z_u_ut_atac <- acc_E$prob_z_u_ut_atac

  rna_E <- cal_E_rna(rlg,replaceg0,phi_rna, w_exp, pi_rna)

  prob_z_rna <- rna_E$prob_z_rna
  prob_z_u_rna <- rna_E$prob_z_u_rna
  prob_z_1_u_rna <- rna_E$prob_z_1_u_rna
  prob_z_u_v_rna <- rna_E$prob_z_u_v_rna
  prob_z_1_u_v_rna <- rna_E$prob_z_1_u_v_rna
  prob_z_v_rna <- rna_E$prob_z_v_rna

  met_E <- cal_E_rna(rdm,replaceh0,phi_met, w_met, pi_met)

  prob_z_met <- met_E$prob_z_rna
  prob_z_u_met <- met_E$prob_z_u_rna
  prob_z_1_u_met <- met_E$prob_z_1_u_rna
  prob_z_u_v_met <- met_E$prob_z_u_v_rna
  prob_z_1_u_v_met <- met_E$prob_z_1_u_v_rna
  prob_z_v_met <- met_E$prob_z_v_rna


  ### update phi
  phi_atac_new<-(1+colSums(prob_z_atac))/(i0+k0)
  phi_rna_new<-(1+colSums(prob_z_rna))/(l0+k0)
  phi_met_new<-(1+colSums(prob_z_met))/(d0+k0)


  ### update qi
  qi_new <- (apply(prob_z_u_ut_atac,1,sum)+alpha_qi-1)/(apply(prob_z_u_atac,1,sum)+alpha_qi+beta_qi-2)


  ### update pi_rna
  pi_rna_new <- pi_rna
  temp1<-apply(prob_z_u_v_rna,1,sum)/apply(prob_z_u_rna,1,sum)
  temp1[which(temp1>=1-10^-5)]=1-10^-5

  tmp1<-apply(prob_z_1_u_v_rna,1,sum)/(apply(prob_z_1_u_rna,1,sum)-1)

  flag1<-which(tmp1<=temp1) # pi_l0 < pi_l1
  pi_rna_new[1,flag1]<-temp1[flag1]
  pi_rna_new[2,flag1]<-tmp1[flag1]


  ### update pi_met
  pi_met_new <- pi_met
  temp2<-apply(prob_z_u_v_met,1,sum)/apply(prob_z_u_met,1,sum)
  tmp2<-(apply(prob_z_1_u_v_met,1,sum)-1)/(apply(prob_z_1_u_met,1,sum)-1)
  tmp2[which(tmp2>=1-10^-5)]=1-10^-5

  flag2<-which(tmp2>=temp2) # pi_d0 > pi_d1
  pi_met_new[1,flag2]<-temp2[flag2]
  pi_met_new[2,flag2]<-tmp2[flag2]


  ### update w_rna linked part via grid search
  w_exp_new <- update_w_linked(colSums(prob_z_rna),apply(prob_z_u_rna,c(2,3),sum),w_acc,w_met)


  ### update w_acc
  mu <- 1/(1+exp(-(eta+gamma*w_exp_new+tau*w_exp_new^2)))
  w_acc_new<-(apply(prob_z_u_atac,c(2,3),sum)+mu*phi_1-1)/(apply(prob_z_atac,2,sum)+phi_1-2)
  w_acc_new[which(w_acc_new>=1-10^-6)]<-1-10^-6
  w_acc_new[which(w_acc_new<=10^-6)]<-10^-6


  ### update w_met
  mu2 <- 1/(1+exp(-(delta + theta*w_exp_new)))
  w_met_new<-(apply(prob_z_u_met,c(2,3),sum)+mu2*phi_2-1)/(apply(prob_z_met,2,sum)+phi_2-2)
  w_met_new[which(w_met_new>=1-10^-6)]<-1-10^-6
  w_met_new[which(w_met_new<=10^-6)]<-10^-6

  # calculate posterior probability pst
  # pst <- cal_post(phi_atac_new,w_acc_new,qi_new,phi_rna_new,w_exp_new,pi_rna_new,
  #                 phi_met_new,w_met_new,pi_met_new)


  return(list(phi_atac=phi_atac_new,phi_rna=phi_rna_new,phi_met=phi_met_new,
              w_exp=w_exp_new,w_acc=w_acc_new,w_met=w_met_new,#w=w_new,
              pi_rna=pi_rna_new,qi=qi_new,pi_met=pi_met_new))#,
              #postprob=pst))


}



### calculate posterior
cal_post <- function(phi_atac_new,w_acc_new,qi_new,
                     phi_rna_new,w_exp_new,pi_rna_new,
                     phi_met_new,w_met_new,pi_met_new){

  ## likelihood for X (i.e. acta data)
  product<-matrix(0,nrow=i0,ncol=k0)
  # replacef0 <- 1-rij
  for(k in c(1:k0)){
    product[,k]<-colSums(log(t(qi_new*rij+(1-qi_new)*replacef0)*w_acc_new[k,]+t(replacef0)*(1-w_acc_new[k,])))
  }
  #product<-apply(log(t1),c(1,2),sum) #sum over p0 (i.e. sum over r) --> a i0*k0 matrix
  # part1<-sum(product[,1])+sum(log(rowSums(t(t(exp(product-product[,1]))*phi_atac_new)))) #sum over k0

  d<-max(product-product[,1])
  part1<-sum(product[,1])+d/2*i0+sum(log(rowSums(t(phi_atac_new*t(exp(product-product[,1]-d/2))))))


  # calculate posterior probability pst
  ## likelihood for Y (i.e. rna data)
  product2<-matrix(0,nrow=l0,ncol=k0)
  #replaceg0 <- matrix(1,nrow=l0,ncol=p0)
  # replaceg0 <- 1-rlg
  for(k in c(1:k0)){
    product2[,k]<-colSums(log(t(rlg*pi_rna_new[1,]+(1-pi_rna_new[1,])*replaceg0)*w_exp_new[k,]+t(rlg*pi_rna_new[2,]+(1-pi_rna_new[2,])*replaceg0)*(1-w_exp_new[k,])))
  }
  #product2<-apply(log(t2),c(1,2),sum)
  a<-max(product2-product2[,1])
  part2<-sum(product2[,1])+a/2*l0+sum(log(rowSums(t(phi_rna_new*t(exp(product2-product2[,1]-a/2))))))


  # calculate posterior probability pst
  ## likelihood for T (i.e. methy. data)
  product3<-matrix(0,nrow=d0,ncol=k0)
  # replaceh0 <- matrix(1,nrow=d0,ncol=p0)
  #replaceh0 <- 1-rdm
  for(k in c(1:k0)){
    product3[,k]<-colSums(log(t(rdm*pi_met_new[1,]+(1-pi_met_new[1,])*replaceh0)*w_met_new[k,]+t(rdm*pi_met_new[2,]+(1-pi_met_new[2,])*replaceh0)*(1-w_met_new[k,])))
  }
  #product3<-apply(log(t3),c(1,2),sum,na.rm=T) ### modified
  b<-max(product3-product3[,1])
  part3<-sum(product3[,1])+b/2*d0+sum(log(rowSums(t(phi_met_new*t(exp(product3-product3[,1]-b/2))))))


  # prior
  part4_phi <- sum(log(phi_atac_new)+log(phi_rna_new)+log(phi_met_new))

  mu_new <- 1/(1+exp(-(eta+gamma*w_exp_new+tau*w_exp_new^2)))
  part4_w_acc_linked<-sum((phi_1*mu_new-1)*log(w_acc_new)+(phi_1-1-phi_1*mu_new)*log(1-w_acc_new)-log(beta(phi_1*mu_new,phi_1-phi_1*mu_new)))

  mu2_new <- 1/(1+exp(-(delta+theta*w_exp_new)))
  part4_w_met_linked<-sum((phi_2*mu2_new-1)*log(w_met_new)+(phi_2-1-phi_2*mu2_new)*log(1-w_met_new)-log(beta(phi_2*mu2_new,phi_2-phi_2*mu2_new)))

  part4_w<-sum((alpha1-1)*log(w_exp_new)+(beta1-1)*log(1-w_exp_new))

  part4_qi <- sum((alpha_qi-1)*log(qi_new)+(beta_qi-1)*log(1-qi_new))

  part4<-part4_phi+part4_w_acc_linked+part4_w+part4_w_met_linked-sum(log(1-pi_rna_new[2,]))-sum(log(pi_met_new[2,]))+part4_qi

  pst<-part1+part2+part3+part4

  return(pst)

}



###----------------------------------------------------------------------###
###----------------------------------------------------------------------###
###----------------------------------------------------------------------###
###----------------------------------------------------------------------###
### functions for seperate estimation
### atac
### M-step
cal_M_step_atac <- function(phi_atac,w_acc,qi){
  # prob_z_atac: n*k
  # prob_z_u_atac: n*k*p

  acc_E <- cal_E_acc(rij, replacef0,phi_atac, w_acc, qi)

  prob_z_atac <- acc_E$prob_z_atac
  prob_z_u_atac <- acc_E$prob_z_u_atac
  prob_z_u_ut_atac <- acc_E$prob_z_u_ut_atac


  ### update phi
  phi_atac_new<-(1+colSums(prob_z_atac))/(i0+k0)

  ### update qi
  qi_new <- (apply(prob_z_u_ut_atac,1,sum)+alpha_qi-1)/(apply(prob_z_u_atac,1,sum)+alpha_qi+beta_qi-2)


  ### update w_acc
  w_acc_new<-(apply(prob_z_u_atac,c(2,3),sum)+alpha1-1)/(apply(prob_z_atac,2,sum)+alpha1+beta1-2)
  w_acc_new[which(w_acc_new>=1-10^-6)]<-1-10^-6
  w_acc_new[which(w_acc_new<=10^-6)]<-10^-6

  return(list(phi_atac=phi_atac_new,w_acc=w_acc_new,qi=qi_new))

  # calculate posterior probability pst
  # pst <- cal_post_atac(phi_atac_new,w_acc_new,qi_new)
  # return(list(phi_atac=phi_atac_new,w_acc=w_acc_new,qi=qi_new,postprob=pst))


}


### calculate posterior
cal_post_atac <- function(phi_atac_new,w_acc_new,qi_new){

  ## likelihood for X (i.e. acta data)
  product<-matrix(0,nrow=i0,ncol=k0)
  # replacef0 <- 1-rij
  for(k in c(1:k0)){
    product[,k]<-colSums(log(t(qi_new*rij+(1-qi_new)*replacef0)*w_acc_new[k,]+t(replacef0)*(1-w_acc_new[k,])))
  }
  #product<-apply(log(t1),c(1,2),sum) #sum over p0 (i.e. sum over r) --> a i0*k0 matrix
  # part1<-sum(product[,1])+sum(log(rowSums(t(t(exp(product-product[,1]))*phi_atac_new)))) #sum over k0

  d<-max(product-product[,1])
  part1<-sum(product[,1])+d/2*i0+sum(log(rowSums(t(phi_atac_new*t(exp(product-product[,1]-d/2))))))



  # prior
  part4_phi <- sum(log(phi_atac_new))

  part4_w<-sum((alpha1-1)*log(w_acc_new)+(beta1-1)*log(1-w_acc_new))

  part4_qi <- sum((alpha_qi-1)*log(qi_new)+(beta_qi-1)*log(1-qi_new))

  part4<-part4_phi+part4_w+part4_qi

  pst<-part1+part4

  return(pst)

}



### rna
### M-step
cal_M_step_rna <- function(phi_rna,w_exp,pi_rna){
  # prob_z_atac: n*k
  # prob_z_u_atac: n*k*p

  rna_E <- cal_E_rna(rlg,replaceg0,phi_rna, w_exp, pi_rna)

  prob_z_rna <- rna_E$prob_z_rna
  prob_z_u_rna <- rna_E$prob_z_u_rna
  prob_z_1_u_rna <- rna_E$prob_z_1_u_rna
  prob_z_u_v_rna <- rna_E$prob_z_u_v_rna
  prob_z_1_u_v_rna <- rna_E$prob_z_1_u_v_rna
  prob_z_v_rna <- rna_E$prob_z_v_rna


  ### update phi
  phi_rna_new<-(1+colSums(prob_z_rna))/(l0+k0)


  ### update pi_rna
  pi_rna_new <- pi_rna
  temp1<-apply(prob_z_u_v_rna,1,sum)/apply(prob_z_u_rna,1,sum)
  temp1[which(temp1>=1-10^-5)]=1-10^-5

  tmp1<-apply(prob_z_1_u_v_rna,1,sum)/(apply(prob_z_1_u_rna,1,sum)-1)

  flag1<-which(tmp1<=temp1) # pi_l0 < pi_l1
  pi_rna_new[1,flag1]<-temp1[flag1]
  pi_rna_new[2,flag1]<-tmp1[flag1]


  ### update w_rna linked part via grid search
  w_exp_new <- (apply(prob_z_u_rna,c(2,3),sum)+alpha1-1)/(apply(prob_z_rna,2,sum)+alpha1+beta1-2)
  w_exp_new[which(w_exp_new>=1-10^-6)]<-1-10^-6
  w_exp_new[which(w_exp_new<=10^-6)]<-10^-6

  return(list(phi_rna=phi_rna_new,w_exp=w_exp_new,pi_rna=pi_rna_new))

  # calculate posterior probability pst
  # pst <- cal_post_rna(phi_rna_new,w_exp_new,pi_rna_new)
  # return(list(phi_rna=phi_rna_new,w_exp=w_exp_new,pi_rna=pi_rna_new,postprob=pst))


}


### calculate posterior
cal_post_rna <- function(phi_rna_new,w_exp_new,pi_rna_new){

  # calculate posterior probability pst
  ## likelihood for Y (i.e. rna data)
  product2<-matrix(0,nrow=l0,ncol=k0)
  #replaceg0 <- matrix(1,nrow=l0,ncol=p0)
  # replaceg0 <- 1-rlg
  for(k in c(1:k0)){
    product2[,k]<-colSums(log(t(rlg*pi_rna_new[1,]+(1-pi_rna_new[1,])*replaceg0)*w_exp_new[k,]+t(rlg*pi_rna_new[2,]+(1-pi_rna_new[2,])*replaceg0)*(1-w_exp_new[k,])))
  }
  #product2<-apply(log(t2),c(1,2),sum)
  a<-max(product2-product2[,1])
  part2<-sum(product2[,1])+a/2*l0+sum(log(rowSums(t(phi_rna_new*t(exp(product2-product2[,1]-a/2))))))



  # prior
  part4_phi <- sum(log(phi_rna_new))

  part4_w<-sum((alpha1-1)*log(w_exp_new)+(beta1-1)*log(1-w_exp_new))


  part4<-part4_phi+part4_w-sum(log(1-pi_rna_new[2,]))

  pst<-part2+part4

  return(pst)

}


### methy.
### M-step
cal_M_step_met <- function(phi_met,w_met,pi_met){
  # prob_z_atac: n*k
  # prob_z_u_atac: n*k*p

  met_E <- cal_E_rna(rdm,replaceh0,phi_met, w_met, pi_met)

  prob_z_met <- met_E$prob_z_rna
  prob_z_u_met <- met_E$prob_z_u_rna
  prob_z_1_u_met <- met_E$prob_z_1_u_rna
  prob_z_u_v_met <- met_E$prob_z_u_v_rna
  prob_z_1_u_v_met <- met_E$prob_z_1_u_v_rna
  prob_z_v_met <- met_E$prob_z_v_rna


  ### update phi
  phi_met_new<-(1+colSums(prob_z_met))/(d0+k0)


  ### update pi_met
  pi_met_new <- pi_met
  temp2<-apply(prob_z_u_v_met,1,sum)/apply(prob_z_u_met,1,sum)
  tmp2<-(apply(prob_z_1_u_v_met,1,sum)-1)/(apply(prob_z_1_u_met,1,sum)-1)
  tmp2[which(tmp2>=1-10^-5)]=1-10^-5

  flag2<-which(tmp2>=temp2) # pi_d0 > pi_d1
  pi_met_new[1,flag2]<-temp2[flag2]
  pi_met_new[2,flag2]<-tmp2[flag2]


  ### update w_rna linked part via grid search
  w_met_new <- (apply(prob_z_u_met,c(2,3),sum)+alpha1-1)/(apply(prob_z_met,2,sum)+alpha1+beta1-2)
  w_met_new[which(w_met_new>=1-10^-6)]<-1-10^-6
  w_met_new[which(w_met_new<=10^-6)]<-10^-6

  return(list(phi_met=phi_met_new,w_met=w_met_new,pi_met=pi_met_new))

  # calculate posterior probability pst
  # pst <- cal_post_met(phi_met_new,w_met_new,pi_met_new)
  # return(list(phi_met=phi_met_new,w_met=w_met_new,pi_met=pi_met_new,postprob=pst))


}



### calculate posterior
cal_post_met <- function(phi_met_new,w_met_new,pi_met_new){

  # calculate posterior probability pst
  ## likelihood for T (i.e. methy. data)
  product3<-matrix(0,nrow=d0,ncol=k0)
  # replaceh0 <- matrix(1,nrow=d0,ncol=p0)
  #replaceh0 <- 1-rdm
  for(k in c(1:k0)){
    product3[,k]<-colSums(log(t(rdm*pi_met_new[1,]+(1-pi_met_new[1,])*replaceh0)*w_met_new[k,]+t(rdm*pi_met_new[2,]+(1-pi_met_new[2,])*replaceh0)*(1-w_met_new[k,])))
  }
  #product3<-apply(log(t3),c(1,2),sum,na.rm=T) ### modified
  b<-max(product3-product3[,1])
  part3<-sum(product3[,1])+b/2*d0+sum(log(rowSums(t(phi_met_new*t(exp(product3-product3[,1]-b/2))))))


  # prior
  part4_phi <- sum(log(phi_met_new))

  part4_w<-sum((alpha1-1)*log(w_met_new)+(beta1-1)*log(1-w_met_new))


  part4<-part4_phi+part4_w-sum(log(pi_met_new[2,]))

  pst<-part3+part4

  return(pst)

}


##-------------------------------------------------------------------------##
##-------------------------------------------------------------------------##
##-------------------------------------------------------------------------##
### Get logit
# use beta regression to estimate phi
get_logit_beta <- function(w_exp_link,w_acc_link,poly){
  if (poly=="linear"){
    logit <- betareg(w_acc_link~w_exp_link,link = "logit")
  } else if (poly=="quadratic"){
    w_exp_link_sq <- w_exp_link^2
    logit <- betareg(w_acc_link~w_exp_link+w_exp_link_sq,link = "logit")
  }

  reg_para <- logit$coefficients$mean
  prec_para <- logit$coefficients$precision

  beta_reg <- list(reg_para=reg_para,prec_para=prec_para)
  return(beta_reg)
}



