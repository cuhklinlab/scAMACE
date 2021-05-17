# A quick guide to scAMACE

#### We present the usage of scAMACE_py in the following sections:
#### [Section 1: Introduction to datasets (Application 1: K562-GM12878 dataset)](#section1)
#### [Section 2: scAMACE](#section2)
#### [Section 3: scAMACE(seperate)](#section3)
#### [Section 4: Beta regression](#section4)
#### [Section 5: Generate simulation data](#section5)


## <a name="section1"></a>1. Datasets
### Application 1: K562-GM12878 dataset
`Feb7_2021_3Types_Data_rna_mean_1000_ratio.RData`: Application 1: K562-GM12878 dataset, including scCAS data, scRNA-Seq data and sc-methylation data.

`Feb7_2021_3Types_Data_rna_mean_1000_ratio_mcmc_ini.RData`: Initialization for Application 1: K562-GM12878 dataset, including initialization data for scCAS data, scRNA-Seq data and sc-methylation data.



## <a name="section2"></a>2. Example of scAMACE
Remarks: We demostrate usage of scAMACE through Application 1.
### 2.1 Load data and prepare for EM algorithm
```{r}
library(scAMACE)

### load the data
data('Feb7_2021_3Types_Data_rna_mean_1000_ratio')
data <- Feb7_2021_3Types_Data_rna_mean_1000_ratio


f10ratio_acc <- data$f10ratio_acc
f1 <- data$f1
f0 <- data$f0
atac_cell_lb <- data$atac_cell

f10ratio_exp <- data$f10ratio_exp
g1 <- data$g1
g0 <- data$g0
rna_cell_lb <- data$rna_cell

f10ratio_met <- data$f10ratio_met
met_data <- data$met_data
met_cell_lb <- data$met_cell


### load the initialization
data("Feb7_2021_3Types_Data_rna_mean_1000_ratio_mcmc_ini")
ini <- Aug4_3Types_Data_ratio_ini


###----------------------------------------------------------------------###
k0 <- 2
i0 <- nrow(f10ratio_acc)
l0 <- nrow(f10ratio_exp)
d0 <- nrow(f10ratio_met)
p0 <- ncol(f10ratio_acc)


niter <- 200


rij <- f1/(f1+f0)
replacef0 <- 1-rij

rlg <- g1/(g1+g0)
replaceg0 <- 1-rlg


ratio <- met_data/(1-met_data) #h1/h0
c <- ratio[which(ratio!=0)]

sf <- quantile(c,probs = 0.5,na.rm=T) # median
sf
rdm <- ratio/sf
rdm <- rdm^2 # power
ind <- which(rdm!="Inf" & rdm!="NA")
rdm[is.na(rdm)] <- mean(rdm[ind])
rdm[which(rdm=="Inf")] <- 3000
rdm[1:5,1:5]

replaceh0 <- matrix(1, nrow = d0, ncol = p0)



alpha1 <- 2
beta1 <- 2


alpha_qi <- 1
beta_qi <- 1


### beta regression
###----------------------------------------------------------------------###
phi_1 <- 2.683904
eta <- -1.190329    
gamma <-  4.376499
tau <- -3.036440


phi_2 <- 3.18635
delta <- 0.1167477
theta <- 0.7307160




### initialization
###----------------------------------------------------------------------###
phi_atac <- rep(1/k0,k0)
phi_rna <- rep(1/k0,k0)
phi_met <- rep(1/k0,k0)

w <- ini$omega_exp
w_acc <- ini$omega_acc
w_acc <- w
w_exp <- ini$omega_exp
w_exp <- w
w_met <- ini$omega_met
w_met <- w

pi_rna <- ini$pi_exp
pi_met <- ini$pi_met
qi <- ini$qi_acc



```


### 2.2 run the EM algorithm
```{r}
print("start EM")
start_time <- Sys.time()
temp<-cal_M_step(phi_atac=phi_atac,w_acc=w_acc,qi=qi,
                  phi_rna=phi_rna,w_exp=w_exp,pi_rna=pi_rna,
                  phi_met=phi_met,w_met=w_met,pi_met=pi_met)

for(i in c(1:niter)){
  temp<-cal_M_step(phi_atac=temp$phi_atac,w_acc=temp$w_acc,qi=temp$qi,
                    phi_rna=temp$phi_rna,w_exp=temp$w_exp,pi_rna=temp$pi_rna,
                    phi_met=temp$phi_met,w_met=temp$w_met,pi_met=temp$pi_met)
  print(i)
  
}

pst <- cal_post(phi_atac=temp$phi_atac,w_acc=temp$w_acc,qi=temp$qi,
                phi_rna=temp$phi_rna,w_exp=temp$w_exp,pi_rna=temp$pi_rna,
                phi_met=temp$phi_met,w_met=temp$w_met,pi_met=temp$pi_met)
pst

end_time <- Sys.time()
end_time - start_time

print("end EM")

```


### 2.3 Summary clustering result (get the cluster assignments)
```{r}
prob_z_atac<-cal_E_acc(rij,replacef0,temp$phi_atac,temp$w_acc,temp$qi)$prob_z_atac
prob_z_rna<-cal_E_rna(rlg,replaceg0,temp$phi_rna,temp$w_exp,temp$pi_rna)$prob_z_rna
prob_z_met<-cal_E_rna(rdm,replaceh0,temp$phi_met,temp$w_met,temp$pi_met)$prob_z_rna

# amp, rmp and mmp are the cluster assignments offered by scAMACE
# scCAS data
amp<-rep(0,i0)
for(i in c(1:i0)){
  amp[i]=which.max(prob_z_atac[i,])
}
amp[1:10]

# scRNA-Seq data
rmp<-rep(0,l0)
for(l in c(1:l0)){
  rmp[l]=which.max(prob_z_rna[l,])
}
rmp[1:10]

# sc-methylation data
mmp<-rep(0,d0)
for(d in c(1:d0)){
  mmp[d]=which.max(prob_z_met[d,])
}
mmp[1:10]


### contingency table
table(atac_cell_lb,amp)

table(rna_cell_lb,rmp)

table(met_cell_lb,mmp)

```


## <a name="section3"></a>3. Example of scAMACE(seperate)
Remarks: We demostrate usage of scAMACE through Application 1.
### 3.1 Load data and prepare for EM algorithm
```{r}
library(scAMACE)

### load the data
data('Feb7_2021_3Types_Data_rna_mean_1000_ratio')
data <- Feb7_2021_3Types_Data_rna_mean_1000_ratio


f10ratio_acc <- data$f10ratio_acc
f1 <- data$f1
f0 <- data$f0
atac_cell_lb <- data$atac_cell

f10ratio_exp <- data$f10ratio_exp
g1 <- data$g1
g0 <- data$g0
rna_cell_lb <- data$rna_cell

f10ratio_met <- data$f10ratio_met
met_data <- data$met_data
met_cell_lb <- data$met_cell


### load the initialization
data("Feb7_2021_3Types_Data_rna_mean_1000_ratio_mcmc_ini")
ini <- Aug4_3Types_Data_ratio_ini


###----------------------------------------------------------------------###
p0 <- ncol(f10ratio_acc)
niter <- 200

```


### 3.2 scCAS data
```{r}
# scCAS
#####
##-------------------------------------------------------------------------##
k0 <- 2 # change to k0=1 if you want to get w_acc for beta regression

i0 <- nrow(f10ratio_acc)

rij <- f1/(f1+f0)
replacef0 <- 1-rij

alpha1 <- 2
beta1 <- 2

alpha_qi <- 1
beta_qi <- 1


# initialization
phi_atac <- rep(1/k0,k0)
w_acc <- ini2$omega_acc
qi <- ini2$qi_acc


print("start ATAC EM")
start_time <- Sys.time()
temp<-cal_M_step_atac(phi_atac=phi_atac,w_acc=w_acc,qi=qi)


for(i in c(1:niter)){
  temp<-cal_M_step_atac(phi_atac=temp$phi_atac,w_acc=temp$w_acc,qi=temp$qi)
  print(i)

  }

pst <- cal_post_atac(phi_atac=temp$phi_atac,w_acc=temp$w_acc,qi=temp$qi)
print(pst)

end_time <- Sys.time()
end_time - start_time
print(end_time - start_time)
print("end ATAC EM")

prob_z_atac<-cal_E_acc(rij,replacef0,temp$phi_atac,temp$w_acc,temp$qi)$prob_z_atac

amp<-rep(0,i0)
for(i in c(1:i0)){
  amp[i]=which.max(prob_z_atac[i,])
}

table(atac_cell_lb,amp)

```


### 3.3 scRNA-Seq data
```{r}
# scRNA-Seq
#####
##-------------------------------------------------------------------------##
k0 <- 2 # change to k0=1 if you want to get w_rna for beta regression

l0 <- nrow(f10ratio_exp)

rlg <- g1/(g1+g0)
replaceg0 <- 1-rlg


alpha1 <- 2
beta1 <- 2


# initialization
phi_rna <- rep(1/k0,k0)
w_exp <- ini$omega_exp
pi_rna <- ini$pi_exp


print("start RNA EM")
start_time <- Sys.time()
temp<-cal_M_step_rna(phi_rna=phi_rna,w_exp=w_exp,pi_rna=pi_rna)

for(i in c(1:niter)){
  temp<-cal_M_step_rna(phi_rna=temp$phi_rna,w_exp=temp$w_exp,pi_rna=temp$pi_rna)
  print(i)

  }

pst <- cal_post_rna(phi_rna=temp$phi_rna,w_exp=temp$w_exp,pi_rna=temp$pi_rna)
print(pst)

end_time <- Sys.time()
end_time - start_time
print(end_time - start_time)
print("end RNA EM")


prob_z_rna<-cal_E_rna(rlg,replaceg0,temp$phi_rna,temp$w_exp,temp$pi_rna)$prob_z_rna

rmp<-rep(0,l0)
for(l in c(1:l0)){
  rmp[l]=which.max(prob_z_rna[l,])
}

table(rna_cell_lb,rmp)

```




### 3.4 sc-methylation data
```{r}
# sc-methylation
#####
##-------------------------------------------------------------------------##
k0 <- 2 # change to k0=1 if you want to get w_met for beta regression

d0 <- nrow(f10ratio_met)

ratio <- met_data/(1-met_data) #h1/h0
c <- ratio[which(ratio!=0)]

sf <- quantile(c,probs = 0.5,na.rm=T) # median
sf
rdm <- ratio/sf

# power
rdm <- rdm^2 # power
ind <- which(rdm!="Inf" & rdm!="NA")
rdm[is.na(rdm)] <- mean(rdm[ind]) ### modified
rdm[which(rdm=="Inf")] <- 3000
rdm[1:5,1:5]

replaceh0 <- matrix(1, nrow = d0, ncol = p0)


alpha1 <- 2
beta1 <- 2


# initialization
phi_met <- rep(1/k0,k0)
w_met <- ini$omega_met
pi_met <- ini$pi_met


print("start methy. EM")
start_time <- Sys.time()
temp<-cal_M_step_met(phi_met=phi_met,w_met=w_met,pi_met=pi_met)

initialpostprob <- temp$postprob

postprob <- c()

for(i in c(1:niter)){
  temp<-cal_M_step_met(phi_met=temp$phi_met,w_met=temp$w_met,pi_met=temp$pi_met)
  postprob <- c(postprob,temp$postprob)
  print(i)
  }

pst <- cal_post_met(phi_met=temp$phi_met,w_met=temp$w_met,pi_met=temp$pi_met)
pst

end_time <- Sys.time()
end_time - start_time
print(end_time - start_time)
print("end methy. EM")


prob_z_met<-cal_E_rna(rdm,replaceh0,temp$phi_met,temp$w_met,temp$pi_met)$prob_z_rna

mmp<-rep(0,d0)
for(d in c(1:d0)){
  mmp[d]=which.max(prob_z_met[d,])
}

table(met_cell_lb,mmp)

```


## <a name="section4"></a>4. Example of beta regression
Remarks: We demostrate usage of beta regression through Application 1.

We first set $K = 1$ and use the model to estimate $\omega_{kg}^{rna}$, $\omega_{kg}^{acc}$ and $\omega_{kg}^{met}$ separately ([Section 3](#section3)), and then fix $\omega_{kg}^{rna}$, $\omega_{kg}^{acc}$ and $\omega_{kg}^{met}$ to estimate \{$\eta, \gamma, \tau, \delta, \theta, \phi^{acc}, \phi^{met}$\} by beta regression.

### 4.1 Load data
```{r}
library(scAMACE)
library(betareg) # for beta regression

# result obtained by setting k0=1 and 'cal_M_step_met'
temp_met <- temp

# result obtained by setting k0=1 and 'cal_M_step_atac'
temp_acc <- temp

# result obtained by setting k0=1 and 'cal_M_step_rna'
temp_rna <- temp


w_acc <- temp_acc$w_acc
w_exp <- temp_rna$w_exp
w_met <- temp_met$w_met

po <- ncol(w_acc)

w_acc_link <- w_acc[,1:po]
w_exp_link <- w_exp[,1:po]
w_met_link <- w_met[,1:po]

ca <- cut(w_acc[,1:po],breaks=seq(0,1,0.1))
cr <- cut(w_exp[,1:po],breaks=seq(0,1,0.1))
cm <- cut(w_met[,1:po],breaks=seq(0,1,0.1))


```


### 4.2 run beta regression
```{r}
###----------------------------------------------------------------------###
## w_acc 
acc_logit <- get_logit_beta(w_exp_link,w_acc_link,"quadratic")
acc_logit$reg_para # eta, gamma, tau
acc_logit$prec_para # phi_1


###----------------------------------------------------------------------###
## w_met
met_logit <- get_logit_beta(w_exp_link,w_met_link,"linear")
met_logit$reg_para # delta,theta
met_logit$prec_para # phi_2


```


## <a name="section5"></a>5. Example of generating simulation data

```{r}
library(scAMACE)
library(magic) # to generate differential w_exp

n1 <- 900
n2 <- 1100
n3 <- 1000
p <- 1000


prop <- 0.05
omega <- 0.8

k1 <- 3
k2 <- 3
k3 <- 3

# you may set some entries in Xprobs to be zero to simulate the unequal number of clusters
Xprobs1 <- c(1/k1, 1/k1, 1/k1)
Xprobs2 <- c(1/k2, 1/k2, 1/k2)
Xprobs3 <- c(1/k3, 1/k3, 1/k3)

shape_rna <- c(7, 1)
scale_rna <- c(0.5, 1)

shape_met1 <- c(0.5, 1)
shape_met2 <- c(0.5, 10)


qi <- rbind(rep(0.2, n1), 0)
pi_rna <- rbind(rep(0.7, n2), 0.3)
pi_met <- rbind(rep(0.4, n3), 0.7)


alpha1 <- beta1 <- 2

phi_1 <- 10
eta <- -1
gamma <-  7
tau <- -2

phi_2 <- 10
delta <- -2
theta <-  5


sim_data <- simData_3data(n1=n1, n2=n2, n3=n3, p=p,
                          Xprobs1=Xprobs1,Xprobs2=Xprobs2,Xprobs3=Xprobs3,
                          shape_rna=shape_rna, scale_rna=scale_rna,
                          shape_met1=shape_met1, shape_met2=shape_met2,
                          qi=qi,pi_rna=pi_rna,pi_met=pi_met,  cutoff=10^-6,
                          alpha1=alpha1,beta1=beta1,phi_1=phi_1,phi_2=phi_2,
                          omega=omega,prop=prop)
ls(sim_data)


```



