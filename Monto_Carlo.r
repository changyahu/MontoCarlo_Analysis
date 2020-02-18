######################################################################################
# Author: Prof. Gordon Cheung (University of Auckland)                               #
# email: gordon.cheung@auckland.ac.nz                                                #
# Purpose: Monte Carlo simulation of conditional indirect effects from Mplus outputs #
# Model: Model 7                                                                     #
# Version: Beta 1.0                                                                  #
# WARNING: This is beta version. Please report all bugs to the author                #
###################################################################################### 


# load the boot package (only needed once per session)
library(MASS)

######## Initial Inputs ########
setwd("/Users/heng/HR/Monto Carlo/")  ## set working directory to "C:/Research/R/Example"
NoSP <- 7   ## Input number of additional parameters NoSP as 7
La1 = 3     ## Specify location of a1 in additional parameters as 3
La3 = 5     ## Specify location of a3 in additional parameters as 5
Lb = 1      ## Specify location of b in additional parameters as 1
LstdW = 7   ## Specify location of stdW in additional parameters as 7
################################

# Calculate number of estimated parameters from ctemp2.txt
est <- scan("./Ctemp2.txt")
NoEP = ((sqrt(1+8*length(est))-1)/2) 

# Calculate starting point for extraction
SPExt <- NoEP-NoSP+1 

# Read estimated parameters from ctemp1.txt
est <- scan("./Ctemp1.txt")
estcoeff <- est[SPExt:NoEP]  # Extract the additional parameters

# Read variance-covariance matrix of parameters from ctemp2.txt
Tech <- matrix(0, NoEP, NoEP)
Tech[outer(1:NoEP, 1:NoEP, '<=')] <- scan('./Ctemp2.txt')
Tech[lower.tri(Tech)] <- t(Tech)[lower.tri(Tech)]
Tech3 <- Tech[SPExt:NoEP, SPExt:NoEP]

# Monte Carlo Simulation of 5,000,000 samples
mcmc <- mvrnorm(n=5000000, mu=estcoeff, Sigma=Tech3, tol = 1e-6)


###### Specific for Model 7 ######
# Define estimated parameters for calculating indirect effects
a1 <- mcmc[,La1]
a3 <- mcmc[,La3]
b <- mcmc[,Lb]
stdW <- mcmc[,LstdW]

# Define Simulated Conditional Indirect Effect as ab (Model 7)
abM2 <- (a1+a3*-2*stdW)*b
abM1 <- (a1+a3*-1*stdW)*b
abM  <- (a1+a3*0*stdW)*b
abP1 <- (a1+a3*1*stdW)*b
abP2 <- (a1+a3*2*stdW)*b
IndexMM <- a3*b

# Capture estimated parameters for calculating indirect effects
Sa1 <- estcoeff[La1]
Sa3 <- estcoeff[La3]
Sb <- estcoeff[Lb]
SstdW <- estcoeff[LstdW]

# Calculate Sample Estimated Conditional Indirect Effect as est (Model 7)
estM2 <- (Sa1+Sa3*-2*SstdW)*Sb
estM1 <- (Sa1+Sa3*-1*SstdW)*Sb
estM <- (Sa1+Sa3*0*SstdW)*Sb
estP1 <- (Sa1+Sa3*1*SstdW)*Sb
estP2 <- (Sa1+Sa3*2*SstdW)*Sb
estIMM <- Sa3*Sb
##################################


#### Confidence Intervals and p-value ####

# Calculate Percentile Probability
if (quantile(abM2,probs=0.5)>0) {
  pM2 = 2*(sum(abM2<0)/5000000)
} else {
  pM2 = 2*(sum(abM2>0)/5000000)
}
if (quantile(abM1,probs=0.5)>0) {
  pM1 = 2*(sum(abM1<0)/5000000)
} else {
  pM1 = 2*(sum(abM1>0)/5000000)
}
if (quantile(abM,probs=0.5)>0) {
  pM = 2*(sum(abM<0)/5000000)
} else {
  pM = 2*(sum(abM>0)/5000000)
}
if (quantile(abP1,probs=0.5)>0) {
  pP1 = 2*(sum(abP1<0)/5000000)
} else {
  pP1 = 2*(sum(abP1>0)/5000000)
}
if (quantile(abP2,probs=0.5)>0) {
  pP2 = 2*(sum(abP2<0)/5000000)
} else {
  pP2 = 2*(sum(abP2>0)/5000000)
}
if (quantile(IndexMM,probs=0.5)>0) {
  pIMM = 2*(sum(IndexMM<0)/5000000)
} else {
  pIMM = 2*(sum(IndexMM>0)/5000000)
}


#### Percentile Confidence Intervals of Conditional Indirect Effects ####

PCI <- matrix(1:48, nrow = 6, dimnames = list(c("Mean-2sd","Mean-1sd","Mean","Mean+1sd","Mean+2sd","Index MM"), c("0.5%","2.5%","5%","Estimate","95%","97.5%","99.5%", "p-value")))

PCI[1,1] <- quantile(abM2,c(0.005))
PCI[1,2] <- quantile(abM2,c(0.025))
PCI[1,3] <- quantile(abM2,c(0.05))
PCI[1,4] <- estM2
PCI[1,5] <- quantile(abM2,c(0.95))
PCI[1,6] <- quantile(abM2,c(0.975))
PCI[1,7] <- quantile(abM2,c(0.995))
PCI[1,8] <- pM2

PCI[2,1] <- quantile(abM1,c(0.005))
PCI[2,2] <- quantile(abM1,c(0.025))
PCI[2,3] <- quantile(abM1,c(0.05))
PCI[2,4] <- estM1
PCI[2,5] <- quantile(abM1,c(0.95))
PCI[2,6] <- quantile(abM1,c(0.975))
PCI[2,7] <- quantile(abM1,c(0.995))
PCI[2,8] <- pM1

PCI[3,1] <- quantile(abM,c(0.005))
PCI[3,2] <- quantile(abM,c(0.025))
PCI[3,3] <- quantile(abM,c(0.05))
PCI[3,4] <- estM
PCI[3,5] <- quantile(abM,c(0.95))
PCI[3,6] <- quantile(abM,c(0.975))
PCI[3,7] <- quantile(abM,c(0.995))
PCI[3,8] <- pM

PCI[4,1] <- quantile(abP1,c(0.005))
PCI[4,2] <- quantile(abP1,c(0.025))
PCI[4,3] <- quantile(abP1,c(0.05))
PCI[4,4] <- estP1
PCI[4,5] <- quantile(abP1,c(0.95))
PCI[4,6] <- quantile(abP1,c(0.975))
PCI[4,7] <- quantile(abP1,c(0.995))
PCI[4,8] <- pP1

PCI[5,1] <- quantile(abP2,c(0.005))
PCI[5,2] <- quantile(abP2,c(0.025))
PCI[5,3] <- quantile(abP2,c(0.05))
PCI[5,4] <- estP2
PCI[5,5] <- quantile(abP2,c(0.95))
PCI[5,6] <- quantile(abP2,c(0.975))
PCI[5,7] <- quantile(abP2,c(0.995))
PCI[5,8] <- pP2

PCI[6,1] <- quantile(IndexMM,c(0.005))
PCI[6,2] <- quantile(IndexMM,c(0.025))
PCI[6,3] <- quantile(IndexMM,c(0.05))
PCI[6,4] <- estIMM
PCI[6,5] <- quantile(IndexMM,c(0.95))
PCI[6,6] <- quantile(IndexMM,c(0.975))
PCI[6,7] <- quantile(IndexMM,c(0.995))
PCI[6,8] <- pIMM

round(PCI, digits=4)


# Bias-Corrected Factor

zM2 = qnorm(sum(abM2<estM2)/5000000)
zM1 = qnorm(sum(abM1<estM1)/5000000)
zM = qnorm(sum(abM<estM)/5000000)
zP1 = qnorm(sum(abP1<estP1)/5000000)
zP2 = qnorm(sum(abP2<estP2)/5000000)
zIMM = qnorm(sum(IndexMM<estIMM)/5000000)

# Calculate Bias-Corrected Probability

if ((estM2>0 & min(abM2)>0) | (estM2<0 & max(abM2)<0)) {
  pbM2 = 0
} else if (qnorm(sum(abM2>0)/5000000)+2*zM2<0) { 
  pbM2 = 2*pnorm((qnorm(sum(abM2>0)/5000000)+2*zM2))
} else {
  pbM2 = 2*pnorm(-1*(qnorm(sum(abM2>0)/5000000)+2*zM2))
}

if ((estM1>0 & min(abM1)>0) | (estM1<0 & max(abM1)<0)) {
  pbM1 = 0
} else if (qnorm(sum(abM1>0)/5000000)+2*zM1<0) { 
  pbM1 = 2*pnorm((qnorm(sum(abM1>0)/5000000)+2*zM1))
} else {
  pbM1 = 2*pnorm(-1*(qnorm(sum(abM1>0)/5000000)+2*zM1))
}

if ((estM>0 & min(abM)>0) | (estM<0 & max(abM)<0)) {
  pbM = 0
} else if (qnorm(sum(abM>0)/5000000)+2*zM<0) { 
  pbM = 2*pnorm((qnorm(sum(abM>0)/5000000)+2*zM))
} else {
  pbM = 2*pnorm(-1*(qnorm(sum(abM>0)/5000000)+2*zM))
}

if ((estP1>0 & min(abP1)>0) | (estP1<0 & max(abP1)<0)) {
  pbP1 = 0
} else if (qnorm(sum(abP1>0)/5000000)+2*zP1<0) { 
  pbP1 = 2*pnorm((qnorm(sum(abP1>0)/5000000)+2*zP1))
} else {
  pbP1 = 2*pnorm(-1*(qnorm(sum(abP1>0)/5000000)+2*zP1))
}

if ((estP2>0 & min(abP2)>0) | (estP2<0 & max(abP2)<0)) {
  pbP2 = 0
} else if (qnorm(sum(abP2>0)/5000000)+2*zP2<0) { 
  pbP2 = 2*pnorm((qnorm(sum(abP2>0)/5000000)+2*zP2))
} else {
  pbP2 = 2*pnorm(-1*(qnorm(sum(abP2>0)/5000000)+2*zP2))
}

if ((estIMM>0 & min(IndexMM)>0) | (estIMM<0 & max(IndexMM)<0)) {
  pbIMM = 0
} else if (qnorm(sum(IndexMM>0)/5000000)+2*zIMM<0) { 
  pbIMM = 2*pnorm((qnorm(sum(IndexMM>0)/5000000)+2*zIMM))
} else {
  pbIMM = 2*pnorm(-1*(qnorm(sum(IndexMM>0)/5000000)+2*zIMM))
}


#### Bias-Corrected Confidence Intervals ####

BCCI <- matrix(1:48, nrow = 6, dimnames = list(c("Mean-2sd","Mean-1sd","Mean","Mean+1sd","Mean+2sd","Index MM"), c("0.5%","2.5%","5%","Estimate","95%","97.5%","99.5%","p-value")))

BCCI[1,1] <- quantile(abM2,probs=pnorm(2*zM2+qnorm(0.005)))
BCCI[1,2] <- quantile(abM2,probs=pnorm(2*zM2+qnorm(0.025)))
BCCI[1,3] <- quantile(abM2,probs=pnorm(2*zM2+qnorm(0.050)))
BCCI[1,4] <- estM2
BCCI[1,5] <- quantile(abM2,probs=pnorm(2*zM2+qnorm(0.950)))
BCCI[1,6] <- quantile(abM2,probs=pnorm(2*zM2+qnorm(0.975)))
BCCI[1,7] <- quantile(abM2,probs=pnorm(2*zM2+qnorm(0.995)))
BCCI[1,8] <- pbM2

BCCI[2,1] <- quantile(abM1,probs=pnorm(2*zM1+qnorm(0.005)))
BCCI[2,2] <- quantile(abM1,probs=pnorm(2*zM1+qnorm(0.025)))
BCCI[2,3] <- quantile(abM1,probs=pnorm(2*zM1+qnorm(0.050)))
BCCI[2,4] <- estM1
BCCI[2,5] <- quantile(abM1,probs=pnorm(2*zM1+qnorm(0.950)))
BCCI[2,6] <- quantile(abM1,probs=pnorm(2*zM1+qnorm(0.975)))
BCCI[2,7] <- quantile(abM1,probs=pnorm(2*zM1+qnorm(0.995)))
BCCI[2,8] <- pbM1

BCCI[3,1] <- quantile(abM,probs=pnorm(2*zM+qnorm(0.005)))
BCCI[3,2] <- quantile(abM,probs=pnorm(2*zM+qnorm(0.025)))
BCCI[3,3] <- quantile(abM,probs=pnorm(2*zM+qnorm(0.050)))
BCCI[3,4] <- estM
BCCI[3,5] <- quantile(abM,probs=pnorm(2*zM+qnorm(0.950)))
BCCI[3,6] <- quantile(abM,probs=pnorm(2*zM+qnorm(0.975)))
BCCI[3,7] <- quantile(abM,probs=pnorm(2*zM+qnorm(0.995)))
BCCI[3,8] <- pbM

BCCI[4,1] <- quantile(abP1,probs=pnorm(2*zP1+qnorm(0.005)))
BCCI[4,2] <- quantile(abP1,probs=pnorm(2*zP1+qnorm(0.025)))
BCCI[4,3] <- quantile(abP1,probs=pnorm(2*zP1+qnorm(0.050)))
BCCI[4,4] <- estP1
BCCI[4,5] <- quantile(abP1,probs=pnorm(2*zP1+qnorm(0.950)))
BCCI[4,6] <- quantile(abP1,probs=pnorm(2*zP1+qnorm(0.975)))
BCCI[4,7] <- quantile(abP1,probs=pnorm(2*zP1+qnorm(0.995)))
BCCI[4,8] <- pbP1

BCCI[5,1] <- quantile(abP2,probs=pnorm(2*zP2+qnorm(0.005)))
BCCI[5,2] <- quantile(abP2,probs=pnorm(2*zP2+qnorm(0.025)))
BCCI[5,3] <- quantile(abP2,probs=pnorm(2*zP2+qnorm(0.050)))
BCCI[5,4] <- estP2
BCCI[5,5] <- quantile(abP2,probs=pnorm(2*zP2+qnorm(0.950)))
BCCI[5,6] <- quantile(abP2,probs=pnorm(2*zP2+qnorm(0.975)))
BCCI[5,7] <- quantile(abP2,probs=pnorm(2*zP2+qnorm(0.995)))
BCCI[5,8] <- pbP2

BCCI[6,1] <- quantile(IndexMM,probs=pnorm(2*zIMM+qnorm(0.005)))
BCCI[6,2] <- quantile(IndexMM,probs=pnorm(2*zIMM+qnorm(0.025)))
BCCI[6,3] <- quantile(IndexMM,probs=pnorm(2*zIMM+qnorm(0.050)))
BCCI[6,4] <- estIMM
BCCI[6,5] <- quantile(IndexMM,probs=pnorm(2*zIMM+qnorm(0.950)))
BCCI[6,6] <- quantile(IndexMM,probs=pnorm(2*zIMM+qnorm(0.975)))
BCCI[6,7] <- quantile(IndexMM,probs=pnorm(2*zIMM+qnorm(0.995)))
BCCI[6,8] <- pbIMM

result <- round(BCCI, digits=4)
