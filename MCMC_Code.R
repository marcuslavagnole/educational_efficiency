# Load libraries
library(BayesCR)
library(gtools)
library(msm)
library(mvtnorm)
library(coda)
library(sn)
library(Matrix)

# Set working directory to file location and load utils
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("Full_Conditionals.R")
load("dados.RData")
load("W_RS.RData")

# Define dependent and independent variables, and covariates that affect the inefficiencies
n <- length(base$ano)/2
t <- length(unique(base$ano))
dados         <- log(base$ideb)
covariaveis_y <- cbind(rep(1,n*t),rep(c(0,1),n),log(cbind(base$profaluno,base$fundebpadef,base$mediaindice)))
covariaveis_u <- log(base$pibpcdef)
numcov_y <-dim(covariaveis_y)[2]
#numcov_u<-dim(covariaveis_u)[2]
numcov_u<-1

# Create a matrix to indicate if the observation i belongs to unit j
m_alpha <- matrix(0,n*t,n)
for(l in 1:n){
  m_alpha[((l-1)*t+1):(l*t),l] <- 1 
}

# Define the number of iterations to run
NN = 100000

# Create auxiliary objects 
beta  <- matrix(NA, NN, numcov_y)
sigma2<- matrix(NA, NN, 1)
u     <- matrix(NA, NN, n*t)
theta <- matrix(NA, NN, numcov_u)
tau2  <- matrix(NA, NN, 1)
eta   <- matrix(NA, NN, t)
alpha <- matrix(NA, NN, n)
psi2  <- matrix(NA, NN, 1)

# Set the initial values
beta[1,]   <- rnorm(numcov_y,0.5,0.1)
sigma2[1,1]<- 1
u[1,]      <- rnorm(n*t,0.2,0.1)
theta[1,]  <- rnorm(numcov_u,0.5,0.1)
tau2[1,1]  <- 0.01
eta[1,]    <- rnorm(t,1,0.1)
alpha[1,]  <- rnorm(n,1,0.1)
psi2[1,1]  <- 1

# MCMC
for(k in 2:NN){
  alpha_aux  <- m_alpha%*%alpha[k-1,]
  beta[k,]   <- atualizarBETA(rep(0,numcov_y),diag(100,numcov_y),dados,covariaveis_y,u[k-1,],sigma2[k-1,1])
  sigma2[k,1]<- atualizarSIGMA2(0.01,0.01,dados,covariaveis_y,u[k-1,],beta[k,],n,t)
  u[k,]      <- atualizarU(dados,covariaveis_y,beta[k,],alpha_aux,covariaveis_u,theta[k-1,],sigma2[k,1],tau2[k-1,1],n,t)
  theta[k,]  <- atualizarTHETA(rep(0,numcov_u),diag(100,numcov_u),u[k,],covariaveis_u,alpha_aux,tau2[k-1,1])
  tau2[k,1]  <- atualizarTAU2(0.01,0.01,u[k,],covariaveis_u,alpha_aux,theta[k,],n,t)
  alpha[k,]  <- atualizarALPHA(psi2[k-1,1],u[k,],covariaveis_u,alpha[k-1,],theta[k,],tau2[k,1],W,t(m_alpha),n,t)
  psi2[k,1]  <- atualizarPSI2(0.1,0.1,W,alpha[k,],n)
}
