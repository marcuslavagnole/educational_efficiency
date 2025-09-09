#' Spatial stochastic frontier models for panel data
#'
#' This function estimates a Spatial stochastic frontier model following the
#' Bayesian paradigm.
#'
#' @param y (nxt)-dimensional vector of responses.
#' @param x Matrix - (nxt)x(p+1) - of predictors (include the intercept).
#' @param x_u Matrix - (nxt)xq - of explanatory variables for the inefficiencies.
#' @param idobs (nxt)-dimensional vector that indicates the n units (categorical 
#' numerical, from 1 to n).
#' @param W Adjacency matrix.
#' @param n Number of units.
#' @param t Number of periods.
#' @param n_mcmc Number of iterations.
#' @param burnin_mcmc Number of initial iterations to be discarded.
#' @param thin_mcmc Thinning parameter.
#' 
#' @return A list with the chains of all parameters of interest.

## Packages
require(mvtnorm); require(msm)

## MCMC
spatial_sfm <- function(y,x,x_u,idobs,W,n,t,n_mcmc,burnin_mcmc,thin_mcmc){
  numcov_y <- ncol(x)
  numcov_u <- ifelse(is.null(ncol(x_u)),1,ncol(x_u))
  resultado <- list()
  # Create auxiliary objects
  m_alpha <- matrix(0,n*t,n)
  for(l in 1:(n*t)){
    m_alpha[l,as.numeric(idobs[l])] <- 1 
  }
  beta  <- matrix(NA, n_mcmc, numcov_y)
  sigma2<- matrix(NA, n_mcmc, 1)
  u     <- matrix(NA, n_mcmc, n*t)
  theta <- matrix(NA, n_mcmc, numcov_u)
  tau2  <- matrix(NA, n_mcmc, 1)
  eta   <- matrix(NA, n_mcmc, t)
  alpha <- matrix(NA, n_mcmc, n)
  psi2  <- matrix(NA, n_mcmc, 1)
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
  for(k in 2:n_mcmc){
    alpha_aux  <- m_alpha%*%alpha[k-1,]
    beta[k,]   <- atualizarBETA(rep(0,numcov_y),diag(100,numcov_y),y,x,u[k-1,],sigma2[k-1,1])
    sigma2[k,1]<- atualizarSIGMA2(0.01,0.01,y,x,u[k-1,],beta[k,],n,t)
    u[k,]      <- atualizarU(y,x,beta[k,],alpha_aux,x_u,theta[k-1,],sigma2[k,1],tau2[k-1,1],n,t)
    theta[k,]  <- atualizarTHETA(rep(0,numcov_u),diag(100,numcov_u),u[k,],x_u,alpha_aux,tau2[k-1,1])
    tau2[k,1]  <- atualizarTAU2(0.01,0.01,u[k,],x_u,alpha_aux,theta[k,],n,t)
    alpha[k,]  <- atualizarALPHA(psi2[k-1,1],u[k,],x_u,alpha[k-1,],theta[k,],tau2[k,1],W,t(m_alpha),n,t)
    psi2[k,1]  <- atualizarPSI2(0.1,0.1,W,alpha[k,],n)
  }
  resultado[['beta']]  <- beta[seq(burnin_mcmc+1,n_mcmc,thin_mcmc),]
  resultado[['sigma2']]<- sigma2[seq(burnin_mcmc+1,n_mcmc,thin_mcmc),1]
  resultado[['u']]     <- u[seq(burnin_mcmc+1,n_mcmc,thin_mcmc),]
  resultado[['theta']] <- theta[seq(burnin_mcmc+1,n_mcmc,thin_mcmc),]
  resultado[['tau2']]  <- tau2[seq(burnin_mcmc+1,n_mcmc,thin_mcmc),1]
  resultado[['alpha']] <- alpha[seq(burnin_mcmc+1,n_mcmc,thin_mcmc),]
  resultado[['psi2']]  <- psi2[seq(burnin_mcmc+1,n_mcmc,thin_mcmc),1]
  
  return(resultado)
}


#######################################################################
## Auxiliary functions: Sampling from full conditional distributions ##
#######################################################################

# Full conditional distribution for beta
atualizarBETA<-function(b,B,y,x,u,sigma2){
  B.inv<- chol2inv(chol(B))
  sigma<- chol2inv(chol(B.inv+(1/sigma2)*(t(x)%*%x)))
  v    <- y+u
  v[is.na(v)] <- 0
  media<- sigma%*%((B.inv%*%b)+(1/sigma2)*(t(x)%*%v))
  beta <- rmvnorm(1,media,sigma)
  return(beta)
}

# Full conditional distribution for sigma2
atualizarSIGMA2<-function(c,C,y,x,u,beta,n,t){
  alpha1<- c + 0.5*n*t
  v <- y-x%*%beta+u
  v[is.na(v)] <- 0
  beta1 <- C + 0.5*(t(v)%*%v)
  sigma2<-1/rgamma(1, alpha1, beta1)
  return(sigma2)
}

# Full conditional distribution for inefficiencies u
atualizarU<-function(y,x,beta,alpha,z,theta,sigma2,tau2,n,t){
  v1   <- y - x%*%beta 
  v2   <- alpha + z*theta 
  media<- (sigma2*v2 - tau2*v1)/(sigma2+tau2)
  sd   <- sqrt(tau2*sigma2/(sigma2+tau2))
  N    <- n*t
  u    <- rtnorm(N,media,sd,lower=0,upper=Inf)
  return(u)
}

# Full conditional distribution for theta
atualizarTHETA<-function(b,B,u,z,alpha,tau2){
  B.inv<- chol2inv(chol(B))
  sigma<- chol2inv(chol(B.inv+(1/tau2)*(t(z)%*%z)))
  v    <- u - alpha
  v[is.na(v)] <- 0
  media<- sigma%*%((B.inv%*%b)+(1/tau2)*(t(z)%*%v))
  theta<- rmvnorm(1,media,sigma)
  return(theta)
}

# Full conditional distribution for tau2
atualizarTAU2<-function(c,C,u,z,alpha,theta,n,t){
  alpha1<- c + 0.5*n*t
  v <- u-alpha-z*theta
  v[is.na(v)] <- 0
  beta1 <- C + 0.5*(t(v)%*%v)
  tau2  <- 1/rgamma(1, alpha1, beta1)
  return(tau2)
}

# Full conditional distribution for the spatial effects alpha
atualizarALPHA<-function(psi2,u,z,alpha,theta,tau2,m_W,m_aux,n,t){
  v      <- u - z*theta 
  v[is.na(v)] <- 0
  n_front<- rowSums(m_W,na.rm=TRUE) 
  alpha.aux <- alpha
  for (j in sample(1:n)){
    sd        <- sqrt(psi2*tau2/(n_front[j]*tau2+t*psi2))
    media     <- (tau2*(m_W[j,]%*%alpha) + psi2*(m_aux[j,]%*%v))/(n_front[j]*tau2+t*psi2)
    alpha.aux[j] <- rnorm(1,media,sd)
  }
  alpha.final <- alpha.aux - mean(alpha.aux)
  return(alpha.final)
}

# Full conditional distribution for psi2
atualizarPSI2<-function(c,C,m_W,alpha,n){
  m1 = (matrix(rep(alpha,n),n,n,byrow=TRUE) - matrix(rep(alpha,n),n,n))^2
  m2 = m_W*(lower.tri(m1, diag = FALSE)*m1)
  alpha1<- c + 0.5*(n-1)
  beta1 <- C + 0.5*sum(m2) 
  psi2  <-1/rgamma(1, alpha1, beta1)
  return(psi2)
}

