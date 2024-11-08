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
  #v2   <- alpha + z%*%theta 
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
  #v <- u-alpha-z%*%theta
  v <- u-alpha-z*theta
  v[is.na(v)] <- 0
  beta1 <- C + 0.5*(t(v)%*%v)
  tau2  <- 1/rgamma(1, alpha1, beta1)
  return(tau2)
}

# Full conditional distribution for the spatial effects alpha
atualizarALPHA<-function(psi2,u,z,alpha,theta,tau2,m_W,m_aux,n,t){
  #v      <- u - z%*%theta 
  v      <- u - z*theta 
  v[is.na(v)] <- 0
  #n_front<- rowSums(m_W) 
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
