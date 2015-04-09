library(LatticeKrig)

fit.Matern.GP <- function(y,X,s,nu){
	
  ## Create a Sequence for Spatial Range
  D <- rdist(s)
  max.dist <- max(D)
  sr.seq <- seq(0.001,max.dist,length=10)
  for(i in 1:length(sr.seq)){
    sr.seq[i] <- 1/Matern.cor.to.range(sr.seq[i],nu=nu,cor.target=0.05)
  }
  
  ## Create a Sequence for %Spatial
  pct.spatial <- seq(0,.95,length=25)
  
  ## log-likelihoods
  ll <- matrix(0,nrow=length(sr.seq),ncol=length(pct.spatial))
  s2.hats <- ll
  beta.hats <- matrix(0,nrow=ncol(X),ncol=length(sr.seq)*length(pct.spatial))
  the.it <- 0
  for(sr in sr.seq){
    for(pct in pct.spatial){
      the.it <- the.it +1
      R <- pct*Matern(D,nu=nu,alpha=sr)+(1-pct)*diag(nrow(s))
      R.chol <- t(chol(R))
      first.piece <- forwardsolve(R.chol,X)
      XpRinvX <- t(first.piece)%*%first.piece
      last.piece <- forwardsolve(R.chol,y)
      beta.hat <- solve(XpRinvX)%*%t(first.piece)%*%last.piece
      beta.hats[,the.it] <- beta.hat
      ss <- forwardsolve(R.chol,y-X%*%beta.hat)
      ss <- t(ss)%*%ss
      s2.hat <- ss/nrow(s)
      s2.hats[sr==sr.seq,pct==pct.spatial] <- s2.hat
      log.like <- -0.5*nrow(s)*log(s2.hat)-sum(log(diag(R.chol)))-0.5*ss/(s2.hat)
      ll[sr==sr.seq,pct==pct.spatial] <- log.like
    }
  }
  
  ## Find maximum likelihood estimate
  the.mle <- which(ll==max(ll),arr.ind=TRUE)
  sr <- sr.seq[the.mle[1,1]]
  pct <- pct.spatial[the.mle[1,2]]
  the.mle <- which(ll==max(ll))
  s2.hat <- s2.hats[the.mle]
  beta.hat <- beta.hats[,the.mle]
  
  ## Return List
  return(list(phi=sr,sigma2=pct*s2.hat,tau2=(1-pct)*s2.hat,beta.hat=beta.hat))
	
}

	