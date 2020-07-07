######################################
##      block minimizer             ##
######################################
soft_threshold = function(x, lambda){
  sign(x) * (abs(x) >= lambda) * (abs(x) - lambda)
}

updateblock = function(eta,tau,k,lambda1,lambda2){
  minimizer = rep(0, k)
  norm_block_soft_thre = sqrt(sum(soft_threshold(tau, lambda1*eta)^2))
  if(norm_block_soft_thre > (lambda2 * eta)){
    minimizer = (1 - lambda2*eta/norm_block_soft_thre) * soft_threshold(tau, lambda1*eta)
  }
  return(minimizer)
}

############################
##Initialisation of E-STEP##
############################

ini.w <- function(k,n,weight=0.9){
  # k: nr. of components
  # n: nr. of observations
  
  EX <- matrix(NA,nrow=n,ncol=k)
  w <- sample(1:k,n,replace=TRUE)
  for (i in 1:n){
    EX[i,w[i]] <- weight
    EX[i,-w[i]] <- 1-weight
    EX[i,] <- EX[i,]/sum(EX[i,])
  }
  EX
}


#######################################################################
## FMRSGL                                                         ##
#######################################################################

fmrsgl<- function(x,y,k,lambda,alpha,ssd.ini,w.ini,
                     term=10^{-6},maxiter=2000,actiter=10,warnings=TRUE){
  
  # ARGUMENTS:
  # x: design matrix (with intercept)
  # y: response vector
  # k: nr of components
  # lambda: 
  # ssd.ini: initialisation for variances (!=0)
  # w.ini: initialisation of aposterior probabilies (n*k matrix)
  # T: termination criterion (see Philip Gill)
  # maxiter: maximal nr of iterations
  # warnings: should warnings (plik not reduced, loglik=-INF) be showed?
  # actiter=10: update active set every 10th EM-iteration

  lambda1 = lambda*alpha
  lambda2 = lambda*(1-alpha)*sqrt(k)
  
  n <- length(y)
  p <- dim(x)[2]
  prob <- rep(1/k,k)
  eta <- matrix(0,nrow=p,ncol=k)
  tau <- 1/rep(ssd.ini,k)
  w <- w.ini
  act <- vector()#to store the active set
  dnregr <- matrix(0,nrow=n,ncol=k)
  
  i <- 0
  err1 <- Inf #convergence of parameters
  err2 <- Inf #convergence of plik
  conv <- FALSE
  plik <- Inf #sets plik to Inf
  theta <- Inf #sets the vector of estimated parameters to Inf 
  warn <- FALSE
  allcoord <- TRUE
  act.iter <- actiter
  
  yy = y*y
  while ((!conv)&(i<maxiter)&!warn){
    #M-STEP: approximate minimization
    #update prob (pi_1, ..., pi_k)
    prob = colSums(w)/n
    #update tau = (tau_1,...,tau_k)
    yhat = x%*%eta
    yyhat = apply(yhat,2,function(v){v*y})
    wyyhat = diag(crossprod(w,yyhat))
    wyy = crossprod(w,yy)
    sj = colSums(w)
    tau = (wyyhat + sqrt(wyyhat^2 + 4*sj*wyy))/(2*wyy)
    ###############################################################print(wyy)
    #update intercept eta_0j = eta[1,]
    tauy = apply(as.matrix(tau), 1, function(a){a*y})
    xeta = x[,-1]%*%eta[-1,]
    H0 = tauy  - xeta
    eta[1,] = diag(crossprod(w,H0))/sj
    
    ###############################################################print(eta[1,])
    #update eta[-1,] = (eta_2. , ..., eta_p.) (exclude intercepts)
    rootw = sqrt(w)
    if (act.iter==actiter){
      act.iter <- 0
      allcoord <- TRUE
    }else {
      allcoord <- FALSE
    }
    
    if (allcoord){
      act = vector()
      for (l in 2:p){
        R = sweep(tauy,2,eta[1,])
        A = rootw*(R - xeta+x[,l,drop=FALSE]%*%eta[l,,drop=FALSE])
        rootwx.l = apply(rootw, 2, function(v){v*x[,l]})
        rootwx.l2 = rootwx.l^2
        colsm = colSums(rootwx.l2)
        # compute jacobian
        term1 = colsm*eta[l,]
        term2 = diag(crossprod(rootwx.l, A))
        jacobian = 2*(term1-term2)
        t = 2*max(colsm)
        eta[l,] = updateblock(2*n/t, eta[l,]-jacobian/t,k,lambda1,lambda2)
        if(sum(abs(eta[l,]) != 0)){
          act = c(act,l)
        }
      }
    }else{
      act.iter <- act.iter+1
      for (l in act){
        R = sweep(tauy,2,eta[1,])
        A = rootw*(R - xeta+x[,l,drop=FALSE]%*%eta[l,,drop=FALSE])
        rootwx.l = apply(rootw, 2, function(v){v*x[,l]})
        rootwx.l2 = rootwx.l^2
        colsm = colSums(rootwx.l2)
        # compute jacobian
        term1 = colsm*eta[l,]
        term2 = diag(crossprod(rootwx.l, A))
        jacobian = 2*(term1-term2)
        t = 2*max(colsm)
        eta[l,] = updateblock(2*n/t, eta[l,]-jacobian/t,k,lambda1,lambda2)
      }
    }  
    
    
    
    #E-STEP
    xetac <- x%*%eta
    for (j in 1:k){
      dnregr[,j] <- dnorm(y,mean=xetac[,j]/tau[j],sd=1/tau[j])
    }
    probdnregr <- t(prob*t(dnregr))
    dmix <- rowSums(probdnregr)
    w <- probdnregr/dmix
    
    #convergence criterion of theta
    thetaold <- theta
    theta <- c(as.vector(eta),as.vector(tau),prob)
    err1 <- max(abs(theta-thetaold)/(1+abs(theta)))#Philip Gill termination criteria
    
    #loglik    
    loglik <- sum(log(dmix))
    if (loglik==-Inf){warning("bad starting value (loglik=-Inf)");warn <- TRUE;break}
    
    #convergence criterion of plik    
    plikold <- plik
    plik <- -loglik/n + lambda1*sum(rowSums(abs(eta[-1,,drop=FALSE]))) + lambda2*sum(sqrt(rowSums((eta[-1,,drop=FALSE])^2)))
    #if(((plik-plikold)>0)&(i>0)) {  #is plik reduced? remark: algorithm is constructed to reduce plik.
      #if (warnings) {
        #break
        #warning("Warning: penalized negative loglik not reduced, lambda is too small")
      #}
    #}
    err2 <- abs(plik-plikold)/(1+abs(plik))
    
    #converged?
    conv <- ((err1<sqrt(term))&(err2<term))
    i <- i+1
  }
  
  # prepare result
  beta= t(apply(eta,1,function(v){v/tau}))
  clust <- apply(-w,1,which.min)
  list(coef = beta, ssd = as.vector(1/tau), prob = prob, niter = i, weight = w, 
       dmix = dmix,probd = probdnregr,component = clust,lambda = lambda,xeta =xetac)
}




fmrsglpath <- function(x,y,k,lambda,alpha,ssd.ini,w.ini,term=10^{-4},maxiter=500,actiter=10){
  
  # ARGUMENTS:
  # x: design matrix (with intercept)
  # y: response vector
  # k: nr of components  
  # lambda:
  # gamma:
  # ssd.ini: initialisation for variances (!=0)
  # ex.ini: initialisation of aposterior probabilies (n*k matrix)
  # T: termination criterion
  # maxiter: maximal nr of iterations
  
  n <- length(y)
  p <- dim(x)[2]-1 
  l.la <- length(lambda)
  
  prob <- matrix(0,nrow=k,ncol=l.la)
  coef <- array(0,dim=c(p+1,k,l.la))
  ssd <- matrix(0,nrow=k,ncol=l.la)
  #plik <- rep(0,l.la)
  weight <- array(0,dim=c(n,k,l.la))
  component <- matrix(0,nrow=n,ncol=l.la)
  niter <- rep(0,l.la)
  
  for (i in 1:l.la){
    
    fit <- fmrsgl(x,y,k,lambda=lambda[i],alpha=alpha,ssd.ini=ssd.ini,w.ini=w.ini,term=term,maxiter=maxiter,actiter=actiter)
    prob[,i] <- fit$prob
    coef[,,i] <- fit$coef
    ssd[,i] <- fit$ssd
    #plik[i] <- fit$plik
    weight[,,i] <- fit$weight
    #component[,i] <- fit$component
    niter[i] <- fit$niter
    #cat('lambda:',lambda[i],'\n')
  }
  
  #prepare results
  dimnames(coef) <- list(coef=0:p,comp=1:k,lambda=lambda)
  dimnames(coef)[[1]][1] <- 'intercept'
  dimnames(prob) <- dimnames(ssd) <- list(comp=1:k,lambda=lambda)
  dimnames(weight) <- list(NULL,comp=1:k,lambda=lambda)
  #dimnames(component) <- list(NULL,lambda=lambda)
  
  res <- list(k=k,lambda=lambda,prob=prob,coef=coef,ssd=ssd,weight=weight,component = component,niter=niter)
  res
  
}

###########
##LOGLOSS##
###########

loglosssgl <- function(coef,ssd,prob,x,y){
  k <- length(prob)
  n <- length(y)
  dnregr <- matrix(0,nrow=n,ncol=k)
  xcoef <- x%*%as.matrix(coef)
  for (j in 1:k){
    dnregr[,j] <- dnorm(y,mean=xcoef[,j],sd=ssd[j])
  }
  probdnregr <- t(prob*t(dnregr))
  dmix <- rowSums(probdnregr)
  -1*sum(log(dmix))
}

#################################################
##PREDICTION ERROR FOR TEST DATA: X.TEST,Y.TEST##
#################################################

predlosssgl <- function(model,x,y){

  coef <- model$coef
  ssd <- model$ssd
  prob <- model$prob
  loss <- loglosssgl(coef,ssd,prob,x,y)
  loss
}

predlosspathsgl <- function(modelseq,x,y){
  
  l <- length(modelseq$lambda)
  coef <- modelseq$coef
  ssd <- modelseq$ssd
  prob <- modelseq$prob
  loss <- rep(0,l)
  for (i in 1:l){
    loss[i] <- loglosssgl(as.matrix(coef[,,i]),ssd[,i],prob[,i],x,y)
  }
  indexopt=which.min(loss)
  optlambda = modelseq$lambda[indexopt]
  optfit = list(coef = as.matrix(coef[,,indexopt]),ssd = as.vector(ssd[,indexopt]), 
                prob= prob[,indexopt],optlambda = optlambda)
  pred <- list(loss=loss,optfit = optfit,optlambda = optlambda)
  pred
}