.construct.arglist <- function(funobj,framepos=-1) {
    
    # evaluates all function arguments at location where construct.arglist is called in that function
    
    # (they might have changed during running the function, or may be in variables)
    
    # construct.arglist gets rid of the variables and so on, and returns the argument list as a list
    
    namedlist=formals(funobj)
    
    argnames=names(namedlist)
    
    for (argn in 1:length(namedlist)) {
        
        if  (exists(argnames[argn],where=sys.frame(framepos))) {
            
            namedlist[[argn]]=get(argnames[argn],envir=sys.frame(framepos))
            
        }
        
    }
    
    return(namedlist)
    
}
MCMC_tvp <- function(Y,X,nburn,nsave,priorbtheta=list(B_1=2,B_2=1,kappa0=1e-07),priorb0=list(a_tau=0.1,c_tau=0.01,d_tau=0.01),priorsig=c(0.01,0.01),priorphi=c(2,2),priormu=c(0,10^2),h0prior="stationary",grid.length=150,thrsh.pct=0.1,thrsh.pct.high=1,sv_on=TRUE,TVS=TRUE,cons.mod=TRUE,nr,thin=0.1,robust=TRUE,a.approx=FALSE, sim.kappa=FALSE,kappa.grid = seq(1e-4, 0.1, 10)){
 # Y <- matrix(rnorm(100),100,1);X <- matrix(rnorm(100*2),100,2);nburn=nsave=100;priorbtheta=list(B_1=2,B_2=1,kappa0=1e-07);priorb0=list(a_tau=0.1,c_tau=0.01,d_tau=0.01);priorsig=c(0.01,0.01);priorphi=c(2,2);priormu=c(0,10^2);h0prior="stationary";grid.length=150;thrsh.pct=0.1;thrsh.pct.high=1;sv_on=TRUE;TVS=TRUE;cons.mod=FALSE;nr=1;thin=0.1;robust=TRUE;a.approx=FALSE; sim.kappa=TRUE
    #-----------------------------------------------------------------------------------------------------------#
    ntot <- nsave+nburn
    a_tau <- priorb0$a_tau
    c_tau <- priorb0$c_tau
    d_tau <- priorb0$d_tau
    
    
    B_gamma1 <- priorbtheta$B_1 #input new
    B_gamma2 <- priorbtheta$B_2 
    kappa0 <- priorbtheta$kappa0 

    a_sig <- priorsig[[1]]
    b_sig <- priorsig[[2]]
    #---------------------------------sourceCpp to source .cpp functions----------------------------------------#
    require(stochvol)
    require(GIGrvg)
    require(snowfall)
    #-----------------------------------------------------------------------------------------------------------#
    class(X) <-"matrix"
    class(Y) <- "matrix"
    T <- nrow(Y)
    K <- ncol(X)
    M <- ncol(Y)
    
    #Create full-data matrices
    if (nr==1) slct <- NULL else slct <- 1:(nr-1)
    Y__ <- Y[,nr]
    X__ <- cbind(Y[,slct],X)
    K_ <- ncol(X__)
    M_ <- M-length(slct)
    #Colnames checking
    if (is.null(colnames(X))) colnames(X) <- seq(1,ncol(X))
    #Initialize everything
    sqrttheta1 <- diag(K_)*0.1
    sqrttheta2 <-diag(K_)*0.01
    #mylm <- lm(Y__ ~ X__)
    
    B.OLS <- solve(crossprod(X__))%*%crossprod(X__,Y__)
    sig.OLS <- crossprod(Y__-X__%*%B.OLS)/(T-K_)
    V.OLS <- as.numeric(sig.OLS)*solve(crossprod(X__))
    
    #Specification of scaling parameter in the lower regime
    sd.OLS <- sqrt(diag(V.OLS))
    if (kappa0<0) kappa0 <- -kappa0 * sd.OLS else kappa0 <- matrix(kappa0,K_,1)
    
    ALPHA0 <- solve(crossprod(X__))%*%crossprod(X__,Y__)
    
    if (a.approx){
      buildCapm <- function(u){
        dlm::dlmModReg(X__, dV = exp(u[1]), dW = exp(u[2:(K_+1)]),addInt = FALSE)
      }
      
      outMLE <- dlm::dlmMLE(Y__, parm = rep(0,K_+1), buildCapm)
      mod <- buildCapm(outMLE$par)
      outS <- dlm::dlmSmooth(Y__, mod)
      states.OLS <- t(matrix(outS$s,T+1,K_))
      
      Achg.OLS <- t(diff(t(states.OLS)))#t(as.numeric(ALPHA0)+ALPHA2)
    }
    
    scale.0 <- 0
    Ytilde1 <- matrix(0,T,1)
    tau2 <- matrix(4,K_,1)
    if (K_>1){        
      V_prior <- diag(as.numeric(tau2))
    }else{
      V_prior <- matrix(as.numeric(tau2))
    }
    D <- matrix(1,2*K_,1)
    d <- d_prior  <- matrix(0,K_,1)
    D_t <- matrix(1,K_,T)
    Omega_t <- matrix(1,K_,T)
    svdraw <- list(para=c(mu=-10,phi=.9,sigma=.2),latent=rep(-3,T))
    hv <- svdraw$latent
    para <- list(mu=-10,phi=.9,sigma=.2)
    H <- matrix(-10,T,M)
    thin_out <- round(seq(nburn,ntot,length.out = thin*nsave))
  
    #storage matrices
    H_store <- matrix(NA,thin*nsave,T)
    ALPHA_store <- array(NA,c(thin*nsave,T,K_))
    V0_store <- array(NA,c(thin*nsave,K_))
    svparms_store <- matrix(NA,thin*nsave,3)
    D_store <- array(NA,c(thin*nsave,K_,T))
    thresholds_store <- array(NA,c(thin*nsave,K_))
    Omega_store <- array(NA,c(thin*nsave,K_,T))
    omega_store <- matrix(NA,thin*nsave,K_)
    sigma2_store <- matrix(NA,thin*nsave,1)
    thrshprior_store <- matrix(NA,thin*nsave,K_)
    kappa_store <- matrix(NA,thin*nsave, 1)
    #Vcov_store <- array(NA,c(thin*nsave,K_,K_,T))
    
    u_ <- matrix(NA,T,1)
    pb <- txtProgressBar(min = 0, max = ntot, style = 3)
    ntot <- nburn+nsave
    ithin <- 0
    accept <- 0
    MaxTrys <- 100
    for (irep in 1:ntot){
        Omega_t <- D_t*diag(sqrttheta1)+(1-D_t)*diag(sqrttheta2)
        #Draw for the time-varying part ALPHA
        ALPHA1 <- try(KF_fast(t(as.matrix(Y__)), X__,as.matrix(exp(hv)),t(Omega_t),K_, 1, T, matrix(0,K_,1), V_prior),silent=TRUE)
        if (is(ALPHA1,"try-error")){
          ALPHA1 <- KF(t(as.matrix(Y__)), X__,as.matrix(exp(hv)),t(Omega_t),K_, 1, T, matrix(0,K_,1), V_prior)
          try0 <- 0
          while (any(abs(ALPHA1$bdraw)>1e+10) && try0<MaxTrys){ #This block resamples if the draw from the state vector is not well behaved
            ALPHA1 <- try(KF(t(as.matrix(Y__)), X__,as.matrix(exp(hv)),t(Omega_t),K_, 1, T, matrix(0,K_,1), V_prior),silent=TRUE)
            try0 <- try0+1
          }
        }
        ALPHA <- ALPHA1$bdraw
        VCOV <- ALPHA1$Vcov
       
        # if (!is(ALPHA1,"try-error")) ALPHA <- ALPHA1

        #sample variances
        Em <- diff(t(ALPHA))
        for (jj in 1:K_){
          if (!cons.mod){
            #First regime
            sig_q <- sqrttheta1[jj,jj]
            if (!a.approx){
              si <- (abs(Em[,jj])>d[jj,1])*1
            }else{
              si <- (abs(Achg.OLS[,jj])>d[jj,1])*1
            }
            Em1 <- Em[si==1,jj]
            scale1 <- sum(Em1^2)
            s_1 <- B_gamma1+sum(si)/2+1/2
            s_2 <- B_gamma2+scale1/2
          
            sig_q  <- 1/rgamma(1,s_1,s_2)
            
            sqrttheta1[jj,jj] <- sig_q
            sqrttheta2[jj,jj] <- kappa0[jj]^2#*sig_q
          }else{
            sqrttheta1[jj,jj] <- sqrttheta2[jj,jj] <- 0
          }
        }
        #Sample hyperparameter lambda from G(a,b)
        lambda2_tau <- rgamma(1,c_tau+a_tau*K_,d_tau+a_tau/2*sum(tau2))
        
        #Sample variances of the time invariant part first
        for (nn in 1:K_){
          tau2_draw <- try(rgig(n=1,lambda=a_tau-0.5,ALPHA[nn,1]^2,a_tau*lambda2_tau),silent=TRUE)
          tau2[nn] <- ifelse(is(tau2_draw,"try-error"),next,tau2_draw)
        }
        
        if (K_>1){
          V_prior <- diag(as.numeric(tau2))
        }else{
          V_prior <- matrix(as.numeric(tau2))
        }
        
        if (TVS){
          #Check whether coefficient is time-varying or constant at each point in time
          Achg<- t(diff(t(ALPHA)))#t(as.numeric(ALPHA0)+ALPHA2)
          Achg <- cbind(matrix(0,K_,1),Achg) #we simply assume that the parameters stayed constant between t=0 and t=1
          if (a.approx) Achg.approx <- Achg.OLS else Achg.approx <- Achg
          grid.mat <- matrix(unlist(lapply(1:K_,function(x) get_grid(Achg[x,],sqrt(sqrttheta1[x,x]),grid.length=grid.length,thrsh.pct=thrsh.pct,thrsh.pct.high=thrsh.pct.high))), ncol = K_)
          probs <- get_threshold(Achg, sqrttheta1, sqrttheta2,grid.mat,Achg.approx)
          
          for (jj in 1:K_){
            post1 <- probs[,jj]
            probs1 <- exp(post1-max(post1))/sum(exp(post1-max(post1)))
            d[jj,] <- sample(grid.mat[,jj],1,prob=probs1)
            if (!a.approx){
              D_t[jj,] <- (abs(Achg[jj,])>d[jj,])*1 #change 2:T usw. here
            }else{
              D_t[jj,] <- (abs(Achg.OLS[jj,])>d[jj,])*1 #change 2:T usw. here
            }
          }
          
          if (sim.kappa){
            grid.kappa <- kappa.grid
            Lik.kappa <- matrix(0,length(grid.kappa),1)
            count <- 0
              for (grid.i in grid.kappa){
                count <- count+1
                sqrttheta.prop <- (grid.i*sd.OLS)^2
                cov.prop <- sqrt(D_t*diag(sqrttheta1)+(1-D_t)*sqrttheta.prop)
                
                Lik.kappa[count,1] <-sum(dnorm(t(Achg),matrix(0,T,2),cov.prop,log=TRUE)) 
              }
          Lik.kappa.norm <- exp(Lik.kappa-max(Lik.kappa))
          probs.kappa <- Lik.kappa.norm/sum(Lik.kappa.norm)
          scale.0 <- sample(grid.kappa,size=1, prob=probs.kappa)
          kappa0 <- scale.0*sd.OLS
          
          }
          
          
        }
        u_ <- matrix(0,T,1)
        for (jj in 1:T) u_[jj,] <- Y__[jj]-X__[jj,]%*%ALPHA[,jj]
        u_[abs(u_)<1e-7] <- 1e-7
        
        if (sv_on){
            #-------------------------if sv_on==TRUE -> use stochastic volatility specification-------------------------#
            svdraw <- svsample2(u_,startpara=para(svdraw),startlatent=latent(svdraw),priorphi=priorphi,priormu=priormu)#,priorlatent0=h0prior
            hv <- latent(svdraw)
            sig_eta <- 0
        }else{
            S_1 <- a_sig+T/2
            S_2 <- b_sig+crossprod(u_)/2
            sig_eta <- 1/rgamma(1,S_1,S_2)
            hv <- matrix(log(sig_eta),T,1)
        }
        if (irep %in% thin_out){
            ithin <- ithin+1
            #Vcov_store[ithin,,,] <- VCOV
            kappa_store[ithin,] <- scale.0
            H_store[ithin,] <- hv
            ALPHA_store[ithin,,] <- t(ALPHA)
            svparms_store[ithin,] <- para(svdraw)
            D_store[ithin,,] <- D_t
            thresholds_store[ithin,] <- d
            omega_store[ithin,] <-diag(sqrttheta1)
            Omega_store[ithin,,] <- Omega_t
            V0_store[ithin,] <- tau2
            sigma2_store[ithin,] <- sig_eta
            thrshprior_store[ithin,] <- d_prior
        }
        setTxtProgressBar(pb, irep)
      #  sfCat(paste("Iteration ",nr, irep), sep="\n")
    }
    
    dimnames(ALPHA_store)[[3]]<-c(replicate(length(slct),"COV"),colnames(X))
    
    return(list(A=ALPHA_store,V0=V0_store,D_dyn=D_store,H=H_store,svparms=svparms_store,threshold.prior=thrshprior_store,
    Omega=Omega_store,sigma2=sigma2_store,omega=omega_store,thresholds=thresholds_store,slct=slct,X=X,Y=Y,sd.OLS=sd.OLS,kappa=kappa_store))
}


VAR_posterior <- function(post_draws,nkeep=thin*save,p=p){
    T <- dim(post_draws[[1]]$A)[[2]]
    K <- dim(post_draws[[1]]$A)[[3]]
    save_o <- dim(post_draws[[1]]$A)[[1]]
    M <- length(post_draws)
    #Storage matrices for the full system
    H_store <- array(0,c(T,M,nkeep))
    ALPHA_store <- array(0,c(T,K,M,nkeep))
    A_store <- array(0,c(T,K,M,nkeep))
    A0_store <- array(0,c(T,M,M,nkeep))
    S_post <- array(0,c(T,M,M,nkeep))
    #Sample from the posterior
    eigs <- matrix(0,T,1)
    stab_ind <- matrix(0,nkeep,1)
    #initialize progress bar
  #  pb <- txtProgressBar(min = 0, max = nkeep, style = 3)
    draw_seq <- round(seq(1,save_o,length.out = nkeep))
    for (ii in 1:nkeep){
        iii <- draw_seq[[ii]]
        for (jj in 1:M){
            slct <- post_draws[[jj]]$slct
            #split and create structural matrix
            A0_store[,jj,slct,ii] <- (-1*post_draws[[jj]]$A[iii,,slct])
            if (jj==1){
                ALPHA_store[,,jj,ii] <- (post_draws[[jj]]$A[iii,,])
            }else{
                ALPHA_store[,,jj,ii] <- (post_draws[[jj]]$A[iii,,-slct])
            }
            H_store[,jj,ii] <- post_draws[[jj]]$H[iii,]
        }
        for (nn in 1:T){
            A0 <- A0_store[nn,,,ii]
            diag(A0) <- 1
            A0inv <- try(solve(A0),silent=TRUE) #just in case A0 becomes numerically unstable
            if (is(A0inv,"try-error")) A0inv <- ginv(A0)
            S_post[nn,,,ii] <- A0inv%*%diag(exp(H_store[nn,,ii]))%*%t(A0inv)
            Atilda <- t(A0inv%*%t(ALPHA_store[nn,,,ii]))
            A_store[nn,,,ii] <- Atilda
            #fmat <- get_companion((Atilda[2:nrow(Atilda),]),c(M,0,p))
            #eigs[nn] <- max(abs(Re(eigen(fmat$MM)$values)))
        }
        #stab_ind[ii,] <- ifelse(max(eigs)>0.99999,0,1)
      #  setTxtProgressBar(pb, ii)
      #  sfCat(paste("Iteration ",nr, irep), sep="\n")
      
    }
    return(list(A_post=A_store,S_post=S_post))
}
mlag <- function(X,lag){
  p <- lag
  X <- as.matrix(X)
  Traw <- nrow(X)
  N <- ncol(X)
  Xlag <- matrix(0,Traw,p*N)
  for (ii in 1:p){
    Xlag[(p+1):Traw,(N*(ii-1)+1):(N*ii)]=X[(p+1-ii):(Traw-ii),(1:N)]
  }
  return(Xlag) 
}

get_companion <- function(Beta_,varndxv){
  # Beta_ <- Atilda
  nn <- varndxv[[1]]
  nd <- varndxv[[2]]
  nl <- varndxv[[3]]
  
  nkk <- nn*nl
  
  Jm <- matrix(0,nkk,nn)
  Jm[1:nn,1:nn] <- diag(nn)
  
  MM <- rbind(t(Beta_),cbind(diag((nl-1)*nn), matrix(0,(nl-1)*nn,nn)))
  
  return(list(MM=MM,Jm=Jm))
}


int_posterior <- function(d=d_prop,Achg1=Achg[jj,],sqrttheta11=sqrttheta1[jj,jj],sqrttheta21=sqrttheta2[jj,jj]){
  D_i <- (abs(Achg1)>d)*1
  
  b_01 <- Achg1[D_i==1]
  b_02 <- Achg1[D_i==0]
  A_1 <- sum(dnorm(b_01,0,sqrt(sqrttheta11),log=TRUE))
  A_2 <- sum(dnorm(b_02,0,sqrt(sqrttheta21),log=TRUE))
  
  post_1 <- A_1+A_2
  #Stability check
  if (is.infinite(post_1) && post_1<0) post_1 <- -10^10
  
  return(post_1)
}

get_grid <- function(Achg,sd.state,grid.length=150,thrsh.pct=0.1,thrsh.pct.high=0.9){
  d_prop <- seq(thrsh.pct*sd.state,thrsh.pct.high*sd.state,length.out=grid.length)
  return(d_prop)
}

MCMC_mix <- function(Y,X,nburn,nsave,priorbtheta=list(B_1=2,B_2=1,kappa0=1e-07),priorb0=list(a_tau=0.1,c_tau=0.01,d_tau=0.01),priorsig=c(0.01,0.01),priorphi=c(2,2),priormu=c(0,10^2),h0prior="stationary",grid.length=150,thrsh.pct=0.1,thrsh.pct.high=1,sv_on=FALSE,TVS=FALSE,cons.mod=FALSE,nr,thin=0.1,robust=TRUE,a.approx=FALSE){
  #Y <- as.matrix(Y1$Y);X <- as.matrix(Y1$X);nburn <- 500;nsave <- 500;priorbtheta=list(B_1=2,B_2=1,kappa0=1e-10);priorb0=list(a_tau=0.1,c_tau=0.01,d_tau=0.01);priorsig=c(0.01,0.01);priorphi=c(2,2);priormu=c(0,10^2);h0prior="stationary";grid.length=150;thrsh.pct=0.1;thrsh.pct.high=1;sv_on=FALSE;TVS=TRUE;cons.mod=FALSE;nr=1;thin=0.1;robust=TRUE;a.approx=FALSE  
  #-----------------------------------------------------------------------------------------------------------#
  ntot <- nsave+nburn
  a_tau <- priorb0$a_tau
  c_tau <- priorb0$c_tau
  d_tau <- priorb0$d_tau
  
  
  B_gamma1 <- priorbtheta$B_1 #input new
  B_gamma2 <- priorbtheta$B_2 
  kappa0 <- priorbtheta$kappa0 
  
  a_sig <- priorsig[[1]]
  b_sig <- priorsig[[2]]
  #---------------------------------sourceCpp to source .cpp functions----------------------------------------#
  require(stochvol)
  require(GIGrvg)
  require(snowfall)
  #-----------------------------------------------------------------------------------------------------------#
  class(X) <-"matrix"
  class(Y) <- "matrix"
  T <- nrow(Y)
  K <- ncol(X)
  M <- ncol(Y)
  
  #Create full-data matrices
  if (nr==1) slct <- NULL else slct <- 1:(nr-1)
  Y__ <- Y[,nr]
  X__ <- cbind(Y[,slct],X)
  K_ <- ncol(X__)
  M_ <- M-length(slct)
  #Colnames checking
  if (is.null(colnames(X))) colnames(X) <- seq(1,ncol(X))
  #Initialize everything
  sqrttheta1 <- diag(K_)*0.1
  sqrttheta2 <-diag(K_)*0.01
  #mylm <- lm(Y__ ~ X__)
  
  kvals = matrix(1,2,1)
  kvals[1,1] = 0
  
  B.OLS <- solve(crossprod(X__))%*%crossprod(X__,Y__)
  sig.OLS <- crossprod(Y__-X__%*%B.OLS)/(T-K_)
  V.OLS <- as.numeric(sig.OLS)*solve(crossprod(X__))
  
  #Specification of scaling parameter in the lower regime
  sd.OLS <- sqrt(diag(V.OLS))
  if (kappa0<0) kappa0 <- -kappa0 * sd.OLS else kappa0 <- matrix(kappa0,K_,1)
  
  ALPHA0 <- solve(crossprod(X__))%*%crossprod(X__,Y__)
  
  
  buildCapm <- function(u){
    dlm::dlmModReg(X__, dV = exp(u[1]), dW = exp(u[2:(K_+1)]),addInt = FALSE)
  }
  
  outMLE <- dlm::dlmMLE(Y__, parm = rep(0,K_+1), buildCapm)
  mod <- buildCapm(outMLE$par)
  outS <- dlm::dlmSmooth(Y__, mod)
  states.OLS <- t(matrix(outS$s,T+1,K_))
  
  Achg.OLS <- t(diff(t(states.OLS)))#t(as.numeric(ALPHA0)+ALPHA2)
  #Achg.OLS  <- cbind(matrix(0,K_,1),Achg.OLS) #we simply assume that the parameters stayed constant between t=0 and t=1
  
  Ytilde1 <- matrix(0,T,1)
  tau2 <- matrix(4,K_,1)
  if (K_>1){        
    V_prior <- diag(as.numeric(tau2))
  }else{
    V_prior <- matrix(as.numeric(tau2))
  }
  D <- matrix(1,2*K_,1)
  d <- d_prior  <- matrix(0,K_,1)
  D_t <- matrix(1,K_,T)
  Omega_t <- matrix(1,K_,T)
  svdraw <- list(para=c(mu=-10,phi=.9,sigma=.2),latent=rep(-3,T))
  hv <- svdraw$latent
  para <- list(mu=-10,phi=.9,sigma=.2)
  H <- matrix(-10,T,M)
  thin_out <- round(seq(nburn,ntot,length.out = thin*nsave))
  kprior <- .5*matrix(1,2,1)
  #storage matrices
  H_store <- matrix(NA,thin*nsave,T)
  ALPHA_store <- array(NA,c(thin*nsave,T,K_))
  V0_store <- array(NA,c(thin*nsave,K_))
  svparms_store <- matrix(NA,thin*nsave,3)
  D_store <- array(NA,c(thin*nsave,K_,T))
  thresholds_store <- array(NA,c(thin*nsave,K_))
  Omega_store <- array(NA,c(thin*nsave,K_,T))
  omega_store <- matrix(NA,thin*nsave,K_)
  sigma2_store <- matrix(NA,thin*nsave,1)
  thrshprior_store <- matrix(NA,thin*nsave,K_)
  #Vcov_store <- array(NA,c(thin*nsave,K_,K_,T))
  
  u_ <- matrix(NA,T,1)
  pb <- txtProgressBar(min = 0, max = ntot, style = 3)
  ntot <- nburn+nsave
  ithin <- 0
  accept <- 0
  MaxTrys <- 100
  for (irep in 1:ntot){
    Omega_t <- D_t*diag(sqrttheta1)+(1-D_t)*diag(sqrttheta2)
    #Draw for the time-varying part ALPHA
    ALPHA1 <- try(KF_fast(t(as.matrix(Y__)), X__,as.matrix(exp(hv)),t(Omega_t),K_, 1, T, matrix(0,K_,1), V_prior),silent=TRUE)
    if (is(ALPHA1,"try-error")){
      ALPHA1 <- KF(t(as.matrix(Y__)), X__,as.matrix(exp(hv)),t(Omega_t),K_, 1, T, matrix(0,K_,1), V_prior)
      try0 <- 0
      while (any(abs(ALPHA1$bdraw)>1e+10) && try0<MaxTrys){ #This block resamples if the draw from the state vector is not well behaved
        ALPHA1 <- try(KF(t(as.matrix(Y__)), X__,as.matrix(exp(hv)),t(Omega_t),K_, 1, T, matrix(0,K_,1), V_prior),silent=TRUE)
        try0 <- try0+1
      }
    }
    ALPHA <- ALPHA1$bdraw
    VCOV <- ALPHA1$Vcov
    
    #sample variances
    Em <- diff(t(ALPHA))
    for (jj in 1:K_){
      if (!cons.mod){
        #First regime
        sig_q <- sqrttheta1[jj,jj]
        
        si <- D_t[jj,2:T]
        
        Em1 <- Em[si==1,jj]
        scale1 <- sum(Em1^2)
        s_1 <- B_gamma1+sum(si)/2+1/2
        s_2 <- B_gamma2+scale1/2
        
        sig_q  <- 1/rgamma(1,s_1,s_2)
        
        sqrttheta1[jj,jj] <- sig_q
        sqrttheta2[jj,jj] <- kappa0[jj]#*sig_q
      }else{
        sqrttheta1[jj,jj] <- sqrttheta2[jj,jj] <- 0
      }
    }
    #Sample hyperparameter lambda from G(a,b)
    lambda2_tau <- rgamma(1,c_tau+a_tau*K_,d_tau+a_tau/2*sum(tau2))
    
    #Sample variances of the time invariant part first
    for (nn in 1:K_){
      tau2_draw <- try(rgig(n=1,lambda=a_tau-0.5,ALPHA[nn,1]^2,a_tau*lambda2_tau),silent=TRUE)
      tau2[nn] <- ifelse(is(tau2_draw,"try-error"),next,tau2_draw)
    }
    
    if (K_>1){
      V_prior <- diag(as.numeric(tau2))
    }else{
      V_prior <- matrix(as.numeric(tau2))
    }
    
    ap = 1 + sum(D_t)
    bp = 1 + T - sum(D_t)
    
    pdrawa = rbeta(1,ap,bp)
    
    kprior[2,1] = pdrawa;
    kprior[1,1] = 1 - kprior[2,1]
    
    
    D_t <- t(gck(Y__,matrix(0,1,T),X__,as.matrix(exp(hv/2)),matrix(0,1,T),kronecker(matrix(1,T,1),diag(1)),sqrt(sqrttheta1),t(D_t),T,matrix(0,1,1),matrix(0,1,1),2,kprior,kvals,1,K))
    u_ <- matrix(0,T,1)
    for (jj in 1:T) u_[jj,] <- Y__[jj]-X__[jj,]%*%ALPHA[,jj]
    u_[abs(u_)<1e-7] <- 1e-7
    
    if (sv_on){
      #-------------------------if sv_on==TRUE -> use stochastic volatility specification-------------------------#
      svdraw <- svsample2(u_,startpara=para(svdraw),startlatent=latent(svdraw),priorphi=priorphi,priormu=priormu)#,priorlatent0=h0prior
      hv <- latent(svdraw)
      sig_eta <- 0
    }else{
      S_1 <- a_sig+T/2
      S_2 <- b_sig+crossprod(u_)/2
      sig_eta <- 1/rgamma(1,S_1,S_2)
      hv <- matrix(log(sig_eta),T,1)
    }
    if (irep %in% thin_out){
      ithin <- ithin+1
      #Vcov_store[ithin,,,] <- VCOV
      H_store[ithin,] <- hv
      ALPHA_store[ithin,,] <- t(ALPHA)
      svparms_store[ithin,] <- para(svdraw)
      D_store[ithin,,] <- D_t
      thresholds_store[ithin,] <- d
      omega_store[ithin,] <-diag(sqrttheta1)
      Omega_store[ithin,,] <- Omega_t
      V0_store[ithin,] <- tau2
      sigma2_store[ithin,] <- sig_eta
      thrshprior_store[ithin,] <- d_prior
    }
    setTxtProgressBar(pb, irep)
    #  sfCat(paste("Iteration ",nr, irep), sep="\n")
  }
  
  dimnames(ALPHA_store)[[3]]<-c(replicate(length(slct),"COV"),colnames(X))
  
  return(list(A=ALPHA_store,V0=V0_store,D_dyn=D_store,H=H_store,svparms=svparms_store,threshold.prior=thrshprior_store,
              Omega=Omega_store,sigma2=sigma2_store,omega=omega_store,thresholds=thresholds_store,slct=slct,X=X,Y=Y,sd.OLS=sd.OLS))
}

gck <- function(yg,gg,hh,capg,f,capf,sigv,kold,t,ex0,vx0,nvalk,kprior,kvals,p,kstate){
  # GCK's Step 1 on page 821
  lpy2n=0;
  mu = matrix(0,t*kstate,1);
  omega = matrix(0,t*kstate,kstate);
  for (i in seq(t-1,1,by=-1)){
    gatplus1 = sigv%*%kold[i+1]
    ftplus1 = capf[(kstate*i+1):(kstate*(i+1)),]
    cgtplus1 = capg[(i*p+1):((i+1)*p),]
    htplus1 = t(hh[(i*p+1):((i+1)*p),])
    
    htt1 <- crossprod(htplus1,gatplus1)
    rtplus1 = tcrossprod(htt1)+tcrossprod(cgtplus1,cgtplus1)
    rtinv = solve(rtplus1)
    btplus1 = tcrossprod(gatplus1)%*%htplus1%*%rtinv
    atplus1 = (diag(kstate)-tcrossprod(btplus1,htplus1))%*%ftplus1
    
    if (kold[i+1] == 0){
      ctplus1 = matrix(0,kstate,kstate)
    }else{
      cct = gatplus1%*%(diag(kstate)-crossprod(gatplus1,htplus1)%*%tcrossprod(rtinv,htplus1)%*%gatplus1)%*%t(gatplus1)
      ctplus1 = t(chol(cct))
    }
    otplus1 = omega[(kstate*i+1):(kstate*(i+1)),]
    
    dtplus1 = crossprod(ctplus1,otplus1)%*%ctplus1+diag(kstate)
    omega[(kstate*(i-1)+1):(kstate*i),] = crossprod(atplus1,(otplus1 - otplus1%*%ctplus1%*%solve(dtplus1)%*%t(ctplus1)%*%otplus1))%*%atplus1+t(ftplus1)%*%htplus1%*%rtinv%*%t(htplus1)%*%ftplus1
    satplus1 = (diag(kstate)-tcrossprod(btplus1,htplus1))%*%(f[,i+1]-btplus1%*%gg[,i+1]) #CHCKCHCKCHKC
    mutplus1 = mu[(kstate*i+1):(kstate*(i+1)),]
    mu[(kstate*(i-1)+1):(kstate*i),] = crossprod(atplus1,(diag(kstate)-otplus1%*%ctplus1%*%solve(dtplus1)%*%t(ctplus1)))%*%(mutplus1-otplus1%*%(satplus1+btplus1%*%yg[i+1]))+t(ftplus1)%*%htplus1%*%rtinv%*%(yg[i+1]-gg[,i+1]-t(htplus1)%*%f[,i+1])
  }
  
  
  # GCKs Step 2 on pages 821-822
  kdraw = kold;
  ht = t(hh[1:p,])
  ft = capf[1:kstate,]
  gat = matrix(0,kstate,kstate)
  # Note: this specification implies no shift in first period -- sensible
  rt = t(ht)%*%ft%*%vx0%*%t(ft)%*%ht + crossprod(ht,gat)%*%crossprod(gat,ht)+ tcrossprod(capg[1:p,])
  rtinv = solve(rt)
  jt = (ft%*%vx0%*%t(ft)%*%ht + tcrossprod(gat)%*%ht)%*%rtinv
  mtm1 = (diag(kstate) - tcrossprod(jt,ht))%*%(f[,1] + ft%*%ex0) + jt%*%(yg[1] - gg[,1])
  vtm1 <- ft%*%tcrossprod(vx0,ft)+tcrossprod(gat)-jt%*%tcrossprod(rt,jt)        
  lprob <- matrix(0,nvalk,1)
  
  for (i in 2:t){
    ht <-  t(hh[((i-1)*p+1):(i*p),])
    ft <- capf[(kstate*(i-1)+1):(kstate*i),]  
    for (j in 1:nvalk){
      gat <- kvals[j,1]%*%sigv
      rt <- crossprod(ht,ft)%*%tcrossprod(vtm1,ft)%*%ht+crossprod(ht,gat)%*%crossprod(gat,ht)+tcrossprod(capg[((i-1)*p+1):(i*p),])
      rtinv <- solve(rt)
      jt <- (ft%*%tcrossprod(vtm1,ft)%*%ht+tcrossprod(gat)%*%ht)%*%rtinv
      mt <- (diag(kstate)-tcrossprod(jt,ht))%*%(f[,i]+ft%*%mtm1)+jt%*%(yg[i]-gg[,i])
      vt <- ft%*%tcrossprod(vtm1,ft)+tcrossprod(gat)-jt%*%tcrossprod(rt,jt)
      
      lpyt = -.5*log(det(rt)) - .5*t(yg[i] - gg[,i] - t(ht)%*%t(f[,i] + ft%*%mtm1))%*%rtinv%*%(yg[i] - gg[,i] - t(ht)%*%(f[,i] + ft%*%mtm1))
      
      if (det(vt)<=0){
        tt <- matrix(0,kstate,kstate)
      }else{
        tt <- t(chol(vt))
      }
      ot = omega[(kstate*(i-1)+1):(kstate*i),]
      mut = mu[(kstate*(i-1)+1):(kstate*i),]
      tempv = diag(kstate) + crossprod(tt,ot)%*%tt
      lpyt1n = -.5*log(det(tempv)) -.5*(crossprod(mt,ot)%*%mt-2*crossprod(mut,mt)-t(mut-ot%*%mt)%*%tt%*%solve(tempv)%*%t(tt)%*%(mut-ot%*%mt))
      lprob[j,1] <- log(kprior[j,1])+lpyt1n+lpyt                       
      if (i==2){
        lpy2n <- lpyt1n+lpyt
      }
    } 
    pprob = exp(lprob-max(lprob))/sum(exp(lprob-max(lprob)))
    tempv = runif(1)
    tempu = 0
    for (j in 1:nvalk){
      tempu <- tempu+pprob[j,1]
      if (tempu> tempv){
        kdraw[i] <- kvals[j,1]
        break
      }
    }
    gat = kdraw[i]%*%sigv
    rt = crossprod(ht,ft)%*%tcrossprod(vtm1,ft)%*%ht+t(ht)%*%tcrossprod(gat)%*%ht+tcrossprod(capg[((i-1)*p+1):(i*p)])
    rtinv = solve(rt)
    jt = (ft%*%tcrossprod(vtm1,ft)%*%ht+tcrossprod(gat)%*%ht)%*%rtinv
    mtm1 <- (diag(kstate)-tcrossprod(jt,ht))%*%(f[,i]+ft%*%mtm1)+jt%*%(yg[i]-gg[,i])  
    vtm1 = ft%*%tcrossprod(vtm1,ft)+tcrossprod(gat)-jt%*%tcrossprod(rt,jt)
  }
  return(kdraw)
}
