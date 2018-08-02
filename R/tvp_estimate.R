estimate_tvp <- function(Y,X,save=15000,burn=15000,priorbtheta=list(B_1=2,B_2=1,kappa0=1e-07),priorb0=list(a_tau=0.1,c_tau=0.01,d_tau=0.01),priorsig=c(0.01,0.01),priorphi=c(2,2),priormu=c(0,10^2),h0prior="stationary",grid.length=150,thrsh.pct=0.1,thrsh.pct.high=1,sv_on=TRUE,TVS=TRUE,cons.mod=FALSE,p=NULL,thin=0.1,CPU=4,approx=FALSE, sim.kappa0=FALSE,kappa0.grid = seq(1e-4, 0.1, 10)){
    require(snowfall)
    #Get dimensions of the data / model
    if (length(ncol(Y))==0) Y <- as.matrix(Y)
    if (is.null(p) && ncol(Y)>1) stop("Specify the number of lags used to create X matrix in the VAR case")
    # priorb0=list(a_tau=0.1,c_tau=0.01,d_tau=0.01);priorsig=c(0.01,0.01);priorphi=c(25,1.5);priormu=c(0,10^2);sv_on=TRUE;TVS=TRUE;nthin=1
    
    K <- ncol(X)
    M <- ncol(Y)
    T <- nrow(Y)
    #--------------------------------------Starts the MCMC algorithm--------------------------------
    if (M>1){
        post_draws <- list()

        if (CPU==1){
          for (i in 1:ncol(Y)) post_draws[[i]] <- MCMC_tvp(Y,X,burn,save,priorbtheta=priorbtheta,priorb0=priorb0,priorsig=priorsig,
                                                             priorphi=priorphi,priormu=priormu,h0prior=h0prior,grid.length=grid.length,thrsh.pct=thrsh.pct,thrsh.pct.high=thrsh.pct.high,sv_on=sv_on,TVS=TVS,cons.mod=cons.mod,nr=i,thin=thin,a.approx=approx, sim.kappa=sim.kappa0,kappa.grid=kappa0.grid)
        }else{
          sfInit(parallel=TRUE,cpus=CPU)#,slaveOutfile="output.txt"
          sfExportAll()
          post_draws <- sfLapply(1:ncol(Y),function(i) MCMC_tvp(Y,X,burn,save,priorbtheta=priorbtheta,priorb0=priorb0,priorsig=priorsig,
                                                                priorphi=priorphi,priormu=priormu,h0prior=h0prior,grid.length=grid.length,thrsh.pct=thrsh.pct,thrsh.pct.high=thrsh.pct.high,sv_on=sv_on,TVS=TVS,cons.mod=cons.mod,nr=i,thin=thin, a.approx=approx, sim.kappa=sim.kappa0,kappa.grid=kappa0.grid))
          sfStop()
        }
        
#       Since we estimate the VAR system equation-by-equation we have to reconstruct the system of equations
        get_var_post <- VAR_posterior(post_draws,thin*save,p=p)
    }else{
        post_draws <- MCMC_tvp(Y,X,burn,save,priorbtheta=priorbtheta,priorb0=priorb0,priorsig=priorsig,
        priorphi=priorphi,priormu=priormu,h0prior=h0prior,grid.length=grid.length,thrsh.pct=thrsh.pct,thrsh.pct.high=thrsh.pct.high,sv_on=sv_on,TVS=TVS,cons.mod=cons.mod,nr=1,thin=thin, a.approx=approx, sim.kappa=sim.kappa0,kappa.grid=kappa0.grid)
        get_var_post <- NULL
    }
    return(list(posterior=post_draws,VAR_coeff=get_var_post,args=.construct.arglist(estimate_tvp)))
}

