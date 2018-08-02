\name{estimate_tvp}
\alias{estimate_tvp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimation of time-varying parameter models with latent thresholds
}
\description{
Estimates a time-varying parameter model with stochastic volatility and a set of latent thresholds that control the time-variation and account for model uncertainty for each parameter.
}
\usage{
estimate_tvp(Y, X, save = 5000, burn = 5000, Q_prmean = 0.01, Qprvar = 30, theta_0 = 1, theta_1 = 1, prior_phi = c(25, 5), prior_mu = c(0, 10), tau0 = 0, sv_on = 0, TVS = 1, alpha1 = 0.01, alpha2 = 0.01, beta0 = 1.5, beta1 = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
Matrix of endogenous variables of dimension T x M. If M>1 then the algorithm estimates a TVP-VAR-SV with latent thresholds.
}
  \item{X}{ Matrix of explanatory variables of dimension T x K. Has to be supplied by the user.
}
  \item{save}{
  Number of saved draws from the posterior. Default is 5000.
}
  \item{burn}{
 Number of burn-ins. Default is 5000.
}
  \item{theta_0}{
A value between 0 and 1. A value of 1 implies that the full range of absolute parameter values is used to construct a prior and the proposal for the underlying MH algorithm for the variable inclusion threshold. Default is 1.
}
  \item{theta_1}{
A value between 0 and 1. A value of 1 implies that the full range of absolute parameter changes is used as a proposal for the MH algorithm for the threshold controlling the time variation. Default is 1.
}
  \item{prior_phi}{
  A numeric vector of length two. Controls the prior moments of the Beta distribution used for the autoregressive parameter of a stochastic volatility model. Defaults to c(25,5)
}
  \item{prior_mu}{
 A numeric vector of length two. Mean and variance for the normally distributed prior on the level of the log-volatily. Defaults to c(0,10)
 }
  \item{sv_on}{
sv_on=1 turns on stochastic volatility as in Kastner & Frühwirth-Schnatter (2013). If sv_on=0 then a homoscedastic model is fitted with a standard inverted gamma prior on the innovation variance of the observation equation.
}
  \item{TVS}{
TVS=1 allows for shrinkage on the time-variation of the parameters. If TVS=0 the model of Nakajima/West (2013, JBES) arises.
}
  \item{alpha1}{
Scalar hyperparameter used for the Inverse Gamma prior on the innovation variance in the observation equation. Default is 0.001.
}
  \item{alpha2}{
Scalar hyperparameter used for the Inverse Gamma prior on the innovation variance in the observation equation. Default is 0.001.
}
  \item{beta0}{
Scalar hyperparameter for the Gamma prior on the threshold controlling the time-variation. Default is 1.5.
}
  \item{beta1}{
Scalar hyperparameter for the Gamma prior on the threshold controlling the time-variation. Default is 1.
}
  \item{kappa_l}{
Scalar hyperparameter. Specifies the prior shape and scale. Defaults to a_kappa=1.
}

  \item{zeta_l}{
Scalar hyperparameter. Specifies the prior shape and scale. Defaults to a_zeta=1.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Florian Huber
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (Y, X, save = 5000, burn = 5000, Q_prmean = 0.01, Qprvar = 30, 
    theta_0 = 1, theta_1 = 1, prior_phi = c(25, 5), prior_mu = c(0, 
        10), tau0 = 0, sv_on = 0, TVS = 1, alpha1 = 0.01, alpha2 = 0.01, 
    beta0 = 1.5, beta1 = 1) 
{
    parm_list <- list("save", "burn", "Q_prmean", "Qprvar", "theta_0", 
        "theta_1", "prior_phi", "prior_mu", "sv_on", "TVS", "tau0", 
        "alpha1", "alpha2", "beta0", "beta1", "nthin")
    if (length(ncol(Y)) == 0) 
        Y <- as.matrix(Y)
    K <- ncol(X)
    M <- ncol(Y)
    T <- nrow(Y)
    if (M > 1) {
        post_draws <- list()
        sfInit(parallel = TRUE, cpus = 4)
        sfExport(list = c(list("Y", "X", "KF", "MCMC_tvpm", "get_variable_incl", 
            "get_moving_incl"), parm_list))
        post_draws <- sfLapply(1:ncol(Y), function(i) MCMC_tvpm(Y, 
            X, nburn = burn, nsave = save, B0_inv = B0inv, Qprmean = Q_prmean, 
            Q_prvar = Qprvar, theta_incl = theta_0, theta_dyn = theta_1, 
            priorphi = prior_phi, priormu = prior_mu, sv_on = sv_on, 
            TVS = TVS, nr = i, tau_0 = tau0, a1 = alpha1, a2 = alpha2, 
            b0 = beta0, b1 = beta1))
        sfStop()
        get_var_post <- VAR_posterior(post_draws, nthin * save)
    }
    else {
        post_draws <- MCMC_tvpm(Y, X, nburn = burn, nsave = save, 
            B0_inv = B0inv, Qprmean = Q_prmean, Q_prvar = Qprvar, 
            theta_incl = theta_0, theta_dyn = theta_1, priorphi = prior_phi, 
            priormu = prior_mu, sv_on = sv_on, TVS = TVS, nr = 1, 
            tau_0 = tau0, a1 = alpha1, a2 = alpha2, b0 = beta0, 
            b1 = beta1)
    }
    return(list(posterior = post_draws, args = .construct.arglist(MCMC_tvpm)))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
