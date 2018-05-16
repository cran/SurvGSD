#' Maximum sample size for a group sequential test
#' under a generalized gamma survival distribution or a log-logistic survival distribution.
#' @description A function obtains maximum sample sizes and associated expected values
#' for a group sequential design under a generalized gamma survival distribution or a log-logistic survival distribution
#' for a given dropout censoring distribution.
#' @import stats flexsurv
#' @importFrom ldbounds bounds drift
#' @importFrom mnormt sadmvn

#' @param t the interim analysis time (vector).
#' @param R the recuritment duration.
#' @param T the study duration.
#' @param FUN.C the cumulative distribution function of dropout censoring. \cr
#' FUN.C = function(y) punif(y,0,h) for a uniform dropout censoring U(0,h); \cr
#' FUN.C = function(y) rep(0,length(y)) for assuming no dropout censoring.
#' @param para0 c(q0,mu0,sigma0), parameters of an assumed generalized gamma distribution for the control arm.
#' A character string q0="LLG" indiactes an assumed log-logistic survival distribution
#' \eqn{F_0(y;\xi,\zeta)=1/(1+(y/\xi)^{-\zeta})} for the control arm, where \eqn{\xi} = mu0 and \eqn{\zeta} = sigma0.
#'
#' @param para1 c(q1,mu1,sigma1), parameters of an assumed generalized gamma distribution for the treatment arm.
#' A character string q1="LLG" indiactes an assumed log-logistic survival distribution
#' \eqn{F_1(y;\xi,\zeta)=1/(1+(y/\xi)^{-\zeta})} for the treatment arm, where \eqn{\xi} = mu1 and \eqn{\zeta} = sigma1.
#'
#' @param haz.r the hazard ratio of the treatment arm to the control arm (numeric or function).
#' @param rho the power in the weight of the Harrington-Fleming statistic.
#' \eqn{\rho=0} for the logrank test; \eqn{\rho=1} for the Wilcoxon test.
#' @param eta,theta parameters of the entry distribution with \eqn{\eta \ge -\theta/R} and \eqn{\eta >0}
#' (\eqn{\theta=0} for the Uniform dropout censoring).
#' @param px the proportion of patients assigned to the treatment arm.
#' The default is px = 0.5 indicating 1:1 allocation.
#' @param spf 1 = O’Brien-Fleming-type; 2 = Pocock-type alpha-spending function. The default is spf = 1.
#' @param alpha the type I error. The default is alpha = 0.05.
#' @param power A desired value of the power. The default is power = 0.8.
#' @return \item{MaxSize}{the maximum sample size.}
#' @return \item{ExpSize}{the expected sample size.}
#' @return \item{ExpEvent}{the expected number of events.}
#' @return \item{A.power}{actual achieved power.}
#' @return \item{Info.fractions}{information fractions at times of all the interim analyses.}
#' @return \item{boundary}{the monitoring boundary values of
#' the standardized Harrington-Fleming statistic at all the interim analyses.}
#' @examples # Assume an exponential (log-logistic) survival distribution
#' # with q0=sigma0=1, mu0=0.367 (xi0=1, zeta0=1.75) for the control arm,
#' # a uniform patient entry (eta=1,theta=0) and a uniform dropout censoring distribution Unif(0,h)
#' # having a 15% censoring probability (lfu=0.15) for a study with R=2, T=3 and the interim
#' # analysis time at t=1,1.5,2,2.5.
#'
#' # To obtain the required h for the uniform dropout censoring distribution.
#' @examples Find.h(lfu=0.15, R=2, T=3, q=1, mu=0.367, sigma=1, eta=1, theta=0) ## exponential
#' @examples Find.h(lfu=0.15, R=2, T=3, q="LLG", mu=1, sigma=1.75, eta=1, theta=0) ## log-logistic
#'
#' # To obtain the maximum sample size for testing a treatment difference of a hazard ratio of 2/3
#' # with a type-I error of 0.05 and a power of 0.8.
#' @examples SSize.FixAlter(t=c(1,1.5,2,2.5), R=2, T=3, FUN.C=function(y) punif(y,0,7.018),
#' para0=c(1,0.367,1), para1=NULL, haz.r=2/3, rho=0, eta=1, theta=0) # exponential
#' @examples SSize.FixAlter(t=c(1,1.5,2,2.5), R=2, T=3, FUN.C=function(y) punif(y,0,7.211),
#' para0=c("LLG",1,1.75), para1=NULL, haz.r=2/3, rho=0, eta=1, theta=0) # log-logistic
#'
#' # To obtain the maximum sample size for testing H_0:F_0=F_1 with a type-I error of 0.05
#' # and a power of 0.8, where F_1 is an exponential (log-logistic) distribution
#' # with the parameter para1=c(1,0.772,1) (para1=c("LLG",1.5,1.75)).
#' @examples SSize.FixAlter(t=c(1,1.5,2,2.5), R=2, T=3, FUN.C=function(y) punif(y,0,7.018),
#' para0=c(1,0.367,1), para1=c(1,0.772,1), haz.r=NULL, rho=0, eta=1, theta=0) # exponential
#' @examples SSize.FixAlter(t=c(1,1.5,2,2.5), R=2, T=3, FUN.C=function(y) punif(y,0,7.211),
#' para0=c("LLG",1,1.75), para1=c("LLG",1.5,1.75), haz.r=NULL, rho=0, eta=1, theta=0) # log-logistic
#'
#' @references Hsu, C.-H., Chen, C.-H, Hsu, K.-N. and Lu, Y.-H. (2018). A useful design utilizing the information fraction
#'                 in a group sequential clinical trial with censored survival data. To appear in Biometrics.
#' @references Azzalini, A. and Genz, A. (2015). The R package `mnormt':
#' The multivariate normal and 't' distributions (version 1.5-3). URL \emph{http://azzalini.stat.unipd.it/SW/Pkg-mnormt}.
#' @references Casper, C. and Perez, O. A. (2014). The R package `ldbounds':
#' Lan-DeMets method for group sequential boundaries (version 1.1-1). URL \emph{https://cran.r-project.org/web/packages/ldbounds/index.html}.
#' @references Jackson, C., Metcalfe, P. and Amdahl, J. (2017). The R package `flexsurv':
#' Flexible Parametric Survival and Multi-State Models (version 1.1). URL \emph{https://github.com/chjackson/flexsurv-dev}.
#' @export
SSize.FixAlter <- function(t, R, T, FUN.C, para0, para1 = NULL,
                              haz.r, rho = 0, eta = 1, theta = 0, px = 0.5, spf = 1, alpha = 0.05, power = 0.8) {

    if (para0[1]=="LLG") {
      q0=para0[1]; xi0=mu0=as.numeric(para0[2]); zeta0=sigma0=as.numeric(para0[3]);
    } else {
      q0=para0[1]; mu0 = para0[2]; sigma0 = para0[3];
    }

    if (is.null(para1) && is.null(haz.r)) {
       print("warning: please enter q1, mu1, sigma1 or a hazard ratio for the PH model")
    } else if (is.numeric(haz.r) && haz.r < 0) {
       print("warning: The hazard ratio needs to be larger than zero")
    } else {
       #### n^{-1/2}S_n(t): Parameters calculation of its Normal approximation ####
       if ( is.null(haz.r) == FALSE ) {
          print("Assume PH model with the specified hazard ratio")
          if (para0[1]=="LLG") {
              para0LLG <- c(xi0, zeta0)
              para <- FUN.Para.FixAlter.LLG(t, R, T, FUN.C, para0LLG, para1=NULL, haz.r, rho, eta, theta, px)
              prob.event <- FUN.Event.Prob.LLG(t, R, T, FUN.C, para0LLG, para1=NULL, haz.r, rho, eta, theta, px)$prob.vec
          } else {
              para <- FUN.Para.FixAlter(t, R, T, FUN.C, para0, para1=NULL, haz.r, rho, eta, theta, px)
              prob.event <- FUN.Event.Prob(t, R, T, FUN.C, para0, para1=NULL, haz.r, rho, eta, theta, px)$prob.vec
          }
       } else {
          print("Assume Non-PH model with the specified parameters for the treatment arm")
          if (para0[1]=="LLG") {
              xi1 = as.numeric(para1[2]); zeta1 = as.numeric(para1[3]);
              para0LLG <- c(xi0, zeta0); para1LLG <- c(xi1, zeta1)
              para <- FUN.Para.FixAlter.LLG(t, R, T, FUN.C, para0LLG, para1LLG, haz.r=NULL, rho, eta, theta, px)
              prob.event <- FUN.Event.Prob.LLG(t, R, T, FUN.C, para0LLG, para1LLG, haz.r=NULL, rho, eta, theta, px)$prob.vec
          } else {
              para <- FUN.Para.FixAlter(t, R, T, FUN.C, para0, para1, haz.r=NULL, rho, eta, theta, px)
              prob.event <- FUN.Event.Prob(t, R, T, FUN.C, para0, para1, haz.r=NULL, rho, eta, theta, px)$prob.vec
          }
       }

       #### Information fractions and boundaries ####
       if (is.null(t)==TRUE) {
          InfT <- 1; upbd <- qnorm(1-alpha[1]); lowbd <- - qnorm(1-alpha[2])
       } else {
          InfT <- General.tau(t, R, T, FUN.C, q0, mu0, sigma0, rho, eta, theta)$Info.fractions  ## information fractions
          InfT <- c(InfT, 1)
          LD <- bounds(InfT, iuse=c(spf,spf), alpha=rep(alpha/2,2))
          upbd  <- LD$upper.bounds  ## bk's
          lowbd <- LD$lower.bounds
       }

       #### establish the covariance matrix of sigma(ti,tj) ####
       cov <- para$sigma.vec
       covm <- matrix(0,length(InfT),length(InfT))
       for (i in length(InfT):1 ) {
           covm[i,] <- cov[i]; covm[,i] <- cov[i]
       }
       varcov <- covm
       upperbd <- sqrt(cov) * upbd
       lowerbd <- sqrt(cov) * lowbd

       #### The required maximum size ####
       aaa <- para$mu.vec/sqrt(para$sigma.vec) / sqrt(InfT)
       SimpleN <- round(drift(lowbd, upbd, t=InfT, pow=power)$drift^2 / aaa[length(InfT)]^2)
       add.n <- round(px/(1-px)) + 1; rem <- SimpleN %% add.n
       SimpleN <- SimpleN - rem + add.n * ceiling(rem/add.n)

       typeIIerror <- 1 - power ## beta
       ini.n <- round(SimpleN/2)
       n <- ini.n - (ini.n %% add.n)
       mean <- sqrt(n) * para$mu.vec
       beta <- sadmvn(lowerbd, upperbd, mean, varcov)
       while(beta>typeIIerror) {
         n <- n + add.n
         mean <- sqrt(n)* para$mu.vec
         beta <- sadmvn(lowerbd, upperbd, mean, varcov)
       }

       #### The required values: d_k, n_k ####
       event.vec <- n * prob.event   ## d_k
       cru.t <- c(t,T); cru.t[cru.t>R] <- R
       size.vec <- n * ( eta*(cru.t)+theta*(cru.t^2/2) )/(eta*R+theta*R^2/2)  ## n_k

       #### P[|T(t_j)|<= b_j, j=1,...,k], k=1,...,K-1 ####
       mean <- sqrt(n)* para$mu.vec
       cum.prob <- rep(0,(length(InfT)-1))
       for (i in 1:(length(InfT)-1)) {
           cum.prob[i] <- sadmvn(lowerbd[1:i], upperbd[1:i], mean[1:i], varcov[1:i,1:i])
       }

       #### The expected values ####
       OriT <- c(t,T)  ## calendar time
       expevent <- event.vec[1]; expsize <- size.vec[1]; exptime <- OriT[1]
       for (i in 2:length(InfT)) {
           expsize <- expsize + (size.vec[i] - size.vec[i-1])*cum.prob[i-1]
           expevent <- expevent + (event.vec[i] - event.vec[i-1])*cum.prob[i-1]
           exptime  <- exptime + (OriT[i] - OriT[i-1])*cum.prob[i-1]
       }

       list(MaxSize=n, ExpSize= round(expsize,1), ExpEvent= round(expevent,1),
            ExpTime = round(exptime,3), A.power=1-beta[1], Info.fractions=InfT, boundary=rbind(upbd,lowbd))

    }

}
