#' To find parameter of a uniform dropout censoring distribution under
#' a generalized gamma survival distribution or a log-logistic survival distribution.
#' @description A function finds parameter h of a uniform dropout censoring distribution U(0,h)
#' with a generalized gamma survival distribution or a log-logistic survival distribution
#' for a given dropout censoring probability.
#' @importFrom flexsurv pgengamma

#' @param lfu the dropout censoring probability.
#' @param R the recuritment duration.
#' @param T the study duration.
#' @param q,mu,sigma shape, location and scale parameters of an assumed generalized gamma distribution for the control arm.
#' A character string q="LLG" indiactes an assumed log-logistic survival distribution
#' \eqn{F_0(y;\xi,\zeta)=1/(1+(y/\xi)^{-\zeta})} for the control arm, where \eqn{\xi} = mu and \eqn{\zeta} = sigma.
#'
#' @param eta,theta parameters of the entry distribution with \eqn{\eta \ge -\theta/R} and \eqn{\eta >0}
#' (\eqn{\theta=0} for the uniform dropout censoring).
#' @return the parameter h of the uniform dropout censoring distribution U(0,h).
#' @examples Find.h(lfu=0.15, R=2, T=3, q=1, mu=0.367, sigma=1, eta=1, theta=0)
#' @examples Find.h(lfu=0.15, R=2, T=3, q="LLG", mu=1, sigma=1.75, eta=1, theta=0)
#' @export
Find.h <- function(lfu, R, T, q, mu, sigma, eta, theta) {

  if (q=="LLG") {
      xi <- mu; zeta <- sigma;
      To.find.h.LLG <- function(h, lfu, R, T, xi, zeta, eta, theta) {
         FUN.LFU.LLG(h=h, R, T, xi, zeta, eta, theta) - lfu
      }
      round(uniroot(To.find.h.LLG, interval=c(0.1,5*T), extendInt="yes", lfu=lfu, R=R, T=T, xi=xi, zeta=zeta, eta=eta, theta=theta)$root, 3)
  } else {
      To.find.h <- function(h, lfu, R, T, q, mu, sigma, eta, theta) {
         FUN.LFU.GGD(h=h, R, T, q, mu, sigma, eta, theta) - lfu
      }
      round(uniroot(To.find.h, interval=c(0.1,5*T), extendInt="yes", lfu=lfu, R=R, T=T, q=q, mu=mu, sigma=sigma, eta=eta, theta=theta)$root, 3)
  }


}

