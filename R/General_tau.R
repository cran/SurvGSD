#' Information fractions under a generalized gamma survival distribution or a log-logistic survival distribution.
#' @description A function calculates information fractions
#' with a generalized gamma survival distribution or a log-logistic survival distribution
#' for a given dropout censoring probability.
#' @importFrom flexsurv pgengamma dgengamma

#' @param t the interim analysis time (vector).
#' @param R the recuritment duration.
#' @param T the study duration.
#' @param FUN.C the cumulative distribution function of dropout censoring. \cr
#' FUN.C = function(y) punif(y,0,h) for a uniform dropout censoring U(0,h); \cr
#' FUN.C = function(y) rep(0,length(y)) for assuming no dropout censoring.
#' @param q,mu,sigma shape, location and scale parameters of an assumed generalized gamma distribution for the control arm.
#' A character string q="LLG" indiactes an assumed log-logistic survival distribution
#' \eqn{F_0(y;\xi,\zeta)=1/(1+(y/\xi)^{-\zeta})} for the control arm, where \eqn{\xi} = mu and \eqn{\zeta} = sigma.
#'
#' @param rho the power in the weight of the Harrington-Fleming statistic.
#' \eqn{\rho=0} for the logrank test; \eqn{\rho=1} for the Wilcoxon test.
#' @param eta,theta parameters of the entry distribution with \eqn{\eta \ge -\theta/R} and \eqn{\eta >0}
#' (\eqn{\theta=0} for the uniform dropout censoring).
#' @return \item{Info.fractions}{information fractions at times of all the interim analyses.}
#' @return \item{Event.prob}{the probability of events accumulated up to T.}
#' @return \item{Total.censor.prob}{the probability of censoring including the dropout and administrative censoring.}
#' @examples General.tau(t=c(1,1.5,2,2.5), R=2, T=3, FUN.C=function(y) punif(y,0,7.018),
#'          q=1, mu=0.367, sigma=1, rho=0, eta=1, theta=0)
#' @examples General.tau(t=c(1,1.5,2,2.5), R=2, T=3, FUN.C=function(y) punif(y,0,7.211),
#'          q="LLG", mu=1, sigma=1.75, rho=0, eta=1, theta=0)
#' @export
General.tau <- function(t, R, T, FUN.C, q, mu, sigma, rho, eta, theta){
    if (q=="LLG") {

      xi <- mu; zeta <- sigma;
      func <- function(t, rho){
        rho.2 <- 2 * rho
        Fun.0 <- function(y) {
          ( 1 - FUN.C(y) ) *
            ( 1/( 1+(y/xi)^zeta ) )^rho.2 *
            (zeta/xi)*(y/xi)^{zeta-1}/( 1+(y/xi)^zeta )^2
        }
        Fun.1 <- function(y) {
          (t-y) * ( 1 - FUN.C(y) ) *
            ( 1/( 1+(y/xi)^zeta ) )^rho.2 *
            (zeta/xi)*(y/xi)^{zeta-1}/( 1+(y/xi)^zeta )^2
        }
        Fun.2 <- function(y) {
          (t-y)^2 * ( 1 - FUN.C(y) ) *
            ( 1/( 1+(y/xi)^zeta ) )^rho.2 *
            (zeta/xi)*(y/xi)^{zeta-1}/( 1+(y/xi)^zeta )^2
        }
        fv1 <- 0; fv2 <- 0
        if ( t > R ) {
          fv1 <- max(     R * integrate(Fun.0, 0, t-R)$value +         integrate(Fun.1, t-R, t)$value, 0)
          fv2 <- max( R^2/2 * integrate(Fun.0, 0, t-R)$value + (1/2) * integrate(Fun.2, t-R, t)$value, 0)
        } else {
          fv1 <- max( integrate(Fun.1, 0, t)$value, 0)
          fv2 <- max( (1/2) * integrate(Fun.2, 0, t)$value, 0)
        }
        (eta/(eta*R + theta*R^2/2)) * fv1 + (theta/(eta*R + theta*R^2/2)) * fv2
      }

    } else {

      func <- function(t, rho){
        rho.2 <- 2 * rho
        Fun.0 <- function(y) {
                  ( 1 - FUN.C(y) ) *
                  pgengamma(y, mu = mu, sigma = sigma, Q = q, lower.tail = FALSE)^rho.2 *
                  dgengamma(y, mu = mu, sigma = sigma, Q = q)
        }
        Fun.1 <- function(y) {
                  (t-y) * ( 1 - FUN.C(y) ) *
                  pgengamma(y, mu = mu, sigma = sigma, Q = q, lower.tail = FALSE)^rho.2 *
                  dgengamma(y, mu = mu, sigma = sigma, Q = q)
        }
        Fun.2 <- function(y) {
                  (t-y)^2 * ( 1 - FUN.C(y) ) *
                  pgengamma(y, mu = mu, sigma = sigma, Q = q, lower.tail = FALSE)^rho.2 *
                  dgengamma(y, mu = mu, sigma = sigma, Q = q)
        }
        fv1 <- 0; fv2 <- 0
        if ( t > R ) {
            fv1 <- max(     R * integrate(Fun.0, 0, t-R)$value +         integrate(Fun.1, t-R, t)$value, 0)
            fv2 <- max( R^2/2 * integrate(Fun.0, 0, t-R)$value + (1/2) * integrate(Fun.2, t-R, t)$value, 0)
        } else {
            fv1 <- max( integrate(Fun.1, 0, t)$value, 0)
            fv2 <- max( (1/2) * integrate(Fun.2, 0, t)$value, 0)
        }
        (eta/(eta*R + theta*R^2/2)) * fv1 + (theta/(eta*R + theta*R^2/2)) * fv2
      }
    }
    max.T <- func(T, rho)
    InfT <- sapply(t, function(x) {
                   round(func(x, rho)/max.T, 5)
            })
    list(Info.fractions=InfT, Event.prob=func(T, rho=0), Total.censor.prob=1-func(T, rho=0))
}

