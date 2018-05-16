### Function for calculating the probability of loss to follow-up (P[W<Y, W<T-V])

FUN.LFU.GGD <- function(h, R, T, q, mu, sigma, eta, theta){
            Fun.0 <- function(y) {
                  pgengamma(y, mu = mu, sigma = sigma, Q = q, lower.tail = FALSE) * (1/h * (y<h))
            }
            Fun.1 <- function(y) {
                  (T-y) * ( pgengamma(y, mu = mu, sigma = sigma, Q = q, lower.tail = FALSE) ) * (1/h * (y<h))
            }
            Fun.2 <- function(y) {
                  (T-y)^2 * ( pgengamma(y, mu = mu, sigma = sigma, Q = q, lower.tail = FALSE) ) * (1/h * (y<h))
            }
            fv1 <- 0; fv2 <- 0

            fv1 <- max(     R * integrate(Fun.0, 0, T-R)$value +         integrate(Fun.1, T-R, T)$value, 0)
            fv2 <- max( R^2/2 * integrate(Fun.0, 0, T-R)$value + (1/2) * integrate(Fun.2, T-R, T)$value, 0)

           (eta/(eta*R + theta*R^2/2)) * fv1 + (theta/(eta*R + theta*R^2/2)) * fv2

}


FUN.LFU.LLG <- function(h, R, T, xi, zeta, eta, theta){
  Fun.0 <- function(y) {
    ( 1/( 1+(y/xi)^zeta ) ) * (1/h * (y<h))
  }
  Fun.1 <- function(y) {
    (T-y) * ( 1/( 1+(y/xi)^zeta ) ) * (1/h * (y<h))
  }
  Fun.2 <- function(y) {
    (T-y)^2 * ( 1/( 1+(y/xi)^zeta ) ) * (1/h * (y<h))
  }
  fv1 <- 0; fv2 <- 0

  fv1 <- max(     R * integrate(Fun.0, 0, T-R)$value +         integrate(Fun.1, T-R, T)$value, 0)
  fv2 <- max( R^2/2 * integrate(Fun.0, 0, T-R)$value + (1/2) * integrate(Fun.2, T-R, T)$value, 0)

  (eta/(eta*R + theta*R^2/2)) * fv1 + (theta/(eta*R + theta*R^2/2)) * fv2
}
