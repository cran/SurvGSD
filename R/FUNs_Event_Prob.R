#### The subroutine for calculating P[Delta_star(t_k)=1] ####

FUN.Event.Prob <- function(t, R, T, FUN.C, para0, para1, haz.r, rho, eta, theta, px) {
       q0 = para0[1]; mu0 = para0[2]; sigma0 = para0[3]
       q1 = para1[1]; mu1 = para1[2]; sigma1 = para1[3]

       if ( is.null(haz.r) == FALSE ) {
          k <- q0^{-2}; alpha <- q0/sigma0; beta <- exp(mu0)*(q0^2)^(sigma0/q0)
          f0 <- function(y) { dgengamma(y, mu0, sigma0, q0) }
          F0.bar <- function(y) { pgengamma(y, mu0, sigma0, q0, lower.tail = FALSE) }
          if ( is.numeric(haz.r) == TRUE ) {
               f1 <- function(y) { haz.r * (F0.bar(y))^(haz.r-1) * f0(y) }
          } else {
            #haz0 <- function(y) {
                     if (q0==1 && sigma0==1) {haz0 <- function(y) rep(1/exp(mu0),length(y))}
                     else if (q0==1 && sigma0!= 1) {haz0 <- function(y) hweibull(y, shape=alpha, scale=beta)}
                     else if (q0==sigma0 && q0!=1) {haz0 <- function(y) hgamma(y, shape=k, rate=1/beta)}
                     else if (q0==0 && sigma0!=0)  {haz0 <- function(y) hlnorm(y, meanlog=mu0, sdlog=sigma0)}
                     else {
                       haz0 <- function(y) f0(y)/F0.bar(y)
                     }
                  #}
            haz1 <- function(t) {
                      haz.r(t) * haz0(t)
                  }
            f1 <- function(y) {
                      sapply(y, function(x) {
                              haz1(x) * exp(-integrate(haz1,0,x)$value)
                      })
                  }
         }
       } else {
         f0 <- function(y) { dgengamma(y, mu0, sigma0, q0)}
         f1 <- function(y) {
                   dgengamma(y, mu=mu1, sigma=sigma1, Q=q1)
               }
       }

       pooled.fY <- function(y) { (1-px) * f0(y)  + px * f1(y) }

       prob.func <- function(t){
                        Fun.0 <- function(y) {
                                     ( 1 - FUN.C(y) ) * pooled.fY(y)
                                 }
                        Fun.1 <- function(y) {
                                     (t-y) * ( 1 - FUN.C(y) ) * pooled.fY(y)
                                 }
                        Fun.2 <- function(y) {
                                     (t-y)^2 * ( 1 - FUN.C(y) ) * pooled.fY(y)
                                 }
                        fv1 <- 0; fv2 <- 0
                        if ( t > R ) {
                             fv1 <-     R * integrate(Fun.0, 0, t-R)$value +         integrate(Fun.1, t-R, t)$value
                             fv2 <- R^2/2 * integrate(Fun.0, 0, t-R)$value + (1/2) * integrate(Fun.2, t-R, t)$value
                        } else {
                             fv1 <-         integrate(Fun.1, 0, t)$value
                             fv2 <- (1/2) * integrate(Fun.2, 0, t)$value
                        }
                        (eta/(eta*R + theta*R^2/2)) * fv1 + (theta/(eta*R + theta*R^2/2)) * fv2
                    }

       tc <- c(t,T)
       prob.vec <- sapply(tc, function(x) {
                          round(prob.func(x), 5)
                   })
       list(prob.vec=prob.vec)
}


FUN.Event.Prob.LLG <- function(t, R, T, FUN.C, para0, para1, haz.r, rho, eta, theta, px){

  xi0 = para0[1]; zeta0 = para0[2]
  xi1 = para1[1]; zeta1 = para1[2]

  if ( is.null(haz.r) == FALSE ) {
    f0 <- function(y) { (zeta0/xi0)*(y/xi0)^{zeta0-1}/( 1+(y/xi0)^zeta0 )^2 }
    F0.bar <- function(y) { 1/( 1+(y/xi0)^zeta0 ) }
    if ( is.numeric(haz.r) == TRUE ) {
      f1 <- function(y) { haz.r * (F0.bar(y))^(haz.r-1) * f0(y) }
    } else {
      haz0 <- function(y) { (zeta0/xi0)*(y/xi0)^{zeta0-1}/( 1+(y/xi0)^zeta0 ) }
      haz1 <- function(t) {
        haz.r(t) * haz0(t)
      }
      f1 <- function(y) {
        sapply(y, function(x) {
          haz1(x) * exp(-integrate(haz1,0,x)$value)
        })
      }
    }
  } else {
    f0 <- function(y) { (zeta0/xi0)*(y/xi0)^{zeta0-1}/( 1+(y/xi0)^zeta0 )^2 }
    f1 <- function(y) { (zeta1/xi1)*(y/xi1)^{zeta1-1}/( 1+(y/xi1)^zeta1 )^2 }
  }

  pooled.fY <- function(y) { (1-px) * f0(y)  + px * f1(y) }

  prob.func <- function(t){
    Fun.0 <- function(y) {
      ( 1 - FUN.C(y) ) * pooled.fY(y)
    }
    Fun.1 <- function(y) {
      (t-y) * ( 1 - FUN.C(y) ) * pooled.fY(y)
    }
    Fun.2 <- function(y) {
      (t-y)^2 * ( 1 - FUN.C(y) ) * pooled.fY(y)
    }
    fv1 <- 0; fv2 <- 0
    if ( t > R ) {
      fv1 <-     R * integrate(Fun.0, 0, t-R)$value +         integrate(Fun.1, t-R, t)$value
      fv2 <- R^2/2 * integrate(Fun.0, 0, t-R)$value + (1/2) * integrate(Fun.2, t-R, t)$value
    } else {
      fv1 <-         integrate(Fun.1, 0, t)$value
      fv2 <- (1/2) * integrate(Fun.2, 0, t)$value
    }
    (eta/(eta*R + theta*R^2/2)) * fv1 + (theta/(eta*R + theta*R^2/2)) * fv2
  }

  tc <- c(t,T)
  prob.vec <- sapply(tc, function(x) {
    round(prob.func(x), 5)
  })
  list(prob.vec=prob.vec)

}
