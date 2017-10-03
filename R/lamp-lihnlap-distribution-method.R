#' Lihn-Laplace process and distribution
#' 
#' Implements some aspects of Lihn-Laplace process
#'
#' @param n numeric, number of observations.
#' @param x numeric, vector of responses.
#' @param s numeric, vector of responses for characteristic function.
#' @param t numeric, the time parameter, of which the variance is t.
#' @param mu numeric, location parameter, default is 0.
#' @param beta numeric, skewness parameter according to skewed lambda distribution, default is 0.
#' @param convo numeric, the convolution number, default is 1, which is without convolution.
#' 
#' @return numeric, standard convention is followed: 
#'         d* returns the density, 
#'         p* returns the distribution function, 
#'         q* returns the quantile function, and 
#'         r* generates random deviates.
#'         The following are our extensions:
#'         k* returns the first 4 cumulants, skewness, and kurtosis,
#'         cf* returns the characteristic function.
#'
#' @keywords Lihn-Laplace
#'
#' @author Stephen H-T. Lihn
#'
#' @importFrom stats runif
#' @importFrom stats rexp
#' @importFrom stats rnorm
#'
#' @export rlihnlap
#' @export dlihnlap
#' @export cflihnlap
#' @export klihnlap
#' @export dlihnlap_poly
#' 
### <======================================================================>
rlihnlap <- function(n, t=1, convo=1, beta=0, mu=0) {
    if (convo==1) {
        if (beta==0) {
            L <- stats::rnorm(n)*sqrt(stats::rexp(n,1))*sqrt(t)
            return(L+mu)
        }
        # skewed Laplace implementation
        B0 <- sqrt(1+beta^2/4)
        B1 <- B0+beta/2 # B plus
        B2 <- B0-beta/2 # B minus
        W <- B2/B0/2
        sgn <- sign(stats::runif(n)-W)
        sigma <- ifelse(sgn > 0, B1, B2)
        L <- stats::rexp(n, 1/sigma)*sgn
        sd <- sqrt(2+beta^2)
        return(L/sd*sqrt(t)+mu)
    }
    # compose convolution recursively
    if (convo < 1) stop("fractional convolution is not supported")
    if (convo!=ceiling(convo)) warning("real-number convolution is only experimental")
    
    L <- rep(0,n)
    frac <- convo+1-ceiling(convo)
    if (frac==0) frac <- 1
    convo1 <- ceiling(convo)
    for (i in 1:convo1) {
        f <- if (i==convo1) frac else 1
        L <- L + rlihnlap(n, t*f, beta=beta)
    }
    return(L/sqrt(convo)+mu)
}
### <---------------------------------------------------------------------->
#' @rdname rlihnlap
dlihnlap <- function(x, t=1, convo=1, beta=0, mu=0) {
    if (convo<=0) stop(paste("ERROR: Lihn-Laplace pdf is not supported for convolution:", convo))

    x <- x-mu # take care of location
    m <- convo
    Sm <- sqrt(t/(2+beta^2)/m)
    B0 <- sqrt(1+beta^2/4)
    w <- abs(B0*x/Sm)

    C <- 1/gamma(m)/sqrt(pi)/Sm 
    K <- (w/2)^(m-1/2) * besselK(w,m-1/2) 
    C0 <- if (convo > 1/2) C*gamma(m-1/2)/2 else NaN
    CK <- ifelse(w==0, C0, C*K) # K(0) is Inf
    Eb <- (B0^-2)^(m-1/2) * exp(beta*x/2/Sm)
    return(CK*Eb)    
}
### <---------------------------------------------------------------------->
#' @rdname rlihnlap
cflihnlap <- function(s, t=1, convo=1, beta=0, mu=0) {

    m <- convo
    Sm <- sqrt(t/(2+beta^2)/m)

    cf_mu <- exp(1i*mu*s)    
    cf1 <- 1 - 1i*beta*s*Sm + s^2*Sm^2 
    return(cf_mu*cf1^(-m))
}
### <---------------------------------------------------------------------->
#' @rdname rlihnlap
klihnlap <- function(t=1, convo=1, beta=0, mu=0) {
    kL <- function(i) {
        if (length(i)>1) return(sapply(i,kL))
        B <- 2+beta^2
        if (i==1) return(mu+beta*sqrt(convo*t/B))
        if (i==2) return(t)
        if (i==3) return(2*t^(3/2)/sqrt(convo)*beta*(3+beta^2)/B^(3/2))
        if (i==4) return(6*t^2/convo*(2+4*beta^2+beta^4)/B^2)
        stop
    }
    skw <- kL(3)/kL(2)^(3/2)
    kur <- kL(4)/kL(2)^2 + 3
    
    k <- c(kL(1:4), skw, kur)
    names(k) <- c("mean", "var", "k3", "k4", "skewness", "kurtosis")
    return(k)
}
### <---------------------------------------------------------------------->
#' @rdname rlihnlap
dlihnlap_poly <- function(x, t=1, convo=1, beta=0, mu=0) {
    x <- x-mu # take care of location
    m <- convo
    Sm <- sqrt(t/(2+beta^2)/m)
    B0 <- sqrt(1+beta^2/4)
    B1 <- B0+beta/2
    B2 <- B0-beta/2
    w <- abs(B0*x/Sm)
    
    B <- ifelse(x>0, B1, B2)
    Ey <- exp(-abs(x/Sm/B))
    
    if (convo==1) return(Ey/Sm/(2*B0))
    if (convo==2) return(Ey/Sm*(w+1)/(4*B0^3))
    if (convo==3) {
        p <- w^2+3*w+3
        return(Ey/Sm*p/(16*B0^5))
    }
    if (convo==4) {
        p <- w^3+6*w^2+15*w+15
        return(Ey/Sm*p/(96*B0^7))
    }
    
    stop(paste("ERROR: Lihn-Laplace pdf is not supported for convolution:", convo))
}
### <---------------------------------------------------------------------->
