#' Stable lambda distribution
#'
#' Implements some aspects of the stable lambda (SL) distribution
#'
#' @param n numeric, number of observations.
#' @param x numeric, vector of responses.
#' @param s numeric, vector of responses for characteristic function.
#' @param nu numeric, vector of nu in the pdf integrand, starting from 0 (not nu0).
#' @param t numeric, the time parameter, where the variance is t, default is 1.
#' @param nu0 numeric, the location parameter, default is 0.
#' @param theta numeric, the scale parameter, default is 1.
#' @param convo numeric, the convolution number, default is 1.
#' @param beta.a numeric, the skewness parameter, default is 0. This number is annualized by sqrt(t).
#' @param mu numeric, the location parameter, default is 0.
#' @param lambda numeric, the shape parameter, default is 4.
#' @param method character, method of characteristic function (CF) calculation. Default is "a".
#'               Method a uses cflihnlap x dstablecnt.
#'               Method b uses dlihnlap x cfstablecnt.
#'               Method c uses direct integration on PDF up to 50 stdev.
#'               They should yield the same result.
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
#' @keywords Modified-Lambda
#'
#' @author Stephen H-T. Lihn
#'
#' @importFrom stats rexp
#' @importFrom stats rnorm
#'
#' @export rqsl
#' @export dqsl
#' @export cfqsl
#' @export kqsl
#' @export rsl
#' @export dsl
#' @export cfsl
#' @export ksl
#' @export qsl_variance_analytic
#' @export qsl_kurtosis_analytic
#' @export qsl_skewness_analytic
#' @export qsl_std_pdf0_analytic
#' @export qsl_pdf_integrand_analytic
#'
### <======================================================================>
rqsl <- function(n, t=1, nu0=0, theta=1, convo=1, beta.a=0, mu=0) {
    rsl(n, t=t, nu0=nu0, theta=theta, convo=convo, beta.a=beta.a, mu=mu, lambda=4)
}
### <---------------------------------------------------------------------->
#' @rdname rqsl
rsl <- function(n, t=1, nu0=0, theta=1, convo=1, beta.a=0, mu=0, lambda=4) {
    L <- rlihnlap(n, t=t, convo=convo, beta=beta.a*sqrt(t))
    N <- rstablecnt(n, alpha=2/lambda, nu0=nu0, theta=theta)
    return(L*N+mu)
}
### <---------------------------------------------------------------------->
#' @rdname rqsl
dsl <- function(x, t=1, nu0=0, theta=1, convo=1, beta.a=0, mu=0, lambda=4) {
    integrand_xp <- function(x, nu) {
        L <- dlihnlap((x-mu)/(nu+nu0), t=t, convo=convo, beta=beta.a*sqrt(t))
        N <- dstablecnt(nu, alpha=2/lambda, nu0=0, theta=theta)
        v <- N*L/(nu+nu0)
        return(ifelse(is.na(v), 0, v)) # bypass the singular points
    }

    pdf <- function(x) {
        f <- function(nu) integrand_xp(x, nu)
        S <- theta*sqrt(t) # scale
        Va <- integrate(f, lower=0.000001*S, upper=0.001*S)$value # very small and problematic
        V0 <- integrate(f, lower=0.001*S, upper=0.1*S)$value
        V1 <- integrate(f, lower=0.1*S, upper=1*S)$value
        V2 <- integrate(f, lower=1*S,   upper=10*S)$value
        V3 <- integrate(f, lower=10*S,  upper=100*S)$value
        V4 <- integrate(f, lower=100*S, upper=1000*S)$value
        return(Va+V0+V1+V2+V3+V4)
    }

    ld <- ecld(4)
    V <- ecd.mcsapply(ld, x, pdf)
    return(V)
}
### <---------------------------------------------------------------------->
#' @rdname rqsl
dqsl <- function(x, t=1, nu0=0, theta=1, convo=1, beta.a=0, mu=0) {
    dsl(x, t=t, nu0=nu0, theta=theta, convo=convo, beta.a=beta.a, mu=mu, lambda=4)
}
### <---------------------------------------------------------------------->
#' @rdname rqsl
kqsl <- function(t=1, nu0=0, theta=1, convo=1, beta.a=0, mu=0) {
    ksl(t=t, nu0=nu0, theta=theta, convo=convo, beta.a=beta.a, mu=mu, lambda=4)
}
### <---------------------------------------------------------------------->
#' @rdname rqsl
ksl <- function(t=1, nu0=0, theta=1, convo=1, beta.a=0, mu=0, lambda=4) {
    kL <- klihnlap(t=t, convo=convo, beta=beta.a*sqrt(t))
    kN <- kstablecnt(alpha=2/lambda, nu0=nu0, theta=theta)
    mC <- k2mnt(kN)*k2mnt(kL)
    kC <- mnt2k(mC)
    skw <- kC[3]/kC[2]^(3/2)
    kur <- kC[4]/kC[2]^2 + 3

    kC[1] <- kC[1]+mu
    k <- c(kC[1:4], skw, kur)
    names(k) <- c("mean", "var", "k3", "k4", "skewness", "kurtosis")
    return(k)
}
### <---------------------------------------------------------------------->
#' @rdname rqsl
qsl_kurtosis_analytic <- function(t=1, nu0=0, theta=1, convo=1, beta.a=0) {
    if (beta.a != 0) stop("qsl_kurtosis_analytic does not support non-zero beta")
    m <- convo
    beta <- beta.a*sqrt(t)
    nu1 <- nu0+6*theta
    A <- nu1^4 + (144+96*m)*theta^2*nu1^2 + 768*(1+m)*theta^3*nu1 + (4032+3456*m)*theta^4
    B <- (nu1^2+24*theta^2)^2
    return(3+(3/m)*A/B)
}
### <---------------------------------------------------------------------->
#' @rdname rqsl
qsl_skewness_analytic <- function(t=1, nu0=0, theta=1, convo=1, beta.a=0) {
    m <- convo
    beta <- beta.a*sqrt(t)
    nu1 <- nu0+6*theta

    P <- 144*beta^2 + 432 + 144*m*(2+beta^2)
    Q <- 192*beta^2*m^2 + 576*(2+beta^2)*m + 384*beta^2 + 1152

    A <- 2*(3+beta^2)*nu1^3 + P*theta^2*nu1 + Q*theta^3
    B <- (2+beta^2)*nu1^2 + 24*theta^2*(2+beta^2*(1+m))
    # print(paste(A,B^(3/2))) # debug this complex beast!
    return(beta/sqrt(m)*A/B^(3/2))
}
### <---------------------------------------------------------------------->
#' @rdname rqsl
qsl_variance_analytic <- function(t=1, nu0=0, theta=1, convo=1, beta.a=0) {
    m <- convo
    beta <- beta.a*sqrt(t)
    nu1 <- nu0+6*theta

    var_annl <- (nu1^2+24*theta^2) + m*beta^2/(2+beta^2)*(24*theta^2)
    return(var_annl*t)
}
### <---------------------------------------------------------------------->
#' @rdname rqsl
qsl_std_pdf0_analytic <- function(t=1, nu0=0, theta=1, convo=1, beta.a=0) {
    if (beta.a != 0) stop("qsl_std_pdf0_analytic does not support non-zero beta")

    if (convo > 100) convo <- ecd.mp1*convo # to ensure gamma function works

    m <- convo
    C <- gamma(m-1/2)*sqrt(2*m)/gamma(m)/8/pi
    n <- nu0/theta
    pdf0 <- 2*sqrt(pi)-pi*sqrt(n)*exp(n/4)*ecd.erfc(sqrt(n/4))
    var_annl <- (n+6)^2+24
    ecd.mp2f(C*pdf0*sqrt(var_annl)) # render numeric always
}
### <---------------------------------------------------------------------->
#' @rdname rqsl
qsl_pdf_integrand_analytic <- function(x, nu, t=1, nu0=0, theta=1, convo=1, beta.a=0, mu=0) {
    if (convo > 100) convo <- ecd.mp1*convo # to ensure gamma function works

    m <- convo
    beta <- beta.a*sqrt(t)
    B0 <- sqrt(1+beta^2/4)
    sm <- sqrt(t/m/(2+beta^2))
    w <- B0*(x-mu)/sm/(nu+nu0)

    C <- 4*pi*gamma(m)*sm*theta^(3/2)
    A <- sqrt(nu)/(nu+nu0)
    K <- abs(w/2/B0^2)^(m-1/2) * besselK(abs(w),m-1/2)
    B <- exp(w*beta/2/B0 - nu/4/theta)
    ecd.mp2f(1/C*A*K*B) # render numeric always
}
### <---------------------------------------------------------------------->
#' @rdname rqsl
cfqsl <- function(s, t=1, nu0=0, theta=1, convo=1, beta.a=0, mu=0, method="a") {
    cfsl(s, t=t, nu0=nu0, theta=theta, convo=convo, beta.a=beta.a, mu=mu, lambda=4, method=method)
}
### <---------------------------------------------------------------------->
#' @rdname rqsl
cfsl <- function(s, t=1, nu0=0, theta=1, convo=1, beta.a=0, mu=0, lambda=4, method="a") {
    # TODO the upper bounds of the integrals are all based on quartic distribution
    # User needs to be careful when extending to other lambda's
    # Future development!!!
    if (method=="a") {
        N <- 100000
        nu <- nu0 + seq(0, 100, length.out=N)*theta
        M <- exp(1i*mu*s)
        L <- cflihnlap(s*nu, t=t, convo=convo, beta=beta.a*sqrt(t))
        N <- dstablecnt(nu, alpha=2/lambda, nu0=nu0, theta=theta)
        dnu <- nu[2]-nu[1]
        return(M*sum(L*N)*dnu)
    }
    if (method=="b") {
        N <- 100000
        x <- seq(-50, 50, length.out=N+1)*sqrt(t)
        dx <- x[2]-x[1]
        M <- exp(1i*mu*s)
        L <- dlihnlap(x, t=t, convo=convo, beta=beta.a*sqrt(t))
        N <- cfstablecnt(s*x, alpha=2/lambda, nu0=nu0, theta=theta)
        return(M*sum(L*N)*dx)
    }
    if (method=="c") {
        # TODO this is a hack for other lambda, but our primary use case is lambda=4
        sd <- qsl_variance_analytic(t=t, nu0=nu0, theta=theta, convo=convo, beta.a=beta.a)^0.5
        M <- exp(1i*mu*s)
        p_re <- function(x) cos(s*x)*dsl(x, t=t, nu0=nu0, theta=theta, convo=convo, beta.a=beta.a)
        p_im <- function(x) sin(s*x)*dsl(x, t=t, nu0=nu0, theta=theta, convo=convo, beta.a=beta.a)
        intg <- function(f,a,b) integrate(f, lower=sd*a, upper=sd*b)$value
        cf_re <- c(intg(p_re, 0, 5), intg(p_re, 5, 10), intg(p_re, 10, 50),
                   intg(p_re, -5, 0), intg(p_re, -10, -5), intg(p_re, -50, -10))

        cf_im <- c(intg(p_im, 0, 5), intg(p_im, 5, 10), intg(p_im, 10, 50),
                   intg(p_im, -5, 0), intg(p_im, -10, -5), intg(p_im, -50, -10))

        # print(paste(sd, cf_re, cf_im))
        return(M*(sum(cf_re) + 1i*sum(cf_im)))
    }

    stop(paste("cfsl: method is not supported", method))
}
### <---------------------------------------------------------------------->
