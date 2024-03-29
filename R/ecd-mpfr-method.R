#' Wrapper to convert numeric to mpfr
#'
#' Convert numeric to mpfr for ecd calculations.
#' \code{ecd.mp1} is the constant 1 wrapped in mpfr class.
#' \code{ecd.mppi} is the function to obtain pi from Rmpfr with an optional precision. This is used to implement \code{ecd.erfq}.
#' \code{ecd.gamma} is a wrapper on \code{ecld.gamma}, which is the incomplete gamma function.
#' \code{ecd.erf} is a wrapper on \code{Rmpfr::erf}.
#' \code{ecd.erfcx} is a wrapper on \code{Rmpfr::erfcx}.
#' \code{ecd.erfc} is a wrapper on \code{Rmpfr::erfc}. This is used to implement \code{ecd.erfq}.
#' \code{ecd.dawson} is a wrapper on \code{gsl::dawson}. Dawson function is used to implement \code{ecd.erfq}.
#' \code{ecd.erfi} is the imaginary scaled error function, which is implemented through \code{ecd.dawson}.
#' \code{ecd.devel} is a developer tool to size down intensive mpfr tests for CRAN. Set \code{ecd_devel} in R options or OS env to change its value.
#'
#' @param x a numeric vector or list. If \code{x} is mpfr class,
#'          it will be passed through.
#' @param precBits an integer for mpfr precBits.
#'                 Default is from \code{getOption("ecd.precBits")}.
#' @param s numeric vector, for the order of incomplete gamma function
#' @param na.stop logical, stop if NaN is generated. The default is \code{TRUE}.
#'
#' @return The \code{mpfr} object
#'
#' @keywords utility
#'
#' @export ecd.mpfr
#' @export ecd.mp1
#' @export ecd.mppi
#' @export ecd.gamma
#' @export ecd.erf
#' @export ecd.erfc
#' @export ecd.erfcx
#' @export ecd.erfi
#' @export ecd.dawson
#' @export ecd.devel
#'
#' @importFrom Rmpfr mpfr
#' @importFrom Rmpfr Const
#' @importFrom Rmpfr erf
#' @importClassesFrom Rmpfr mpfr
#' @importFrom gsl dawson
#'
#' @examples
#' x <- ecd.mpfr(1)
#' y <- ecd.mpfr(c(1,2,3))
#' z <- ecd.mp1
#' p <- ecd.mppi()
### <======================================================================>
"ecd.mpfr" <- function(x, precBits=getOption("ecd.precBits"))
{
	if (is.null(precBits)) {
	    precBits <- 120L
	} else if (is.na(precBits) | is.nan(precBits) | precBits < 2L) {
	    precBits <- 120L
	}
	
    if (is(x,"numeric")) return(Rmpfr::mpfr(x, precBits))
    if (is(x,"mpfr")) return(x)
    
    # the list is usually caused by mclapply (the apply family)
    if (is(x,"list")) return(new("mpfr", unlist(x)))
    
    
    stop(paste("Unknown class to convert to mpfr:", c))
}
### <---------------------------------------------------------------------->
#' @rdname ecd.mpfr
ecd.mp1 <- ecd.mpfr(1)
### <---------------------------------------------------------------------->
#' @rdname ecd.mpfr
ecd.mppi <- function(precBits=getOption("ecd.precBits"))
{
    if (is.null(precBits)) {
        precBits <- 120L
    } else if (is.na(precBits) | is.nan(precBits) | precBits < 2L) {
        precBits <- 120L
    }
    return(Rmpfr::Const("pi", precBits))
}
### <---------------------------------------------------------------------->
#' @rdname ecd.mpfr
ecd.gamma <- function(s,x,na.stop=TRUE) ecld.gamma(s,x,na.stop)
### <---------------------------------------------------------------------->
#' @rdname ecd.mpfr
ecd.erf <- function(x) {
    if (is(x,"complex")) stop("unable to handle complex input, RcppFaddeeva no longer available")
    Rmpfr::erf(x)
}
### <---------------------------------------------------------------------->
#' @rdname ecd.mpfr
ecd.erfc <- function(x) {
    if (is(x,"complex")) stop("unable to handle complex input, RcppFaddeeva no longer available")
    Rmpfr::erfc(x)
}
### <---------------------------------------------------------------------->
#' @rdname ecd.mpfr
ecd.erfcx <- function(x) {
    if (is(x,"complex")) stop("unable to handle complex input, RcppFaddeeva no longer available")
    Rmpfr::erfc(x)*exp(x^2)
}
### <---------------------------------------------------------------------->
#' @rdname ecd.mpfr
ecd.dawson <- function(x) {
    if (is(x,"complex")) stop("unable to handle complex input, RcppFaddeeva no longer available")
    one <- x*0.0 + 1 # this is to preserve the MPFR type
    y <- gsl::dawson(ecd.mp2f(x))
    return(y*one)
}
### <---------------------------------------------------------------------->
#' @rdname ecd.mpfr
ecd.erfi <- function(x) {
    if (is(x,"complex")) stop("unable to handle complex input, RcppFaddeeva no longer available")
    y <- ecd.dawson(x)
    return(y*exp(x^2)*2/sqrt(pi))
}
### <---------------------------------------------------------------------->
#' @rdname ecd.mpfr
# export ecd_devel=1 to reduce the scope when not in "devel" mode
# CRAN has a 20-minute limit on tests
ecd.devel <- function() {
    dev <- getOption("ecd_devel")
    if (! is.null(dev)) { if (dev == TRUE) return(TRUE) }
    if (nchar(Sys.getenv("ecd_devel"))>0) return(TRUE)
    return(FALSE)
}
### <---------------------------------------------------------------------->


