#' Read sample data config
#' 
#' Read sample data config. This is used in ecd.data() and related tests.
#'
#' @param symbol Optional, if provided, only return config related to it.
#'               default = NULL
#' @param extdata_dir optionally specify user's own extdata folder which can
#'               also comes with user's data_conf.yml
#'
#' @return The data.frame object for the config
#'
#' @keywords sample data
#'
#' @examples
#' c <- .ecd.data_config()
#' c <- .ecd.data_config("dji")
#' @noRd
### <======================================================================>
".ecd.data_config" <- function(symbol=NULL, extdata_dir=NULL)
{
    filename <- "data_conf.yml" # the data config should be in this file
    conf_file <- .ecd.locate_file(filename, extdata_dir=extdata_dir)
    # read the conf file
    if (! file.exists(conf_file)) {
        stop(paste("conf_file does not exist:", conf_file))
    }
    conf_all <- yaml.load_file(conf_file)
    
    m0 <- function() {
        cols = c("symbol", "cols", "col_in", "date_format", "test_date", "test_val")
        matrix(vector(), 0, 6, dimnames=list(c(), cols)) 
    }
    add <- function(m, c) {
        if (is.null(symbol) || (symbol == c[1])) {
            m2 <- rbind(m, c)
            rownames(m2) <- NULL
            m2
        }else{
            m
        }
    }
    
    mx <- m0()
    for (s in names(conf_all)) {
        cf1 <- conf_all[[s]]
        mx <- add(mx, c(s, cf1$cols, cf1$col_in, cf1$date_format, cf1$test_date, cf1$test_val))
    }
    df <- data.frame(mx, stringsAsFactors=FALSE)
    rownames(df) <- df$symbol
    return(df)

    # ----------------------------------------------------------
    # obsolete, keep them in comment for future reference
    #dt <- "Close"
    #dta <- "Adjusted"
    #fmt <- "%m/%d/%Y"

    #df <- m0()
    #df <- add(df, c("dji",   10, dta, fmt, "19350103", 105.14))
    #df <- add(df, c("dji28", 10, dta, fmt, "19281002", 238.14))
    #df <- add(df, c("ge",     3, dta, fmt, "19890301", 2.283))
    #df <- add(df, c("gox",    7, dta, fmt, "19960501", 129.44))
    #df <- add(df, c("gold",   2, dt, fmt, "19720301", 47.8))
    #df <- add(df, c("hkhsi",  2, dt, fmt, "19870302", 2894.3))
    #df <- add(df, c("ndq",    2, dt, fmt, "19710301", 101.78))
    #df <- add(df, c("nikkei", 2, dt, fmt, "19840301", 9920))
    #df <- add(df, c("r10y",   2, dt, fmt, "19620301", 3.98))
    #df <- add(df, c("spxrv",  2, dt, fmt, "20000131", 0.00019816))
    #df <- add(df, c("rus", 2, dt, fmt, "19871201", 111.27))
    #df <- add(df, c("spx", 7, dt, fmt, "19500501", 18.22))
    #df <- add(df, c("chf", 2, dt, fmt, "20080505", 1.0539))
    #df <- add(df, c("eur", 2, dt, fmt, "19750131", 0.7423))
    #df <- add(df, c("vix", 7, dt, fmt, "19900110", 22.44))
    #df <- add(df, c("wti", 2, dt, fmt, "19860106", 26.53))
    
    #dto  = "LAST"
    #fmto = "%Y%m%d"
    #df <- add(df, c("spxoption",  201109, dto, fmto, NA, NA)) # CBOE market data
    #df <- add(df, c("spxoption2", 201506, dto, fmto, NA, NA)) # CBOE market data
    #df <- add(df, c("spxoption2f", 201506, dto, fmto, NA, NA)) # CBOE market data
    #df <- add(df, c("spxoption3", 201508, dto, fmto, NA, NA)) # CBOE market data
    #df <- add(df, c("gldoption",  201109, dto, fmto, NA, NA)) # CBOE market data
    #df <- add(df, c("gldoption3", 201508, dto, fmto, NA, NA)) # CBOE market data
    
    #df <- data.frame(df, stringsAsFactors=FALSE)
    #rownames(df) <- df$symbol
    #df
}
### <---------------------------------------------------------------------->
