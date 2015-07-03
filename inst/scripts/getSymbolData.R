rm(list = ls())

library(quantmod)

#' Function to download closing symbol data. This function must be used to download
#' and store data in the data folder when symbol list is updated or dates change. 
#' After executing the function store the symdata object in the data folder.
#' 
#' @details
#' This method takes in a list of symbols and downloads them for specific dates
#' @param symbols list of symbols in string format
#' @param from start date in yyyy-mm-dd format
#' @param to   end   date in yyyy-mm-dd format if missing then defaults to current date
#' @param returns boolean to indicate if prices must be converted to log returns
#' @param frequency string, "monthly"/"daily" to indicate the type of time series.
#'        defaults to monthly
#' @author Rohit Arora
#' 
get.closing.symdata <- function(symbols, from, to=Sys.Date(), returns=TRUE,
                                        frequency = c("monthly", "daily"))
  {
  
  if (length(symbols) == 0) stop("No symbols to download")
  
  d <- try( as.Date( from, format= "%Y-%m-%d" ) )
  if( class( d ) == "try-error" || is.na( d ) ) stop( "Invalid start date" )
  
  d <- try( as.Date( to, format= "%Y-%m-%d" ) )
  if( class( d ) == "try-error" || is.na( d ) ) stop( "Invalid end date" )  
  

  getSymbols(symbols, adjust=TRUE, from = from, to = to)
  
  frequency = frequency[1]
  
  for(symbol in symbols) {
    temp0  <- temp <- adjustOHLC(get(symbol), symbol.name=symbol)
    if(frequency == "monthly") {
      temp    <- to.monthly(temp0, indexAt='endof', drop.time=FALSE)
      colnames(temp) <- colnames(temp0)      
    }
    assign(x=symbol, value=temp)
  }
  
  symbolData <- do.call(merge, lapply(symbols, function(x) Cl(get(x))))
  colnames(symbolData) <- symbols
  if (returns) symbolData <- diff(log(symbolData))[-1,]
  symbolData
}

symbols <-  c('BABA', 'TWTR', 'LNKD', 'YHOO', 'GE', 'LAZ', 'V')
symdata <- get.closing.symdata(symbols, from = "2007-04-01")

remove(symbols); remove(get.closing.symdata)