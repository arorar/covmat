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

menu <- function() {
  print("1 . Download symbols for factor data")  
  print("2 . Download symbols for filtering using RMT")  
  print("3 . Exit")  
}

while(TRUE) {
  menu()
  opt <- readline("Select any option : ")
  
  switch(opt, 
         "1" = {
           symbols <-  c('BABA', 'TWTR', 'LNKD', 'YHOO', 'GE', 'LAZ', 'V')
           
           symdata <- get.closing.symdata(symbols, from = "2007-04-01")
           l <- ls(); rm(list = c(l[l != "symdata"],"l")) 
           break
         }, 
         "2" = {
           symbols <- c('ABT','ANF','ACN','ACE','ACT','ADBE','AMD','AES',
                        'AET','AFL','A','GAS','APD','ARG','AKAM','AA','ALXN',
                        'ATI','AGN','ALL','ALTR','AMZN','AEE','AEP','AXP',
                        'AIG','AMT','AMP','ABC','AMGN','APH','APC','ADI','AON',
                        'APA','AIV','APOL','AAPL','AMAT','ADM','AIZ','T','ADSK',
                        'ADP','AN','AZO','AVB','AVY','AVP','BHI','BLL','BAC',
                        'BCR','BAX','BBT','BDX','BBBY','BMS','BBY','BIIB',
                        'BLK','HRB','BA','BWA','BXP','BSX','BMY','BRCM',
                        'CA','CVC','COG','CAM','CPB','COF','CAH')
           
           largesymdata <- get.closing.symdata(symbols, from = "2007-04-01")
           largesymdata <- na.omit(largesymdata)
           l <- ls(); rm(list = c(l[l != "largesymdata"],"l"))
           break
         },
         "3" = {
            print('Exiting...')
           break
         },
         {
            print("Invalid Option. Please re-select")  
           cat("\n")
         } 
  )
}

