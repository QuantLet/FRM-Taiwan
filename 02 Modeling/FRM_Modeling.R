################################################################################

setwd("/Users/johnnyjheng/Documents/研究所/論文/20240619 David JHE TEN David Ten Hae FRM/David/20240608 FRM-Taiwan/03 Modeling")
library(quantreg)
source("FRM_qrL1.r")
source("quantilelasso.r")

################################################################################


# variables

companies <- c('Cap', 'Fin', 'Ele', 'FinEle')
tau = 0.05                    # quantile level
ws  = 63                      # window size

for (company in companies) {

  input_path <- sprintf("../01 Raw Data/%s", company)
  output_path <- sprintf("%s_%s/", company, ws)
  
  macro <- read.csv(file = '../01 Raw Data/Macro/macro.csv')
  return <- read.csv(file = sprintf('%s/Stock returns.csv', input_path))
  data <- return[return$Date %in% macro$Date, ]
  data <- data[order(data$Date), ]
  row.names(data) <- NULL
  
  Rank <- read.csv(file = sprintf('%s/Top rank company.csv', input_path))
  rank <- Rank[Rank$Date %in% macro$Date, ]
  row.names(rank) <- NULL
  
  # generate a data matrix with standardized return
  xx0 <- data[,-1]
  xx0 <- scale(xx0) 
  
  # generate a data matrix with standardized macroeconomic variables
  m <- macro[,-1]
  m <- scale(m)
  
  # generate a data matrix with standardized top rank firm list
  r <- rank[,-1]
  
  
  
  # start the linear quantile lasso estimation for each coins
  
  for (k in 1:ncol(xx0)) {
    
    cat("Coin:", k, "\n")
    
    
    y   = as.matrix(xx0[, k])     # log return of firm k
    xx1 = as.matrix(xx0[, -k])    # log returns of firms except firm k
    rr  = as.matrix(r)            # top rank of all firms
    n   = nrow(xx1)               # number of rows of log return
    p   = ncol(r) + ncol(m) - 1   # number of covariates
    
                                  # also try the ws = 24 and compare
    kk  = gsub("^X(\\d.*)$", "\\1", colnames(xx0)[k]) # firm k's stock code
    
    # lambda calculated from linear quantile lasso
    lambda_l       = matrix(0, (n - ws), 1)
    beta_l       = matrix(0, (n - ws), p)
    
    ############################### start the moving window estimation
    
    for (i in 1:(n - ws)) {
      cat("time:",i,"\n")
  
      rw  = rr[(i + ws), ]        # top rank of firms in rolling window
      if (!(kk %in% rw)) {next}   # if firm k is not top cap in rolling window, then skip
      
      yw  = y[i:(i + ws)]         # log return of firm k in rolling window
      xx  = xx1[i:(i + ws), ]     # log return of firms except firm k in rolling window
      mw  = m[i:(i + ws), ]       # log return of macro variables in rolling window
      xx  = xx1[i:(i + ws), colnames(xx1)[gsub("^X(\\d.*)$", "\\1", colnames(xx1)) %in% rw]]
      # all the independent variables
      
      # xxw: merge xx and mw
      xxw         = cbind(xx, mw)
    
      fit         = linear(yw, xxw, tau, i, k)
      
      lambda_l[i] = fit$lambda.in
      beta_l[i,] = fit$beta.in
      
      cat("lambda:",lambda_l[i],"\n")
      cat("beta:",beta_l[i,],"\n")
      
      
    }
    
    
    write.csv(lambda_l, file = sprintf('%s/lambda_%s.csv', output_path, k))
    write.csv(beta_l, file = sprintf('%s/beta_%s.csv', output_path, k))
    
  } 

  
}
