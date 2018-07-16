#!/usr/bin/env Rscript
# Functions for plotting the data used in this study
# Benjamin Sanchez


## @knitr plotLM
# Linear model plotting function
LMplot <- function(x,y,scaling,allInOne) {
  #Make linear model with fixed slope = 1
  lmodel    <- lm(y - x ~ 1)
  intercept <- lmodel$coefficients[1]
  if(scaling == 0) {
    pos <- 0
  } else {
    pos <- 1
    intercept <- log10(scaling)  #Also fix intercept
  }
  col_opt <- ifelse(scaling == 0,'red','green')
  abline(a = intercept, b = 1, col = col_opt)
  #Compute and display R2 adjusted:
  n       <- length(x)
  yp      <- x + intercept
  R2      <- 1 - (sum((y - yp)^2)/sum((y - mean(y))^2))
  R2adj   <- round(1-((1-R2)*(n-1)/(n-1-1)),2)
  if(!allInOne) {
    text(min(x, na.rm = TRUE),max(y, na.rm = TRUE)-pos-0.5, bquote('R'^2 ~ '=' ~ .(R2adj)), pos = 4, col = col_opt)
  }
}

## @knitr plotES
# ES plotting function:
ESplot <- function(ESdata,scaling,name,allInOne,first) {
  x    <- log10(ESdata[[paste0('iBAQ.L.T4h_',tolower(name))]])
  y    <- log10(ESdata$amount.pg)
  y    <- y[!is.na(x)]
  x    <- x[!is.na(x)]
  if(length(scaling) > 1) {
    scaling <- scaling[grep(name,names(scaling))]
  }
  if(allInOne) {
    if(first) {
      plot(x,y, col = 'blue', xaxt = 'n', yaxt = 'n')
    } else {
      points(x,y, col = 'blue', xaxt = 'n', yaxt = 'n')
    }
  } else {
    plot(x,y, col = 'blue', xaxt = 'n', yaxt = 'n', main = name)
  }
  LMplot(x,y,0,allInOne)
  LMplot(x,y,scaling,allInOne)
}
# plot all ES plots
ESplots <- function(ESdata,scaling,allInOne) {
  if(allInOne) {
    num_plots <- c(1,1)
  } else {
    num_plots <- c(2,3)
  }
  par(mfrow = num_plots, mar = c(0, 0, 1, 0) + 0.5, cex = 1)
  ESplot(ESdata,scaling,'top5_batch1',allInOne,TRUE)
  ESplot(ESdata,scaling,'top5_batch2',allInOne,FALSE)
  ESplot(ESdata,scaling,'top5_batch3',allInOne,FALSE)
  ESplot(ESdata,scaling,'top10_batch1',allInOne,FALSE)
  ESplot(ESdata,scaling,'top10_batch2',allInOne,FALSE)
  ESplot(ESdata,scaling,'top10_batch3',allInOne,FALSE)
}


## @knitr plotTotalProt
protPlot <- function(data,pattern) {
  abundances <- data[,grep(pattern,names(data))]
  names(abundances) <- gsub(pattern,'',names(abundances))
  totProt <- colSums(abundances, na.rm = TRUE)/1e6  #ug in sample
  par(mfcol = c(1,1), mar = c(2.5,2.5,1,1), cex = 1)
  show_label <- ifelse(length(totProt) > 6,'n','s')
  barplot(totProt, col = factor(names(abundances)), ylab = 'Total protein content in sample [ug]', 
          cex.names = 0.8, mgp = c(1.5, 0.5, 0), xaxt = show_label)
}

