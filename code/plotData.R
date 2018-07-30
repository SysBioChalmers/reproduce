#!/usr/bin/env Rscript
# Functions for plotting the data used in this study
# Benjamin Sanchez


## @knitr plotLM
plotLM <- function(x,y,scaling,allInOne) {
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
  #Compute and display adjusted R2:
  n       <- length(x)
  yp      <- x + intercept
  R2      <- 1 - (sum((y - yp)^2)/sum((y - mean(y))^2))
  R2adj   <- round(1-((1-R2)*(n-1)/(n-1-1)),2)
  if(!allInOne) {
    text(min(x, na.rm = TRUE),max(y, na.rm = TRUE)-pos-0.5,
         bquote('R'^2 ~ '=' ~ .(R2adj)), pos = 4, col = col_opt)
  }
}


## @knitr plotES
plotES <- function(ESdata,vars,scaling,name,allInOne,first,CVm) {
  #Plot data:
  x    <- log10(ESdata[[paste0(vars,tolower(name))]])
  y    <- log10(ESdata$amount.pg)
  y    <- y[!is.na(x)]
  x    <- x[!is.na(x)]
  if(length(scaling) > 1) { scaling <- scaling[grep(name,names(scaling))] }
  if(allInOne) {
    if(first) {
      min_x <- floor(min(x, na.rm = TRUE))
      max_x <- ceiling(max(x, na.rm = TRUE))
      min_y <- floor(min(y, na.rm = TRUE))
      max_y <- ceiling(max(y, na.rm = TRUE))
      x_lab <- paste0('log10(',gsub('.L.T4h_','',vars),' values)')
      plot(x,y, col = 'blue', xaxs = 'i', yaxs = 'i', xaxt = 'n', yaxt = 'n',
           xlim = c(min_x, max_x), ylim = c(min_y, max_y), asp = 1,
           xlab = x_lab, ylab = 'log10(abundance)')
      text(min_x, max_y-0.5, bquote('CV'['m'] ~ '=' ~ .(CVm) ~ '%'), pos = 4)
      axis(side=1, at = seq(min_x, max_x, by = 1), labels = min_x:max_x, tck = 0.015)
      axis(side=2, at = seq(min_y, max_y, by = 1), labels = min_y:max_y, tck = 0.015)
      axis(side=3, at = seq(min_x, max_x, by = 1), labels = min_x:max_x, tck = 0.015)
      axis(side=4, at = seq(min_y, max_y, by = 1), labels = min_y:max_y, tck = 0.015)
    } else {
      points(x,y, col = 'blue')
    }
  } else {
    plot(x,y, col = 'blue', xaxt = 'n', yaxt = 'n', main = name)
  }
  #Get linear fits:
  plotLM(x,y,0,allInOne)
  plotLM(x,y,scaling,allInOne)
}


## @knitr plotAllES
plotAllES <- function(ESdata,vars,scaling,allInOne) {
  #Compute the mean coefficient of variation:
  data <- ESdata[,grep(vars,names(ESdata))]
  data[data == 0] <- NA
  data <- log10(data)
  SD   <- apply(data, 1, sd, na.rm = TRUE)    #Standard deviation for each protein
  mu   <- apply(data, 1, mean, na.rm = TRUE)  #Mean for each protein
  CV   <- SD/mu                               #Coefficient of variation for each protein
  CVm  <- mean(CV, na.rm = TRUE)*100          #Mean coefficient of variation (as percentage)
  CVm  <- round(CVm, digits = 1)
  #Plot data:
  if(allInOne) { par(mfrow = c(1,1), mar = c(4,4,1,1), pty = "s", cex = 1) }
  else         { par(mfrow = c(2,3), mar = c(0, 0, 1, 0) + 0.5, cex = 1) }
  plotES(ESdata,vars,scaling,'top5_batch1',allInOne,TRUE,CVm)
  plotES(ESdata,vars,scaling,'top5_batch2',allInOne,FALSE,CVm)
  plotES(ESdata,vars,scaling,'top5_batch3',allInOne,FALSE,CVm)
  plotES(ESdata,vars,scaling,'top10_batch1',allInOne,FALSE,CVm)
  plotES(ESdata,vars,scaling,'top10_batch2',allInOne,FALSE,CVm)
  plotES(ESdata,vars,scaling,'top10_batch3',allInOne,FALSE,CVm)
}


## @knitr plotTotalProt
plotTotalProt <- function(data,pattern) {
  abundances <- data[,grep(pattern,names(data))]
  names(abundances) <- gsub(pattern,'',names(abundances))
  totProt <- colSums(abundances, na.rm = TRUE)/1e6  #ug in sample
  if(length(totProt) > 6) {
    show_label <- 'n'
    margins    <- c(1,2.5,1,1)
    } else {
    show_label <- 's'
    margins    <- c(2.5,2.5,1,1)
    }
  par(mfcol = c(1,1), mar = margins, cex = 1)
  barplot(totProt, col = factor(names(abundances)), cex.names = 0.8, mgp = c(1.5, 0.5, 0),
          ylab = 'Total detected protein in sample [ug]', xaxt = show_label)
  return(mean(totProt)*1e6) #average total pg in sample
}


## @knitr plotScatter
plotScatter <- function(data1,data2,title) {
  # Remove NA values - zeros - Infs:
  no_na  <- is.na(data1) + is.na(data2) == 0
  data1  <- data1[no_na]
  data2  <- data2[no_na]
  check1 <- (data1 > 0)*(!is.infinite((data1))) == 1
  check2 <- (data2 > 0)*(!is.infinite((data2))) == 1
  data1  <- data1[check1*check2 == 1]
  data2  <- data2[check1*check2 == 1]
  # Define colors:
  FC         <- data2/data1
  FCm        <- round(median(10^abs(log10(FC))), digits = 2)
  col_scheme <- ifelse(abs(log10(FC))>log10(2),
                       ifelse(abs(log10(FC))>1, '#EB2426', '#F3B70D'), '#0F8141')
  # Plot data:
  data1   <- log10(data1)
  data2   <- log10(data2)
  min_val <- floor(min(c(data1,data2)))
  max_val <- ceiling(max(c(data1,data2)))
  plot(data1, data2, pch = 1, col = col_scheme, xlab = '', ylab = '',
       xaxs = 'i', yaxs = 'i', xaxt = 'n', yaxt = 'n', main = title,
       xlim = c(min_val,max_val), ylim = c(min_val,max_val))
  axis(side=1, at = seq(min_val, max_val, by = 1), labels = FALSE, tck = 0.015)
  axis(side=2, at = seq(min_val, max_val, by = 1), labels = FALSE, tck = 0.015)
  axis(side=3, at = seq(min_val, max_val, by = 1), labels = FALSE, tck = 0.015)
  axis(side=4, at = seq(min_val, max_val, by = 1), labels = FALSE, tck = 0.015)
  abline(0,1,col='black')
  # Show fit to a linear model:
  lmodel <- lm(data2 ~ data1)
  lsum   <- summary(lmodel)
  R2     <- round(lsum$adj.r.squared,2)
  text(min_val,max_val-0.5, bquote('R'^2 ~ '=' ~ .(R2)), pos = 4)
  text(min_val,max_val-1.5, bquote('FC'['m'] ~ '=' ~ .(FCm)), pos = 4)
}


## @knitr plotVariability
plotVariability <- function(data,groupNames,title) {
  # Erase distinction from name:
  for(i in 1:length(groupNames)) {
    names(data) <- gsub(groupNames[i],'',names(data))
  }
  # Create data with all possible combinations:
  data1 <- NULL
  data2 <- NULL
  for(i in 2:(length(names(data))-1)) {
    for(j in (i+1):length(names(data))) {
      if(names(data)[i] == names(data)[j]) {
        data1 <- c(data1,data[,i])
        data2 <- c(data2,data[,j])
      }
    }
  }
  # Plot all combinations:
  plotScatter(data1,data2,title)
}


## @knitr plotPCA
plotPCA <- function(data,title){
  # Take log from data:
  log_data                        <- log10(as.matrix(data))
  log_data[is.infinite(log_data)] <- NA
  log_data[is.nan(log_data)]      <- NA
  log_data                        <- log_data[!is.na(rowSums(log_data)),]
  log_data                        <- t(log_data)
  # Do PCA:
  pca  <- prcomp(log_data)
  var  <- pca$sdev/sum(pca$sdev)*100
  var1 <- round(var[1], digits = 1)
  var2 <- round(var[2], digits = 1)
  # Plotting options:
  pch_opt <- NULL
  col_opt <- NULL
  for(i in 2:length(names(data))) {
    # Shape by bio. rep.
    if(length(grep('R1.1',names(data)[i])) == 1)      { pch_opt[i] <- 1 } 
    else if(length(grep('R2.1',names(data)[i])) == 1) { pch_opt[i] <- 2 }
    else if(length(grep('R3.1',names(data)[i])) == 1) { pch_opt[i] <- 3 }
    # Color by tech. rep.
    if(length(grep('Batch1',names(data)[i])) == 1)      { col_opt[i] <- '#EB2426' }
    else if(length(grep('Batch2',names(data)[i])) == 1) { col_opt[i] <- '#3953A3' }
    else if(length(grep('Batch3',names(data)[i])) == 1) { col_opt[i] <- '#0F8141' }
  }
  # Plot PCA:
  deltax <- max(pca$x[,1]) - min(pca$x[,1])
  deltay <- max(pca$x[,2]) - min(pca$x[,2])
  xmin   <- min(pca$x[,1]) - deltax/8
  xmax   <- max(pca$x[,1]) + deltax/8
  ymin   <- min(pca$x[,2]) - deltay/8
  ymax   <- max(pca$x[,2]) + deltay/8
  plot(pca$x, col = col_opt, pch = pch_opt, cex = 1.5, xaxt = 'n', yaxt = 'n',
       xlim = c(xmin, xmax), ylim = c(ymin, ymax), main = title)
  mtext(bquote('PC1 = ' ~ .(var1) ~ '%'), side = 1, line = -1)
  mtext(bquote('PC2 = ' ~ .(var2) ~ '%'), side = 2, line = -1)
}

