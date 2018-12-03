#!/usr/bin/env Rscript
# Functions for plotting the data used in this study
# Benjamin Sanchez

## @knitr getColors
getColors <- function(n) {
  colormap <- read.table('../data/raw_external/colormaps/cividis.txt', header = FALSE)
  colormap <- as.matrix(colormap)
  N        <- length(colormap[,1])
  pos      <- seq(from = 1, to = N, by = floor((N-1)/(n-1)))
  cols <- NULL
  for(i in 1:n) {
    cols[i] <- rgb(red = colormap[pos[i],1],
                   green = colormap[pos[i],2],
                   blue = colormap[pos[i],3])
  }
  return(cols)
}

## @knitr plotScatter
plotScatter <- function(data1,data2,title,labelx,labely) {
  # Remove NA values - zeros - Infs:
  no_na  <- is.na(data1) + is.na(data2) == 0
  data1  <- data1[no_na]
  data2  <- data2[no_na]
  check1 <- (data1 > 0)*(!is.infinite((data1))) == 1
  check2 <- (data2 > 0)*(!is.infinite((data2))) == 1
  data1  <- data1[check1*check2 == 1]
  data2  <- data2[check1*check2 == 1]
  # Define colors:
  FC         <- abs(log10(data2/data1))
  FCm        <- round(median(10^FC), digits = 2)
  cols       <- getColors(6)
  col_scheme <- ifelse(10^FC>2, ifelse(10^FC>10, cols[4], cols[6]), cols[2])
  # Plot data:
  data1   <- log10(data1)
  data2   <- log10(data2)
  min_val <- round(min(c(data1,data2)))-1
  max_val <- round(max(c(data1,data2)))
  plot(data1, data2, col = col_scheme, xaxt = 'n', yaxt = 'n',
       xaxs = 'i', yaxs = 'i', main = title,  xlab = '', ylab = '',
       xlim = c(min_val,max_val), ylim = c(min_val,max_val))
  show_labels <- ifelse(labelx == '', FALSE, TRUE)
  axis(side=1, at = seq(min_val, max_val, by = 1), labels = show_labels, tck = 0.015)
  axis(side=2, at = seq(min_val, max_val, by = 1), labels = show_labels, tck = 0.015)
  axis(side=3, at = seq(min_val, max_val, by = 1), labels = FALSE, tck = 0.015)
  axis(side=4, at = seq(min_val, max_val, by = 1), labels = FALSE, tck = 0.015)
  if(labelx != '') {
    title(xlab=labelx, line=2.5)
    title(ylab=labely, line=2.5)
  }
  abline(0,1,col='black')
  # Show fit to y = x:
  R2 <- round(1 - (sum((data1 - data2)^2)/sum((data1 - mean(data1))^2)),2)
  text(max_val, min_val+1.5, bquote('R'^2 ~ '=' ~ .(R2)), pos = 2)
  text(max_val, min_val+0.5, bquote('FC'['m'] ~ '=' ~ .(FCm)), pos = 2)
  return(FC)
}


## @knitr plotVariability
plotVariability <- function(data,groupNames,title,labelx='',labely='',repeatData=TRUE) {
  # Get data:
  data <- getReplicateData(data,groupNames,1,repeatData)
  # Plot all combinations:
  tmp <- plotScatter(data[,1],data[,2],title,labelx,labely)
}


## @knitr plotVsLength
plotVsLength <- function(data,varNames,titleNames) {
  for(i in 1:length(varNames)) {
    if(length(varNames) == 1) {
      x <- data$Sequence.length
      y <- data[[varNames]]
    } else {
      name   <- varNames[i]
      data_i <- data[,grep(name,names(data))]
      x <- NULL
      y <- NULL
      for(j in 1:length(names(data_i))) {
        x <- c(x,data$Sequence.length)
        y <- c(y,log10(data_i[,j]))
      }
    }
    y[is.infinite(y)] <- NA
    x <- x[!is.na(y)]
    y <- y[!is.na(y)]
    col_opt <- rgb(red = 0, green = 0, blue = 0, alpha = 0.03)
    plot(x, y, col = col_opt, main = titleNames[i], xlab = '', ylab = '')
    lmodel <- lm(y ~ x)
    a  <- lmodel$coefficients[1]
    b  <- lmodel$coefficients[2]
    R2 <- round(summary(lmodel)$r.squared,2)
    abline(a, b, col = 'red')
    text(round(max(x, na.rm = TRUE)),floor(max(y, na.rm = TRUE)),
         bquote('R'^2 ~ '=' ~ .(R2)), pos = 2)
  }
}

## @knitr plotESdata
plotESdata <- function(ESdata,method,title) {
  pattern  <- paste0('Abundance.',method,'.ES')
  dataExp  <- NULL
  dataPred <- NULL
  for(i in 1:length(names(ESdata))) {
    if(grepl(pattern,names(ESdata)[i])) {
      dataExp  <- c(dataExp,ESdata$amount.fmoles)
      dataPred <- c(dataPred,ESdata[,i])
    }
  }
  FC <- plotScatter(dataExp,dataPred,title,bquote('log'['10'] ~ '(measured)'),
                    bquote('log'['10'] ~ '(predicted)'))
  return(FC)
}


## @knitr plotRPdata
plotRPdata <- function(RPdata,title) {
  # Remove zeros, compute median and error, and then take log:
  data <- as.matrix(RPdata[,-1])
  data[data == 0] <- NA
  mead_val <- median(data, na.rm = TRUE)
  FC       <- as.vector(abs(log10(data/mead_val)))
  FC       <- FC[!is.na(FC)]
  FCm      <- round(median(10^FC),2)
  data     <- log10(data)
  mead_val <- log10(mead_val)
  # Color by tech. rep:
  col_opt <- NULL
  cols    <- getColors(6)
  for(i in 2:length(names(RPdata))) {
    if(length(grep('batch1',names(RPdata)[i])) == 1)      { col_opt[i] <- cols[1] }
    else if(length(grep('batch2',names(RPdata)[i])) == 1) { col_opt[i] <- cols[3] }
    else if(length(grep('batch3',names(RPdata)[i])) == 1) { col_opt[i] <- cols[5] }
  }
  # Plot data:
  matplot(data, pch = 1, xaxt = 'n', col = col_opt, main = title,
          ylab = bquote('log'['10'] ~ '(abundance)'))
  axis(side=1, at = 1:length(RPdata[,1]), labels = RPdata[,1], las=2, cex.axis = 0.7)
  max_x <- length(RPdata[,1])
  lines(c(1,max_x),c(mead_val,mead_val), col = 'black', lwd = 2, lty = 2)
  text(max_x, mead_val-1, bquote('FC'['m'] ~ '=' ~ .(FCm)), pos = 2)
  return(FC)
}


## @knitr plotCumulativeDistrib
plotCumulativeDistrib <- function(FCs,varName){
  # Assign names to dataframe:
  names(FCs) <- c('iBAQ','iBAQrescaled','TPA','TPAnorm')
  # Compute differences between distributions:
  htest1 <- ks.test(FCs$iBAQ, FCs$iBAQrescaled)
  htest2 <- ks.test(FCs$iBAQ, FCs$TPA)
  htest3 <- ks.test(FCs$iBAQ, FCs$TPAnorm)
  htest4 <- ks.test(FCs$iBAQrescaled, FCs$TPA)
  htest5 <- ks.test(FCs$iBAQrescaled, FCs$TPAnorm)
  htest6 <- ks.test(FCs$TPA, FCs$TPAnorm)
  print(paste(varName, '- number of FC compared =',length(FCs$iBAQ)))
  print(paste(varName, '- number of FC compared =',length(FCs$iBAQrescaled)))
  print(paste(varName, '- number of FC compared =',length(FCs$TPA)))
  print(paste(varName, '- number of FC compared =',length(FCs$TPAnorm)))
  print(paste(varName, 'of iBAQ = iBAQrescaled: p-val =',round(htest1$p.value,4)))
  print(paste(varName, 'of iBAQ = TPA: p-val =',round(htest2$p.value,4)))
  print(paste(varName, 'of iBAQ = TPAnorm: p-val =',round(htest3$p.value,4)))
  print(paste(varName, 'of iBAQrescaled = TPA: p-val =',round(htest4$p.value,4)))
  print(paste(varName, 'of iBAQrescaled = TPAnorm: p-val =',round(htest5$p.value,4)))
  print(paste(varName, 'of TPA = TPAnorm: p-val =',round(htest6$p.value,4)))
  # Plot data:
  min_x <- 0
  max_x <- 1
  if(nchar(varName) > 4) {
    title <- varName
  } else {
    title <- ''
  }
  plot(1,1, xaxs = 'i', yaxs = 'i', xaxt = 'n', yaxt = 'n', type='n',
       xlim = c(min_x, max_x), ylim = c(0, 1), main = title, xlab = '', ylab = '')
  axis(side=1, at = seq(min_x, max_x, by = 0.2), labels = TRUE,  tck = 0.015)
  axis(side=2, at = seq(0, 1, by = 0.2),         labels = TRUE,  tck = 0.015)
  axis(side=3, at = seq(min_x, max_x, by = 0.2), labels = FALSE, tck = 0.015)
  axis(side=4, at = seq(0, 1, by = 0.2),         labels = FALSE, tck = 0.015)
  title(xlab=bquote('abs(log'['10'] ~ '(FC))'), line=2.5)
  title(ylab='Cumulative Distribution', line=2.5)
  lines(c(log10(2),log10(2)),c(0,1), col = 'black', lwd = 2, lty = 2)
  # Plot fold changes as a cdf:
  N    <- length(names(FCs))
  cols <- getColors(N)
  for(i in 1:N) {
    FC   <- sort(FCs[[i]])
    step <- 1/(length(FC)-1)
    cdf  <- seq(0, 1, by = step)
    if(startsWith(varName,'Tech') && i == 4) {
      lines(FC,cdf, col = cols[i], lwd = 2, lty = 2)
    } else {
      lines(FC,cdf, col = cols[i], lwd = 2)
    }
  }
  # Plot values for 2-fold position:
  for(i in 1:N) {
    FC  <- sort(FCs[[i]])
    pos <- which.min(abs(FC - log10(2)))
    points(log10(2),cdf[pos[1]], pch = 21, col = 'black', bg = cols[i])
  }
}


## @knitr plotTotalProt
plotTotalProt <- function(data,pattern,titleName) {
  pos <- grep(pattern,names(data))
  # Display number of proteins detected:
  meanVals <- rowMeans(data[,pos], na.rm = TRUE)
  coverage <- sum(!is.na(meanVals))
  print(paste(titleName,'-',coverage,'proteins detected'))
  # Get total protein detected:
  data[,pos]     <- data[,pos]*data$Mol..weight..kDa.      #fmol/sample*kDa = pg/sample
  totProt        <- colSums(data[,pos], na.rm = TRUE)/1e6  #ug/sample
  meanProt       <- round(mean(totProt))+1
  names(totProt) <- gsub(pattern,'',names(totProt))
  factors <- factor(names(totProt))
  cols    <- getColors(nlevels(factors))
  barplot(totProt, col = cols[factors], mgp = c(1.5, 0.5, 0), cex.names = 0.8,
          ylab = 'Total detected protein [ug]', xaxt = 'n', ylim = c(0,meanProt))
  title(main = titleName, cex.main = 0.9)
}


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
  cols    <- getColors(3)
  col_opt <- ifelse(scaling == 0,cols[1],cols[2])
  abline(a = intercept, b = 1, col = col_opt)
  #Compute and display R2:
  yp <- x + intercept
  R2 <- round(1 - (sum((y - yp)^2)/sum((y - mean(y))^2)),2)
  if(!allInOne) {
    text(round(min(x, na.rm = TRUE))-1,max(y, na.rm = TRUE)-pos,
         bquote('R'^2 ~ '=' ~ .(R2)), pos = 4, col = col_opt)
  }
}


## @knitr plotES
plotES <- function(ESdata,pattern,scaling,name,allInOne,first,CVm) {
  #Plot data:
  x <- log10(ESdata[[paste0(pattern,name)]])
  y <- log10(ESdata$amount.fmoles)
  s <- rep(0:18,length.out = length(x))
  y <- y[!is.na(x)]
  s <- s[!is.na(x)]
  x <- x[!is.na(x)]
  min_x <- round(min(x, na.rm = TRUE))-1
  max_x <- ceiling(max(x, na.rm = TRUE))
  min_y <- floor(min(y, na.rm = TRUE))
  max_y <- ceiling(max(y, na.rm = TRUE))
  cols  <- getColors(3)
  if(length(scaling) > 1) { scaling <- scaling[grep(name,names(scaling))] }
  if(allInOne) {
    interSize <- 0.5
    if(first) {
      plot(x,y, col = cols[3], xaxs = 'i', yaxs = 'i', xaxt = 'n', yaxt = 'n',
           xlim = c(min_x, max_x), ylim = c(min_y, max_y), asp = 1,
           xlab = bquote('log'['10'] ~ '(intensity value)'),
           ylab = bquote('log'['10'] ~ '(abundance)'), mgp = c(2, 0.5, 0))
      text(min_x, max_y-0.5, bquote('CV'['m'] ~ '=' ~ .(CVm) ~ '%'), pos = 4)
    } else {
      points(x,y, col = cols[3])
    }
  } else {
    interSize <- 0.3
    plot(x,y, pch = s,col = cols[3], xaxs = 'i', yaxs = 'i', xaxt = 'n', yaxt = 'n',
         asp = 1, xlim = c(min_x, max_x), ylim = c(min_y, max_y), main = name)
  }
  if(!allInOne || first) {
    axis(side=1, at = seq(min_x, max_x, by = 1), labels = min_x:max_x,
         tck = 0.015, mgp = c(2, interSize*2/3, 0))
    axis(side=2, at = seq(min_y, max_y, by = 1), labels = min_y:max_y,
         tck = 0.015, mgp = c(2, interSize, 0))
    axis(side=3, at = seq(min_x, max_x, by = 1), labels = FALSE, tck = 0.015)
    axis(side=4, at = seq(min_y, max_y, by = 1), labels = FALSE, tck = 0.015)
  }
  #Get linear fits:
  plotLM(x,y,0,allInOne)
  plotLM(x,y,scaling,allInOne)
}


## @knitr plotAllES
plotAllES <- function(ESdata,pattern,scaling,allInOne) {
  #Compute the mean coefficient of variation:
  data <- ESdata[,grep(pattern,names(ESdata))]
  SD   <- apply(data, 1, sd, na.rm = TRUE)    #Standard deviation for each protein
  mu   <- apply(data, 1, mean, na.rm = TRUE)  #Mean for each protein
  CV   <- SD/mu                               #Coefficient of variation for each protein
  CVm  <- mean(CV, na.rm = TRUE)*100          #Mean coefficient of variation [%]
  CVm  <- round(CVm, digits = 1)
  #Plot data:
  if(allInOne) { par(mfrow = c(1,1), mar = c(3,3,0,0), pty = "s", cex = 1) }
  else         { par(mfrow = c(2,3), mar = c(3, 3, 2, 1), cex = 0.5) }
  ESdata <- ESdata[order(ESdata$amount.fmoles),]
  plotES(ESdata,pattern,scaling,'top5_batch1',allInOne,TRUE,CVm)
  plotES(ESdata,pattern,scaling,'top5_batch2',allInOne,FALSE,CVm)
  plotES(ESdata,pattern,scaling,'top5_batch3',allInOne,FALSE,CVm)
  plotES(ESdata,pattern,scaling,'top10_batch1',allInOne,FALSE,CVm)
  plotES(ESdata,pattern,scaling,'top10_batch2',allInOne,FALSE,CVm)
  plotES(ESdata,pattern,scaling,'top10_batch3',allInOne,FALSE,CVm)
}


## @knitr plotPCA
plotPCA <- function(data,title,outside = FALSE){
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
  cols    <- getColors(6)
  for(i in 1:length(names(data))) {
    # Shape by bio. rep.
    if(length(grep('R1.1',names(data)[i])) == 1)      { pch_opt[i] <- 1 } 
    else if(length(grep('R2.1',names(data)[i])) == 1) { pch_opt[i] <- 2 }
    else if(length(grep('R3.1',names(data)[i])) == 1) { pch_opt[i] <- 3 }
    # Color by tech. rep.
    if(length(grep('batch1',names(data)[i])) == 1)      { col_opt[i] <- cols[1] }
    else if(length(grep('batch2',names(data)[i])) == 1) { col_opt[i] <- cols[3] }
    else if(length(grep('batch3',names(data)[i])) == 1) { col_opt[i] <- cols[5] }
  }
  # Plot PCA:
  deltax <- max(pca$x[,1]) - min(pca$x[,1])
  deltay <- max(pca$x[,2]) - min(pca$x[,2])
  xmin   <- min(pca$x[,1]) - deltax/8
  xmax   <- max(pca$x[,1]) + deltax/8
  ymin   <- min(pca$x[,2]) - deltay/8
  ymax   <- max(pca$x[,2]) + deltay/8
  plot(pca$x, col = col_opt, pch = pch_opt, cex = 1.5, xaxs = 'i', yaxs = 'i',
       xaxt = 'n', yaxt = 'n', xlab = '', main = title, ylab = '',
       xlim = c(xmin, xmax), ylim = c(ymin, ymax))
  if(outside) {
    linePos <- 1.5
    axis(side=1, at = seq(xmin, xmax, by = deltax/8), labels = FALSE, tck = 0.015)
    axis(side=2, at = seq(ymin, ymax, by = deltay/8), labels = FALSE, tck = 0.015)
    axis(side=3, at = seq(xmin, xmax, by = deltax/8), labels = FALSE, tck = 0.015)
    axis(side=4, at = seq(ymin, ymax, by = deltay/8), labels = FALSE, tck = 0.015)
  } else {
    linePos <- -1
  }
  mtext(bquote('PC'['1'] ~ ' = ' ~ .(var1) ~ '%'), side = 1, line = linePos)
  mtext(bquote('PC'['2'] ~ ' = ' ~ .(var2) ~ '%'), side = 2, line = linePos)
}


## @knitr plotAllVariability
plotAllVariability <- function(abundance,showTitle) {
  if(showTitle) {
    titleNames = c('Biological Variability','Technical Variability','PCA')
  } else {
    titleNames = c('','','')
  }
  plotVariability(abundance[,-1],c('.R1.1','.R2.1','.R3.1'),titleNames[1])
  plotVariability(abundance[,-1],c('_batch1','_batch2','_batch3'),titleNames[2])
  plotPCA(abundance[,-1],titleNames[3])
}


## @knitr plotFCvsAbundance
plotFCvsAbundance <- function(sampleData,ESdata,pattern,titleName,allInOne){
  # Options depending on type of plot:
  cols <- getColors(4)
  rgbs <- col2rgb(cols)/255
  if(grepl('iBAQ.R',pattern)) {
    color_data <- rgb(red = rgbs[1,1], green = rgbs[2,1], blue = rgbs[3,1], alpha = 0.03)
  } else if(grepl('iBAQrescaled.R',pattern)) {
    color_data <- rgb(red = rgbs[1,2], green = rgbs[2,2], blue = rgbs[3,2], alpha = 0.03)
  } else if(grepl('TPA.R',pattern)) {
    color_data <- rgb(red = rgbs[1,3], green = rgbs[2,3], blue = rgbs[3,3], alpha = 0.04)
  } else if(grepl('TPAnorm.R',pattern)) {
    color_data <- rgb(red = rgbs[1,4], green = rgbs[2,4], blue = rgbs[3,4], alpha = 0.05)
  }
  if(allInOne) {
    color_data <- rgb(red = 0, green = 0, blue = 0, alpha = 0.02)
  }
  # Get FC values for sample data:
  abundance  <- sampleData[,grep(pattern,names(sampleData))]
  sampleData <- getReplicateData(abundance,c('_batch1','_batch2','_batch3'),2)
  # Get FC values for ES data:
  ESdata_pattern <- gsub('.R..1_','.ES',pattern)
  ESabundance    <- ESdata[,grep(ESdata_pattern,names(ESdata))]
  ESdata         <- getReplicateData(ESabundance,c('_batch1','_batch2','_batch3'),2)
  # Print number of yeast proteins below minimum detected UPS2 protein
  min_es <- min(ESdata[,1], na.rm = TRUE)
  max_es <- max(ESdata[,1], na.rm = TRUE)
  if(!allInOne) {
    belowThreshold <- sum(sampleData[,1] < min_es, na.rm = TRUE)/length(names(abundance))
    print(paste(titleName,'->',round(belowThreshold),'proteins below UPS2 detection range'))
  }
  # Plot FC of abundanceData
  min_x <- round(min(sampleData[,1])) + 1
  max_x <- round(max(sampleData[,1]))
  min_y <- 0
  max_y <- round(max(sampleData[,2])) - 1
  if(!allInOne || grepl('iBAQ.R',pattern)) {
    plot(sampleData[,1],sampleData[,2], xaxs = 'i', yaxs = 'i',
         xaxt = 'n', yaxt = 'n', col = color_data,
         xlim = c(min_x, max_x), ylim = c(min_y, max_y),
         xlab = bquote('log'['10'] ~ '(abundance [fmol/sample])'), ylab = '')
    title(ylab = bquote('abs(log'['10'] ~ '(FC))'), line=2.5)
    axis(side=1, at = seq(min_x, max_x, by = 1), labels = min_x:max_x, tck = 0.015)
    axis(side=2, at = seq(min_y, max_y, by = 1), labels = min_y:max_y, tck = 0.015)
    axis(side=3, at = seq(min_x, max_x, by = 1), labels = FALSE, tck = 0.015)
    axis(side=4, at = seq(min_y, max_y, by = 1), labels = FALSE, tck = 0.015)
  } else {
    points(sampleData[,1],sampleData[,2], col = color_data)
  }
  # Plot UPS2 window:
  if(!allInOne || grepl('TPAnorm.R',pattern)) {
    polygon(c(min_es,min_es,max_es,max_es,min_es),c(min_y,max_y,max_y,min_y,min_y),
            col = rgb(red = 0, green = 0, blue = 0, alpha = 0.1), border = NA)
  }
  # Plot UPS2 points:
  if(!allInOne) {
    title(main = titleName)
    points(ESdata[,1],ESdata[,2], col = 'black')
    # Compute fraction in window:
    in_window <- (sampleData[,1] > min_es)*(sampleData[,1] < max_es)
    fraction  <- sum(in_window)/length(in_window)*100
    fraction  <- round(fraction, digits = 1)
    text(min_es-1, max_y-0.5, bquote('Fraction in window =' ~ .(fraction) ~ '%'), pos = 4)
  }
  return(sampleData)
}


## @knitr plotSplines
plotSplines <- function(SILACdata,ESdata){
  #Plot all data:
  sample1 <- plotFCvsAbundance(SILACdata,ESdata,'Abundance.iBAQ.R..1_','',TRUE)
  sample2 <- plotFCvsAbundance(SILACdata,ESdata,'Abundance.iBAQrescaled.R..1_','',TRUE)
  sample3 <- plotFCvsAbundance(SILACdata,ESdata,'Abundance.TPA.R..1_','',TRUE)
  sample4 <- plotFCvsAbundance(SILACdata,ESdata,'Abundance.TPAnorm.R..1_','',TRUE)
  # Create smoothing splines:
  ss1 <- smooth.spline(sample1[,1], sample1[,2], df = 10)
  ss2 <- smooth.spline(sample2[,1], sample2[,2], df = 10)
  ss3 <- smooth.spline(sample3[,1], sample3[,2], df = 10)
  ss4 <- smooth.spline(sample4[,1], sample4[,2], df = 10)
  cols <- getColors(4)
  lines(ss1, lwd = 2, col = cols[1])
  lines(ss2, lwd = 2, col = cols[2])
  lines(ss3, lwd = 2, col = cols[3])
  lines(ss4, lwd = 2, col = cols[4])
}

