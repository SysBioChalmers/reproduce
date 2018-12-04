#!/usr/bin/env Rscript
# Functions for processing different datasets used in this study
# Benjamin Sanchez


## @knitr splitIBAQdata
# ES data is inside the iBAQ data (among other measurements)
ESdata <- iBAQdata[grep('ups', iBAQdata$Protein.IDs),]
ISdata <- iBAQdata[grep('ups', iBAQdata$Protein.IDs, invert = TRUE),]
# Leave only the proteins from UPS2 (ends in "ups")
for(i in 1:length(ESdata$Protein.IDs)) {
  group <- strsplit(as.character(ESdata$Protein.IDs[i]),';')[[1]]
  levels(ESdata$Protein.IDs) <- c(levels(ESdata$Protein.IDs), group[grep('ups',group)])
  ESdata$Protein.IDs[i]      <- group[grep('ups',group)]
}
# Merge the ES data with the known abundances (fmol/sample) from UPS2:
ESdata <- merge(UPS2, ESdata, by = 'Protein.IDs', all.x = FALSE, all.y = FALSE)
# Take out zero values:
ESdata[ESdata == 0] <- NA


## @knitr addNtheoPeptides
addNtheoPeptides <- function(data,NTPdata,fraction) {
  intPos  <- names(NTPdata) == paste0('Intensity',fraction)
  iBAQpos <- names(NTPdata) == paste0('iBAQ',fraction)
  NTP     <- NTPdata[,intPos]/NTPdata[,iBAQpos]
  NTP     <- data.frame(Protein.IDs = NTPdata$Protein.IDs, theo.peptides = NTP)
  data    <- merge(NTP, data, by = 'Protein.IDs', all.x = FALSE, all.y = TRUE)
  return(data)
}


## @knitr getSampleAbundance
getSampleAbundance <- function(SILACdata,ISdata,method) {
  # Merge abundance data from IS into SILAC dataset:
  ISpattern  <- paste0('Abundance.',method,'.IS.')
  abundances <- ISdata[,c(1,grep(ISpattern,names(ISdata)))]
  SILACdata  <- merge(SILACdata,abundances, by = 'Protein.IDs',
                      all.x = TRUE, all.y = FALSE)
  # Compute sample abundances:
  LHNratios <- grep('Ratio.H.L.normalized.R',names(SILACdata))
  for(LHNratio in LHNratios) {
    # Define name of relevant variables:
    LHNratio_name <- names(SILACdata)[LHNratio]
    root_name     <- gsub('^.*?_','',LHNratio_name)
    IS_name       <- paste0(ISpattern,root_name)
    # Compute abundance for sample and rescale:
    abundance <- SILACdata[,LHNratio]^-1        #Ratios are stored as H/L
    abundance <- abundance*SILACdata[[IS_name]] #(L/H)*abundance [fmol/sample]
    # Add abundances to dataset:
    new_name <- gsub('Ratio.H.L.normalized','',LHNratio_name)
    new_name <- paste0('Abundance.',method,new_name)  #name for abundance of sample
    SILACdata[[new_name]] <- abundance
  }
  return(SILACdata)
}


## @knitr rescaleData
rescaleData <- function(data,pattern,name,totProt) {
  pos <- grep(pattern,names(data))   #Values in the data to rescale
  for(i in pos) {
    # Compute new abundance:
    abundance <- data[,i]*data$Mol..weight..kDa.
    abundance <- abundance/sum(abundance, na.rm = TRUE) #g/g in sample
    abundance <- abundance*totProt                      #pg in sample
    abundance <- abundance/data$Mol..weight..kDa.       #fmol in sample
    # Add abundances to dataset:
    new_name <- gsub(pattern,paste0('Abundance.',name),names(data)[i])
    data[[new_name]] <- abundance
  }
  return(data)
}


## @knitr normalizeIntensities
normalizeIntensities <- function(data,varName) {
  # Create a copy of data and normalize intensities:
  Ndata <- data
  MSpos <- grep('Intensity',names(Ndata))
  Ndata[,MSpos] <- Ndata[,MSpos]/Ndata[[varName]]
  #Change variable names:
  varName <- gsub('Sequence.length','length',varName)
  varName <- gsub('theo.peptides','Ntheo',varName)
  names(Ndata) <- gsub('Intensity',paste0('normInt.',varName),names(Ndata))
  #Merge back:
  data <- merge(data, Ndata)
  return(data)
}


## @knitr interpolateAbundance
interpolateAbundance <- function(data,pattern,abundance_pattern) {
  pos <- grep(pattern,names(data))   #All values to interpolate
  for(i in pos) {
    # Define name of relevant variables:
    name_i <- names(data)[i]
    Lname  <- gsub('.H.','.L.',name_i)  #name of ES intensity
    # Build linear model and apply to get abundance of H:
    UPS2abundances <- log10(ESdata$amount.fmoles)
    UPS2values     <- log10(ESdata[[Lname]])
    lmodel         <- lm(UPS2abundances - UPS2values ~ 1) #ES curve with fixed slope = 1
    intercept      <- lmodel[1]$coefficients[1]           #intercept
    abundance      <- 10^(log10(data[,i]) + intercept)    #LT for log(data) [fmol/sample]
    # Add abundances to dataset:
    name_i         <- gsub(pattern,paste0('Abundance.',abundance_pattern),name_i)
    data[[name_i]] <- abundance
  }
  return(data)
}


## @knitr getRPdata
getRPdata <- function(data,RP) {
  data <- merge(RP,data, by = 'Protein.IDs', all.x = FALSE, all.y = FALSE)
  data$Protein.IDs <- NULL
  # Remove ending "A" or "B" of paralogs:
  data$wikigene_name <- gsub('[AB]$','',data$wikigene_name)
  # Add up paralogs:
  data <- ddply(data, .(wikigene_name), numcolwise(sum))
  # Sort dataframe:
  letters <- substr(data$wikigene_name, 0, 3)
  numbers <- as.numeric(substr(data$wikigene_name, 4, nchar(data$wikigene_name)))
  data    <- data[order(letters,numbers),]
}


## @knitr getReplicateData
getReplicateData <- function(data,groupNames,option,repeatData=TRUE){
  # Erase distinction from name:
  for(i in 1:length(groupNames)) {
    names(data) <- gsub(groupNames[i],'',names(data))
  }
  data1 <- NULL
  data2 <- NULL
  # Define if data will or will not be repeated:
  if(repeatData) { iend <- length(names(data))   }
  else           { iend <- length(names(data))-1 }
  # Go through the dataset and compute all combinations abundance-FC:
  for(i in 1:iend) {
    if(repeatData) { j1 <- 1   }
    else           { j1 <- i+1 }
    for(j in j1:length(names(data))) {
      if((names(data)[i] == names(data)[j]) && (i!=j)) {
        data1 <- c(data1,data[,i])
        data2 <- c(data2,data[,j])
      }
    }
  }
  if(option == 1) {
    data <- cbind(data1,data2)
  } else if(option == 2) {
    # Remove NA, NaN or Inf:
    x  <- log10(data1)
    FC <- abs(log10(data2/data1))
    x[is.infinite(x)]   <- NA
    FC[is.infinite(FC)] <- NA
    FC   <- FC[!is.na(x)]
    x    <- x[!is.na(x)]
    x    <- x[!is.na(FC)]
    FC   <- FC[!is.na(FC)]
    data <- cbind(x,FC)
  }
  return(data)
}


## @knitr FCbreakdownRep
FCbreakdownRep <- function(data,groupNames) {
  for(i in 1:length(groupNames)) {
    names(data) <- gsub(groupNames[i],'',names(data))
  }
  FCs <- NULL
  for(i in 1:(length(names(data))-1)) {
    for(j in (i+1):length(names(data))) {
      if(names(data)[i] == names(data)[j]) {
        FCs <- cbind(FCs,10^abs(log10(data[,i]/data[,j])))
      }
    }
  }
  FCs[is.infinite(FCs)] <- NA
  maxFC <- apply(FCs, 1, max, na.rm = TRUE)
  x     <- NULL
  x[1]  <- sum(maxFC < 2, na.rm = TRUE)/length(maxFC)
  x[2]  <- sum((maxFC >= 2)*(maxFC <= 10), na.rm = TRUE)/length(maxFC)
  x[3]  <- sum(maxFC > 10, na.rm = TRUE)/length(maxFC)
  x     <- paste0(round(x*100,digits = 1),"%")
  return(x)
}


## @knitr FCbreakdown
FCbreakdown <- function(data) {
  x      <- NULL
  x[1:3] <- FCbreakdownRep(data,c('.R1.1','.R2.1','.R3.1'))
  x[4:6] <- FCbreakdownRep(data,c('_batch1','_batch2','_batch3'))
  x[7:9] <- FCbreakdownRep(data,c('.R1.1','.R2.1','.R3.1','_batch1','_batch2','_batch3'))
  return(x)
}


## @knitr getCVs
getCVs <- function(data,groupNames,method){
  # Erase distinction from name:
  for(i in 1:length(groupNames)) {
    names(data) <- gsub(groupNames[i],'',names(data))
  }
  #Compute all possible CVs (within groups and within proteins):
  uniqueGroups <- unique(names(data))
  CVs <- NULL
  for(group in uniqueGroups) {
    groupPos  <- grep(group,names(data))
    for(i in 1:length(data[,1])) {
      data_i <- as.numeric(data[i,groupPos])
      if(sum(is.na(data_i)) == 0) {
        if(sum(data_i) > 0) {
          mean_i <- mean(data_i)
          std_i  <- sd(data_i)
          CV_i   <- log10(std_i/mean_i*100)
          CVs    <- c(CVs,CV_i)
        }
      }
    }
  }
  return(CVs)
}

