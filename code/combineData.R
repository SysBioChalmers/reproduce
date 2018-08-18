#!/usr/bin/env Rscript
# Functions for combining different datasets used in this study
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


## @knitr normalizeIntensities
normalizeIntensities <- function(data) {
  MSpos <- grep('Intensity',names(data))
  data[,MSpos] <- data[,MSpos]/data$Sequence.length
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
    abundance      <- 10^(log10(data[,i]) + intercept)    #L.T. for log(data) [fmol/sample]
    # Add abundances to dataset:
    name_i         <- gsub(pattern,paste0('Abundance.',abundance_pattern),name_i)
    data[[name_i]] <- abundance
  }
  return(data)
}


## @knitr getSamplesAbundance
getSampleAbundance <- function(SILACdata,ISdata,pattern) {
  # Merge abundance data from IS into SILAC dataset:
  pattern    <- paste0('Abundance.',pattern)
  abundances <- ISdata[,c(1,grep(paste0(pattern,'.IS.'),names(ISdata)))]
  SILACdata  <- merge(SILACdata,abundances, by = 'Protein.IDs', all.x = TRUE, all.y = FALSE)
  # Compute sample abundances:
  LHNratios <- grep('Ratio.H.L.normalized.R',names(SILACdata))
  for(LHNratio in LHNratios) {
    # Define name of relevant variables:
    LHNratio_name <- names(SILACdata)[LHNratio]
    root_name     <- tolower(gsub('^.*?_','',LHNratio_name))
    IS_name       <- paste0(pattern,'.IS.',root_name)
    # Compute abundance for sample and rescale:
    abundance <- SILACdata[,LHNratio]^-1        #Ratios are stored as H/L
    abundance <- abundance*SILACdata[[IS_name]] #(L/H)*abundance [fmol/sample]
    # Add abundances to dataset:
    new_name <- gsub('Ratio.H.L.normalized','',LHNratio_name)
    new_name <- paste0(pattern,new_name)            #name for abundance of sample
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


## @knitr getFCvsAbundance
getFCvsAbundance <- function(abundanceData){
  x  <- NULL
  FC <- NULL
  # Go through the dataset and compute all combinations abundance-FC:
  for(i in 1:length(names(abundanceData))) {
    for(j in 1:length(names(abundanceData))) {
      if(i!=j) {
        x  <- c(x,log10(abundanceData[,i]))
        FC <- c(FC,abs(log10(abundanceData[,j]/abundanceData[,i])))
      }
    }
  }
  # Remove NA, NaN or Inf:
  x[is.infinite(x)]   <- NA
  FC[is.infinite(FC)] <- NA
  FC   <- FC[!is.na(x)]
  x    <- x[!is.na(x)]
  x    <- x[!is.na(FC)]
  FC   <- FC[!is.na(FC)]
  data <- cbind(x,FC)
  return(data)
}

