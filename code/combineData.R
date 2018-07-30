#!/usr/bin/env Rscript
# Functions for combining different datasets used in this study
# Benjamin Sanchez


## @knitr getESdata
# ES data is inside the iBAQ data (among other measurements)
ESdata   <- iBAQdata[grep('ups', iBAQdata$Protein.IDs),]
iBAQdata <- iBAQdata[grep('ups', iBAQdata$Protein.IDs, invert = TRUE),]
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
# Get mass abundances (pg/sample) from molar abundance (fmol/sample) and molecular weight (kDA = g/mmol)
ESdata$amount.pg <- ESdata$amount.fmoles*ESdata$Mol..weight..kDa.


## @knitr getISabundance
getISabundance <- function(iBAQdata,pattern) {
  Hposs <- grep(pattern,names(iBAQdata))   #All 6 values from the IS (H fraction)
  for(Hpos in Hposs) {
    # Define name of relevant variables:
    Hname <- names(iBAQdata)[Hpos]
    Lname <- gsub('.H.','.L.',Hname)  #name of ES intensity
    # Build linear model and apply to iBAQ(H) to get abundance of H:
    UPS2abundances <- log10(ESdata$amount.pg)                 
    UPS2values     <- log10(ESdata[[Lname]])
    lmodel         <- lm(UPS2abundances - UPS2values ~ 1) #ES curve with fixed slope = 1
    intercept      <- lmodel[1]$coefficients[1]           #intercept
    iBAQH          <- iBAQdata[Hpos]                      #iBAQ(H)
    ISabundance    <- 10^(log10(iBAQH) + intercept)       #L.T. for log(iBAQ(H)) [pg in sample]
    # Add abundances to dataset:
    new_name <- gsub('^.*?_','',Hname)
    if(pattern == 'iBAQ.H.T4h') { new_name <- paste0('AbundanceIS.',new_name) }
    else if(pattern == 'Intensity.H.T4h') { new_name <- paste0('Abundance2IS.',new_name) }
    iBAQdata[[new_name]] <- ISabundance
  }
  return(iBAQdata)
}


## @knitr getSamplesAbundance
getSampleAbundance <- function(SILACdata,iBAQdata,pattern) {
  # Merge abundance data from iBAQ into SILAC dataset:
  abundances <- iBAQdata[,c(1,grep(paste0(pattern,'IS.'),names(iBAQdata)))]
  SILACdata  <- merge(SILACdata,abundances, by = 'Protein.IDs', all.x = FALSE, all.y = FALSE)
  # Compute sample abundances:
  LHNratios <- grep('Ratio.H.L.normalized.R',names(SILACdata))
  for(LHNratio in LHNratios) {
    # Define name of relevant variables:
    LHNratio_name <- names(SILACdata)[LHNratio]
    root_name     <- tolower(gsub('^.*?_','',LHNratio_name))
    IS_name       <- paste0(pattern,'IS.',root_name)
    # Compute abundance for sample and rescale:
    abundance <- SILACdata[,LHNratio]^-1        #Ratios are stored as H/L
    abundance <- abundance*SILACdata[[IS_name]] #(L/H)*abundance [pg in iBAQ sample]
    abundance <- abundance/12*15                #Rescale to new size [pg in SILAC sample]
    # Add abundances to dataset:
    new_name <- gsub('Ratio.H.L.normalized','',LHNratio_name)
    new_name <- paste0(pattern,new_name)            #name for abundance of sample
    SILACdata[[new_name]] <- abundance
  }
  return(SILACdata)
}


## @knitr rescaleIS
rescaleIS <- function(data,pattern,mean_totProt) {
  Hposs <- grep(pattern,names(data))   #Values in the IS (H fraction)
  for(Hpos in Hposs) {
    # Compute new abundance:
    abundance <- data[,Hpos]
    abundance <- abundance/sum(abundance, na.rm = TRUE)   #g/g in sample
    abundance <- abundance*mean_totProt                   #pg in sample
    # Add abundances to dataset:
    if(pattern == 'iBAQ.H.T4h') {
      new_name <- gsub('^.*?_','AbundanceRescaledIS.',names(data)[Hpos])
    } else if(pattern == 'Intensity.H.T4h') {
      new_name <- gsub('^.*?_','Abundance2RescaledIS.',names(data)[Hpos])
    }
    data[[new_name]] <- abundance
  }
  return(data)
}

