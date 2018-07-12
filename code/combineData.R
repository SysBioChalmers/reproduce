#!/usr/bin/env Rscript
# Functions for combining different datasets used in this study
# Benjamin Sanchez


## @knitr getESdata
# ES data is inside the IS data (among other measurements)
ESdata <- ISdata[grep('ups',ISdata$Protein.IDs),]
# Leave only the protein from UPS2 (ends in "ups")
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
Hpeaks <- grep('iBAQ.H.T4h',names(ISdata))   #All 6 absolute iBAQ peaks from the IS (H fraction)
for(Hpeak in Hpeaks) {
  #Define name of relevant variables:
  Hpeak_name <- names(ISdata)[Hpeak]
  Lpeak_name <- gsub('.H.','.L.',Hpeak_name)  #name of ES peak
  #Build linear model and apply to iBAQ(H) to get abundance of H:
  UPS2abundances <- log10(ESdata$amount.pg)                 
  UPS2peaks      <- log10(ESdata[[Lpeak_name]])
  lmodel         <- lm(UPS2abundances ~ UPS2peaks)  #ES curve
  coeff1         <- lmodel[1]$coefficients[1]       #slope
  coeff2         <- lmodel[1]$coefficients[2]       #intercept
  iBAQH          <- ISdata[Hpeak]                   #iBAQ(H)
  ISabundance    <- 10^(coeff1 + coeff2*log10(iBAQH))   #L.T. for log(iBAQ(H)) [pg in sample]
  #Add abundances to dataset:
  new_name <- gsub('^.*?_','',Hpeak_name)
  new_name <- paste0('Abundance.',new_name)     #name for abundance of sample
  ISdata[[new_name]] <- ISabundance
}

## @knitr skippUPS2
Hpeaks <- grep('iBAQ.H.T4h',names(ISdata))   #All 6 absolute iBAQ peaks from the IS (H fraction)
for(Hpeak in Hpeaks) {
  abundance <- ISdata[,Hpeak]
  abundance <- abundance/sum(abundance, na.rm = TRUE)   #g/g in sample
  abundance <- abundance*12*1e-6                        #g in sample
  abundance <- abundance*1e12                           #pg in sample
  #Add abundances to dataset:
  new_name <- gsub('^.*?_','',names(ISdata)[Hpeak])
  new_name <- paste0('AbundanceNoUPS2.',new_name)     #name for abundance of sample
  ISdata[[new_name]] <- abundance
}

