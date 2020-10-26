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


## @knitr standardizePeptideData
standardizePeptideData <- function(data,fraction) {
  # The 8 hour gradient setting will not be used in this study:
  if(length(grep('T8h_',data$Experiment))>0) {
    data <- data[-grep('T8h_',data$Experiment),]
  }
  # Other format fixes:
  data$Experiment <- gsub('T4h_','',data$Experiment)
  data$Experiment <- gsub('_Batch','_batch',data$Experiment)
  data$Experiment <- gsub('-1_top','.1_top',data$Experiment)
  # Create new dataframe:
  variables <- c(paste('Intensity',fraction,sep='.'),'Charge','Retention.time')
  experiments <- unique(data$Experiment)
  varNames <- expand.grid(variables,experiments)
  varNames <- paste(varNames[,1],varNames[,2],sep='.')
  varNames <- paste(varNames,collapse=',')
  varNames <- paste('Sequence',varNames,sep=',')
  standardData <- read.csv(text=varNames)
  sequences <- unique(data$Sequence)
  standardData[1:length(sequences),1] <- sequences
  for(i in 1:length(sequences)) {
    sequence <- data$Sequence[i]
    for(variable in variables) {
      rowPos <- grep(sequence,standardData$Sequence)
      colName <- paste(variable,data$Experiment[i],sep='.')
      standardData[[colName]][rowPos] <- data[[variable]][i]
    }
  }
  return(standardData)
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
rescaleData <- function(data,pattern,name) {
  pos <- grep(pattern,names(data))   #Values in the data to rescale
  for(i in pos) {
    # Compute new abundance:
    abundance <- data[,i]*data$Mol..weight..kDa.
    abundance <- abundance/sum(abundance, na.rm = TRUE) #g/g protein
    abundance <- abundance/data$Mol..weight..kDa.       #mmol/g protein
    abundance <- abundance*1e12/1e6                     #fmol/ug protein
    # Add abundances to dataset:
    new_name <- gsub(pattern,paste0('Abundance.',name),names(data)[i])
    data[[new_name]] <- abundance
  }
  return(data)
}


## @knitr normalizeIntensities
normalizeIntensities <- function(data,pattern,varName) {
  # Find subset of intensities and normalize them:
  Ndata <- data[,grep(pattern,names(data))]/data[[varName]]
  #Change variable names:
  varName <- gsub('Sequence.length','length',varName)
  varName <- gsub('theo.peptides','Ntheo',varName)
  names(Ndata) <- gsub('Intensity',paste0('normInt.',varName),names(Ndata))
  #Combine both dataframes:
  data <- cbind(data, Ndata)
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
  for(groupName in groupNames) {
    names(data) <- gsub(groupName,'',names(data))
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


## @knitr getReplicateGroups
getReplicateGroups <- function(replicateType) {
  if(replicateType == 'Bio'){
    groupNames <- c('.R1.1','.R2.1','.R3.1')
  } else if(replicateType == 'Tech') {
    groupNames <- c('_batch1','_batch2','_batch3')
  } else if(replicateType == 'All') {
    groupNames <- c('.R1.1','.R2.1','.R3.1','_batch1','_batch2','_batch3')
  }
  return(groupNames)
}


## @knitr FCbreakdownRep
FCbreakdownRep <- function(data,replicateType) {
  for(groupName in getReplicateGroups(replicateType)) {
    names(data) <- gsub(groupName,'',names(data))
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
  x[1:3] <- FCbreakdownRep(data,'Bio')
  x[4:6] <- FCbreakdownRep(data,'Tech')
  x[7:9] <- FCbreakdownRep(data,'All')
  return(x)
}


## @knitr peptideDifferences
peptideDifferences <- function(data,unionPeptides,replicateType) {
  for(groupName in getReplicateGroups(replicateType)) {
    names(data) <- gsub(groupName,'',names(data))
  }
  peptideDiffs <- NULL
  groups <- unique(names(data))
  for(group in groups) {
    peptidesDiff <- t(t(as.integer(unionPeptides))) - data[,grep(group,names(data))]
    peptidesDiff <- apply(peptidesDiff, 1, mean, na.rm = TRUE)
    peptidesDiff <- na.omit(peptidesDiff)
    peptideDiffs <- c(peptideDiffs,peptidesDiff)
  }
  print(paste(replicateType,'median peptide difference:',round(mean(peptideDiffs),digits = 2)))
  peptideDiffSummary    <- NULL
  peptideDiffSummary[1] <- sum(peptideDiffs <= 2)/length(peptideDiffs)
  peptideDiffSummary[2] <- sum((peptideDiffs > 2)*(peptideDiffs < 5))/length(peptideDiffs)
  peptideDiffSummary[3] <- sum(peptideDiffs >= 5)/length(peptideDiffs)
  peptideDiffSummary    <- paste0(round(peptideDiffSummary*100,digits = 1),"%")
  return(peptideDiffSummary)
}

