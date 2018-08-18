#!/usr/bin/env Rscript
# Functions for loading and minor formatting of all different datasets used in this study
# Benjamin Sanchez


## @knitr loadSwissprot
swissprot <- read.csv(file = '../data/raw_external/swissprot.tab', sep = '\t', header = TRUE)


## @knitr loadUPS2
UPS2 <- read.csv(file = '../data/raw_external/ups2.csv', sep = ',', header = TRUE)
# Convert fmol/kit to fmol in sample:
UPS2$amount.fmoles <- UPS2$amount.fmoles/10.6   #fmol/ug (each kit has 10.6 ug of protein)
UPS2$amount.fmoles <- UPS2$amount.fmoles*1.1    #fmol in sample (only 1.1 ug are injected)


## @knitr loadIBAQdata
iBAQdata <- read.csv(file = '../data/raw_internal/mq14_int_std_peaklength-2-points_proteinGroups.txt', sep = '\t', header = TRUE)
iBAQdata$Protein.IDs <- gsub(';P99999ups','',iBAQdata$Protein.IDs)  #P00044 was erroneously match as P99999 (the human homolog)
names(iBAQdata) <- gsub('iBAQ.L.T4h_','Abundance.MaxQuant.ES.',names(iBAQdata))
names(iBAQdata) <- gsub('iBAQ.H.T4h_','Abundance.MaxQuant.IS.',names(iBAQdata))
names(iBAQdata) <- gsub('_Batch','_batch',names(iBAQdata))

## @knitr loadSILACdata
SILACdata <- read.csv(file = '../data/raw_internal/1710_mq14_sample_default_settings_proteinGroups.txt', sep = '\t', header = TRUE)
names(SILACdata) <- gsub('_Batch','_batch',names(SILACdata))
