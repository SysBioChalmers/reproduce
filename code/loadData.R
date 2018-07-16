#!/usr/bin/env Rscript
# Functions for loading and minor formatting of all different datasets used in this study
# Benjamin Sanchez


## @knitr loadSwissprot
swissprot <- read.csv(file = '../data/raw_external/swissprot.tab', sep = '\t', header = TRUE)


## @knitr loadUPS2
UPS2 <- read.csv(file = '../data/raw_external/ups2.csv', sep = ',', header = TRUE)
# Convert fmol/kit to fmol in sample:
UPS2$amount.fmoles <- UPS2$amount.fmoles/10.6   #each kit has 10.6 ug of protein
UPS2$amount.fmoles <- UPS2$amount.fmoles*2.2    #fmol in sample (only 2.2 ug were used)


## @knitr loadIBAQdata
iBAQdata <- read.csv(file = '../data/raw_internal/mq14_int_std_peaklength-2-points_proteinGroups.txt', sep = '\t', header = TRUE)


## @knitr loadSILACdata
SILACdata <- read.csv(file = '../data/raw_internal/1710_mq14_sample_default_settings_proteinGroups.txt', sep = '\t', header = TRUE)
