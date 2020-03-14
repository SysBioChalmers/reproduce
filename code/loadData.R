#!/usr/bin/env Rscript
# Functions for loading and minor formatting of all different datasets used in this study
# Benjamin Sanchez


## @knitr loadUPS2
UPS2 <- read.csv(file = '../data/raw_external/ups2.csv', sep = ',', header = TRUE)
# Convert fmol/kit to fmol in sample:
UPS2$amount.fmoles <- UPS2$amount.fmoles/10.6   #fmol/ug (each kit has 10.6 ug of protein)
UPS2$amount.fmoles <- UPS2$amount.fmoles*1.1    #fmol in sample (only 1.1 ug are injected)


## @knitr loadIBAQdata
fileName <- '../data/raw_internal/mq14_int_std_peaklength-2-points_proteinGroups.txt'
iBAQdata <- read.csv(file = fileName, sep = '\t', header = TRUE)
#P00044 was erroneously match as P99999 (the human homolog):
iBAQdata$Protein.IDs <- gsub(';P99999ups','',iBAQdata$Protein.IDs)  
names(iBAQdata) <- gsub('iBAQ.L.T4h_','Abundance.iBAQ.ES.',names(iBAQdata))
names(iBAQdata) <- gsub('iBAQ.H.T4h_','Abundance.iBAQ.IS.',names(iBAQdata))
names(iBAQdata) <- gsub('_Batch','_batch',names(iBAQdata))
#The 8 hour gradient setting will not be used in this study
iBAQdata[,grep('T8h_',names(iBAQdata))] <- NULL


## @knitr loadSILACdata
fileName  <- '../data/raw_internal/1710_mq14_sample_default_settings_proteinGroups.txt'
SILACdata <- read.csv(file = fileName, sep = '\t', header = TRUE)
names(SILACdata) <- gsub('_Batch','_batch',names(SILACdata))


## @knitr loadNtheoPeptides
fileName    <- '../data/raw_internal/int_std_proteinGroups_theoretical_peptides.txt'
NTPdata     <- read.csv(file = fileName, sep = '\t', header = TRUE)
fileName    <- '../data/raw_internal/181203_UPS2_theoretical_peptide_nrs.txt'
NTPdataUPS2 <-read.csv(file = fileName, sep = '\t', header = TRUE)
NTPdataUPS2$Protein.IDs <- gsub(';.*','',NTPdataUPS2$Protein.IDs)


## @knitr loadRibProteins
RP <- read.csv('../data/raw_external/RP-list.csv')
# Add uniprot IDs:
ensembl <- useMart('ENSEMBL_MART_ENSEMBL', dataset='scerevisiae_gene_ensembl',
                   host='http://dec2016.archive.ensembl.org/')
RP      <- getBM(attributes=c('wikigene_name','uniprot_swissprot'),
                 filters = 'wikigene_name', values = RP$wikigene_name, mart = ensembl)
# Remove duplicated protein entries:
RP <- RP[!duplicated(RP$uniprot_swissprot),]
# Change variable name:
names(RP) <- gsub('uniprot_swissprot','Protein.IDs',names(RP))
