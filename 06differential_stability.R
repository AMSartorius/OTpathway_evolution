######################################
##### SET UP WORKING ENVIRONMENT #####
######################################

# rm(list = ls()) # USE WITH CAUTION, delete all objects in the workspace, USE WITH CAUTION
# gc(reset = T)   # USE WITH CAUTION, resest memory (especially useful when working with large data sets), USE WITH CAUTION

options(stringsAsFactors = F) # disables automatic conversion of char. strings into factors
Sys.setenv(LANG = "en")

### .......... loading packages
library(readr)           # load .csv files
library(tidyverse)       # for everything and pipes
library(writexl)         # export .xlsx files 



## ================================= Differential stability analysis


OTpw.genes = c("OXT", "OXTR", "GNAQ", "HRAS", "KRAS", "NRAS", "RAF1", "MAP2K1", "MAP2K2", "MAPK1",
               "MAPK3" ,"PLA2G4E", "PLA2G4A", "JMJD7-PLA2G4B", "PLA2G4B", "PLA2G4C", "PLA2G4D", "PLA2G4F", "PTGS2", "MAP2K5",
               "MAPK7", "JUN", "FOS", "MEF2C", "CCND1", "ELK1", "RYR1", "RYR2", "RYR3", "CD38",         
               "TRPM2", "KCNJ2", "KCNJ12", "KCNJ18", "KCNJ4", "KCNJ14", "PLCB1", "PLCB2", "PLCB3", "PLCB4",        
               "PRKCA", "PRKCB", "PRKCG", "EEF2K", "EEF2", "CACNA1C", "CACNA1D", "CACNA1F", "CACNA1S", "CACNB1" ,      
               "CACNB2", "CACNB3", "CACNB4", "CACNA2D1", "CACNA2D2", "CACNA2D3", "CACNA2D4", "CACNG1", "CACNG2", "CACNG3",       
               "CACNG4", "CACNG5", "CACNG6", "CACNG7", "CACNG8", "ITPR1", "ITPR2", "ITPR3", "CALML3", "CALM2",
               "CALM3", "CALM1", "CALML6", "CALML5", "CALML4", "PPP3CA", "PPP3CB", "PPP3CC", "PPP3R1", "PPP3R2",       
               "NFATC1", "NFATC2", "NFATC3", "NFATC4", "RGS2", "RCAN1", "CAMKK2", "PRKAA1", "PRKAA2", "PRKAB1",       
               "PRKAB2", "PRKAG1", "PRKAG3", "PRKAG2", "CAMK1D", "CAMK1G", "CAMK1", "CAMK2A", "CAMK2D",  "CAMK2B",       
               "CAMK2G", "CAMK4", "NOS3", "GUCY1A2", "GUCY1A1", "GUCY1B1", "NPPA", "NPR1", "NPR2", "MYLK",        
               "MYLK2", "MYLK3", "MYLK4", "MYL6B", "MYL6", "MYL9", "ACTG1", "ACTB", "RHOA", "ROCK1",       
               "ROCK2", "PPP1CA", "PPP1CB", "PPP1CC", "PPP1R12A", "PPP1R12B", "PPP1R12C", "GNAS", "ADCY1", "ADCY2",       
               "ADCY3", "ADCY4", "ADCY5", "ADCY6", "ADCY7", "ADCY8", "ADCY9", "PRKACA", "PRKACB", "PRKACG",      
               "GNAI1", "GNAI3", "GNAI2", "GNAO1", "PIK3CG", "PIK3R5", "PIK3R6", "SRC", "KCNJ3", "KCNJ6",        
               "KCNJ9", "KCNJ5", "EGFR", "CDKN1A")     
# load gene list
OTpthwyDF <- data.frame(OTpw.genes) # Create OT pathway genes list


# load data generated with abagen
DSabagen      <- read_csv("data/processed/abagen_dk/DS_values_abagen.csv")   # each row corresponds to a gene column from the expression matrix, i.e., rows are genes here
gnames        <- read_csv("data/processed/abagen_dk/AHBAdk_220513_9861.csv")
genames       <- data.frame(colnames(gnames))     # get gene names
gene.names    <- data.frame(genames[-1,])
DSabagen.anno <- cbind(DSabagen, gene.names)   # annotate w/ gene names
15633/2                                          # estimate top 50% cut-off value
# 7816.5. Cut off: 7817, top 50% are N = 7816

# order by DS value
DSabagen.anno.desc <- DSabagen.anno %>% 
  arrange(desc(`0`))

DStop50   <- DSabagen.anno.desc[1:7816, ]             # extract top 50%
DSlower50 <- DSabagen.anno.desc[7817:15633, ]         # extract lowest 50%

DSOTtop50   <- dplyr::inner_join(OTpthwyDF, DStop50,   by = c("OTpw.genes" = "genames..1..."))       # match with *all* OT genes
DSOTlower50 <- dplyr::inner_join(OTpthwyDF, DSlower50, by = c("OTpw.genes" = "genames..1..."))     # match with OT *all* genes

OTgenes.out <- dplyr::anti_join(OTpthwyDF, DSabagen.anno.desc, by = c("OTpw.genes" = "genames..1..."))    # check whether some genes were not available

# no data available for 18 genes
# out of 136 genes, 105 were among the top 50% genes with highest DS values, 31 were among the genes with lower DS values
# out of the genes for which data was available (n = 136), 77.21% have high DS values and are thus consistently and reliably expressed across brain regions independent of donor characteristics. 


# prepare supp mat excel files
DSabagen.anno.desc$OTpathwaygene <- "n"
DSabagen.anno.desc[DSabagen.anno.desc$genames..1... %in% OTpw.genes, "OTpathwaygene"] <- "YES"   # this does not work, why not?
DSabagen.anno.desc$division <- "bottom50"
DSabagen.anno.desc$division[1:7816] <- "top50"

names(DSabagen.anno.desc)[names(DSabagen.anno.desc) == "X1"] <- "index"
names(DSabagen.anno.desc)[names(DSabagen.anno.desc) == "0"] <- "DS"
names(DSabagen.anno.desc)[names(DSabagen.anno.desc) == "genames..1..."] <- "GeneLabel"

which(DSabagen.anno.desc$OTpathwaygene == "YES" & DSabagen.anno.desc$division == "top50")
# ---> 105 genes in the OT pathway (out of 136, 18 not available) are in the top50 DS 

# for supplementary materials
write_xlsx(DSabagen.anno.desc, "output/supp_mat_new/supp_mat_09.xlsx")



######################################
########### END OF SCRIPT ############
######################################



