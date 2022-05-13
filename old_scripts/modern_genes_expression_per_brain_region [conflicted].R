
library(psych)
library(powerAnalysis)


BASE="C:/01Dokumente/07Erwerbstaetigkeit/01PhD_UiO/projects/first_author/oxytocinHealth/RQ4/"


genesOT = c('GUCY1A2', 'GUCY1A1', 'CACNA2D2', 'CACNA2D1', 'CACNA2D4',  'CACNA2D3','PRKAB1', 'PRKAB2', 'PPP1CA', 'PPP1CB',
            'MYL6', 'MYL9', 'CALM2','CALM3','RAF1', 'CALM1','PRKAA2', 'ROCK1','ADCY1','ROCK2',
            'CALML6', 'ADCY5', 'CALML3', 'ADCY4','CALML4', 'ADCY3','ADCY2','CALML5', 'ADCY9','EGFR',   
            'PRKAA1', 'ADCY8','ADCY7','ADCY6','PPP1CC', 'NRAS', 'ACTB',  'CCND1', 'MAP2K5', 'PRKCG',
            'MAP2K2', 'PLA2G4E','PLA2G4F','MAP2K1', 'PLA2G4C','PLA2G4D','PLA2G4A','PLA2G4B','PRKCB','EEF2',
            'PRKCA',  'CACNB2', 'RCAN1',  'CACNB3', 'CACNB4', 'GNAQ', 'GNAS','CACNB1','OXTR', 'GUCY1B1',
            'CACNA1C','ACTG1','CACNA1D','ELK1',  'CACNA1F', 'MYL6B','MYLK', 'TRPM2','PPP3R1', 'PRKACG',
            'EEF2K','PPP3R2', 'RGS2', 'NPPA', 'CACNA1S','PRKACB', 'PRKACA', 'KCNJ5','RHOA', 'KCNJ6',
            'PPP1R12A', 'PPP1R12B', 'JUN', 'KCNJ12', 'KCNJ9', 'KCNJ14', 'NFATC4', 'NFATC3', 'NFATC2', 'NFATC1',
            'GNAO1', 'CAMK4', 'CAMK1', 'PPP1R12C', 'NPR1', 'NPR2', 'PIK3R5', 'ITPR2', 'ITPR3', 'OXT',
            'ITPR1', 'PIK3R6',  'PPP3CA', 'PPP3CB', 'PPP3CC', 'CD38', 'CAMK1D', 'MEF2C', 'NOS3', 'FOS',
            'PLCB4', 'JMJD7-PLA2G4B', 'KRAS', 'PLCB2', 'CAMK1G', 'PLCB3', 'PLCB1', 'MYLK2', 'RYR2', 'RYR3',
            'CAMK2D', 'MYLK3', 'SRC', 'RYR1', 'CDKN1A', 'CAMK2A', 'CAMK2B', 'PRKAG2', 'PRKAG3', 'PRKAG1',
            'GNAI1', 'GNAI2', 'CAMKK2', 'MYLK4', 'PTGS2', 'GNAI3', 'PIK3CG', 'CACNG7', 'CACNG8', 'MAPK3',
            'CACNG1', 'MAPK1', 'HRAS', 'CACNG2', 'CACNG3', 'MAPK7', 'CAMK2G', 'CACNG4', 'KCNJ2', 'KCNJ3',
            'CACNG5', 'CACNG6', 'KCNJ4')

get_data <-
  function(region) 
  {
    exp <- get_expression(structure_ids=c(region),
                          gene_ids=kop, 
                          dataset='5_stages') # Get expression values
    
    
    exp_df <- as.data.frame(do.call(rbind, exp)) 
    rownames(exp_df) <- NULL
    rownames(exp_df) <- c("Prenatal", "Infant", "Child", "Adolescent", "Adult")
    exp_t <- as.data.frame(t(exp_df))
    exp_t <- rownames_to_column(exp_t, var = "gene")                         
    
    ens_gene_dev <- ensembldb::select(EnsDb.Hsapiens.v79, 
                                      keys= exp_t$gene, 
                                      keytype = "SYMBOL", 
                                      columns = c("SYMBOL","GENEID"))
    
    ens_gene_dev <- ens_gene_dev %>% distinct(SYMBOL, .keep_all = TRUE)
    
    names(ens_gene_dev)[names(ens_gene_dev) == "SYMBOL"] <- "gene"
    
    ex_t_f <- dplyr::full_join(ens_gene_dev, exp_t, by = "gene")
    
    ex_t_f <- dplyr::select(ex_t_f, -gene)
    
    names(ex_t_f)[names(ex_t_f) == "GENEID"] <- "Gene_ID"
    
    HomoSapiens.PhyloMap <- read_excel(paste0(BASE, "RProject/OT-bd-analyses_loc/MBE_2008_Homo_Sapiens_PhyloMap.xls"), 
                                       sheet = 1, skip = 1)
    
    HomoSapiens.PhyloMap <- HomoSapiens.PhyloMap[ , c("Phylostratum","Gene_ID")]
    
    mm <- MatchMap(HomoSapiens.PhyloMap, ex_t_f)
    
    mm$GeneID <- toupper(mm$GeneID)
    
    mm_ids <- merge(ens_gene_dev, mm, by.x="GENEID", by.y = "GeneID", all= F)
    
    #mm$GeneName <- genesOT
    
    value <- list(
      mm = mm,
      mm_ids = mm_ids
    ) # Create a list of output objects
    attr(value, "class") <- "get_data"
    value
  }

get_name(10163)
M1C_dat <- get_data("Allen:10163")
M1C_mm <- M1C_dat$mm_ids

get_name(10173)
DFC_dat <- get_data("Allen:10173")
DFC_mm <- DFC_dat$mm_ids

get_name(10185)
VFC_dat <- get_data("Allen:10185")
VFC_mm <- VFC_dat$mm_ids

get_name(10194)
OFC_dat <- get_data("Allen:10194")
OFC_mm <- OFC_dat$mm_ids

get_name(10209)
S1C_dat <- get_data("Allen:10209")
S1C_mm <- S1C_dat$mm_ids

get_name(10225)
IPC_dat <- get_data("Allen:10225")
IPC_mm <- IPC_dat$mm_ids

get_name(10236)
A1C_dat <- get_data("Allen:10236")
A1C_mm <- A1C_dat$mm_ids

get_name(10243)
STC_dat <- get_data("Allen:10243")
STC_mm <- STC_dat$mm_ids

get_name(10252)
ITC_dat <- get_data("Allen:10252")
ITC_mm <- ITC_dat$mm_ids

get_name(10269)
V1C_dat <- get_data("Allen:10269")
V1C_mm <- V1C_dat$mm_ids

get_name(10278)
MFC_dat <- get_data("Allen:10278")
MFC_mm <- MFC_dat$mm_ids

get_name(10294)
HIP_dat <- get_data("Allen:10294")
HIP_mm <- HIP_dat$mm_ids

get_name(10333)
STR_dat <- get_data("Allen:10333")
STR_mm <- STR_dat$mm_ids

get_name(10361)
AMY_dat <- get_data("Allen:10361")
AMY_mm <- AMY_dat$mm_ids

get_name(10398)
MD_dat <- get_data("Allen:10398")
MD_mm <- MD_dat$mm_ids

get_name(10657)
CBC_dat <- get_data("Allen:10657")
CBC_mm <- CBC_dat$mm_ids


# get gene names for column names



temp <- M1C_mm[ ! M1C_mm$Phylostratum %in% c(1,2,3), ] # remove old/ancient phylostrata
col_n <- temp$gene                                     # use this object for all further formatting
temp <- data.frame(t(temp))
colnames(temp) <- NULL
colnames(temp) <- col_n
temp <- temp[-c(1:3), ]
region <- rep("M1C", 5)
temp <- cbind(region, temp)
M1C_fin <- rownames_to_column(temp, var = "stage") 

temp <- DFC_mm[ ! DFC_mm$Phylostratum %in% c(1,2,3), ]
temp <- data.frame(t(temp))
colnames(temp) <- NULL
colnames(temp) <- col_n
temp <- temp[-c(1:3), ]
region <- rep("DFC", 5)
temp <- cbind(region, temp)
DFC_fin <- rownames_to_column(temp, var = "stage") 

temp <- VFC_mm[ ! VFC_mm$Phylostratum %in% c(1,2,3), ]
temp <- data.frame(t(temp))
colnames(temp) <- NULL
colnames(temp) <- col_n
temp <- temp[-c(1:3), ]
region <- rep("VFC", 5)
temp <- cbind(region, temp)
VFC_fin <- rownames_to_column(temp, var = "stage") 

temp <- OFC_mm[ ! OFC_mm$Phylostratum %in% c(1,2,3), ]
temp <- data.frame(t(temp))
colnames(temp) <- NULL
colnames(temp) <- col_n
temp <- temp[-c(1:3), ]
region <- rep("OFC", 5)
temp <- cbind(region, temp)
OFC_fin <- rownames_to_column(temp, var = "stage") 

temp <- S1C_mm[ ! S1C_mm$Phylostratum %in% c(1,2,3), ]
temp <- data.frame(t(temp))
colnames(temp) <- NULL
colnames(temp) <- col_n
temp <- temp[-c(1:3), ]
region <- rep("S1C", 5)
temp <- cbind(region, temp)
S1C_fin <- rownames_to_column(temp, var = "stage") 

temp <- IPC_mm[ ! IPC_mm$Phylostratum %in% c(1,2,3), ]
temp <- data.frame(t(temp))
colnames(temp) <- NULL
colnames(temp) <- col_n
temp <- temp[-c(1:3), ]
region <- rep("IPC", 5)
temp <- cbind(region, temp)
IPC_fin <- rownames_to_column(temp, var = "stage") 

temp <- A1C_mm[ ! A1C_mm$Phylostratum %in% c(1,2,3), ]
temp <- data.frame(t(temp))
colnames(temp) <- NULL
colnames(temp) <- col_n
temp <- temp[-c(1:3), ]
region <- rep("A1C", 5)
temp <- cbind(region, temp)
A1C_fin <- rownames_to_column(temp, var = "stage") 

temp <- STC_mm[ ! STC_mm$Phylostratum %in% c(1,2,3), ]
temp <- data.frame(t(temp))
colnames(temp) <- NULL
colnames(temp) <- col_n
temp <- temp[-c(1:3), ]
region <- rep("STC", 5)
temp <- cbind(region, temp)
STC_fin <- rownames_to_column(temp, var = "stage") 

temp <- ITC_mm[ ! ITC_mm$Phylostratum %in% c(1,2,3), ]
temp <- data.frame(t(temp))
colnames(temp) <- NULL
colnames(temp) <- col_n
temp <- temp[-c(1:3), ]
region <- rep("ITC", 5)
temp <- cbind(region, temp)
ITC_fin <- rownames_to_column(temp, var = "stage") 

temp <- V1C_mm[ ! V1C_mm$Phylostratum %in% c(1,2,3), ]
temp <- data.frame(t(temp))
colnames(temp) <- NULL
colnames(temp) <- col_n
temp <- temp[-c(1:3), ]
region <- rep("V1C", 5)
temp <- cbind(region, temp)
V1C_fin <- rownames_to_column(temp, var = "stage") 

temp <- MFC_mm[ ! MFC_mm$Phylostratum %in% c(1,2,3), ]
temp <- data.frame(t(temp))
colnames(temp) <- NULL
colnames(temp) <- col_n
temp <- temp[-c(1:3), ]
region <- rep("MFC", 5)
temp <- cbind(region, temp)
MFC_fin <- rownames_to_column(temp, var = "stage") 

temp <- HIP_mm[ ! HIP_mm$Phylostratum %in% c(1,2,3), ]
temp <- data.frame(t(temp))
colnames(temp) <- NULL
colnames(temp) <- col_n
temp <- temp[-c(1:3), ]
region <- rep("HIP", 5)
temp <- cbind(region, temp)
HIP_fin <- rownames_to_column(temp, var = "stage") 

temp <- STR_mm[ ! STR_mm$Phylostratum %in% c(1,2,3), ]
temp <- data.frame(t(temp))
colnames(temp) <- NULL
colnames(temp) <- col_n
temp <- temp[-c(1:3), ]
region <- rep("STR", 5)
temp <- cbind(region, temp)
STR_fin <- rownames_to_column(temp, var = "stage") 

temp <- AMY_mm[ ! AMY_mm$Phylostratum %in% c(1,2,3), ]
temp <- data.frame(t(temp))
colnames(temp) <- NULL
colnames(temp) <- col_n
temp <- temp[-c(1:3), ]
region <- rep("AMY", 5)
temp <- cbind(region, temp)
AMY_fin <- rownames_to_column(temp, var = "stage") 

temp <- MD_mm[ ! MD_mm$Phylostratum %in% c(1,2,3), ]
temp <- data.frame(t(temp))
colnames(temp) <- NULL
colnames(temp) <- col_n
temp <- temp[-c(1:3), ]
region <- rep("MD", 5)
temp <- cbind(region, temp)
MD_fin <- rownames_to_column(temp, var = "stage") 


temp <- CBC_mm[ ! CBC_mm$Phylostratum %in% c(1,2,3), ]
temp <- data.frame(t(temp))
colnames(temp) <- NULL
colnames(temp) <- col_n
temp <- temp[-c(1:3), ]
region <- rep("CBC", 5)
temp <- cbind(region, temp)
CBC_fin <- rownames_to_column(temp, var = "stage") 


temp2 <- M1C_fin %>%
  bind_rows(DFC_fin) %>%
  bind_rows(VFC_fin) %>%
  bind_rows(OFC_fin) %>%
  bind_rows(S1C_fin) %>%
  bind_rows(IPC_fin) %>%
  bind_rows(A1C_fin) %>%
  bind_rows(STC_fin) %>%
  bind_rows(ITC_fin) %>%
  bind_rows(V1C_fin) %>%
  bind_rows(MFC_fin) %>%
  bind_rows(HIP_fin) %>%
  bind_rows(STR_fin) %>%
  bind_rows(AMY_fin) %>%
  bind_rows(MD_fin) %>%
  bind_rows(CBC_fin) 

# ----------------------------- CD38

temp2$CD38 <- as.numeric(temp2$CD38)
CD38_mean <- mean(temp2$CD38)
CD38_regions <- unique(temp2$region)
CD38_pvals <- numeric()
CD38_t <- numeric()
CD38_cohd <- numeric()

for (i in 1:length(CD38_regions)){
  Temp <- t.test(temp2[temp2$region==CD38_regions[i],"CD38"], 
                 mu = CD38_mean, alternative = "two.sided")
  CohD <- ES.t.one(m=Temp$estimate,
                   sd=sd(temp2[temp2$region==CD38_regions[i],"CD38"]),
                   mu=CD38_mean,
                   alternative = "two.sided")
  CD38_pvals <- c(CD38_pvals,Temp$p.value)
  CD38_t <- c(CD38_t,Temp$statistic)
  CD38_cohd <- c(CD38_cohd,CohD$d)
}

CD38_results <- data.frame(Region=CD38_regions,P_unadjusted=CD38_pvals,
                           t_stat=CD38_t,CohD=CD38_cohd)
CD38_results$P_fdr <- p.adjust(CD38_results$P_unadjusted, method="fdr")

CD38_results$stat_sign <- ""
CD38_results$stat_sign[CD38_results$P_fdr <= 0.05]  <- "*"

CD38_results 

## generating output file
fileout=paste0(BASE, "RProject/Ot-bd-analyses/output/CD38_results.txt")
write.table(CD38_results, fileout, quote=F, sep="\t", col.names = T, row.names = F)


# ----------------------------- CACNG3

temp2$CACNG3 <- as.numeric(temp2$CACNG3)
CACNG3_mean <- mean(temp2$CACNG3)
CACNG3_regions <- unique(temp2$region)
CACNG3_pvals <- numeric()
CACNG3_t <- numeric()
CACNG3_cohd <- numeric()

for (i in 1:length(CACNG3_regions)){
  Temp <- t.test(temp2[temp2$region==CACNG3_regions[i],"CACNG3"], 
                 mu = CACNG3_mean, alternative = "two.sided")
  CohD <- ES.t.one(m=Temp$estimate,
                   sd=sd(temp2[temp2$region==CACNG3_regions[i],"CACNG3"]),
                   mu=CACNG3_mean,
                   alternative = "two.sided")
  CACNG3_pvals <- c(CACNG3_pvals,Temp$p.value)
  CACNG3_t <- c(CACNG3_t,Temp$statistic)
  CACNG3_cohd <- c(CACNG3_cohd,CohD$d)
}

CACNG3_results <- data.frame(Region=CACNG3_regions,P_unadjusted=CACNG3_pvals,
                           t_stat=CACNG3_t,CohD=CACNG3_cohd)
CACNG3_results$P_fdr <- p.adjust(CACNG3_results$P_unadjusted, method="fdr")

CACNG3_results$stat_sign <- ""
CACNG3_results$stat_sign[CACNG3_results$P_fdr <= 0.05]  <- "*"

CACNG3_results 

## generating output file
fileout=paste0(BASE, "RProject/Ot-bd-analyses/output/CACNG3_results.txt")
write.table(CACNG3_results, fileout, quote=F, sep="\t", col.names = T, row.names = F)


# ----------------------------- CACNB1

temp2$CACNB1 <- as.numeric(temp2$CACNB1)
CACNB1_mean <- mean(temp2$CACNB1)
CACNB1_regions <- unique(temp2$region)
CACNB1_pvals <- numeric()
CACNB1_t <- numeric()
CACNB1_cohd <- numeric()

for (i in 1:length(CACNB1_regions)){
  Temp <- t.test(temp2[temp2$region==CACNB1_regions[i],"CACNB1"], 
                 mu = CACNB1_mean, alternative = "two.sided")
  CohD <- ES.t.one(m=Temp$estimate,
                   sd=sd(temp2[temp2$region==CACNB1_regions[i],"CACNB1"]),
                   mu=CACNB1_mean,
                   alternative = "two.sided")
  CACNB1_pvals <- c(CACNB1_pvals,Temp$p.value)
  CACNB1_t <- c(CACNB1_t,Temp$statistic)
  CACNB1_cohd <- c(CACNB1_cohd,CohD$d)
}

CACNB1_results <- data.frame(Region=CACNB1_regions,P_unadjusted=CACNB1_pvals,
                             t_stat=CACNB1_t,CohD=CACNB1_cohd)
CACNB1_results$P_fdr <- p.adjust(CACNB1_results$P_unadjusted, method="fdr")

CACNB1_results$stat_sign <- ""
CACNB1_results$stat_sign[CACNB1_results$P_fdr <= 0.05]  <- "*"

CACNB1_results 

## generating output file
fileout=paste0(BASE, "RProject/Ot-bd-analyses/output/CACNB1_results.txt")
write.table(CACNB1_results, fileout, quote=F, sep="\t", col.names = T, row.names = F)


# ----------------------------- NFATC3

temp2$NFATC3 <- as.numeric(temp2$NFATC3)
NFATC3_mean <- mean(temp2$NFATC3)
NFATC3_regions <- unique(temp2$region)
NFATC3_pvals <- numeric()
NFATC3_t <- numeric()
NFATC3_cohd <- numeric()

for (i in 1:length(NFATC3_regions)){
  Temp <- t.test(temp2[temp2$region==NFATC3_regions[i],"NFATC3"], 
                 mu = NFATC3_mean, alternative = "two.sided")
  CohD <- ES.t.one(m=Temp$estimate,
                   sd=sd(temp2[temp2$region==NFATC3_regions[i],"NFATC3"]),
                   mu=NFATC3_mean,
                   alternative = "two.sided")
  NFATC3_pvals <- c(NFATC3_pvals,Temp$p.value)
  NFATC3_t <- c(NFATC3_t,Temp$statistic)
  NFATC3_cohd <- c(NFATC3_cohd,CohD$d)
}

NFATC3_results <- data.frame(Region=NFATC3_regions,P_unadjusted=NFATC3_pvals,
                             t_stat=NFATC3_t,CohD=NFATC3_cohd)
NFATC3_results$P_fdr <- p.adjust(NFATC3_results$P_unadjusted, method="fdr")

NFATC3_results$stat_sign <- ""
NFATC3_results$stat_sign[NFATC3_results$P_fdr <= 0.05]  <- "*"

NFATC3_results 

## generating output file
fileout=paste0(BASE, "RProject/Ot-bd-analyses/output/NFATC3_results.txt")
write.table(NFATC3_results, fileout, quote=F, sep="\t", col.names = T, row.names = F)


# ----------------------------- NFATC4

temp2$NFATC4 <- as.numeric(temp2$NFATC4)
NFATC4_mean <- mean(temp2$NFATC4)
NFATC4_regions <- unique(temp2$region)
NFATC4_pvals <- numeric()
NFATC4_t <- numeric()
NFATC4_cohd <- numeric()

for (i in 1:length(NFATC4_regions)){
  Temp <- t.test(temp2[temp2$region==NFATC4_regions[i],"NFATC4"], 
                 mu = NFATC4_mean, alternative = "two.sided")
  CohD <- ES.t.one(m=Temp$estimate,
                   sd=sd(temp2[temp2$region==NFATC4_regions[i],"NFATC4"]),
                   mu=NFATC4_mean,
                   alternative = "two.sided")
  NFATC4_pvals <- c(NFATC4_pvals,Temp$p.value)
  NFATC4_t <- c(NFATC4_t,Temp$statistic)
  NFATC4_cohd <- c(NFATC4_cohd,CohD$d)
}

NFATC4_results <- data.frame(Region=NFATC4_regions,P_unadjusted=NFATC4_pvals,
                             t_stat=NFATC4_t,CohD=NFATC4_cohd)
NFATC4_results$P_fdr <- p.adjust(NFATC4_results$P_unadjusted, method="fdr")

NFATC4_results$stat_sign <- ""
NFATC4_results$stat_sign[NFATC4_results$P_fdr <= 0.05]  <- "*"

NFATC4_results 

## generating output file
fileout=paste0(BASE, "RProject/Ot-bd-analyses/output/NFATC4_results.txt")
write.table(NFATC4_results, fileout, quote=F, sep="\t", col.names = T, row.names = F)


# ----------------------------- NFATC2

temp2$NFATC2 <- as.numeric(temp2$NFATC2)
NFATC2_mean <- mean(temp2$NFATC2)
NFATC2_regions <- unique(temp2$region)
NFATC2_pvals <- numeric()
NFATC2_t <- numeric()
NFATC2_cohd <- numeric()

for (i in 1:length(NFATC2_regions)){
  Temp <- t.test(temp2[temp2$region==NFATC2_regions[i],"NFATC2"], 
                 mu = NFATC2_mean, alternative = "two.sided")
  CohD <- ES.t.one(m=Temp$estimate,
                   sd=sd(temp2[temp2$region==NFATC2_regions[i],"NFATC2"]),
                   mu=NFATC2_mean,
                   alternative = "two.sided")
  NFATC2_pvals <- c(NFATC2_pvals,Temp$p.value)
  NFATC2_t <- c(NFATC2_t,Temp$statistic)
  NFATC2_cohd <- c(NFATC2_cohd,CohD$d)
}

NFATC2_results <- data.frame(Region=NFATC2_regions,P_unadjusted=NFATC2_pvals,
                             t_stat=NFATC2_t,CohD=NFATC2_cohd)
NFATC2_results$P_fdr <- p.adjust(NFATC2_results$P_unadjusted, method="fdr")

NFATC2_results$stat_sign <- ""
NFATC2_results$stat_sign[NFATC2_results$P_fdr <= 0.05]  <- "*"

NFATC2_results 

## generating output file
fileout=paste0(BASE, "RProject/Ot-bd-analyses/output/NFATC2_results.txt")
write.table(NFATC2_results, fileout, quote=F, sep="\t", col.names = T, row.names = F)


# ----------------------------- OXT

temp2$OXT <- as.numeric(temp2$OXT)
OXT_mean <- mean(temp2$OXT)
OXT_regions <- unique(temp2$region)
OXT_pvals <- numeric()
OXT_t <- numeric()
OXT_cohd <- numeric()

for (i in 1:length(OXT_regions)){
  Temp <- t.test(temp2[temp2$region==OXT_regions[i],"OXT"], 
                 mu = OXT_mean, alternative = "two.sided")
  CohD <- ES.t.one(m=Temp$estimate,
                   sd=sd(temp2[temp2$region==OXT_regions[i],"OXT"]),
                   mu=OXT_mean,
                   alternative = "two.sided")
  OXT_pvals <- c(OXT_pvals,Temp$p.value)
  OXT_t <- c(OXT_t,Temp$statistic)
  OXT_cohd <- c(OXT_cohd,CohD$d)
}

OXT_results <- data.frame(Region=OXT_regions,P_unadjusted=OXT_pvals,
                             t_stat=OXT_t,CohD=OXT_cohd)
OXT_results$P_fdr <- p.adjust(OXT_results$P_unadjusted, method="fdr")

OXT_results$stat_sign <- ""
OXT_results$stat_sign[OXT_results$P_fdr <= 0.05]  <- "*"

OXT_results 

## generating output file
fileout=paste0(BASE, "RProject/Ot-bd-analyses/output/OXT_results.txt")
write.table(OXT_results, fileout, quote=F, sep="\t", col.names = T, row.names = F)


# ----------------------------- CACNG7

temp2$CACNG7 <- as.numeric(temp2$CACNG7)
CACNG7_mean <- mean(temp2$CACNG7)
CACNG7_regions <- unique(temp2$region)
CACNG7_pvals <- numeric()
CACNG7_t <- numeric()
CACNG7_cohd <- numeric()

for (i in 1:length(CACNG7_regions)){
  Temp <- t.test(temp2[temp2$region==CACNG7_regions[i],"CACNG7"], 
                 mu = CACNG7_mean, alternative = "two.sided")
  CohD <- ES.t.one(m=Temp$estimate,
                   sd=sd(temp2[temp2$region==CACNG7_regions[i],"CACNG7"]),
                   mu=CACNG7_mean,
                   alternative = "two.sided")
  CACNG7_pvals <- c(CACNG7_pvals,Temp$p.value)
  CACNG7_t <- c(CACNG7_t,Temp$statistic)
  CACNG7_cohd <- c(CACNG7_cohd,CohD$d)
}

CACNG7_results <- data.frame(Region=CACNG7_regions,P_unadjusted=CACNG7_pvals,
                          t_stat=CACNG7_t,CohD=CACNG7_cohd)
CACNG7_results$P_fdr <- p.adjust(CACNG7_results$P_unadjusted, method="fdr")

CACNG7_results$stat_sign <- ""
CACNG7_results$stat_sign[CACNG7_results$P_fdr <= 0.05]  <- "*"

CACNG7_results 

## generating output file
fileout=paste0(BASE, "RProject/Ot-bd-analyses/output/CACNG7_results.txt")
write.table(CACNG7_results, fileout, quote=F, sep="\t", col.names = T, row.names = F)


# ----------------------------- CACNG1

temp2$CACNG1 <- as.numeric(temp2$CACNG1)
CACNG1_mean <- mean(temp2$CACNG1)
CACNG1_regions <- unique(temp2$region)
CACNG1_pvals <- numeric()
CACNG1_t <- numeric()
CACNG1_cohd <- numeric()

for (i in 1:length(CACNG1_regions)){
  Temp <- t.test(temp2[temp2$region==CACNG1_regions[i],"CACNG1"], 
                 mu = CACNG1_mean, alternative = "two.sided")
  CohD <- ES.t.one(m=Temp$estimate,
                   sd=sd(temp2[temp2$region==CACNG1_regions[i],"CACNG1"]),
                   mu=CACNG1_mean,
                   alternative = "two.sided")
  CACNG1_pvals <- c(CACNG1_pvals,Temp$p.value)
  CACNG1_t <- c(CACNG1_t,Temp$statistic)
  CACNG1_cohd <- c(CACNG1_cohd,CohD$d)
}

CACNG1_results <- data.frame(Region=CACNG1_regions,P_unadjusted=CACNG1_pvals,
                             t_stat=CACNG1_t,CohD=CACNG1_cohd)
CACNG1_results$P_fdr <- p.adjust(CACNG1_results$P_unadjusted, method="fdr")

CACNG1_results$stat_sign <- ""
CACNG1_results$stat_sign[CACNG1_results$P_fdr <= 0.05]  <- "*"

CACNG1_results 

## generating output file
fileout=paste0(BASE, "RProject/Ot-bd-analyses/output/CACNG1_results.txt")
write.table(CACNG1_results, fileout, quote=F, sep="\t", col.names = T, row.names = F)


# ----------------------------- CDKN1A

temp2$CDKN1A <- as.numeric(temp2$CDKN1A)
CDKN1A_mean <- mean(temp2$CDKN1A)
CDKN1A_regions <- unique(temp2$region)
CDKN1A_pvals <- numeric()
CDKN1A_t <- numeric()
CDKN1A_cohd <- numeric()

for (i in 1:length(CDKN1A_regions)){
  Temp <- t.test(temp2[temp2$region==CDKN1A_regions[i],"CDKN1A"], 
                 mu = CDKN1A_mean, alternative = "two.sided")
  CohD <- ES.t.one(m=Temp$estimate,
                   sd=sd(temp2[temp2$region==CDKN1A_regions[i],"CDKN1A"]),
                   mu=CDKN1A_mean,
                   alternative = "two.sided")
  CDKN1A_pvals <- c(CDKN1A_pvals,Temp$p.value)
  CDKN1A_t <- c(CDKN1A_t,Temp$statistic)
  CDKN1A_cohd <- c(CDKN1A_cohd,CohD$d)
}

CDKN1A_results <- data.frame(Region=CDKN1A_regions,P_unadjusted=CDKN1A_pvals,
                             t_stat=CDKN1A_t,CohD=CDKN1A_cohd)
CDKN1A_results$P_fdr <- p.adjust(CDKN1A_results$P_unadjusted, method="fdr")

CDKN1A_results$stat_sign <- ""
CDKN1A_results$stat_sign[CDKN1A_results$P_fdr <= 0.05]  <- "*"

CDKN1A_results 

## generating output file
fileout=paste0(BASE, "RProject/Ot-bd-analyses/output/CDKN1A_results.txt")
write.table(CDKN1A_results, fileout, quote=F, sep="\t", col.names = T, row.names = F)


# ----------------------------- ELK1

temp2$ELK1 <- as.numeric(temp2$ELK1)
ELK1_mean <- mean(temp2$ELK1)
ELK1_regions <- unique(temp2$region)
ELK1_pvals <- numeric()
ELK1_t <- numeric()
ELK1_cohd <- numeric()

for (i in 1:length(ELK1_regions)){
  Temp <- t.test(temp2[temp2$region==ELK1_regions[i],"ELK1"], 
                 mu = ELK1_mean, alternative = "two.sided")
  CohD <- ES.t.one(m=Temp$estimate,
                   sd=sd(temp2[temp2$region==ELK1_regions[i],"ELK1"]),
                   mu=ELK1_mean,
                   alternative = "two.sided")
  ELK1_pvals <- c(ELK1_pvals,Temp$p.value)
  ELK1_t <- c(ELK1_t,Temp$statistic)
  ELK1_cohd <- c(ELK1_cohd,CohD$d)
}

ELK1_results <- data.frame(Region=ELK1_regions,P_unadjusted=ELK1_pvals,
                             t_stat=ELK1_t,CohD=ELK1_cohd)
ELK1_results$P_fdr <- p.adjust(ELK1_results$P_unadjusted, method="fdr")

ELK1_results$stat_sign <- ""
ELK1_results$stat_sign[ELK1_results$P_fdr <= 0.05]  <- "*"

ELK1_results 

## generating output file
fileout=paste0(BASE, "RProject/Ot-bd-analyses/output/ELK1_results.txt")
write.table(ELK1_results, fileout, quote=F, sep="\t", col.names = T, row.names = F)


# ----------------------------- CACNG6

temp2$CACNG6 <- as.numeric(temp2$CACNG6)
CACNG6_mean <- mean(temp2$CACNG6)
CACNG6_regions <- unique(temp2$region)
CACNG6_pvals <- numeric()
CACNG6_t <- numeric()
CACNG6_cohd <- numeric()

for (i in 1:length(CACNG6_regions)){
  Temp <- t.test(temp2[temp2$region==CACNG6_regions[i],"CACNG6"], 
                 mu = CACNG6_mean, alternative = "two.sided")
  CohD <- ES.t.one(m=Temp$estimate,
                   sd=sd(temp2[temp2$region==CACNG6_regions[i],"CACNG6"]),
                   mu=CACNG6_mean,
                   alternative = "two.sided")
  CACNG6_pvals <- c(CACNG6_pvals,Temp$p.value)
  CACNG6_t <- c(CACNG6_t,Temp$statistic)
  CACNG6_cohd <- c(CACNG6_cohd,CohD$d)
}

CACNG6_results <- data.frame(Region=CACNG6_regions,P_unadjusted=CACNG6_pvals,
                           t_stat=CACNG6_t,CohD=CACNG6_cohd)
CACNG6_results$P_fdr <- p.adjust(CACNG6_results$P_unadjusted, method="fdr")

CACNG6_results$stat_sign <- ""
CACNG6_results$stat_sign[CACNG6_results$P_fdr <= 0.05]  <- "*"

CACNG6_results 

## generating output file
fileout=paste0(BASE, "RProject/Ot-bd-analyses/output/CACNG6_results.txt")
write.table(CACNG6_results, fileout, quote=F, sep="\t", col.names = T, row.names = F)


# ----------------------------- PIK3R5

temp2$PIK3R5 <- as.numeric(temp2$PIK3R5)
PIK3R5_mean <- mean(temp2$PIK3R5)
PIK3R5_regions <- unique(temp2$region)
PIK3R5_pvals <- numeric()
PIK3R5_t <- numeric()
PIK3R5_cohd <- numeric()

for (i in 1:length(PIK3R5_regions)){
  Temp <- t.test(temp2[temp2$region==PIK3R5_regions[i],"PIK3R5"], 
                 mu = PIK3R5_mean, alternative = "two.sided")
  CohD <- ES.t.one(m=Temp$estimate,
                   sd=sd(temp2[temp2$region==PIK3R5_regions[i],"PIK3R5"]),
                   mu=PIK3R5_mean,
                   alternative = "two.sided")
  PIK3R5_pvals <- c(PIK3R5_pvals,Temp$p.value)
  PIK3R5_t <- c(PIK3R5_t,Temp$statistic)
  PIK3R5_cohd <- c(PIK3R5_cohd,CohD$d)
}

PIK3R5_results <- data.frame(Region=PIK3R5_regions,P_unadjusted=PIK3R5_pvals,
                             t_stat=PIK3R5_t,CohD=PIK3R5_cohd)
PIK3R5_results$P_fdr <- p.adjust(PIK3R5_results$P_unadjusted, method="fdr")

PIK3R5_results$stat_sign <- ""
PIK3R5_results$stat_sign[PIK3R5_results$P_fdr <= 0.05]  <- "*"

PIK3R5_results 

## generating output file
fileout=paste0(BASE, "RProject/Ot-bd-analyses/output/PIK3R5_results.txt")
write.table(PIK3R5_results, fileout, quote=F, sep="\t", col.names = T, row.names = F)


# ----------------------------- CACNB2

temp2$CACNB2 <- as.numeric(temp2$CACNB2)
CACNB2_mean <- mean(temp2$CACNB2)
CACNB2_regions <- unique(temp2$region)
CACNB2_pvals <- numeric()
CACNB2_t <- numeric()
CACNB2_cohd <- numeric()

for (i in 1:length(CACNB2_regions)){
  Temp <- t.test(temp2[temp2$region==CACNB2_regions[i],"CACNB2"], 
                 mu = CACNB2_mean, alternative = "two.sided")
  CohD <- ES.t.one(m=Temp$estimate,
                   sd=sd(temp2[temp2$region==CACNB2_regions[i],"CACNB2"]),
                   mu=CACNB2_mean,
                   alternative = "two.sided")
  CACNB2_pvals <- c(CACNB2_pvals,Temp$p.value)
  CACNB2_t <- c(CACNB2_t,Temp$statistic)
  CACNB2_cohd <- c(CACNB2_cohd,CohD$d)
}

CACNB2_results <- data.frame(Region=CACNB2_regions,P_unadjusted=CACNB2_pvals,
                             t_stat=CACNB2_t,CohD=CACNB2_cohd)
CACNB2_results$P_fdr <- p.adjust(CACNB2_results$P_unadjusted, method="fdr")

CACNB2_results$stat_sign <- ""
CACNB2_results$stat_sign[CACNB2_results$P_fdr <= 0.05]  <- "*"

CACNB2_results 

## generating output file
fileout=paste0(BASE, "RProject/Ot-bd-analyses/output/CACNB2_results.txt")
write.table(CACNB2_results, fileout, quote=F, sep="\t", col.names = T, row.names = F)


# ----------------------------- CACNG2

temp2$CACNG2 <- as.numeric(temp2$CACNG2)
CACNG2_mean <- mean(temp2$CACNG2)
CACNG2_regions <- unique(temp2$region)
CACNG2_pvals <- numeric()
CACNG2_t <- numeric()
CACNG2_cohd <- numeric()

for (i in 1:length(CACNG2_regions)){
  Temp <- t.test(temp2[temp2$region==CACNG2_regions[i],"CACNG2"], 
                 mu = CACNG2_mean, alternative = "two.sided")
  CohD <- ES.t.one(m=Temp$estimate,
                   sd=sd(temp2[temp2$region==CACNG2_regions[i],"CACNG2"]),
                   mu=CACNG2_mean,
                   alternative = "two.sided")
  CACNG2_pvals <- c(CACNG2_pvals,Temp$p.value)
  CACNG2_t <- c(CACNG2_t,Temp$statistic)
  CACNG2_cohd <- c(CACNG2_cohd,CohD$d)
}

CACNG2_results <- data.frame(Region=CACNG2_regions,P_unadjusted=CACNG2_pvals,
                             t_stat=CACNG2_t,CohD=CACNG2_cohd)
CACNG2_results$P_fdr <- p.adjust(CACNG2_results$P_unadjusted, method="fdr")

CACNG2_results$stat_sign <- ""
CACNG2_results$stat_sign[CACNG2_results$P_fdr <= 0.05]  <- "*"

CACNG2_results 

## generating output file
fileout=paste0(BASE, "RProject/Ot-bd-analyses/output/CACNG2_results.txt")
write.table(CACNG2_results, fileout, quote=F, sep="\t", col.names = T, row.names = F)


# ----------------------------- CACNB3

temp2$CACNB3 <- as.numeric(temp2$CACNB3)
CACNB3_mean <- mean(temp2$CACNB3)
CACNB3_regions <- unique(temp2$region)
CACNB3_pvals <- numeric()
CACNB3_t <- numeric()
CACNB3_cohd <- numeric()

for (i in 1:length(CACNB3_regions)){
  Temp <- t.test(temp2[temp2$region==CACNB3_regions[i],"CACNB3"], 
                 mu = CACNB3_mean, alternative = "two.sided")
  CohD <- ES.t.one(m=Temp$estimate,
                   sd=sd(temp2[temp2$region==CACNB3_regions[i],"CACNB3"]),
                   mu=CACNB3_mean,
                   alternative = "two.sided")
  CACNB3_pvals <- c(CACNB3_pvals,Temp$p.value)
  CACNB3_t <- c(CACNB3_t,Temp$statistic)
  CACNB3_cohd <- c(CACNB3_cohd,CohD$d)
}

CACNB3_results <- data.frame(Region=CACNB3_regions,P_unadjusted=CACNB3_pvals,
                             t_stat=CACNB3_t,CohD=CACNB3_cohd)
CACNB3_results$P_fdr <- p.adjust(CACNB3_results$P_unadjusted, method="fdr")

CACNB3_results$stat_sign <- ""
CACNB3_results$stat_sign[CACNB3_results$P_fdr <= 0.05]  <- "*"

CACNB3_results 

## generating output file
fileout=paste0(BASE, "RProject/Ot-bd-analyses/output/CACNB3_results.txt")
write.table(CACNB3_results, fileout, quote=F, sep="\t", col.names = T, row.names = F)


# ----------------------------- FOS

temp2$FOS <- as.numeric(temp2$FOS)
FOS_mean <- mean(temp2$FOS)
FOS_regions <- unique(temp2$region)
FOS_pvals <- numeric()
FOS_t <- numeric()
FOS_cohd <- numeric()

for (i in 1:length(FOS_regions)){
  Temp <- t.test(temp2[temp2$region==FOS_regions[i],"FOS"], 
                 mu = FOS_mean, alternative = "two.sided")
  CohD <- ES.t.one(m=Temp$estimate,
                   sd=sd(temp2[temp2$region==FOS_regions[i],"FOS"]),
                   mu=FOS_mean,
                   alternative = "two.sided")
  FOS_pvals <- c(FOS_pvals,Temp$p.value)
  FOS_t <- c(FOS_t,Temp$statistic)
  FOS_cohd <- c(FOS_cohd,CohD$d)
}

FOS_results <- data.frame(Region=FOS_regions,P_unadjusted=FOS_pvals,
                             t_stat=FOS_t,CohD=FOS_cohd)
FOS_results$P_fdr <- p.adjust(FOS_results$P_unadjusted, method="fdr")

FOS_results$stat_sign <- ""
FOS_results$stat_sign[FOS_results$P_fdr <= 0.05]  <- "*"

FOS_results 

## generating output file
fileout=paste0(BASE, "RProject/Ot-bd-analyses/output/FOS_results.txt")
write.table(FOS_results, fileout, quote=F, sep="\t", col.names = T, row.names = F)


# ----------------------------- NPPA

temp2$NPPA <- as.numeric(temp2$NPPA)
NPPA_mean <- mean(temp2$NPPA)
NPPA_regions <- unique(temp2$region)
NPPA_pvals <- numeric()
NPPA_t <- numeric()
NPPA_cohd <- numeric()

for (i in 1:length(NPPA_regions)){
  Temp <- t.test(temp2[temp2$region==NPPA_regions[i],"NPPA"], 
                 mu = NPPA_mean, alternative = "two.sided")
  CohD <- ES.t.one(m=Temp$estimate,
                   sd=sd(temp2[temp2$region==NPPA_regions[i],"NPPA"]),
                   mu=NPPA_mean,
                   alternative = "two.sided")
  NPPA_pvals <- c(NPPA_pvals,Temp$p.value)
  NPPA_t <- c(NPPA_t,Temp$statistic)
  NPPA_cohd <- c(NPPA_cohd,CohD$d)
}

NPPA_results <- data.frame(Region=NPPA_regions,P_unadjusted=NPPA_pvals,
                          t_stat=NPPA_t,CohD=NPPA_cohd)
NPPA_results$P_fdr <- p.adjust(NPPA_results$P_unadjusted, method="fdr")

NPPA_results$stat_sign <- ""
NPPA_results$stat_sign[NPPA_results$P_fdr <= 0.05]  <- "*"

NPPA_results 

## generating output file
fileout=paste0(BASE, "RProject/Ot-bd-analyses/output/NPPA_results.txt")
write.table(NPPA_results, fileout, quote=F, sep="\t", col.names = T, row.names = F)


# ----------------------------- OXTR

temp2$OXTR <- as.numeric(temp2$OXTR)
OXTR_mean <- mean(temp2$OXTR)
OXTR_regions <- unique(temp2$region)
OXTR_pvals <- numeric()
OXTR_t <- numeric()
OXTR_cohd <- numeric()

for (i in 1:length(OXTR_regions)){
  Temp <- t.test(temp2[temp2$region==OXTR_regions[i],"OXTR"], 
                 mu = OXTR_mean, alternative = "two.sided")
  CohD <- ES.t.one(m=Temp$estimate,
                   sd=sd(temp2[temp2$region==OXTR_regions[i],"OXTR"]),
                   mu=OXTR_mean,
                   alternative = "two.sided")
  OXTR_pvals <- c(OXTR_pvals,Temp$p.value)
  OXTR_t <- c(OXTR_t,Temp$statistic)
  OXTR_cohd <- c(OXTR_cohd,CohD$d)
}

OXTR_results <- data.frame(Region=OXTR_regions,P_unadjusted=OXTR_pvals,
                           t_stat=OXTR_t,CohD=OXTR_cohd)
OXTR_results$P_fdr <- p.adjust(OXTR_results$P_unadjusted, method="fdr")

OXTR_results$stat_sign <- ""
OXTR_results$stat_sign[OXTR_results$P_fdr <= 0.05]  <- "*"

OXTR_results 

## generating output file
fileout=paste0(BASE, "RProject/Ot-bd-analyses/output/OXTR_results.txt")
write.table(OXTR_results, fileout, quote=F, sep="\t", col.names = T, row.names = F)


# ----------------------------- CACNB4

temp2$CACNB4 <- as.numeric(temp2$CACNB4)
CACNB4_mean <- mean(temp2$CACNB4)
CACNB4_regions <- unique(temp2$region)
CACNB4_pvals <- numeric()
CACNB4_t <- numeric()
CACNB4_cohd <- numeric()

for (i in 1:length(CACNB4_regions)){
  Temp <- t.test(temp2[temp2$region==CACNB4_regions[i],"CACNB4"], 
                 mu = CACNB4_mean, alternative = "two.sided")
  CohD <- ES.t.one(m=Temp$estimate,
                   sd=sd(temp2[temp2$region==CACNB4_regions[i],"CACNB4"]),
                   mu=CACNB4_mean,
                   alternative = "two.sided")
  CACNB4_pvals <- c(CACNB4_pvals,Temp$p.value)
  CACNB4_t <- c(CACNB4_t,Temp$statistic)
  CACNB4_cohd <- c(CACNB4_cohd,CohD$d)
}

CACNB4_results <- data.frame(Region=CACNB4_regions,P_unadjusted=CACNB4_pvals,
                           t_stat=CACNB4_t,CohD=CACNB4_cohd)
CACNB4_results$P_fdr <- p.adjust(CACNB4_results$P_unadjusted, method="fdr")

CACNB4_results$stat_sign <- ""
CACNB4_results$stat_sign[CACNB4_results$P_fdr <= 0.05]  <- "*"

CACNB4_results 

## generating output file
fileout=paste0(OUT_DIR, "RProject/Ot-bd-analyses/output/CACNB4_results.txt")
write.table(CACNB4_results, fileout, quote=F, sep="\t", col.names = T, row.names = F)




################



