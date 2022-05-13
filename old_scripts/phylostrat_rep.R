### 210409, amsartorius ###

#########################
# INFO
#########################

# Project: Enrichment of oxytocin pathway genes in human evolution

# R script template: xx

# Correspondence to: 
# Dr. Daniel S. Quintana, NORMENT, University of Oslo, daniel.quintana@medisin.uio.no
# Alina M. Sartorius, Department of Psychology, NORMENT, University of Oslo, a.m.sartorius@psykologi.uio.no

#########################
# setting up working environment 
#########################

rm(list=ls()) # delete all objects in the workspace
gc(reset=T) # resest memory (especially useful when working with large data sets)

options(stringsAsFactors=F) # disables automatic conversion of char. strings into factors

BASE_PATH="D:/eigene_dateien/02Dokumente/07Erwerbstaetigkeit/01PhD_UiO/projects/first_author/oxytocinHealth/RQ4/analyses/"     # base directory, adjust to local dir structure
OUT_DIR=paste0(BASE_PATH,"data/processed/out/phylostratigraphy/")                                                             # output directory, adjust to local dir structure


#########################
# (install and) load packages 
#########################

#library(dplyr)
library(readr)
library(readxl) 
#library(powerAnalysis)   # for effect size (Cohen's D) calculation     
#library(psych)
#library(tidyverse)


#########################
# load and prepare data
#########################

ex_m <- read_csv(paste0(BASE_PATH, "data/raw/brainspanatlas/expression_matrix.csv"), col_names=F) 
ex_m <- as.data.frame(ex_m)
rows_meta <- read.csv(paste0(BASE_PATH, "data/raw/brainspanatlas/rows_metadata.csv"))
cols_meta <- read.csv(paste0(BASE_PATH, "data/raw/brainspanatlas/columns_metadata.csv"))
phylo <- read_excel(paste0(BASE_PATH, "data/raw/phylo.xlsx"))
OT_pathway_genes <- read_excel(paste0(BASE_PATH, "data/raw/OT_pathway_geneset_153genes_R.xlsx"))


# ontogenic stages: 4-7 pcw Embryonic
#                   8-38 pcw prenatal
#                   Birth-18 months infancy (1.5 yrs)
#                   19 months-11 yrs childhood
#                   12-19 yrs Adolescence
#                   20-60+ yrs Adulthood


# regions of interest: M1C, DFC, VFC, OFC, S1C, IPC, A1C, STC, ITC, V1C, MFC, HIP, STR, AMY, MD, CBC


# restructure and annotate dataframe

ex_m_t2 <- t(ex_m)
ex_m_t2 <- ex_m_t2[-c(1), ]
ex_m_t2 <- cbind(col_num = c(1:524), ex_m_t2)
ex_m_col <- merge(cols_meta, ex_m_t2, by.x="column_num", by.y="col_num", all=T)
ex_m_col$age<- gsub('\\s+', '', ex_m_col$age)                              # removes all spaces in the column 'age'


# create separate DFs for each brain region

structures = NULL
str = c("M1C", "DFC", "VFC", "OFC", "S1C", "IPC", "A1C", "STC", "ITC", "V1C", "MFC", "HIP", 
        "STR","AMY", "MD", "CBC")
for (i in str) {
  str_temp <- ex_m_col[which(ex_m_col$structure_acronym == i), ]
  str_temp <- as.data.frame(str_temp)
  assign(paste0("df_",i), str_temp)                                 #make new df
}


# proceed with each brain region individually

# -----------------------------------------------A1C

hd_A1C <- df_A1C$age
df_A1C <- df_A1C[,-c(1:8)]
A1C_t <- t(df_A1C)
A1C_t <- as.data.frame((A1C_t))
colnames(A1C_t) <- hd_A1C


prenatal_A1C <- rowMeans(A1C_t[, c(1:14)])
infant_A1C <- rowMeans(A1C_t[, c(15:17)])
child_A1C <- rowMeans(A1C_t[, c(18:22)])
adols_A1C <- rowMeans(A1C_t[, c(23:25)])
adult_A1C <- rowMeans(A1C_t[, c(26:31)])

A1C_fin <- Reduce(function(...) merge(..., all=TRUE), 
                  list(data.frame(prenatal_A1C, infant_A1C, child_A1C, adols_A1C, adult_A1C)))


# -----------------------------------------------DFC

hd_DFC <- df_DFC$age
df_DFC <- df_DFC[,-c(1:8)]
DFC_t <- t(df_DFC)
DFC_t <- as.data.frame((DFC_t))
colnames(DFC_t) <- hd_DFC


prenatal_DFC <- rowMeans(DFC_t[, c(1:17)])
infant_DFC <- rowMeans(DFC_t[, c(18:21)])
child_DFC <- rowMeans(DFC_t[, c(22:27)])
adols_DFC <- rowMeans(DFC_t[, c(28:30)])
adult_DFC <- rowMeans(DFC_t[, c(31:35)])

DFC_fin <- Reduce(function(...) merge(..., all=TRUE), 
                  list(data.frame(prenatal_DFC, infant_DFC, child_DFC, adols_DFC, adult_DFC)))


# -----------------------------------------------VFC

hd_VFC <- df_VFC$age
df_VFC <- df_VFC[,-c(1:8)]
VFC_t <- t(df_VFC)
VFC_t <- as.data.frame((VFC_t))
colnames(VFC_t) <- hd_VFC


prenatal_VFC <- rowMeans(VFC_t[, c(1:16)])
infant_VFC <- rowMeans(VFC_t[, c(17:19)])
child_VFC <- rowMeans(VFC_t[, c(20:26)])
adols_VFC <- rowMeans(VFC_t[, c(27:29)])
adult_VFC <- rowMeans(VFC_t[, c(30:35)])

VFC_fin <- Reduce(function(...) merge(..., all=TRUE), 
                  list(data.frame(prenatal_VFC, infant_VFC, child_VFC, adols_VFC, adult_VFC)))

View(VFC_t)
head(hd_VFC)


# ..... tbc



#########################
# time for the analyses
#########################












