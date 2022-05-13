######################################
##### SET UP WORKING ENVIRONMENT #####
######################################

rm(list = ls()) # delete all objects in the workspace
gc(reset = T) # resest memory (especially useful when working with large data sets)

options(stringsAsFactors = F) # disables automatic conversion of char. strings into factors

Sys.setenv(LANG = "en")

# set working directory, only necessary when not working with R project
#setwd("D:/eigene_dateien/02Dokumente/07Erwerbstaetigkeit/01PhD_UiO/projects/first_author/oxytocinHealth/RQ4/analyses/phylogeny_microsynteny/BLASTp/")

# load packages
library(data.table) # for loading data
library(tidyr)      # always useful
library(ggplot2)    # plotting
library(stringr)    # replacing strings
library(dplyr)
library(readr)

### load example BLASTp runs
#cacnb1.ex <- read_csv("data/raw/ref_BLASTp/BLASTp_cox1_csavignyi_ccrescentus_refBurgeri.csv")

###################################


path = "D:/eigene_dateien/02Dokumente/07Erwerbstaetigkeit/01PhD_UiO/projects/first_author/OTpathway_evolution/RProject/OT-bd-analyses_loc/"
temp <- list.files(path = paste0(path, "data/raw/BLASTp/all/"), full.names = T)
temp3 <- list.files(path = paste0(path, "data/raw/BLASTp/all/"), full.names = F)


L1 <- list()
for (i in 1:length(temp)) {
  L1[[i]] = data.frame(read_csv(temp[i])) # read in csvs into list
  names(L1)[i] <- paste(temp3[i])
  names(L1[[i]]) <- make.names(names(L1[[i]]), unique = TRUE)
  #L1 <- tools::file_path_sans_ext(L1) 
}

# list2env(L1 ,.GlobalEnv)   # unlists L1 list with dataframes and 

pat = data.frame(c("Ciona savignyi", "Ciona intestinalis", "Branchiostoma lanceolatum", "Branchiostoma floridae", "Strongylocentrotus purpuratus", 
                   "Asterias rubens", "Saccoglossus kowalevskii", "Xenoturbella bocki", "Symsagittifera roscoffensis", "Caenorhabditis elegans", 
                   "Priapulus caudatus", "Octopus vulgaris", "Hydra vulgaris", "Acropora millepora", "Clytia hemisphaerica", "Ephydatia muelleri", 
                   "Amphimedon queenslandica", "Salpingoeca rosetta", "Monosiga brevicollis", "Fusarium culmorum", "Aspergillus nidulans", 
                   "Saccharomyces cerevisiae", "Dictyostelium discoideum", "Spironucleus salmonicida", "Methanolobus zinderi", "Escherichia coli",
                   "Caulobacter crescentus", "Caulobacter vibrioides"))
pat$c..Ciona.savignyi....Ciona.intestinalis....Branchiostoma.lanceolatum... <- as.character(pat$c..Ciona.savignyi....Ciona.intestinalis....Branchiostoma.lanceolatum...)

names <- data.frame(names(L1))
View(names)


# loops through ONE ref gene csv file and returns list with max scores, accession lengths, 
# protein length of ref gene for each species, and results


### ----------------------------------

# -------------------------------------------- template
df.temp <- data.frame(1:nrow(pat))
L.temp <- list()

for (k in 1:nrow(pat)) {
  L.temp[[k]] <- data.frame(L1[[0000]][grep(pat[k,], L1[[0000]][["Scientific.Name"]]), ])
  df.temp[k,] <- max((data.frame(L1[[0000]][grep(pat[k,], L1[[0000]][["Scientific.Name"]]), ]))$Max.Score)
  df.temp[k, 2] <- first(L.temp[[k]][["Acc..Len"]][L.temp[[k]][["Max.Score"]] == max(L.temp[[k]][["Max.Score"]])])
  df.temp[, 3] <- rep(0000, 28)
  df.temp[, 4] <- df.temp[, 1] / ((df.temp[, 3]) + df.temp[, 2])
}

### ----------------------------------



# OXT
df.001oxt <- data.frame(1:nrow(pat))
L.001oxt  <- list()

for (k in 1:nrow(pat)) {
  L.001oxt[[k]] <- data.frame(L1[[1]][grep(pat[k,], L1[[1]][["Scientific.Name"]]), ])
  df.001oxt[k,] <- max((data.frame(L1[[1]][grep(pat[k,], L1[[1]][["Scientific.Name"]]), ]))$Max.Score)
  df.001oxt[k, 2] <- first(L.001oxt[[k]][["Acc..Len"]][L.001oxt[[k]][["Max.Score"]] == max(L.001oxt[[k]][["Max.Score"]])])
  df.001oxt[, 3] <- rep(126, 28)
  df.001oxt[, 4] <- df.001oxt[, 1] / ((df.001oxt[, 3]) + df.001oxt[, 2])
  colnames(df.001oxt) <- c("maxscore", "p.length", "ref.p.length", "OXT.maxscore.scaled")
}


# OXTR
df.002oxtr <- data.frame(1:nrow(pat))
L.002oxtr  <- list()

for (k in 1:nrow(pat)) {
  L.002oxtr[[k]] <- data.frame(L1[[2]][grep(pat[k,], L1[[2]][["Scientific.Name"]]), ])
  df.002oxtr[k,] <- max((data.frame(L1[[2]][grep(pat[k,], L1[[2]][["Scientific.Name"]]), ]))$Max.Score)
  df.002oxtr[k, 2] <- first(L.002oxtr[[k]][["Acc..Len"]][L.002oxtr[[k]][["Max.Score"]] == max(L.002oxtr[[k]][["Max.Score"]])])
  df.002oxtr[, 3] <- rep(431, 28)
  df.002oxtr[, 4] <- df.002oxtr[, 1] / ((df.002oxtr[, 3]) + df.002oxtr[, 2])
  colnames(df.002oxtr) <- c("maxscore", "p.length", "ref.p.length", "OXTR.maxscore.scaled")
}


# GNAQ
df.003gnaq <- data.frame(1:nrow(pat))
L.003gnaq  <- list()

for (k in 1:nrow(pat)) {
  L.003gnaq[[k]] <- data.frame(L1[[3]][grep(pat[k,], L1[[3]][["Scientific.Name"]]), ])
  df.003gnaq[k,] <- max((data.frame(L1[[3]][grep(pat[k,], L1[[3]][["Scientific.Name"]]), ]))$Max.Score)
  df.003gnaq[k, 2] <- first(L.003gnaq[[k]][["Acc..Len"]][L.003gnaq[[k]][["Max.Score"]] == max(L.003gnaq[[k]][["Max.Score"]])])
  df.003gnaq[, 3] <- rep(359, 28)
  df.003gnaq[, 4] <- df.003gnaq[, 1] / ((df.003gnaq[, 3]) + df.003gnaq[, 2])
  colnames(df.003gnaq) <- c("maxscore", "p.length", "ref.p.length", "GNAQ.maxscore.scaled")
}


# HRAS
df.004hras <- data.frame(1:nrow(pat))
L.004hras <- list()

for (k in 1:nrow(pat)) {
  L.004hras[[k]] <- data.frame(L1[[4]][grep(pat[k,], L1[[4]][["Scientific.Name"]]), ])
  df.004hras[k,] <- max((data.frame(L1[[4]][grep(pat[k,], L1[[4]][["Scientific.Name"]]), ]))$Max.Score)
  df.004hras[k, 2] <- first(L.004hras[[k]][["Acc..Len"]][L.004hras[[k]][["Max.Score"]] == max(L.004hras[[k]][["Max.Score"]])])
  df.004hras[, 3] <- rep(188, 28)
  df.004hras[, 4] <- df.004hras[, 1] / ((df.004hras[, 3]) + df.004hras[, 2])
  colnames(df.004hras) <- c("maxscore", "p.length", "ref.p.length", "HRAS.maxscore.scaled")
}

# KRAS
df.005kras <- data.frame(1:nrow(pat))
L.005kras <- list()

for (k in 1:nrow(pat)) {
  L.005kras[[k]] <- data.frame(L1[[5]][grep(pat[k,], L1[[5]][["Scientific.Name"]]), ])
  df.005kras[k,] <- max((data.frame(L1[[5]][grep(pat[k,], L1[[5]][["Scientific.Name"]]), ]))$Max.Score)
  df.005kras[k, 2] <- first(L.005kras[[k]][["Acc..Len"]][L.005kras[[k]][["Max.Score"]] == max(L.005kras[[k]][["Max.Score"]])])
  df.005kras[, 3] <- rep(207, 28)
  df.005kras[, 4] <- df.005kras[, 1] / ((df.005kras[, 3]) + df.005kras[, 2])
  colnames(df.005kras) <- c("maxscore", "p.length", "ref.p.length", "KRAS.maxscore.scaled")
}


# NRAS
df.006nras <- data.frame(1:nrow(pat))
L.006nras <- list()

for (k in 1:nrow(pat)) {
  L.006nras[[k]] <- data.frame(L1[[6]][grep(pat[k,], L1[[6]][["Scientific.Name"]]), ])
  df.006nras[k,] <- max((data.frame(L1[[6]][grep(pat[k,], L1[[6]][["Scientific.Name"]]), ]))$Max.Score)
  df.006nras[k, 2] <- first(L.006nras[[k]][["Acc..Len"]][L.006nras[[k]][["Max.Score"]] == max(L.006nras[[k]][["Max.Score"]])])
  df.006nras[, 3] <- rep(189, 28)
  df.006nras[, 4] <- df.006nras[, 1] / ((df.006nras[, 3]) + df.006nras[, 2])
  colnames(df.006nras) <- c("maxscore", "p.length", "ref.p.length", "NRAS.maxscore.scaled")
}


# RAF1
df.007raf1 <- data.frame(1:nrow(pat))
L.007raf1 <- list()

for (k in 1:nrow(pat)) {
  L.007raf1[[k]] <- data.frame(L1[[7]][grep(pat[k,], L1[[7]][["Scientific.Name"]]), ])
  df.007raf1[k,] <- max((data.frame(L1[[7]][grep(pat[k,], L1[[7]][["Scientific.Name"]]), ]))$Max.Score)
  df.007raf1[k, 2] <- first(L.007raf1[[k]][["Acc..Len"]][L.007raf1[[k]][["Max.Score"]] == max(L.007raf1[[k]][["Max.Score"]])])
  df.007raf1[, 3] <- rep(650, 28)
  df.007raf1[, 4] <- df.007raf1[, 1] / ((df.007raf1[, 3]) + df.007raf1[, 2])
  colnames(df.007raf1) <- c("maxscore", "p.length", "ref.p.length", "RAF1.maxscore.scaled")
}


# MAP2K1
df.008map2k1 <- data.frame(1:nrow(pat))
L.008map2k1 <- list()

for (k in 1:nrow(pat)) {
  L.008map2k1[[k]] <- data.frame(L1[[8]][grep(pat[k,], L1[[8]][["Scientific.Name"]]), ])
  df.008map2k1[k,] <- max((data.frame(L1[[8]][grep(pat[k,], L1[[8]][["Scientific.Name"]]), ]))$Max.Score)
  df.008map2k1[k, 2] <- first(L.008map2k1[[k]][["Acc..Len"]][L.008map2k1[[k]][["Max.Score"]] == max(L.008map2k1[[k]][["Max.Score"]])])
  df.008map2k1[, 3] <- rep(393, 28)
  df.008map2k1[, 4] <- df.008map2k1[, 1] / ((df.008map2k1[, 3]) + df.008map2k1[, 2])
  colnames(df.008map2k1) <- c("maxscore", "p.length", "ref.p.length",  "MAP2K1.maxscore.scaled")
}



# MAP2K2
df.009map2k2 <- data.frame(1:nrow(pat))
L.009map2k2 <- list()

for (k in 1:nrow(pat)) {
  L.009map2k2[[k]] <- data.frame(L1[[9]][grep(pat[k,], L1[[9]][["Scientific.Name"]]), ])
  df.009map2k2[k,] <- max((data.frame(L1[[9]][grep(pat[k,], L1[[9]][["Scientific.Name"]]), ]))$Max.Score)
  df.009map2k2[k, 2] <- first(L.009map2k2[[k]][["Acc..Len"]][L.009map2k2[[k]][["Max.Score"]] == max(L.009map2k2[[k]][["Max.Score"]])])
  df.009map2k2[, 3] <- rep(400, 28)
  df.009map2k2[, 4] <- df.009map2k2[, 1] / ((df.009map2k2[, 3]) + df.009map2k2[, 2])
  colnames(df.009map2k2) <- c("maxscore", "p.length", "ref.p.length",  "MAP2K2.maxscore.scaled")
}


# MAPK1
df.010mapk1 <- data.frame(1:nrow(pat))
L.010mapk1 <- list()

for (k in 1:nrow(pat)) {
  L.010mapk1[[k]] <- data.frame(L1[[10]][grep(pat[k,], L1[[10]][["Scientific.Name"]]), ])
  df.010mapk1[k,] <- max((data.frame(L1[[10]][grep(pat[k,], L1[[10]][["Scientific.Name"]]), ]))$Max.Score)
  df.010mapk1[k, 2] <- first(L.010mapk1[[k]][["Acc..Len"]][L.010mapk1[[k]][["Max.Score"]] == max(L.010mapk1[[k]][["Max.Score"]])])
  df.010mapk1[, 3] <- rep(364, 28)
  df.010mapk1[, 4] <- df.010mapk1[, 1] / ((df.010mapk1[, 3]) + df.010mapk1[, 2])
  colnames(df.010mapk1) <- c("maxscore", "p.length", "ref.p.length",  "MAPK1.maxscore.scaled")
}


# PLA2G4A
df.013pla2g4a <- data.frame(1:nrow(pat))
L.013pla2g4a <- list()

for (k in 1:nrow(pat)) {
  L.013pla2g4a[[k]] <- data.frame(L1[[11]][grep(pat[k,], L1[[11]][["Scientific.Name"]]), ])
  df.013pla2g4a[k,] <- max((data.frame(L1[[11]][grep(pat[k,], L1[[11]][["Scientific.Name"]]), ]))$Max.Score)
  df.013pla2g4a[k, 2] <- first(L.013pla2g4a[[k]][["Acc..Len"]][L.013pla2g4a[[k]][["Max.Score"]] == max(L.013pla2g4a[[k]][["Max.Score"]])])
  df.013pla2g4a[, 3] <- rep(780, 28)
  df.013pla2g4a[, 4] <- df.013pla2g4a[, 1] / ((df.013pla2g4a[, 3]) + df.013pla2g4a[, 2])
  colnames(df.013pla2g4a) <- c("maxscore", "p.length", "ref.p.length",  "PLA2G4A.maxscore.scaled")
}


# PLA2G4F
df.018pla2g4f <- data.frame(1:nrow(pat))
L.018pla2g4f <- list()

for (k in 1:nrow(pat)) {
  L.018pla2g4f[[k]] <- data.frame(L1[[12]][grep(pat[k,], L1[[12]][["Scientific.Name"]]), ])
  df.018pla2g4f[k,] <- max((data.frame(L1[[12]][grep(pat[k,], L1[[12]][["Scientific.Name"]]), ]))$Max.Score)
  df.018pla2g4f[k, 2] <- first(L.018pla2g4f[[k]][["Acc..Len"]][L.018pla2g4f[[k]][["Max.Score"]] == max(L.018pla2g4f[[k]][["Max.Score"]])])
  df.018pla2g4f[, 3] <- rep(856, 28)
  df.018pla2g4f[, 4] <- df.018pla2g4f[, 1] / ((df.018pla2g4f[, 3]) + df.018pla2g4f[, 2])
  colnames(df.018pla2g4f) <- c("maxscore", "p.length", "ref.p.length",  "PLA2G4F.maxscore.scaled")
}


# PTGS2
df.019ptgs2 <- data.frame(1:nrow(pat))
L.019ptgs2 <- list()

for (k in 1:nrow(pat)) {
  L.019ptgs2[[k]] <- data.frame(L1[[13]][grep(pat[k,], L1[[13]][["Scientific.Name"]]), ])
  df.019ptgs2[k,] <- max((data.frame(L1[[13]][grep(pat[k,], L1[[13]][["Scientific.Name"]]), ]))$Max.Score)
  df.019ptgs2[k, 2] <- first(L.019ptgs2[[k]][["Acc..Len"]][L.019ptgs2[[k]][["Max.Score"]] == max(L.019ptgs2[[k]][["Max.Score"]])])
  df.019ptgs2[, 3] <- rep(731, 28)
  df.019ptgs2[, 4] <- df.019ptgs2[, 1] / ((df.019ptgs2[, 3]) + df.019ptgs2[, 2])
  colnames(df.019ptgs2) <- c("maxscore", "p.length", "ref.p.length",  "PTGS2.maxscore.scaled")
}


# MAP2K5
df.020map2k5 <- data.frame(1:nrow(pat))
L.020map2k5 <- list()

for (k in 1:nrow(pat)) {
  L.020map2k5[[k]] <- data.frame(L1[[14]][grep(pat[k,], L1[[14]][["Scientific.Name"]]), ])
  df.020map2k5[k,] <- max((data.frame(L1[[14]][grep(pat[k,], L1[[14]][["Scientific.Name"]]), ]))$Max.Score)
  df.020map2k5[k, 2] <- first(L.020map2k5[[k]][["Acc..Len"]][L.020map2k5[[k]][["Max.Score"]] == max(L.020map2k5[[k]][["Max.Score"]])])
  df.020map2k5[, 3] <- rep(444, 28)
  df.020map2k5[, 4] <- df.020map2k5[, 1] / ((df.020map2k5[, 3]) + df.020map2k5[, 2])
  colnames(df.020map2k5) <- c("maxscore", "p.length", "ref.p.length",  "MAP2K5.maxscore.scaled")
}


# JUN
df.022jun <- data.frame(1:nrow(pat))
L.022jun <- list()

for (k in 1:nrow(pat)) {
  L.022jun[[k]] <- data.frame(L1[[15]][grep(pat[k,], L1[[15]][["Scientific.Name"]]), ])
  df.022jun[k,] <- max((data.frame(L1[[15]][grep(pat[k,], L1[[15]][["Scientific.Name"]]), ]))$Max.Score)
  df.022jun[k, 2] <- first(L.022jun[[k]][["Acc..Len"]][L.022jun[[k]][["Max.Score"]] == max(L.022jun[[k]][["Max.Score"]])])
  df.022jun[, 3] <- rep(318, 28)
  df.022jun[, 4] <- df.022jun[, 1] / ((df.022jun[, 3]) + df.022jun[, 2])
  colnames(df.022jun) <- c("maxscore", "p.length", "ref.p.length",  "JUN.maxscore.scaled")
}


# FOS
df.023fos <- data.frame(1:nrow(pat))
L.023fos <- list()

for (k in 1:nrow(pat)) {
  L.023fos[[k]] <- data.frame(L1[[16]][grep(pat[k,], L1[[16]][["Scientific.Name"]]), ])
  df.023fos[k,] <- max((data.frame(L1[[16]][grep(pat[k,], L1[[16]][["Scientific.Name"]]), ]))$Max.Score)
  df.023fos[k, 2] <- first(L.023fos[[k]][["Acc..Len"]][L.023fos[[k]][["Max.Score"]] == max(L.023fos[[k]][["Max.Score"]])])
  df.023fos[, 3] <- rep(301, 28)
  df.023fos[, 4] <- df.023fos[, 1] / ((df.023fos[, 3]) + df.023fos[, 2])
  colnames(df.023fos) <- c("maxscore", "p.length", "ref.p.length",  "FOS.maxscore.scaled")
}


# 024MEF2C
df.024mef2c <- data.frame(1:nrow(pat))
L.024mef2c <- list()

for (k in 1:nrow(pat)) {
  L.024mef2c[[k]] <- data.frame(L1[[17]][grep(pat[k,], L1[[17]][["Scientific.Name"]]), ])
  df.024mef2c[k,] <- max((data.frame(L1[[17]][grep(pat[k,], L1[[17]][["Scientific.Name"]]), ]))$Max.Score)
  df.024mef2c[k, 2] <- first(L.024mef2c[[k]][["Acc..Len"]][L.024mef2c[[k]][["Max.Score"]] == max(L.024mef2c[[k]][["Max.Score"]])])
  df.024mef2c[, 3] <- rep(465, 28)
  df.024mef2c[, 4] <- df.024mef2c[, 1] / ((df.024mef2c[, 3]) + df.024mef2c[, 2])
  colnames(df.024mef2c) <- c("maxscore", "p.length", "ref.p.length",  "MEF2C.maxscore.scaled")
}


# 025CCND1
df.025ccnd1 <- data.frame(1:nrow(pat))
L.025ccnd1 <- list()

for (k in 1:nrow(pat)) {
  L.025ccnd1[[k]] <- data.frame(L1[[18]][grep(pat[k,], L1[[18]][["Scientific.Name"]]), ])
  df.025ccnd1[k,] <- max((data.frame(L1[[18]][grep(pat[k,], L1[[18]][["Scientific.Name"]]), ]))$Max.Score)
  df.025ccnd1[k, 2] <- first(L.025ccnd1[[k]][["Acc..Len"]][L.025ccnd1[[k]][["Max.Score"]] == max(L.025ccnd1[[k]][["Max.Score"]])])
  df.025ccnd1[, 3] <- rep(292, 28)
  df.025ccnd1[, 4] <- df.025ccnd1[, 1] / ((df.025ccnd1[, 3]) + df.025ccnd1[, 2])
  colnames(df.025ccnd1) <- c("maxscore", "p.length", "ref.p.length",  "CCND1.maxscore.scaled")
}


# 026ELK1
df.026elk1 <- data.frame(1:nrow(pat))
L.026elk1 <- list()

for (k in 1:nrow(pat)) {
  L.026elk1[[k]] <- data.frame(L1[[19]][grep(pat[k,], L1[[19]][["Scientific.Name"]]), ])
  df.026elk1[k,] <- max((data.frame(L1[[19]][grep(pat[k,], L1[[19]][["Scientific.Name"]]), ]))$Max.Score)
  df.026elk1[k, 2] <- first(L.026elk1[[k]][["Acc..Len"]][L.026elk1[[k]][["Max.Score"]] == max(L.026elk1[[k]][["Max.Score"]])])
  df.026elk1[, 3] <- rep(372, 28)
  df.026elk1[, 4] <- df.026elk1[, 1] / ((df.026elk1[, 3]) + df.026elk1[, 2])
  colnames(df.026elk1) <- c("maxscore", "p.length", "ref.p.length",  "ELK1.maxscore.scaled")
}


# 028RYR2
df.028ryr2 <- data.frame(1:nrow(pat))
L.028ryr2 <- list()

for (k in 1:nrow(pat)) {
  L.028ryr2[[k]] <- data.frame(L1[[20]][grep(pat[k,], L1[[20]][["Scientific.Name"]]), ])
  df.028ryr2[k,] <- max((data.frame(L1[[20]][grep(pat[k,], L1[[20]][["Scientific.Name"]]), ]))$Max.Score)
  df.028ryr2[k, 2] <- first(L.028ryr2[[k]][["Acc..Len"]][L.028ryr2[[k]][["Max.Score"]] == max(L.028ryr2[[k]][["Max.Score"]])])
  df.028ryr2[, 3] <- rep(4953, 28)
  df.028ryr2[, 4] <- df.028ryr2[, 1] / ((df.028ryr2[, 3]) + df.028ryr2[, 2])
  colnames(df.028ryr2) <- c("maxscore", "p.length", "ref.p.length",  "RYR2.maxscore.scaled")
}


# 029RYR3
df.029RYR3 <- data.frame(1:nrow(pat))
L.029RYR3 <- list()

for (k in 1:nrow(pat)) {
  L.029RYR3[[k]] <- data.frame(L1[[21]][grep(pat[k,], L1[[21]][["Scientific.Name"]]), ])
  df.029RYR3[k,] <- max((data.frame(L1[[21]][grep(pat[k,], L1[[21]][["Scientific.Name"]]), ]))$Max.Score)
  df.029RYR3[k, 2] <- first(L.029RYR3[[k]][["Acc..Len"]][L.029RYR3[[k]][["Max.Score"]] == max(L.029RYR3[[k]][["Max.Score"]])])
  df.029RYR3[, 3] <- rep(4881, 28)
  df.029RYR3[, 4] <- df.029RYR3[, 1] / ((df.029RYR3[, 3]) + df.029RYR3[, 2])
  colnames(df.029RYR3) <- c("maxscore", "p.length", "ref.p.length",  "RYR3.maxscore.scaled")
}


# -------------------------------------------- 030CD38
df.030CD38 <- data.frame(1:nrow(pat))
L.030CD38 <- list()

for (k in 1:nrow(pat)) {
  L.030CD38[[k]] <- data.frame(L1[[22]][grep(pat[k,], L1[[22]][["Scientific.Name"]]), ])
  df.030CD38[k,] <- max((data.frame(L1[[22]][grep(pat[k,], L1[[22]][["Scientific.Name"]]), ]))$Max.Score)
  df.030CD38[k, 2] <- first(L.030CD38[[k]][["Acc..Len"]][L.030CD38[[k]][["Max.Score"]] == max(L.030CD38[[k]][["Max.Score"]])])
  df.030CD38[, 3] <- rep(291, 28)
  df.030CD38[, 4] <- df.030CD38[, 1] / ((df.030CD38[, 3]) + df.030CD38[, 2])
  colnames(df.030CD38) <- c("maxscore", "p.length", "ref.p.length",  "CD38.maxscore.scaled")
}


# -------------------------------------------- 032KCNJ2
df.032KCNJ2 <- data.frame(1:nrow(pat))
L.032KCNJ2 <- list()

for (k in 1:nrow(pat)) {
  L.032KCNJ2[[k]] <- data.frame(L1[[23]][grep(pat[k,], L1[[23]][["Scientific.Name"]]), ])
  df.032KCNJ2[k,] <- max((data.frame(L1[[23]][grep(pat[k,], L1[[23]][["Scientific.Name"]]), ]))$Max.Score)
  df.032KCNJ2[k, 2] <- first(L.032KCNJ2[[k]][["Acc..Len"]][L.032KCNJ2[[k]][["Max.Score"]] == max(L.032KCNJ2[[k]][["Max.Score"]])])
  df.032KCNJ2[, 3] <- rep(426, 28)
  df.032KCNJ2[, 4] <- df.032KCNJ2[, 1] / ((df.032KCNJ2[, 3]) + df.032KCNJ2[, 2])
  colnames(df.032KCNJ2) <- c("maxscore", "p.length", "ref.p.length",  "KCNJ2.maxscore.scaled")
}


# -------------------------------------------- 035KCNJ4
df.035KCNJ4 <- data.frame(1:nrow(pat))
L.035KCNJ4 <- list()

for (k in 1:nrow(pat)) {
  L.035KCNJ4[[k]] <- data.frame(L1[[24]][grep(pat[k,], L1[[24]][["Scientific.Name"]]), ])
  df.035KCNJ4[k,] <- max((data.frame(L1[[24]][grep(pat[k,], L1[[24]][["Scientific.Name"]]), ]))$Max.Score)
  df.035KCNJ4[k, 2] <- first(L.035KCNJ4[[k]][["Acc..Len"]][L.035KCNJ4[[k]][["Max.Score"]] == max(L.035KCNJ4[[k]][["Max.Score"]])])
  df.035KCNJ4[, 3] <- rep(426, 28)
  df.035KCNJ4[, 4] <- df.035KCNJ4[, 1] / ((df.035KCNJ4[, 3]) + df.035KCNJ4[, 2])
  colnames(df.035KCNJ4) <- c("maxscore", "p.length", "ref.p.length",  "KCNJ4.maxscore.scaled")
}


# -------------------------------------------- 037PLCB1
df.037PLCB1 <- data.frame(1:nrow(pat))
L.037PLCB1 <- list()

for (k in 1:nrow(pat)) {
  L.037PLCB1[[k]] <- data.frame(L1[[25]][grep(pat[k,], L1[[25]][["Scientific.Name"]]), ])
  df.037PLCB1[k,] <- max((data.frame(L1[[25]][grep(pat[k,], L1[[25]][["Scientific.Name"]]), ]))$Max.Score)
  df.037PLCB1[k, 2] <- first(L.037PLCB1[[k]][["Acc..Len"]][L.037PLCB1[[k]][["Max.Score"]] == max(L.037PLCB1[[k]][["Max.Score"]])])
  df.037PLCB1[, 3] <- rep(1260, 28)
  df.037PLCB1[, 4] <- df.037PLCB1[, 1] / ((df.037PLCB1[, 3]) + df.037PLCB1[, 2])
  colnames(df.037PLCB1) <- c("maxscore", "p.length", "ref.p.length",  "PLCB1.maxscore.scaled")
}


# -------------------------------------------- 040PLCB4
df.040PLCB4 <- data.frame(1:nrow(pat))
L.040PLCB4 <- list()

for (k in 1:nrow(pat)) {
  L.040PLCB4[[k]] <- data.frame(L1[[26]][grep(pat[k,], L1[[26]][["Scientific.Name"]]), ])
  df.040PLCB4[k,] <- max((data.frame(L1[[26]][grep(pat[k,], L1[[26]][["Scientific.Name"]]), ]))$Max.Score)
  df.040PLCB4[k, 2] <- first(L.040PLCB4[[k]][["Acc..Len"]][L.040PLCB4[[k]][["Max.Score"]] == max(L.040PLCB4[[k]][["Max.Score"]])])
  df.040PLCB4[, 3] <- rep(1207, 28)
  df.040PLCB4[, 4] <- df.040PLCB4[, 1] / ((df.040PLCB4[, 3]) + df.040PLCB4[, 2])
  colnames(df.040PLCB4) <- c("maxscore", "p.length", "ref.p.length",  "PLCB4.maxscore.scaled")
}


# -------------------------------------------- 041PRKCA
df.041PRKCA <- data.frame(1:nrow(pat))
L.041PRKCA <- list()

for (k in 1:nrow(pat)) {
  L.041PRKCA[[k]] <- data.frame(L1[[27]][grep(pat[k,], L1[[27]][["Scientific.Name"]]), ])
  df.041PRKCA[k,] <- max((data.frame(L1[[27]][grep(pat[k,], L1[[27]][["Scientific.Name"]]), ]))$Max.Score)
  df.041PRKCA[k, 2] <- first(L.041PRKCA[[k]][["Acc..Len"]][L.041PRKCA[[k]][["Max.Score"]] == max(L.041PRKCA[[k]][["Max.Score"]])])
  df.041PRKCA[, 3] <- rep(664, 28)
  df.041PRKCA[, 4] <- df.041PRKCA[, 1] / ((df.041PRKCA[, 3]) + df.041PRKCA[, 2])
  colnames(df.041PRKCA) <- c("maxscore", "p.length", "ref.p.length",  "PRKCA.maxscore.scaled")
}


# -------------------------------------------- 042PRKCB
df.042PRKCB <- data.frame(1:nrow(pat))
L.042PRKCB <- list()

for (k in 1:nrow(pat)) {
  L.042PRKCB[[k]] <- data.frame(L1[[28]][grep(pat[k,], L1[[28]][["Scientific.Name"]]), ])
  df.042PRKCB[k,] <- max((data.frame(L1[[28]][grep(pat[k,], L1[[28]][["Scientific.Name"]]), ]))$Max.Score)
  df.042PRKCB[k, 2] <- first(L.042PRKCB[[k]][["Acc..Len"]][L.042PRKCB[[k]][["Max.Score"]] == max(L.042PRKCB[[k]][["Max.Score"]])])
  df.042PRKCB[, 3] <- rep(671, 28)
  df.042PRKCB[, 4] <- df.042PRKCB[, 1] / ((df.042PRKCB[, 3]) + df.042PRKCB[, 2])
  colnames(df.042PRKCB) <- c("maxscore", "p.length", "ref.p.length",  "PRKCB.maxscore.scaled")
}


# -------------------------------------------- 044EEF2K
df.044EEF2K <- data.frame(1:nrow(pat))
L.044EEF2K <- list()

for (k in 1:nrow(pat)) {
  L.044EEF2K[[k]] <- data.frame(L1[[29]][grep(pat[k,], L1[[29]][["Scientific.Name"]]), ])
  df.044EEF2K[k,] <- max((data.frame(L1[[29]][grep(pat[k,], L1[[29]][["Scientific.Name"]]), ]))$Max.Score)
  df.044EEF2K[k, 2] <- first(L.044EEF2K[[k]][["Acc..Len"]][L.044EEF2K[[k]][["Max.Score"]] == max(L.044EEF2K[[k]][["Max.Score"]])])
  df.044EEF2K[, 3] <- rep(715, 28)
  df.044EEF2K[, 4] <- df.044EEF2K[, 1] / ((df.044EEF2K[, 3]) + df.044EEF2K[, 2])
  colnames(df.044EEF2K) <- c("maxscore", "p.length", "ref.p.length",  "EEF2K.maxscore.scaled")
}


# -------------------------------------------- 045EEF2
df.045EEF2 <- data.frame(1:nrow(pat))
L.045EEF2 <- list()

for (k in 1:nrow(pat)) {
  L.045EEF2[[k]] <- data.frame(L1[[30]][grep(pat[k,], L1[[30]][["Scientific.Name"]]), ])
  df.045EEF2[k,] <- max((data.frame(L1[[30]][grep(pat[k,], L1[[30]][["Scientific.Name"]]), ]))$Max.Score)
  df.045EEF2[k, 2] <- first(L.045EEF2[[k]][["Acc..Len"]][L.045EEF2[[k]][["Max.Score"]] == max(L.045EEF2[[k]][["Max.Score"]])])
  df.045EEF2[, 3] <- rep(858, 28)
  df.045EEF2[, 4] <- df.045EEF2[, 1] / ((df.045EEF2[, 3]) + df.045EEF2[, 2])
  colnames(df.045EEF2) <- c("maxscore", "p.length", "ref.p.length",  "EEF2.maxscore.scaled")
}


# -------------------------------------------- 046CACNA1C
df.046CACNA1C <- data.frame(1:nrow(pat))
L.046CACNA1C <- list()

for (k in 1:nrow(pat)) {
  L.046CACNA1C[[k]] <- data.frame(L1[[31]][grep(pat[k,], L1[[31]][["Scientific.Name"]]), ])
  df.046CACNA1C[k,] <- max((data.frame(L1[[31]][grep(pat[k,], L1[[31]][["Scientific.Name"]]), ]))$Max.Score)
  df.046CACNA1C[k, 2] <- first(L.046CACNA1C[[k]][["Acc..Len"]][L.046CACNA1C[[k]][["Max.Score"]] == max(L.046CACNA1C[[k]][["Max.Score"]])])
  df.046CACNA1C[, 3] <- rep(2200, 28)
  df.046CACNA1C[, 4] <- df.046CACNA1C[, 1] / ((df.046CACNA1C[, 3]) + df.046CACNA1C[, 2])
  colnames(df.046CACNA1C) <- c("maxscore", "p.length", "ref.p.length",  "CACNA1C.maxscore.scaled")
}


# -------------------------------------------- 047CACNA1D
df.047CACNA1D <- data.frame(1:nrow(pat))
L.047CACNA1D <- list()

for (k in 1:nrow(pat)) {
  L.047CACNA1D[[k]] <- data.frame(L1[[32]][grep(pat[k,], L1[[32]][["Scientific.Name"]]), ])
  df.047CACNA1D[k,] <- max((data.frame(L1[[32]][grep(pat[k,], L1[[32]][["Scientific.Name"]]), ]))$Max.Score)
  df.047CACNA1D[k, 2] <- first(L.047CACNA1D[[k]][["Acc..Len"]][L.047CACNA1D[[k]][["Max.Score"]] == max(L.047CACNA1D[[k]][["Max.Score"]])])
  df.047CACNA1D[, 3] <- rep(2130, 28)
  df.047CACNA1D[, 4] <- df.047CACNA1D[, 1] / ((df.047CACNA1D[, 3]) + df.047CACNA1D[, 2])
  colnames(df.047CACNA1D) <- c("maxscore", "p.length", "ref.p.length",  "CACNA1D.maxscore.scaled")
}


# -------------------------------------------- 050CACNB1
df.050CACNB1 <- data.frame(1:nrow(pat))
L.050CACNB1 <- list()

for (k in 1:nrow(pat)) {
  L.050CACNB1[[k]] <- data.frame(L1[[33]][grep(pat[k,], L1[[33]][["Scientific.Name"]]), ])
  df.050CACNB1[k,] <- max((data.frame(L1[[33]][grep(pat[k,], L1[[33]][["Scientific.Name"]]), ]))$Max.Score)
  df.050CACNB1[k, 2] <- first(L.050CACNB1[[k]][["Acc..Len"]][L.050CACNB1[[k]][["Max.Score"]] == max(L.050CACNB1[[k]][["Max.Score"]])])
  df.050CACNB1[, 3] <- rep(800, 28)
  df.050CACNB1[, 4] <- df.050CACNB1[, 1] / ((df.050CACNB1[, 3]) + df.050CACNB1[, 2])
  colnames(df.050CACNB1) <- c("maxscore", "p.length", "ref.p.length",  "CACNB1.maxscore.scaled")
}


# -------------------------------------------- 051CACNB2
df.051CACNB2 <- data.frame(1:nrow(pat))
L.051CACNB2 <- list()

for (k in 1:nrow(pat)) {
  L.051CACNB2[[k]] <- data.frame(L1[[34]][grep(pat[k,], L1[[34]][["Scientific.Name"]]), ])
  df.051CACNB2[k,] <- max((data.frame(L1[[34]][grep(pat[k,], L1[[34]][["Scientific.Name"]]), ]))$Max.Score)
  df.051CACNB2[k, 2] <- first(L.051CACNB2[[k]][["Acc..Len"]][L.051CACNB2[[k]][["Max.Score"]] == max(L.051CACNB2[[k]][["Max.Score"]])])
  df.051CACNB2[, 3] <- rep(613, 28)
  df.051CACNB2[, 4] <- df.051CACNB2[, 1] / ((df.051CACNB2[, 3]) + df.051CACNB2[, 2])
  colnames(df.051CACNB2) <- c("maxscore", "p.length", "ref.p.length",  "CACNB2.maxscore.scaled")
}


# -------------------------------------------- 053CACNB4
df.053CACNB4 <- data.frame(1:nrow(pat))
L.053CACNB4 <- list()

for (k in 1:nrow(pat)) {
  L.053CACNB4[[k]] <- data.frame(L1[[35]][grep(pat[k,], L1[[35]][["Scientific.Name"]]), ])
  df.053CACNB4[k,] <- max((data.frame(L1[[35]][grep(pat[k,], L1[[35]][["Scientific.Name"]]), ]))$Max.Score)
  df.053CACNB4[k, 2] <- first(L.053CACNB4[[k]][["Acc..Len"]][L.053CACNB4[[k]][["Max.Score"]] == max(L.053CACNB4[[k]][["Max.Score"]])])
  df.053CACNB4[, 3] <- rep(476, 28)
  df.053CACNB4[, 4] <- df.053CACNB4[, 1] / ((df.053CACNB4[, 3]) + df.053CACNB4[, 2])
  colnames(df.053CACNB4) <- c("maxscore", "p.length", "ref.p.length",  "CACNB4.maxscore.scaled")
}


# -------------------------------------------- 054CACNA2D1
df.054CACNA2D1 <- data.frame(1:nrow(pat))
L.054CACNA2D1 <- list()

for (k in 1:nrow(pat)) {
  L.054CACNA2D1[[k]] <- data.frame(L1[[36]][grep(pat[k,], L1[[36]][["Scientific.Name"]]), ])
  df.054CACNA2D1[k,] <- max((data.frame(L1[[36]][grep(pat[k,], L1[[36]][["Scientific.Name"]]), ]))$Max.Score)
  df.054CACNA2D1[k, 2] <- first(L.054CACNA2D1[[k]][["Acc..Len"]][L.054CACNA2D1[[k]][["Max.Score"]] == max(L.054CACNA2D1[[k]][["Max.Score"]])])
  df.054CACNA2D1[, 3] <- rep(1086, 28)
  df.054CACNA2D1[, 4] <- df.054CACNA2D1[, 1] / ((df.054CACNA2D1[, 3]) + df.054CACNA2D1[, 2])
  colnames(df.054CACNA2D1) <- c("maxscore", "p.length", "ref.p.length",  "CACNA2D1.maxscore.scaled")
}


# -------------------------------------------- 056CACNA2D3
df.056CACNA2D3 <- data.frame(1:nrow(pat))
L.056CACNA2D3 <- list()

for (k in 1:nrow(pat)) {
  L.056CACNA2D3[[k]] <- data.frame(L1[[37]][grep(pat[k,], L1[[37]][["Scientific.Name"]]), ])
  df.056CACNA2D3[k,] <- max((data.frame(L1[[37]][grep(pat[k,], L1[[37]][["Scientific.Name"]]), ]))$Max.Score)
  df.056CACNA2D3[k, 2] <- first(L.056CACNA2D3[[k]][["Acc..Len"]][L.056CACNA2D3[[k]][["Max.Score"]] == max(L.056CACNA2D3[[k]][["Max.Score"]])])
  df.056CACNA2D3[, 3] <- rep(1085, 28)
  df.056CACNA2D3[, 4] <- df.056CACNA2D3[, 1] / ((df.056CACNA2D3[, 3]) + df.056CACNA2D3[, 2])
  colnames(df.056CACNA2D3) <- c("maxscore", "p.length", "ref.p.length",  "CACNA2D3.maxscore.scaled")
}


# -------------------------------------------- 057CACNA2D4
df.057CACNA2D4 <- data.frame(1:nrow(pat))
L.057CACNA2D4 <- list()

for (k in 1:nrow(pat)) {
  L.057CACNA2D4[[k]] <- data.frame(L1[[38]][grep(pat[k,], L1[[38]][["Scientific.Name"]]), ])
  df.057CACNA2D4[k,] <- max((data.frame(L1[[38]][grep(pat[k,], L1[[38]][["Scientific.Name"]]), ]))$Max.Score)
  df.057CACNA2D4[k, 2] <- first(L.057CACNA2D4[[k]][["Acc..Len"]][L.057CACNA2D4[[k]][["Max.Score"]] == max(L.057CACNA2D4[[k]][["Max.Score"]])])
  df.057CACNA2D4[, 3] <- rep(1195, 28)
  df.057CACNA2D4[, 4] <- df.057CACNA2D4[, 1] / ((df.057CACNA2D4[, 3]) + df.057CACNA2D4[, 2])
  colnames(df.057CACNA2D4) <- c("maxscore", "p.length", "ref.p.length",  "CACNA2D4.maxscore.scaled")
}


# -------------------------------------------- 058CACNG1
df.058CACNG1 <- data.frame(1:nrow(pat))
L.058CACNG1 <- list()

for (k in 1:nrow(pat)) {
  L.058CACNG1[[k]] <- data.frame(L1[[39]][grep(pat[k,], L1[[39]][["Scientific.Name"]]), ])
  df.058CACNG1[k,] <- max((data.frame(L1[[39]][grep(pat[k,], L1[[39]][["Scientific.Name"]]), ]))$Max.Score)
  df.058CACNG1[k, 2] <- first(L.058CACNG1[[k]][["Acc..Len"]][L.058CACNG1[[k]][["Max.Score"]] == max(L.058CACNG1[[k]][["Max.Score"]])])
  df.058CACNG1[, 3] <- rep(227, 28)
  df.058CACNG1[, 4] <- df.058CACNG1[, 1] / ((df.058CACNG1[, 3]) + df.058CACNG1[, 2])
  colnames(df.058CACNG1) <- c("maxscore", "p.length", "ref.p.length",  "CACNG1.maxscore.scaled")
}


# -------------------------------------------- 059CACNG2
df.059CACNG2 <- data.frame(1:nrow(pat))
L.059CACNG2 <- list()

for (k in 1:nrow(pat)) {
  L.059CACNG2[[k]] <- data.frame(L1[[40]][grep(pat[k,], L1[[40]][["Scientific.Name"]]), ])
  df.059CACNG2[k,] <- max((data.frame(L1[[40]][grep(pat[k,], L1[[40]][["Scientific.Name"]]), ]))$Max.Score)
  df.059CACNG2[k, 2] <- first(L.059CACNG2[[k]][["Acc..Len"]][L.059CACNG2[[k]][["Max.Score"]] == max(L.059CACNG2[[k]][["Max.Score"]])])
  df.059CACNG2[, 3] <- rep(322, 28)
  df.059CACNG2[, 4] <- df.059CACNG2[, 1] / ((df.059CACNG2[, 3]) + df.059CACNG2[, 2])
  colnames(df.059CACNG2) <- c("maxscore", "p.length", "ref.p.length",  "CACNG2.maxscore.scaled")
}


# -------------------------------------------- 060CACNG3
df.060CACNG3 <- data.frame(1:nrow(pat))
L.060CACNG3 <- list()

for (k in 1:nrow(pat)) {
  L.060CACNG3[[k]] <- data.frame(L1[[41]][grep(pat[k,], L1[[41]][["Scientific.Name"]]), ])
  df.060CACNG3[k,] <- max((data.frame(L1[[41]][grep(pat[k,], L1[[41]][["Scientific.Name"]]), ]))$Max.Score)
  df.060CACNG3[k, 2] <- first(L.060CACNG3[[k]][["Acc..Len"]][L.060CACNG3[[k]][["Max.Score"]] == max(L.060CACNG3[[k]][["Max.Score"]])])
  df.060CACNG3[, 3] <- rep(317, 28)
  df.060CACNG3[, 4] <- df.060CACNG3[, 1] / ((df.060CACNG3[, 3]) + df.060CACNG3[, 2])
  colnames(df.060CACNG3) <- c("maxscore", "p.length", "ref.p.length",  "CACNG3.maxscore.scaled")
}


# -------------------------------------------- 061CACNG4
df.061CACNG4 <- data.frame(1:nrow(pat))
L.061CACNG4 <- list()

for (k in 1:nrow(pat)) {
  L.061CACNG4[[k]] <- data.frame(L1[[42]][grep(pat[k,], L1[[42]][["Scientific.Name"]]), ])
  df.061CACNG4[k,] <- max((data.frame(L1[[42]][grep(pat[k,], L1[[42]][["Scientific.Name"]]), ]))$Max.Score)
  df.061CACNG4[k, 2] <- first(L.061CACNG4[[k]][["Acc..Len"]][L.061CACNG4[[k]][["Max.Score"]] == max(L.061CACNG4[[k]][["Max.Score"]])])
  df.061CACNG4[, 3] <- rep(328, 28)
  df.061CACNG4[, 4] <- df.061CACNG4[, 1] / ((df.061CACNG4[, 3]) + df.061CACNG4[, 2])
  colnames(df.061CACNG4) <- c("maxscore", "p.length", "ref.p.length",  "CACNG4.maxscore.scaled")
}


# -------------------------------------------- 062CACNG5
df.062CACNG5 <- data.frame(1:nrow(pat))
L.062CACNG5 <- list()

for (k in 1:nrow(pat)) {
  L.062CACNG5[[k]] <- data.frame(L1[[43]][grep(pat[k,], L1[[43]][["Scientific.Name"]]), ])
  df.062CACNG5[k,] <- max((data.frame(L1[[43]][grep(pat[k,], L1[[43]][["Scientific.Name"]]), ]))$Max.Score)
  df.062CACNG5[k, 2] <- first(L.062CACNG5[[k]][["Acc..Len"]][L.062CACNG5[[k]][["Max.Score"]] == max(L.062CACNG5[[k]][["Max.Score"]])])
  df.062CACNG5[, 3] <- rep(279, 28)
  df.062CACNG5[, 4] <- df.062CACNG5[, 1] / ((df.062CACNG5[, 3]) + df.062CACNG5[, 2])
  colnames(df.062CACNG5) <- c("maxscore", "p.length", "ref.p.length",  "CACNG5.maxscore.scaled")
}


# -------------------------------------------- 066ITPR1
df.066ITPR1 <- data.frame(1:nrow(pat))
L.066ITPR1 <- list()

for (k in 1:nrow(pat)) {
  L.066ITPR1[[k]] <- data.frame(L1[[44]][grep(pat[k,], L1[[44]][["Scientific.Name"]]), ])
  df.066ITPR1[k,] <- max((data.frame(L1[[44]][grep(pat[k,], L1[[44]][["Scientific.Name"]]), ]))$Max.Score)
  df.066ITPR1[k, 2] <- first(L.066ITPR1[[k]][["Acc..Len"]][L.066ITPR1[[k]][["Max.Score"]] == max(L.066ITPR1[[k]][["Max.Score"]])])
  df.066ITPR1[, 3] <- rep(2757, 28)
  df.066ITPR1[, 4] <- df.066ITPR1[, 1] / ((df.066ITPR1[, 3]) + df.066ITPR1[, 2])
  colnames(df.066ITPR1) <- c("maxscore", "p.length", "ref.p.length",  "ITPR1.maxscore.scaled")
}


# -------------------------------------------- 067ITPR2
df.067ITPR2 <- data.frame(1:nrow(pat))
L.067ITPR2 <- list()

for (k in 1:nrow(pat)) {
  L.067ITPR2[[k]] <- data.frame(L1[[45]][grep(pat[k,], L1[[45]][["Scientific.Name"]]), ])
  df.067ITPR2[k,] <- max((data.frame(L1[[45]][grep(pat[k,], L1[[45]][["Scientific.Name"]]), ]))$Max.Score)
  df.067ITPR2[k, 2] <- first(L.067ITPR2[[k]][["Acc..Len"]][L.067ITPR2[[k]][["Max.Score"]] == max(L.067ITPR2[[k]][["Max.Score"]])])
  df.067ITPR2[, 3] <- rep(2691, 28)
  df.067ITPR2[, 4] <- df.067ITPR2[, 1] / ((df.067ITPR2[, 3]) + df.067ITPR2[, 2])
  colnames(df.067ITPR2) <- c("maxscore", "p.length", "ref.p.length",  "ITPR2.maxscore.scaled")
}


# -------------------------------------------- 070CALM2
df.070CALM2 <- data.frame(1:nrow(pat))
L.070CALM2 <- list()

for (k in 1:nrow(pat)) {
  L.070CALM2[[k]] <- data.frame(L1[[46]][grep(pat[k,], L1[[46]][["Scientific.Name"]]), ])
  df.070CALM2[k,] <- max((data.frame(L1[[46]][grep(pat[k,], L1[[46]][["Scientific.Name"]]), ]))$Max.Score)
  df.070CALM2[k, 2] <- first(L.070CALM2[[k]][["Acc..Len"]][L.070CALM2[[k]][["Max.Score"]] == max(L.070CALM2[[k]][["Max.Score"]])])
  df.070CALM2[, 3] <- rep(149, 28)
  df.070CALM2[, 4] <- df.070CALM2[, 1] / ((df.070CALM2[, 3]) + df.070CALM2[, 2])
  colnames(df.070CALM2) <- c("maxscore", "p.length", "ref.p.length",  "CALM2.maxscore.scaled")
}


# -------------------------------------------- 072CALM1
df.072CALM1 <- data.frame(1:nrow(pat))
L.072CALM1 <- list()

for (k in 1:nrow(pat)) {
  L.072CALM1[[k]] <- data.frame(L1[[47]][grep(pat[k,], L1[[47]][["Scientific.Name"]]), ])
  df.072CALM1[k,] <- max((data.frame(L1[[47]][grep(pat[k,], L1[[47]][["Scientific.Name"]]), ]))$Max.Score)
  df.072CALM1[k, 2] <- first(L.072CALM1[[k]][["Acc..Len"]][L.072CALM1[[k]][["Max.Score"]] == max(L.072CALM1[[k]][["Max.Score"]])])
  df.072CALM1[, 3] <- rep(149, 28)
  df.072CALM1[, 4] <- df.072CALM1[, 1] / ((df.072CALM1[, 3]) + df.072CALM1[, 2])
  colnames(df.072CALM1) <- c("maxscore", "p.length", "ref.p.length",  "CALM1.maxscore.scaled")
}



# -------------------------------------------- 073CALML6
df.073CALML6 <- data.frame(1:nrow(pat))
L.073CALML6 <- list()

for (k in 1:nrow(pat)) {
  L.073CALML6[[k]] <- data.frame(L1[[48]][grep(pat[k,], L1[[48]][["Scientific.Name"]]), ])
  df.073CALML6[k,] <- max((data.frame(L1[[48]][grep(pat[k,], L1[[48]][["Scientific.Name"]]), ]))$Max.Score)
  df.073CALML6[k, 2] <- first(L.073CALML6[[k]][["Acc..Len"]][L.073CALML6[[k]][["Max.Score"]] == max(L.073CALML6[[k]][["Max.Score"]])])
  df.073CALML6[, 3] <- rep(226, 28)
  df.073CALML6[, 4] <- df.073CALML6[, 1] / ((df.073CALML6[, 3]) + df.073CALML6[, 2])
  colnames(df.073CALML6) <- c("maxscore", "p.length", "ref.p.length",  "CALML6.maxscore.scaled")
}


# -------------------------------------------- 075CALML4
df.075CALML4 <- data.frame(1:nrow(pat))
L.075CALML4 <- list()

for (k in 1:nrow(pat)) {
  L.075CALML4[[k]] <- data.frame(L1[[49]][grep(pat[k,], L1[[49]][["Scientific.Name"]]), ])
  df.075CALML4[k,] <- max((data.frame(L1[[49]][grep(pat[k,], L1[[49]][["Scientific.Name"]]), ]))$Max.Score)
  df.075CALML4[k, 2] <- first(L.075CALML4[[k]][["Acc..Len"]][L.075CALML4[[k]][["Max.Score"]] == max(L.075CALML4[[k]][["Max.Score"]])])
  df.075CALML4[, 3] <- rep(153, 28)
  df.075CALML4[, 4] <- df.075CALML4[, 1] / ((df.075CALML4[, 3]) + df.075CALML4[, 2])
  colnames(df.075CALML4) <- c("maxscore", "p.length", "ref.p.length",  "CALML4.maxscore.scaled")
}


# -------------------------------------------- 076PPP3CA
df.076PPP3CA <- data.frame(1:nrow(pat))
L.076PPP3CA <- list()

for (k in 1:nrow(pat)) {
  L.076PPP3CA[[k]] <- data.frame(L1[[50]][grep(pat[k,], L1[[50]][["Scientific.Name"]]), ])
  df.076PPP3CA[k,] <- max((data.frame(L1[[50]][grep(pat[k,], L1[[50]][["Scientific.Name"]]), ]))$Max.Score)
  df.076PPP3CA[k, 2] <- first(L.076PPP3CA[[k]][["Acc..Len"]][L.076PPP3CA[[k]][["Max.Score"]] == max(L.076PPP3CA[[k]][["Max.Score"]])])
  df.076PPP3CA[, 3] <- rep(515, 28)
  df.076PPP3CA[, 4] <- df.076PPP3CA[, 1] / ((df.076PPP3CA[, 3]) + df.076PPP3CA[, 2])
  colnames(df.076PPP3CA) <- c("maxscore", "p.length", "ref.p.length",  "PPP3CA.maxscore.scaled")
}


# -------------------------------------------- 077PPP3CB
df.077PPP3CB <- data.frame(1:nrow(pat))
L.077PPP3CB <- list()

for (k in 1:nrow(pat)) {
  L.077PPP3CB[[k]] <- data.frame(L1[[51]][grep(pat[k,], L1[[51]][["Scientific.Name"]]), ])
  df.077PPP3CB[k,] <- max((data.frame(L1[[51]][grep(pat[k,], L1[[51]][["Scientific.Name"]]), ]))$Max.Score)
  df.077PPP3CB[k, 2] <- first(L.077PPP3CB[[k]][["Acc..Len"]][L.077PPP3CB[[k]][["Max.Score"]] == max(L.077PPP3CB[[k]][["Max.Score"]])])
  df.077PPP3CB[, 3] <- rep(515, 28)
  df.077PPP3CB[, 4] <- df.077PPP3CB[, 1] / ((df.077PPP3CB[, 3]) + df.077PPP3CB[, 2])
  colnames(df.077PPP3CB) <- c("maxscore", "p.length", "ref.p.length",  "PPP3CB.maxscore.scaled")
}


# -------------------------------------------- 078PPP3CC
df.078PPP3CC <- data.frame(1:nrow(pat))
L.078PPP3CC <- list()

for (k in 1:nrow(pat)) {
  L.078PPP3CC[[k]] <- data.frame(L1[[52]][grep(pat[k,], L1[[52]][["Scientific.Name"]]), ])
  df.078PPP3CC[k,] <- max((data.frame(L1[[52]][grep(pat[k,], L1[[52]][["Scientific.Name"]]), ]))$Max.Score)
  df.078PPP3CC[k, 2] <- first(L.078PPP3CC[[k]][["Acc..Len"]][L.078PPP3CC[[k]][["Max.Score"]] == max(L.078PPP3CC[[k]][["Max.Score"]])])
  df.078PPP3CC[, 3] <- rep(581, 28)
  df.078PPP3CC[, 4] <- df.078PPP3CC[, 1] / ((df.078PPP3CC[, 3]) + df.078PPP3CC[, 2])
  colnames(df.078PPP3CC) <- c("maxscore", "p.length", "ref.p.length",  "PPP3CC.maxscore.scaled")
}


# -------------------------------------------- 079PPP3R1
df.079PPP3R1 <- data.frame(1:nrow(pat))
L.079PPP3R1 <- list()

for (k in 1:nrow(pat)) {
  L.079PPP3R1[[k]] <- data.frame(L1[[53]][grep(pat[k,], L1[[53]][["Scientific.Name"]]), ])
  df.079PPP3R1[k,] <- max((data.frame(L1[[53]][grep(pat[k,], L1[[53]][["Scientific.Name"]]), ]))$Max.Score)
  df.079PPP3R1[k, 2] <- first(L.079PPP3R1[[k]][["Acc..Len"]][L.079PPP3R1[[k]][["Max.Score"]] == max(L.079PPP3R1[[k]][["Max.Score"]])])
  df.079PPP3R1[, 3] <- rep(177, 28)
  df.079PPP3R1[, 4] <- df.079PPP3R1[, 1] / ((df.079PPP3R1[, 3]) + df.079PPP3R1[, 2])
  colnames(df.079PPP3R1) <- c("maxscore", "p.length", "ref.p.length",  "PPP3R1.maxscore.scaled")
}


# -------------------------------------------- 081NFATC1
df.081NFATC1 <- data.frame(1:nrow(pat))
L.081NFATC1 <- list()

for (k in 1:nrow(pat)) {
  L.081NFATC1[[k]] <- data.frame(L1[[54]][grep(pat[k,], L1[[54]][["Scientific.Name"]]), ])
  df.081NFATC1[k,] <- max((data.frame(L1[[54]][grep(pat[k,], L1[[54]][["Scientific.Name"]]), ]))$Max.Score)
  df.081NFATC1[k, 2] <- first(L.081NFATC1[[k]][["Acc..Len"]][L.081NFATC1[[k]][["Max.Score"]] == max(L.081NFATC1[[k]][["Max.Score"]])])
  df.081NFATC1[, 3] <- rep(922, 28)
  df.081NFATC1[, 4] <- df.081NFATC1[, 1] / ((df.081NFATC1[, 3]) + df.081NFATC1[, 2])
  colnames(df.081NFATC1) <- c("maxscore", "p.length", "ref.p.length",  "NFATC1.maxscore.scaled")
}


# -------------------------------------------- 082NFATC2
df.082NFATC2 <- data.frame(1:nrow(pat))
L.082NFATC2 <- list()

for (k in 1:nrow(pat)) {
  L.082NFATC2[[k]] <- data.frame(L1[[55]][grep(pat[k,], L1[[55]][["Scientific.Name"]]), ])
  df.082NFATC2[k,] <- max((data.frame(L1[[55]][grep(pat[k,], L1[[55]][["Scientific.Name"]]), ]))$Max.Score)
  df.082NFATC2[k, 2] <- first(L.082NFATC2[[k]][["Acc..Len"]][L.082NFATC2[[k]][["Max.Score"]] == max(L.082NFATC2[[k]][["Max.Score"]])])
  df.082NFATC2[, 3] <- rep(960, 28)
  df.082NFATC2[, 4] <- df.082NFATC2[, 1] / ((df.082NFATC2[, 3]) + df.082NFATC2[, 2])
  colnames(df.082NFATC2) <- c("maxscore", "p.length", "ref.p.length",  "NFATC2.maxscore.scaled")
}


# -------------------------------------------- 083NFATC3
df.083NFATC3 <- data.frame(1:nrow(pat))
L.083NFATC3 <- list()

for (k in 1:nrow(pat)) {
  L.083NFATC3[[k]] <- data.frame(L1[[56]][grep(pat[k,], L1[[56]][["Scientific.Name"]]), ])
  df.083NFATC3[k,] <- max((data.frame(L1[[56]][grep(pat[k,], L1[[56]][["Scientific.Name"]]), ]))$Max.Score)
  df.083NFATC3[k, 2] <- first(L.083NFATC3[[k]][["Acc..Len"]][L.083NFATC3[[k]][["Max.Score"]] == max(L.083NFATC3[[k]][["Max.Score"]])])
  df.083NFATC3[, 3] <- rep(981, 28)
  df.083NFATC3[, 4] <- df.083NFATC3[, 1] / ((df.083NFATC3[, 3]) + df.083NFATC3[, 2])
  colnames(df.083NFATC3) <- c("maxscore", "p.length", "ref.p.length",  "NFATC3.maxscore.scaled")
}


# -------------------------------------------- 087CAMK2
df.087CAMK2 <- data.frame(1:nrow(pat))
L.087CAMK2 <- list()

for (k in 1:nrow(pat)) {
  L.087CAMK2[[k]] <- data.frame(L1[[57]][grep(pat[k,], L1[[57]][["Scientific.Name"]]), ])
  df.087CAMK2[k,] <- max((data.frame(L1[[57]][grep(pat[k,], L1[[57]][["Scientific.Name"]]), ]))$Max.Score)
  df.087CAMK2[k, 2] <- first(L.087CAMK2[[k]][["Acc..Len"]][L.087CAMK2[[k]][["Max.Score"]] == max(L.087CAMK2[[k]][["Max.Score"]])])
  df.087CAMK2[, 3] <- rep(557, 28)
  df.087CAMK2[, 4] <- df.087CAMK2[, 1] / ((df.087CAMK2[, 3]) + df.087CAMK2[, 2])
  colnames(df.087CAMK2) <- c("maxscore", "p.length", "ref.p.length",  "CAMK2.maxscore.scaled")
}


# -------------------------------------------- 088PRKAA1
df.088PRKAA1 <- data.frame(1:nrow(pat))
L.088PRKAA1 <- list()

for (k in 1:nrow(pat)) {
  L.088PRKAA1[[k]] <- data.frame(L1[[58]][grep(pat[k,], L1[[58]][["Scientific.Name"]]), ])
  df.088PRKAA1[k,] <- max((data.frame(L1[[58]][grep(pat[k,], L1[[58]][["Scientific.Name"]]), ]))$Max.Score)
  df.088PRKAA1[k, 2] <- first(L.088PRKAA1[[k]][["Acc..Len"]][L.088PRKAA1[[k]][["Max.Score"]] == max(L.088PRKAA1[[k]][["Max.Score"]])])
  df.088PRKAA1[, 3] <- rep(558, 28)
  df.088PRKAA1[, 4] <- df.088PRKAA1[, 1] / ((df.088PRKAA1[, 3]) + df.088PRKAA1[, 2])
  colnames(df.088PRKAA1) <- c("maxscore", "p.length", "ref.p.length",  "PRKAA1.maxscore.scaled")
}


# -------------------------------------------- 089PRKAA2
df.089PRKAA2 <- data.frame(1:nrow(pat))
L.089PRKAA2 <- list()

for (k in 1:nrow(pat)) {
  L.089PRKAA2[[k]] <- data.frame(L1[[59]][grep(pat[k,], L1[[59]][["Scientific.Name"]]), ])
  df.089PRKAA2[k,] <- max((data.frame(L1[[59]][grep(pat[k,], L1[[59]][["Scientific.Name"]]), ]))$Max.Score)
  df.089PRKAA2[k, 2] <- first(L.089PRKAA2[[k]][["Acc..Len"]][L.089PRKAA2[[k]][["Max.Score"]] == max(L.089PRKAA2[[k]][["Max.Score"]])])
  df.089PRKAA2[, 3] <- rep(553, 28)
  df.089PRKAA2[, 4] <- df.089PRKAA2[, 1] / ((df.089PRKAA2[, 3]) + df.089PRKAA2[, 2])
  colnames(df.089PRKAA2) <- c("maxscore", "p.length", "ref.p.length",  "PRKAA2.maxscore.scaled")
}


# -------------------------------------------- 090PRKAB1
df.090PRKAB1 <- data.frame(1:nrow(pat))
L.090PRKAB1 <- list()

for (k in 1:nrow(pat)) {
  L.090PRKAB1[[k]] <- data.frame(L1[[60]][grep(pat[k,], L1[[60]][["Scientific.Name"]]), ])
  df.090PRKAB1[k,] <- max((data.frame(L1[[60]][grep(pat[k,], L1[[60]][["Scientific.Name"]]), ]))$Max.Score)
  df.090PRKAB1[k, 2] <- first(L.090PRKAB1[[k]][["Acc..Len"]][L.090PRKAB1[[k]][["Max.Score"]] == max(L.090PRKAB1[[k]][["Max.Score"]])])
  df.090PRKAB1[, 3] <- rep(274, 28)
  df.090PRKAB1[, 4] <- df.090PRKAB1[, 1] / ((df.090PRKAB1[, 3]) + df.090PRKAB1[, 2])
  colnames(df.090PRKAB1) <- c("maxscore", "p.length", "ref.p.length",  "PRKAB1.maxscore.scaled")
}


# -------------------------------------------- 092PRKAG1
df.092PRKAG1 <- data.frame(1:nrow(pat))
L.092PRKAG1 <- list()

for (k in 1:nrow(pat)) {
  L.092PRKAG1[[k]] <- data.frame(L1[[61]][grep(pat[k,], L1[[61]][["Scientific.Name"]]), ])
  df.092PRKAG1[k,] <- max((data.frame(L1[[61]][grep(pat[k,], L1[[61]][["Scientific.Name"]]), ]))$Max.Score)
  df.092PRKAG1[k, 2] <- first(L.092PRKAG1[[k]][["Acc..Len"]][L.092PRKAG1[[k]][["Max.Score"]] == max(L.092PRKAG1[[k]][["Max.Score"]])])
  df.092PRKAG1[, 3] <- rep(304, 28)
  df.092PRKAG1[, 4] <- df.092PRKAG1[, 1] / ((df.092PRKAG1[, 3]) + df.092PRKAG1[, 2])
  colnames(df.092PRKAG1) <- c("maxscore", "p.length", "ref.p.length",  "PRKAG1.maxscore.scaled")
}


# -------------------------------------------- 094PRKAG2
df.094PRKAG2 <- data.frame(1:nrow(pat))
L.094PRKAG2 <- list()

for (k in 1:nrow(pat)) {
  L.094PRKAG2[[k]] <- data.frame(L1[[62]][grep(pat[k,], L1[[62]][["Scientific.Name"]]), ])
  df.094PRKAG2[k,] <- max((data.frame(L1[[62]][grep(pat[k,], L1[[62]][["Scientific.Name"]]), ]))$Max.Score)
  df.094PRKAG2[k, 2] <- first(L.094PRKAG2[[k]][["Acc..Len"]][L.094PRKAG2[[k]][["Max.Score"]] == max(L.094PRKAG2[[k]][["Max.Score"]])])
  df.094PRKAG2[, 3] <- rep(551, 28)
  df.094PRKAG2[, 4] <- df.094PRKAG2[, 1] / ((df.094PRKAG2[, 3]) + df.094PRKAG2[, 2])
  colnames(df.094PRKAG2) <- c("maxscore", "p.length", "ref.p.length",  "PRKAG2.maxscore.scaled")
}


# -------------------------------------------- 095CAMK1D
df.095CAMK1D <- data.frame(1:nrow(pat))
L.095CAMK1D <- list()

for (k in 1:nrow(pat)) {
  L.095CAMK1D[[k]] <- data.frame(L1[[63]][grep(pat[k,], L1[[63]][["Scientific.Name"]]), ])
  df.095CAMK1D[k,] <- max((data.frame(L1[[63]][grep(pat[k,], L1[[63]][["Scientific.Name"]]), ]))$Max.Score)
  df.095CAMK1D[k, 2] <- first(L.095CAMK1D[[k]][["Acc..Len"]][L.095CAMK1D[[k]][["Max.Score"]] == max(L.095CAMK1D[[k]][["Max.Score"]])])
  df.095CAMK1D[, 3] <- rep(392, 28)
  df.095CAMK1D[, 4] <- df.095CAMK1D[, 1] / ((df.095CAMK1D[, 3]) + df.095CAMK1D[, 2])
  colnames(df.095CAMK1D) <- c("maxscore", "p.length", "ref.p.length",  "CAMK1D.maxscore.scaled")
}


# -------------------------------------------- 096CAMK1G
df.096CAMK1G <- data.frame(1:nrow(pat))
L.096CAMK1G <- list()

for (k in 1:nrow(pat)) {
  L.096CAMK1G[[k]] <- data.frame(L1[[64]][grep(pat[k,], L1[[64]][["Scientific.Name"]]), ])
  df.096CAMK1G[k,] <- max((data.frame(L1[[64]][grep(pat[k,], L1[[64]][["Scientific.Name"]]), ]))$Max.Score)
  df.096CAMK1G[k, 2] <- first(L.096CAMK1G[[k]][["Acc..Len"]][L.096CAMK1G[[k]][["Max.Score"]] == max(L.096CAMK1G[[k]][["Max.Score"]])])
  df.096CAMK1G[, 3] <- rep(468, 28)
  df.096CAMK1G[, 4] <- df.096CAMK1G[, 1] / ((df.096CAMK1G[, 3]) + df.096CAMK1G[, 2])
  colnames(df.096CAMK1G) <- c("maxscore", "p.length", "ref.p.length",  "CAMK1G.maxscore.scaled")
}


# -------------------------------------------- 097CAMK1
df.097CAMK1 <- data.frame(1:nrow(pat))
L.097CAMK1 <- list()

for (k in 1:nrow(pat)) {
  L.097CAMK1[[k]] <- data.frame(L1[[65]][grep(pat[k,], L1[[65]][["Scientific.Name"]]), ])
  df.097CAMK1[k,] <- max((data.frame(L1[[65]][grep(pat[k,], L1[[65]][["Scientific.Name"]]), ]))$Max.Score)
  df.097CAMK1[k, 2] <- first(L.097CAMK1[[k]][["Acc..Len"]][L.097CAMK1[[k]][["Max.Score"]] == max(L.097CAMK1[[k]][["Max.Score"]])])
  df.097CAMK1[, 3] <- rep(386, 28)
  df.097CAMK1[, 4] <- df.097CAMK1[, 1] / ((df.097CAMK1[, 3]) + df.097CAMK1[, 2])
  colnames(df.097CAMK1) <- c("maxscore", "p.length", "ref.p.length",  "CAMK1.maxscore.scaled")
}


# -------------------------------------------- 098CAMK2A
df.098CAMK2A <- data.frame(1:nrow(pat))
L.098CAMK2A <- list()

for (k in 1:nrow(pat)) {
  L.098CAMK2A[[k]] <- data.frame(L1[[66]][grep(pat[k,], L1[[66]][["Scientific.Name"]]), ])
  df.098CAMK2A[k,] <- max((data.frame(L1[[66]][grep(pat[k,], L1[[66]][["Scientific.Name"]]), ]))$Max.Score)
  df.098CAMK2A[k, 2] <- first(L.098CAMK2A[[k]][["Acc..Len"]][L.098CAMK2A[[k]][["Max.Score"]] == max(L.098CAMK2A[[k]][["Max.Score"]])])
  df.098CAMK2A[, 3] <- rep(478, 28)
  df.098CAMK2A[, 4] <- df.098CAMK2A[, 1] / ((df.098CAMK2A[, 3]) + df.098CAMK2A[, 2])
  colnames(df.098CAMK2A) <- c("maxscore", "p.length", "ref.p.length",  "CAMK2A.maxscore.scaled")
}


# -------------------------------------------- 099CAMK2D
df.099CAMK2D <- data.frame(1:nrow(pat))
L.099CAMK2D <- list()

for (k in 1:nrow(pat)) {
  L.099CAMK2D[[k]] <- data.frame(L1[[67]][grep(pat[k,], L1[[67]][["Scientific.Name"]]), ])
  df.099CAMK2D[k,] <- max((data.frame(L1[[67]][grep(pat[k,], L1[[67]][["Scientific.Name"]]), ]))$Max.Score)
  df.099CAMK2D[k, 2] <- first(L.099CAMK2D[[k]][["Acc..Len"]][L.099CAMK2D[[k]][["Max.Score"]] == max(L.099CAMK2D[[k]][["Max.Score"]])])
  df.099CAMK2D[, 3] <- rep(494, 28)
  df.099CAMK2D[, 4] <- df.099CAMK2D[, 1] / ((df.099CAMK2D[, 3]) + df.099CAMK2D[, 2])
  colnames(df.099CAMK2D) <- c("maxscore", "p.length", "ref.p.length",  "CAMK2D.maxscore.scaled")
}


# -------------------------------------------- 101CAMK2G
df.101CAMK2G <- data.frame(1:nrow(pat))
L.101CAMK2G <- list()

for (k in 1:nrow(pat)) {
  L.101CAMK2G[[k]] <- data.frame(L1[[68]][grep(pat[k,], L1[[68]][["Scientific.Name"]]), ])
  df.101CAMK2G[k,] <- max((data.frame(L1[[68]][grep(pat[k,], L1[[68]][["Scientific.Name"]]), ]))$Max.Score)
  df.101CAMK2G[k, 2] <- first(L.101CAMK2G[[k]][["Acc..Len"]][L.101CAMK2G[[k]][["Max.Score"]] == max(L.101CAMK2G[[k]][["Max.Score"]])])
  df.101CAMK2G[, 3] <- rep(546, 28)
  df.101CAMK2G[, 4] <- df.101CAMK2G[, 1] / ((df.101CAMK2G[, 3]) + df.101CAMK2G[, 2])
  colnames(df.101CAMK2G) <- c("maxscore", "p.length", "ref.p.length",  "CAMK2G.maxscore.scaled")
}


# -------------------------------------------- 102CAMK4
df.102CAMK4 <- data.frame(1:nrow(pat))
L.102CAMK4 <- list()

for (k in 1:nrow(pat)) {
  L.102CAMK4[[k]] <- data.frame(L1[[69]][grep(pat[k,], L1[[69]][["Scientific.Name"]]), ])
  df.102CAMK4[k,] <- max((data.frame(L1[[69]][grep(pat[k,], L1[[69]][["Scientific.Name"]]), ]))$Max.Score)
  df.102CAMK4[k, 2] <- first(L.102CAMK4[[k]][["Acc..Len"]][L.102CAMK4[[k]][["Max.Score"]] == max(L.102CAMK4[[k]][["Max.Score"]])])
  df.102CAMK4[, 3] <- rep(428, 28)
  df.102CAMK4[, 4] <- df.102CAMK4[, 1] / ((df.102CAMK4[, 3]) + df.102CAMK4[, 2])
  colnames(df.102CAMK4) <- c("maxscore", "p.length", "ref.p.length",  "CAMK4.maxscore.scaled")
}


# -------------------------------------------- 104GUCY1A2
df.104GUCY1A2 <- data.frame(1:nrow(pat))
L.104GUCY1A2 <- list()

for (k in 1:nrow(pat)) {
  L.104GUCY1A2[[k]] <- data.frame(L1[[70]][grep(pat[k,], L1[[70]][["Scientific.Name"]]), ])
  df.104GUCY1A2[k,] <- max((data.frame(L1[[70]][grep(pat[k,], L1[[70]][["Scientific.Name"]]), ]))$Max.Score)
  df.104GUCY1A2[k, 2] <- first(L.104GUCY1A2[[k]][["Acc..Len"]][L.104GUCY1A2[[k]][["Max.Score"]] == max(L.104GUCY1A2[[k]][["Max.Score"]])])
  df.104GUCY1A2[, 3] <- rep(702, 28)
  df.104GUCY1A2[, 4] <- df.104GUCY1A2[, 1] / ((df.104GUCY1A2[, 3]) + df.104GUCY1A2[, 2])
  colnames(df.104GUCY1A2) <- c("maxscore", "p.length", "ref.p.length",  "GUCY1A2.maxscore.scaled")
}


# -------------------------------------------- 106GUCY1B1
df.106GUCY1B1 <- data.frame(1:nrow(pat))
L.106GUCY1B1 <- list()

for (k in 1:nrow(pat)) {
  L.106GUCY1B1[[k]] <- data.frame(L1[[71]][grep(pat[k,], L1[[71]][["Scientific.Name"]]), ])
  df.106GUCY1B1[k,] <- max((data.frame(L1[[71]][grep(pat[k,], L1[[71]][["Scientific.Name"]]), ]))$Max.Score)
  df.106GUCY1B1[k, 2] <- first(L.106GUCY1B1[[k]][["Acc..Len"]][L.106GUCY1B1[[k]][["Max.Score"]] == max(L.106GUCY1B1[[k]][["Max.Score"]])])
  df.106GUCY1B1[, 3] <- rep(502, 28)
  df.106GUCY1B1[, 4] <- df.106GUCY1B1[, 1] / ((df.106GUCY1B1[, 3]) + df.106GUCY1B1[, 2])
  colnames(df.106GUCY1B1) <- c("maxscore", "p.length", "ref.p.length",  "GUCY1B1.maxscore.scaled")
}


# -------------------------------------------- 108NPR1
df.108NPR1 <- data.frame(1:nrow(pat))
L.108NPR1 <- list()

for (k in 1:nrow(pat)) {
  L.108NPR1[[k]] <- data.frame(L1[[72]][grep(pat[k,], L1[[72]][["Scientific.Name"]]), ])
  df.108NPR1[k,] <- max((data.frame(L1[[72]][grep(pat[k,], L1[[72]][["Scientific.Name"]]), ]))$Max.Score)
  df.108NPR1[k, 2] <- first(L.108NPR1[[k]][["Acc..Len"]][L.108NPR1[[k]][["Max.Score"]] == max(L.108NPR1[[k]][["Max.Score"]])])
  df.108NPR1[, 3] <- rep(1054, 28)
  df.108NPR1[, 4] <- df.108NPR1[, 1] / ((df.108NPR1[, 3]) + df.108NPR1[, 2])
  colnames(df.108NPR1) <- c("maxscore", "p.length", "ref.p.length",  "NPR1.maxscore.scaled")
}


# -------------------------------------------- 112MYLK3
df.112MYLK3 <- data.frame(1:nrow(pat))
L.112MYLK3 <- list()

for (k in 1:nrow(pat)) {
  L.112MYLK3[[k]] <- data.frame(L1[[73]][grep(pat[k,], L1[[73]][["Scientific.Name"]]), ])
  df.112MYLK3[k,] <- max((data.frame(L1[[73]][grep(pat[k,], L1[[73]][["Scientific.Name"]]), ]))$Max.Score)
  df.112MYLK3[k, 2] <- first(L.112MYLK3[[k]][["Acc..Len"]][L.112MYLK3[[k]][["Max.Score"]] == max(L.112MYLK3[[k]][["Max.Score"]])])
  df.112MYLK3[, 3] <- rep(895, 28)
  df.112MYLK3[, 4] <- df.112MYLK3[, 1] / ((df.112MYLK3[, 3]) + df.112MYLK3[, 2])
  colnames(df.112MYLK3) <- c("maxscore", "p.length", "ref.p.length",  "MYLK3.maxscore.scaled")
}


# -------------------------------------------- 116MYL9
df.116MYL9 <- data.frame(1:nrow(pat))
L.116MYL9 <- list()

for (k in 1:nrow(pat)) {
  L.116MYL9[[k]] <- data.frame(L1[[74]][grep(pat[k,], L1[[74]][["Scientific.Name"]]), ])
  df.116MYL9[k,] <- max((data.frame(L1[[74]][grep(pat[k,], L1[[74]][["Scientific.Name"]]), ]))$Max.Score)
  df.116MYL9[k, 2] <- first(L.116MYL9[[k]][["Acc..Len"]][L.116MYL9[[k]][["Max.Score"]] == max(L.116MYL9[[k]][["Max.Score"]])])
  df.116MYL9[, 3] <- rep(172, 28)
  df.116MYL9[, 4] <- df.116MYL9[, 1] / ((df.116MYL9[, 3]) + df.116MYL9[, 2])
  colnames(df.116MYL9) <- c("maxscore", "p.length", "ref.p.length",  "MYL9.maxscore.scaled")
}


# -------------------------------------------- 119RHOA
df.119RHOA <- data.frame(1:nrow(pat))
L.119RHOA <- list()

for (k in 1:nrow(pat)) {
  L.119RHOA[[k]] <- data.frame(L1[[75]][grep(pat[k,], L1[[75]][["Scientific.Name"]]), ])
  df.119RHOA[k,] <- max((data.frame(L1[[75]][grep(pat[k,], L1[[75]][["Scientific.Name"]]), ]))$Max.Score)
  df.119RHOA[k, 2] <- first(L.119RHOA[[k]][["Acc..Len"]][L.119RHOA[[k]][["Max.Score"]] == max(L.119RHOA[[k]][["Max.Score"]])])
  df.119RHOA[, 3] <- rep(194, 28)
  df.119RHOA[, 4] <- df.119RHOA[, 1] / ((df.119RHOA[, 3]) + df.119RHOA[, 2])
  colnames(df.119RHOA) <- c("maxscore", "p.length", "ref.p.length",  "RHOA.maxscore.scaled")
}


# -------------------------------------------- 120ROCK1
df.120ROCK1 <- data.frame(1:nrow(pat))
L.120ROCK1 <- list()

for (k in 1:nrow(pat)) {
  L.120ROCK1[[k]] <- data.frame(L1[[76]][grep(pat[k,], L1[[76]][["Scientific.Name"]]), ])
  df.120ROCK1[k,] <- max((data.frame(L1[[76]][grep(pat[k,], L1[[76]][["Scientific.Name"]]), ]))$Max.Score)
  df.120ROCK1[k, 2] <- first(L.120ROCK1[[k]][["Acc..Len"]][L.120ROCK1[[k]][["Max.Score"]] == max(L.120ROCK1[[k]][["Max.Score"]])])
  df.120ROCK1[, 3] <- rep(1364, 28)
  df.120ROCK1[, 4] <- df.120ROCK1[, 1] / ((df.120ROCK1[, 3]) + df.120ROCK1[, 2])
  colnames(df.120ROCK1) <- c("maxscore", "p.length", "ref.p.length",  "ROCK1.maxscore.scaled")
}


# -------------------------------------------- 121ROCK2
df.121ROCK2 <- data.frame(1:nrow(pat))
L.121ROCK2 <- list()

for (k in 1:nrow(pat)) {
  L.121ROCK2[[k]] <- data.frame(L1[[77]][grep(pat[k,], L1[[77]][["Scientific.Name"]]), ])
  df.121ROCK2[k,] <- max((data.frame(L1[[77]][grep(pat[k,], L1[[77]][["Scientific.Name"]]), ]))$Max.Score)
  df.121ROCK2[k, 2] <- first(L.121ROCK2[[k]][["Acc..Len"]][L.121ROCK2[[k]][["Max.Score"]] == max(L.121ROCK2[[k]][["Max.Score"]])])
  df.121ROCK2[, 3] <- rep(1403, 28)
  df.121ROCK2[, 4] <- df.121ROCK2[, 1] / ((df.121ROCK2[, 3]) + df.121ROCK2[, 2])
  colnames(df.121ROCK2) <- c("maxscore", "p.length", "ref.p.length",  "ROCK2.maxscore.scaled")
}


# -------------------------------------------- 123PPP1CB
df.123PPP1CB <- data.frame(1:nrow(pat))
L.123PPP1CB <- list()

for (k in 1:nrow(pat)) {
  L.123PPP1CB[[k]] <- data.frame(L1[[78]][grep(pat[k,], L1[[78]][["Scientific.Name"]]), ])
  df.123PPP1CB[k,] <- max((data.frame(L1[[78]][grep(pat[k,], L1[[78]][["Scientific.Name"]]), ]))$Max.Score)
  df.123PPP1CB[k, 2] <- first(L.123PPP1CB[[k]][["Acc..Len"]][L.123PPP1CB[[k]][["Max.Score"]] == max(L.123PPP1CB[[k]][["Max.Score"]])])
  df.123PPP1CB[, 3] <- rep(327, 28)
  df.123PPP1CB[, 4] <- df.123PPP1CB[, 1] / ((df.123PPP1CB[, 3]) + df.123PPP1CB[, 2])
  colnames(df.123PPP1CB) <- c("maxscore", "p.length", "ref.p.length",  "PPP1CB.maxscore.scaled")
}


# -------------------------------------------- 124PPP1CC
df.124PPP1CC <- data.frame(1:nrow(pat))
L.124PPP1CC <- list()

for (k in 1:nrow(pat)) {
  L.124PPP1CC[[k]] <- data.frame(L1[[79]][grep(pat[k,], L1[[79]][["Scientific.Name"]]), ])
  df.124PPP1CC[k,] <- max((data.frame(L1[[79]][grep(pat[k,], L1[[79]][["Scientific.Name"]]), ]))$Max.Score)
  df.124PPP1CC[k, 2] <- first(L.124PPP1CC[[k]][["Acc..Len"]][L.124PPP1CC[[k]][["Max.Score"]] == max(L.124PPP1CC[[k]][["Max.Score"]])])
  df.124PPP1CC[, 3] <- rep(325, 28)
  df.124PPP1CC[, 4] <- df.124PPP1CC[, 1] / ((df.124PPP1CC[, 3]) + df.124PPP1CC[, 2])
  colnames(df.124PPP1CC) <- c("maxscore", "p.length", "ref.p.length",  "PPP1CC.maxscore.scaled")
}


# -------------------------------------------- 125PPP1R12A
df.125PPP1R12A <- data.frame(1:nrow(pat))
L.125PPP1R12A <- list()

for (k in 1:nrow(pat)) {
  L.125PPP1R12A[[k]] <- data.frame(L1[[80]][grep(pat[k,], L1[[80]][["Scientific.Name"]]), ])
  df.125PPP1R12A[k,] <- max((data.frame(L1[[80]][grep(pat[k,], L1[[80]][["Scientific.Name"]]), ]))$Max.Score)
  df.125PPP1R12A[k, 2] <- first(L.125PPP1R12A[[k]][["Acc..Len"]][L.125PPP1R12A[[k]][["Max.Score"]] == max(L.125PPP1R12A[[k]][["Max.Score"]])])
  df.125PPP1R12A[, 3] <- rep(1117, 28)
  df.125PPP1R12A[, 4] <- df.125PPP1R12A[, 1] / ((df.125PPP1R12A[, 3]) + df.125PPP1R12A[, 2])
  colnames(df.125PPP1R12A) <- c("maxscore", "p.length", "ref.p.length",  "PPP1R12A.maxscore.scaled")
}


# -------------------------------------------- 129ADCY1
df.129ADCY1 <- data.frame(1:nrow(pat))
L.129ADCY1 <- list()

for (k in 1:nrow(pat)) {
  L.129ADCY1[[k]] <- data.frame(L1[[81]][grep(pat[k,], L1[[81]][["Scientific.Name"]]), ])
  df.129ADCY1[k,] <- max((data.frame(L1[[81]][grep(pat[k,], L1[[81]][["Scientific.Name"]]), ]))$Max.Score)
  df.129ADCY1[k, 2] <- first(L.129ADCY1[[k]][["Acc..Len"]][L.129ADCY1[[k]][["Max.Score"]] == max(L.129ADCY1[[k]][["Max.Score"]])])
  df.129ADCY1[, 3] <- rep(1112, 28)
  df.129ADCY1[, 4] <- df.129ADCY1[, 1] / ((df.129ADCY1[, 3]) + df.129ADCY1[, 2])
  colnames(df.129ADCY1) <- c("maxscore", "p.length", "ref.p.length",  "ADCY1.maxscore.scaled")
}


# -------------------------------------------- 130ADCY2
df.130ADCY2 <- data.frame(1:nrow(pat))
L.130ADCY2 <- list()

for (k in 1:nrow(pat)) {
  L.130ADCY2[[k]] <- data.frame(L1[[82]][grep(pat[k,], L1[[82]][["Scientific.Name"]]), ])
  df.130ADCY2[k,] <- max((data.frame(L1[[82]][grep(pat[k,], L1[[82]][["Scientific.Name"]]), ]))$Max.Score)
  df.130ADCY2[k, 2] <- first(L.130ADCY2[[k]][["Acc..Len"]][L.130ADCY2[[k]][["Max.Score"]] == max(L.130ADCY2[[k]][["Max.Score"]])])
  df.130ADCY2[, 3] <- rep(1093, 28)
  df.130ADCY2[, 4] <- df.130ADCY2[, 1] / ((df.130ADCY2[, 3]) + df.130ADCY2[, 2])
  colnames(df.130ADCY2) <- c("maxscore", "p.length", "ref.p.length",  "ADCY2.maxscore.scaled")
}


# -------------------------------------------- 133ADCY5
df.133ADCY5 <- data.frame(1:nrow(pat))
L.133ADCY5 <- list()

for (k in 1:nrow(pat)) {
  L.133ADCY5[[k]] <- data.frame(L1[[83]][grep(pat[k,], L1[[83]][["Scientific.Name"]]), ])
  df.133ADCY5[k,] <- max((data.frame(L1[[83]][grep(pat[k,], L1[[83]][["Scientific.Name"]]), ]))$Max.Score)
  df.133ADCY5[k, 2] <- first(L.133ADCY5[[k]][["Acc..Len"]][L.133ADCY5[[k]][["Max.Score"]] == max(L.133ADCY5[[k]][["Max.Score"]])])
  df.133ADCY5[, 3] <- rep(1189, 28)
  df.133ADCY5[, 4] <- df.133ADCY5[, 1] / ((df.133ADCY5[, 3]) + df.133ADCY5[, 2])
  colnames(df.133ADCY5) <- c("maxscore", "p.length", "ref.p.length",  "ADCY5.maxscore.scaled")
}


# -------------------------------------------- 135ADCY7
df.135ADCY7 <- data.frame(1:nrow(pat))
L.135ADCY7 <- list()

for (k in 1:nrow(pat)) {
  L.135ADCY7[[k]] <- data.frame(L1[[84]][grep(pat[k,], L1[[84]][["Scientific.Name"]]), ])
  df.135ADCY7[k,] <- max((data.frame(L1[[84]][grep(pat[k,], L1[[84]][["Scientific.Name"]]), ]))$Max.Score)
  df.135ADCY7[k, 2] <- first(L.135ADCY7[[k]][["Acc..Len"]][L.135ADCY7[[k]][["Max.Score"]] == max(L.135ADCY7[[k]][["Max.Score"]])])
  df.135ADCY7[, 3] <- rep(1062, 28)
  df.135ADCY7[, 4] <- df.135ADCY7[, 1] / ((df.135ADCY7[, 3]) + df.135ADCY7[, 2])
  colnames(df.135ADCY7) <- c("maxscore", "p.length", "ref.p.length",  "ADCY7.maxscore.scaled")
}


# -------------------------------------------- 136ADCY8
df.136ADCY8 <- data.frame(1:nrow(pat))
L.136ADCY8 <- list()

for (k in 1:nrow(pat)) {
  L.136ADCY8[[k]] <- data.frame(L1[[85]][grep(pat[k,], L1[[85]][["Scientific.Name"]]), ])
  df.136ADCY8[k,] <- max((data.frame(L1[[85]][grep(pat[k,], L1[[85]][["Scientific.Name"]]), ]))$Max.Score)
  df.136ADCY8[k, 2] <- first(L.136ADCY8[[k]][["Acc..Len"]][L.136ADCY8[[k]][["Max.Score"]] == max(L.136ADCY8[[k]][["Max.Score"]])])
  df.136ADCY8[, 3] <- rep(1217, 28)
  df.136ADCY8[, 4] <- df.136ADCY8[, 1] / ((df.136ADCY8[, 3]) + df.136ADCY8[, 2])
  colnames(df.136ADCY8) <- c("maxscore", "p.length", "ref.p.length",  "ADCY8.maxscore.scaled")
}


# -------------------------------------------- 137ADCY9
df.137ADCY9 <- data.frame(1:nrow(pat))
L.137ADCY9 <- list()

for (k in 1:nrow(pat)) {
  L.137ADCY9[[k]] <- data.frame(L1[[86]][grep(pat[k,], L1[[86]][["Scientific.Name"]]), ])
  df.137ADCY9[k,] <- max((data.frame(L1[[86]][grep(pat[k,], L1[[86]][["Scientific.Name"]]), ]))$Max.Score)
  df.137ADCY9[k, 2] <- first(L.137ADCY9[[k]][["Acc..Len"]][L.137ADCY9[[k]][["Max.Score"]] == max(L.137ADCY9[[k]][["Max.Score"]])])
  df.137ADCY9[, 3] <- rep(1552, 28)
  df.137ADCY9[, 4] <- df.137ADCY9[, 1] / ((df.137ADCY9[, 3]) + df.137ADCY9[, 2])
  colnames(df.137ADCY9) <- c("maxscore", "p.length", "ref.p.length",  "ADCY9.maxscore.scaled")
}


# -------------------------------------------- 138PRKACA
df.138PRKACA <- data.frame(1:nrow(pat))
L.138PRKACA <- list()

for (k in 1:nrow(pat)) {
  L.138PRKACA[[k]] <- data.frame(L1[[87]][grep(pat[k,], L1[[87]][["Scientific.Name"]]), ])
  df.138PRKACA[k,] <- max((data.frame(L1[[87]][grep(pat[k,], L1[[87]][["Scientific.Name"]]), ]))$Max.Score)
  df.138PRKACA[k, 2] <- first(L.138PRKACA[[k]][["Acc..Len"]][L.138PRKACA[[k]][["Max.Score"]] == max(L.138PRKACA[[k]][["Max.Score"]])])
  df.138PRKACA[, 3] <- rep(403, 28)
  df.138PRKACA[, 4] <- df.138PRKACA[, 1] / ((df.138PRKACA[, 3]) + df.138PRKACA[, 2])
  colnames(df.138PRKACA) <- c("maxscore", "p.length", "ref.p.length",  "PRKACA.maxscore.scaled")
}


# -------------------------------------------- 139PRKACB
df.139PRKACB <- data.frame(1:nrow(pat))
L.139PRKACB <- list()

for (k in 1:nrow(pat)) {
  L.139PRKACB[[k]] <- data.frame(L1[[88]][grep(pat[k,], L1[[88]][["Scientific.Name"]]), ])
  df.139PRKACB[k,] <- max((data.frame(L1[[88]][grep(pat[k,], L1[[88]][["Scientific.Name"]]), ]))$Max.Score)
  df.139PRKACB[k, 2] <- first(L.139PRKACB[[k]][["Acc..Len"]][L.139PRKACB[[k]][["Max.Score"]] == max(L.139PRKACB[[k]][["Max.Score"]])])
  df.139PRKACB[, 3] <- rep(403, 28)
  df.139PRKACB[, 4] <- df.139PRKACB[, 1] / ((df.139PRKACB[, 3]) + df.139PRKACB[, 2])
  colnames(df.139PRKACB) <- c("maxscore", "p.length", "ref.p.length",  "PRKACB.maxscore.scaled")
}


# -------------------------------------------- 141GNAI1
df.141GNAI1 <- data.frame(1:nrow(pat))
L.141GNAI1 <- list()

for (k in 1:nrow(pat)) {
  L.141GNAI1[[k]] <- data.frame(L1[[89]][grep(pat[k,], L1[[89]][["Scientific.Name"]]), ])
  df.141GNAI1[k,] <- max((data.frame(L1[[89]][grep(pat[k,], L1[[89]][["Scientific.Name"]]), ]))$Max.Score)
  df.141GNAI1[k, 2] <- first(L.141GNAI1[[k]][["Acc..Len"]][L.141GNAI1[[k]][["Max.Score"]] == max(L.141GNAI1[[k]][["Max.Score"]])])
  df.141GNAI1[, 3] <- rep(354, 28)
  df.141GNAI1[, 4] <- df.141GNAI1[, 1] / ((df.141GNAI1[, 3]) + df.141GNAI1[, 2])
  colnames(df.141GNAI1) <- c("maxscore", "p.length", "ref.p.length",  "GNAI1.maxscore.scaled")
}


# -------------------------------------------- 146PIK3R5
df.146PIK3R5 <- data.frame(1:nrow(pat))
L.146PIK3R5 <- list()

for (k in 1:nrow(pat)) {
  L.146PIK3R5[[k]] <- data.frame(L1[[90]][grep(pat[k,], L1[[90]][["Scientific.Name"]]), ])
  df.146PIK3R5[k,] <- max((data.frame(L1[[90]][grep(pat[k,], L1[[90]][["Scientific.Name"]]), ]))$Max.Score)
  df.146PIK3R5[k, 2] <- first(L.146PIK3R5[[k]][["Acc..Len"]][L.146PIK3R5[[k]][["Max.Score"]] == max(L.146PIK3R5[[k]][["Max.Score"]])])
  df.146PIK3R5[, 3] <- rep(904, 28)
  df.146PIK3R5[, 4] <- df.146PIK3R5[, 1] / ((df.146PIK3R5[, 3]) + df.146PIK3R5[, 2])
  colnames(df.146PIK3R5) <- c("maxscore", "p.length", "ref.p.length",  "PIK3R5.maxscore.scaled")
}


# -------------------------------------------- 148SRC
df.148SRC <- data.frame(1:nrow(pat))
L.148SRC <- list()

for (k in 1:nrow(pat)) {
  L.148SRC[[k]] <- data.frame(L1[[91]][grep(pat[k,], L1[[91]][["Scientific.Name"]]), ])
  df.148SRC[k,] <- max((data.frame(L1[[91]][grep(pat[k,], L1[[91]][["Scientific.Name"]]), ]))$Max.Score)
  df.148SRC[k, 2] <- first(L.148SRC[[k]][["Acc..Len"]][L.148SRC[[k]][["Max.Score"]] == max(L.148SRC[[k]][["Max.Score"]])])
  df.148SRC[, 3] <- rep(540, 28)
  df.148SRC[, 4] <- df.148SRC[, 1] / ((df.148SRC[, 3]) + df.148SRC[, 2])
  colnames(df.148SRC) <- c("maxscore", "p.length", "ref.p.length",  "SRC.maxscore.scaled")
}


# -------------------------------------------- 149KCNJ3
df.149KCNJ3 <- data.frame(1:nrow(pat))
L.149KCNJ3 <- list()

for (k in 1:nrow(pat)) {
  L.149KCNJ3[[k]] <- data.frame(L1[[92]][grep(pat[k,], L1[[92]][["Scientific.Name"]]), ])
  df.149KCNJ3[k,] <- max((data.frame(L1[[92]][grep(pat[k,], L1[[92]][["Scientific.Name"]]), ]))$Max.Score)
  df.149KCNJ3[k, 2] <- first(L.149KCNJ3[[k]][["Acc..Len"]][L.149KCNJ3[[k]][["Max.Score"]] == max(L.149KCNJ3[[k]][["Max.Score"]])])
  df.149KCNJ3[, 3] <- rep(493, 28)
  df.149KCNJ3[, 4] <- df.149KCNJ3[, 1] / ((df.149KCNJ3[, 3]) + df.149KCNJ3[, 2])
  colnames(df.149KCNJ3) <- c("maxscore", "p.length", "ref.p.length",  "KCNJ3.maxscore.scaled")
}


# -------------------------------------------- 150KCNJ6
df.150KCNJ6 <- data.frame(1:nrow(pat))
L.150KCNJ6 <- list()

for (k in 1:nrow(pat)) {
  L.150KCNJ6[[k]] <- data.frame(L1[[93]][grep(pat[k,], L1[[93]][["Scientific.Name"]]), ])
  df.150KCNJ6[k,] <- max((data.frame(L1[[93]][grep(pat[k,], L1[[93]][["Scientific.Name"]]), ]))$Max.Score)
  df.150KCNJ6[k, 2] <- first(L.150KCNJ6[[k]][["Acc..Len"]][L.150KCNJ6[[k]][["Max.Score"]] == max(L.150KCNJ6[[k]][["Max.Score"]])])
  df.150KCNJ6[, 3] <- rep(407, 28)
  df.150KCNJ6[, 4] <- df.150KCNJ6[, 1] / ((df.150KCNJ6[, 3]) + df.150KCNJ6[, 2])
  colnames(df.150KCNJ6) <- c("maxscore", "p.length", "ref.p.length",  "KCNJ6.maxscore.scaled")
}


# -------------------------------------------- 152KCNJ5
df.152KCNJ5 <- data.frame(1:nrow(pat))
L.152KCNJ5 <- list()

for (k in 1:nrow(pat)) {
  L.152KCNJ5[[k]] <- data.frame(L1[[94]][grep(pat[k,], L1[[94]][["Scientific.Name"]]), ])
  df.152KCNJ5[k,] <- max((data.frame(L1[[94]][grep(pat[k,], L1[[94]][["Scientific.Name"]]), ]))$Max.Score)
  df.152KCNJ5[k, 2] <- first(L.152KCNJ5[[k]][["Acc..Len"]][L.152KCNJ5[[k]][["Max.Score"]] == max(L.152KCNJ5[[k]][["Max.Score"]])])
  df.152KCNJ5[, 3] <- rep(444, 28)
  df.152KCNJ5[, 4] <- df.152KCNJ5[, 1] / ((df.152KCNJ5[, 3]) + df.152KCNJ5[, 2])
  colnames(df.152KCNJ5) <- c("maxscore", "p.length", "ref.p.length",  "KCNJ5.maxscore.scaled")
}


# -------------------------------------------- 153EGFR
df.153EGFR <- data.frame(1:nrow(pat))
L.153EGFR <- list()

for (k in 1:nrow(pat)) {
  L.153EGFR[[k]] <- data.frame(L1[[95]][grep(pat[k,], L1[[95]][["Scientific.Name"]]), ])
  df.153EGFR[k,] <- max((data.frame(L1[[95]][grep(pat[k,], L1[[95]][["Scientific.Name"]]), ]))$Max.Score)
  df.153EGFR[k, 2] <- first(L.153EGFR[[k]][["Acc..Len"]][L.153EGFR[[k]][["Max.Score"]] == max(L.153EGFR[[k]][["Max.Score"]])])
  df.153EGFR[, 3] <- rep(1216, 28)
  df.153EGFR[, 4] <- df.153EGFR[, 1] / ((df.153EGFR[, 3]) + df.153EGFR[, 2])
  colnames(df.153EGFR) <- c("maxscore", "p.length", "ref.p.length", "EGFR.maxscore.scaled")
}




res <- cbind(df.001oxt, df.002oxtr, df.003gnaq, df.004hras, df.005kras, df.006nras, df.007raf1, df.008map2k1, df.009map2k2, df.010mapk1,
             df.013pla2g4a, df.018pla2g4f, df.019ptgs2, df.020map2k5, df.022jun, df.023fos, df.024mef2c, df.025ccnd1, df.026elk1,
             df.028ryr2, df.029RYR3, df.030CD38, df.032KCNJ2, df.035KCNJ4, df.037PLCB1, df.040PLCB4, df.041PRKCA, df.042PRKCB,
             df.044EEF2K, df.045EEF2, df.046CACNA1C, df.047CACNA1D, df.050CACNB1, df.051CACNB2, df.053CACNB4, df.054CACNA2D1, 
             df.056CACNA2D3, df.057CACNA2D4, df.058CACNG1, df.059CACNG2, df.060CACNG3, df.061CACNG4, df.062CACNG5, df.066ITPR1, 
             df.067ITPR2, df.070CALM2, df.072CALM1, df.073CALML6, df.075CALML4, df.076PPP3CA, df.077PPP3CB, df.078PPP3CC, df.079PPP3R1,
             df.081NFATC1, df.082NFATC2, df.083NFATC3, df.087CAMK2, df.088PRKAA1, df.089PRKAA2, df.090PRKAB1, df.092PRKAG1, 
             df.094PRKAG2, df.095CAMK1D, df.096CAMK1G, df.097CAMK1, df.098CAMK2A, df.099CAMK2D, df.101CAMK2G, df.102CAMK4, df.104GUCY1A2,
             df.106GUCY1B1, df.108NPR1, df.112MYLK3, df.116MYL9, df.119RHOA, df.120ROCK1, df.121ROCK2, df.123PPP1CB, df.124PPP1CC,
             df.125PPP1R12A, df.129ADCY1, df.130ADCY2, df.133ADCY5, df.135ADCY7, df.136ADCY8, df.137ADCY9, df.138PRKACA, df.139PRKACB,
             df.141GNAI1, df.146PIK3R5, df.148SRC, df.149KCNJ3, df.150KCNJ6, df.152KCNJ5, df.153EGFR)
  
res.c <- res %>% 
  select(-starts_with("max"), -p.length, -ref.p.length)

res.c$rowM <- rowMeans(res.c, na.rm = T)

res.c <- relocate(res.c, rowM, .before = "OXT.maxscore.scaled")

res.c.ex <- res.c %>% 
  mutate(pat, .before = "rowM") %>% 
  rename(species.name = c..Ciona.savignyi....Ciona.intestinalis....Branchiostoma.lanceolatum...)

library("writexl")
write_xlsx(res.c.ex, paste0(path, "output/BLASTp/scaled_max_scores.xlsx"))



L.hom <- list()
for (i in 1:ncol(res.c)) {
  L.hom[i] <- data.frame(res.c[[i]] > res.c[["rowM"]])
  names(L.hom)[i] <- paste(colnames(res.c[i]))
}

gene.ages <- as.data.frame(L.hom) %>% 
  mutate(pat, .before = "rowM") %>% 
  rename(species.name = c..Ciona.savignyi....Ciona.intestinalis....Branchiostoma.lanceolatum...) 

write_xlsx(gene.ages, paste0(path, "output/BLASTp/homology_boolean.xlsx"))

### --------------------------

library(readxl)
OTpthwy.res <- read_excel(paste0(path, "data/processed/gene_ages.xlsx"), sheet = 1)

for (i in 1:10) {
OTpthwy.res$Ortholog[OTpthwy.res$Phylostratum == i] <- "n"
}
for (i in 11:20) {
  OTpthwy.res$Ortholog[OTpthwy.res$Phylostratum == i] <- "y"
}



for (i in 1:10) {
  OTpthwy.res$Homolog[OTpthwy.res$Phylostratum == i] <- "y"
}
for (i in 11:20) {
  OTpthwy.res$Homolog[OTpthwy.res$Phylostratum == i] <- "n"
}


OTpthwy.res$age.2.levels <- NA
for (i in 1:3) {
  OTpthwy.res$age.2.levels[OTpthwy.res$Phylostratum == i] <- "ancient"
}
for (i in 4:20) {
  OTpthwy.res$age.2.levels[OTpthwy.res$Phylostratum == i] <- "modern"
}


OTpthwy.res$age.3.levels <- NA
for (i in 1:3) {
  OTpthwy.res$age.3.levels[OTpthwy.res$Phylostratum == i] <- "ancient"
}
for (i in 4:10) {
  OTpthwy.res$age.3.levels[OTpthwy.res$Phylostratum == i] <- "mod.non.vertebrate"
}
for (i in 11:20) {
  OTpthwy.res$age.3.levels[OTpthwy.res$Phylostratum == i] <- "mod.vertebrate"
}

OTpthwy.res %>% 
  count(age.3.levels)

write_xlsx(OTpthwy.res, paste0(path, "output/BLASTp/OTpthwy_res.xlsx"))





#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################
#############################################################################

### load 154 datasets
# OXT
oxt.hs <- read_csv("data/raw/BLASTp/001OXT/OXT-h-sapiens.csv")
oxt.cc <- read_csv("data/raw/BLASTp/001OXT/OXT-c-carcharias.csv")

# OXTR
oxtr.hs <- read_csv("data/raw/BLASTp/002OXTR/OXTR-h-sapiens.csv")
oxtr.pm <- read_csv("data/raw/BLASTp/002OXTR/OXTR-p-marinus.csv")

# GNAQ
gnaq.hs <- read_csv("data/raw/BLASTp/003GNAQ/GNAQ-h-sapiens.csv")
gnaq.cc <- read_csv("data/raw/BLASTp/003GNAQ/GNAQ-c-carcharias.csv")

# HRAS
hras.hs <- read_csv("data/raw/BLASTp/004HRAS/HRAS-h-sapiens.csv")
hras.pm <- read_csv("data/raw/BLASTp/004HRAS/HRAS-p-marinus.csv")

# KRAS
kras.hs <- read_csv("data/raw/BLASTp/005KRAS/KRAS-h-sapiens.csv")
kras.cc <- read_csv("data/raw/BLASTp/005KRAS/KRAS-c-carcharias.csv")

# NRAS
nras.hs <- read_csv("data/raw/BLASTp/006NRAS/NRAS-h-sapiens.csv")
nras.cc <- read_csv("data/raw/BLASTp/006NRAS/NRAS-c-carcharias.csv")

# RAF1
raf1.hs <- read_csv("data/raw/BLASTp/007RAF1/RAF1-h-sapiens.csv")
raf1.cc <- read_csv("data/raw/BLASTp/007RAF1/RAF1-c-carcharias.csv")

# MAP2K1
map2k1.hs <- read_csv("data/raw/BLASTp/008MAP2K1/MAP2K1-h-sapiens.csv")
map2k1.cc <- read_csv("data/raw/BLASTp/008MAP2K1/MAP2K1-c-carcharias.csv")

# MAP2K2
map2k2.hs <- read_csv("data/raw/BLASTp/009MAP2K2/MAP2K2-h-sapiens.csv")
map2k2.cc <- read_csv("data/raw/BLASTp/009MAP2K2/MAP2K2-c-carcharias.csv")

# MAPK1
mapk1.hs <- read_csv("data/raw/BLASTp/010MAPK1/MAPK1-h-sapiens.csv")
mapk1.cc <- read_csv("data/raw/BLASTp/010MAPK1/MAPK1-c-carcharias.csv")

# MAPK3
# no hits accoridng to microsynteny in c. carcharias or p. marinus

# PLA2G4E
# no hits accoridng to microsynteny in c. carcharias or p. marinus

# PLA2G4A
pla2g4a.hs <- read_csv("data/raw/BLASTp/013PLA2G4A/PLA2G4A-h-sapiens.csv")
pla2g4a.pm <- read_csv("data/raw/BLASTp/013PLA2G4A/PLA2G4A-p-marinus.csv")

# JMJD7-PLA2G4B
# no hits accoridng to microsynteny in c. carcharias or p. marinus

# PLA2G4B
# no hits accoridng to microsynteny in c. carcharias or p. marinus

# PLA2G4C
# no hits accoridng to microsynteny in c. carcharias or p. marinus

# PLA2G4D
# no hits accoridng to microsynteny in c. carcharias or p. marinus

# PLA2G4F
pla2g4f.hs <- read_csv("data/raw/BLASTp/018PLA2G4F/PLA2G4F-h-sapiens.csv")
pla2g4f.cc <- read_csv("data/raw/BLASTp/018PLA2G4F/PLA2G4F-c-carcharias.csv")

# PTGS2
ptgs2.hs <- read_csv("data/raw/BLASTp/019PTGS2/PTGS2-h-sapiens.csv")
ptgs2.pm <- read_csv("data/raw/BLASTp/019PTGS2/PTGS2-p-marinus.csv")

# MAP2K5
map2k5.hs <- read_csv("data/raw/BLASTp/020MAP2K5/MAP2K5-h-sapiens.csv")
map2k5.pm <- read_csv("data/raw/BLASTp/020MAP2K5/MAP2K5-p-marinus.csv")

# MAPK7
# no hits accoridng to microsynteny in c. carcharias or p. marinus

# JUN
jun.hs <- read_csv("data/raw/BLASTp/022JUN/JUN-h-sapiens.csv")
jun.cc <- read_csv("data/raw/BLASTp/022JUN/JUN-c-carcharias.csv")

# FOS
fos.hs <- read_csv("data/raw/BLASTp/023FOS/FOS-h-sapiens.csv")
fos.cc <- read_csv("data/raw/BLASTp/023FOS/FOS-c-carcharias.csv")

# MEF2C
mef2c.hs <- read_csv("data/raw/BLASTp/024MEF2C/MEF2C-h-sapiens.csv")
mef2c.cc <- read_csv("data/raw/BLASTp/024MEF2C/MEF2C-c-carcharias.csv")

# CCND1
ccnd1.hs <- read_csv("data/raw/BLASTp/025CCND1/CCND1-h-sapiens.csv")
ccnd1.cc <- read_csv("data/raw/BLASTp/025CCND1/CCND1-c-carcharias.csv")

# ELK1
elk1.hs <- read_csv("data/raw/BLASTp/026ELK1/ELK1-h-sapiens.csv")
elk1.cc <- read_csv("data/raw/BLASTp/026ELK1/ELK1-c-carcharias.csv")

# RYR1
# no reliable hits according to microsynteny in c. carcharias or p. marinus

# RYR2
ryr2.hs <- read_csv("data/raw/BLASTp/028RYR2/RYR2-h-sapiens.csv")
ryr2.cc <- read_csv("data/raw/BLASTp/028RYR2/RYR2-c-carcharias.csv")

# RYR3
ryr3.hs <- read_csv("data/raw/BLASTp/029RYR3/RYR3-h-sapiens.csv")
ryr3.cc <- read_csv("data/raw/BLASTp/029RYR3/RYR3-c-carcharias.csv")

# CD38
cd38.hs <- read_csv("data/raw/BLASTp/030CD38/CD38-h-sapiens.csv")
cd38.cc <- read_csv("data/raw/BLASTp/030CD38/CD38-c-carcharias.csv")

# TRPM2
# no reliable hits according to microsynteny in c. carcharias or p. marinus

# KCNJ2
kcnj2.hs <- read_csv("data/raw/BLASTp/032KCNJ2/KCNJ2-h-sapiens.csv")
kcnj2.cc <- read_csv("data/raw/BLASTp/032KCNJ2/KCNJ2-c-carcharias.csv")

# KCNJ12 
# no reliable hits according to microsynteny in c. carcharias or p. marinus

# KCNJ18
# no reliable hits according to microsynteny in c. carcharias or p. marinus

# KCNJ4
kcnj4.hs <- read_csv("data/raw/BLASTp/035KCNJ4/KCNJ4-h-sapiens.csv")
kcnj4.cc <- read_csv("data/raw/BLASTp/035KCNJ4/KCNJ4-c-carcharias.csv")

# KCNJ14
# no reliable hits according to microsynteny in c. carcharias or p. marinus

# PLCB1
plcb1.hs <- read_csv("data/raw/BLASTp/037PLCB1/PLCB1-h-sapiens.csv")
plcb1.cc <- read_csv("data/raw/BLASTp/037PLCB1/PLCB1-c-carcharias.csv")

# PLCB2
# no reliable hits according to microsynteny in c. carcharias or p. marinus

# PLCB3
# no reliable hits according to microsynteny in c. carcharias or p. marinus

# PLCB4
plcb4.hs <- read_csv("data/raw/BLASTp/040PLCB4/PLCB4-h-sapiens.csv")
plcb4.cc <- read_csv("data/raw/BLASTp/040PLCB4/PLCB4-c-carcharias.csv")

# PRKCA
prkca.hs <- read_csv("data/raw/BLASTp/041PRKCA/PRKCA-h-sapiens.csv")
prkca.cc <- read_csv("data/raw/BLASTp/041PRKCA/PRKCA-c-carcharias.csv")

# PRKCB
prkcb.hs <- read_csv("data/raw/BLASTp/042PRKCB/PRKCB-h-sapiens.csv")
prkcb.cc <- read_csv("data/raw/BLASTp/042PRKCB/PRKCB-c-carcharias.csv")

# PRKCG
# no reliable hits according to microsynteny in c. carcharias or p. marinus

# EEF2K
eef2k.hs <- read_csv("data/raw/BLASTp/044EEF2K/EEF2K-h-sapiens.csv")
eef2k.cc <- read_csv("data/raw/BLASTp/044EEF2K/EEF2K-c-carcharias.csv")

# EEF2
eef2.hs <- read_csv("data/raw/BLASTp/045EEF2/EEF2-h-sapiens.csv")
eef2.cc <- read_csv("data/raw/BLASTp/045EEF2/EEF2-c-carcharias.csv")

# CACNA1C
cacna1c.hs <- read_csv("data/raw/BLASTp/046CACNA1C/CACNA1C-h-sapiens.csv")
cacna1c.cc <- read_csv("data/raw/BLASTp/046CACNA1C/CACNA1C-c-carcharias.csv")

# CACNA1D
cacna1d.hs <- read_csv("data/raw/BLASTp/047CACNA1D/CACNA1D-h-sapiens.csv")
cacna1d.cc <- read_csv("data/raw/BLASTp/047CACNA1D/CACNA1D-c-carcharias.csv")

# CACNA1F
# no reliable hits according to microsynteny in c. carcharias or p. marinus

# CACNA1S
# no reliable hits according to microsynteny in c. carcharias or p. marinus

# CACNB1
cacnb1.hs <- read_csv("data/raw/BLASTp/050CACNB1/CACNB1-h-sapiens.csv")
cacnb1.pm <- read_csv("data/raw/BLASTp/050CACNB1/CACNB1-p-marinus.csv")

# CACNB2
cacnb2.hs <- read_csv("data/raw/BLASTp/051CACNB2/CACNB2-h-sapiens.csv")
cacnb2.cc <- read_csv("data/raw/BLASTp/051CACNB2/CACNB2-c-carcharias.csv")

# CACNB3
# no reliable hits according to microsynteny in c. carcharias or p. marinus

# CACNB4
cacnb4.hs <- read_csv("data/raw/BLASTp/053CACNB4/CACNB4-h-sapiens.csv")
cacnb4.cc <- read_csv("data/raw/BLASTp/053CACNB4/CACNB4-c-carcharias.csv")

# CACNA2D1
cacna2d1.hs <- read_csv("data/raw/BLASTp/054CACNA2D1/CACNA2D1-h-sapiens.csv")
cacna2d1.cc <- read_csv("data/raw/BLASTp/054CACNA2D1/CACNA2D1-c-carcharias.csv")

# CACNA2D2
# no reliable hits according to microsynteny in c. carcharias or p. marinus

# CACNA2D3
cacna2d3.hs <- read_csv("data/raw/BLASTp/056CACNA2D3/CACNA2D3-h-sapiens.csv")
cacna2d3.pm <- read_csv("data/raw/BLASTp/056CACNA2D3/CACNA2D3-p-marinus.csv")

# CACNA2D4
cacna2d4.hs <- read_csv("data/raw/BLASTp/057CACNA2D4/CACNA2D4-h-sapiens.csv")
cacna2d4.pm <- read_csv("data/raw/BLASTp/057CACNA2D4/CACNA2D4-p-marinus.csv")

# CACNG1
cacng1.hs <- read_csv("data/raw/BLASTp/058CACNG1/CACNG1-h-sapiens.csv")
cacng1.cc <- read_csv("data/raw/BLASTp/058CACNG1/CACNG1-c-carcharias.csv")

# CACNG2
cacng2.hs <- read_csv("data/raw/BLASTp/059CACNG2/CACNG2-h-sapiens.csv")
cacng2.cc <- read_csv("data/raw/BLASTp/059CACNG2/CACNG2-c-carcharias.csv")

# CACNG3
cacng3.hs <- read_csv("data/raw/BLASTp/060CACNG3/CACNG3-h-sapiens.csv")
cacng3.cc <- read_csv("data/raw/BLASTp/060CACNG3/CACNG3-c-carcharias.csv")

# CACNG4
cacng4.hs <- read_csv("data/raw/BLASTp/061CACNG4/CACNG4-h-sapiens.csv")
cacng4.cc <- read_csv("data/raw/BLASTp/061CACNG4/CACNG4-c-carcharias.csv")

# CACNG5
cacng5.hs <- read_csv("data/raw/BLASTp/062CACNG5/CACNG5-h-sapiens.csv")
cacng5.cc <- read_csv("data/raw/BLASTp/062CACNG5/CACNG5-c-carcharias.csv")

# CACNG6
# no reliable hits according to microsynteny in c. carcharias or p. marinus

# CACNG7
# no reliable hits according to microsynteny in c. carcharias or p. marinus

# CACNG8
# no reliable hits according to microsynteny in c. carcharias or p. marinus

# ITPR1
itpr1.hs <- read_csv("data/raw/BLASTp/066ITPR1/ITPR1-h-sapiens.csv")
itpr1.cc <- read_csv("data/raw/BLASTp/066ITPR1/ITPR1-c-carcharias.csv")

# ITPR2
itpr2.hs <- read_csv("data/raw/BLASTp/067ITPR2/ITPR2-h-sapiens.csv")
itpr2.cc <- read_csv("data/raw/BLASTp/067ITPR2/ITPR2-c-carcharias.csv")

# ITPR3
# no reliable hits according to microsynteny in c. carcharias or p. marinus

# CALML3
# no reliable hits according to microsynteny in c. carcharias or p. marinus

# CALM2
calm2.hs <- read_csv("data/raw/BLASTp/070CALM2/CALM2-h-sapiens.csv")
calm2.cc <- read_csv("data/raw/BLASTp/070CALM2/CALM2-c-carcharias.csv")

# CALM3
# no reliable hits according to microsynteny in c. carcharias or p. marinus

# CALM1
calm1.hs <- read_csv("data/raw/BLASTp/072CALM1/CALM1-h-sapiens.csv")
calm1.cc <- read_csv("data/raw/BLASTp/072CALM1/CALM1-c-carcharias.csv")

# CALML6
calml6.hs <- read_csv("data/raw/BLASTp/073CALML6/CALML6-h-sapiens.csv")
calml6.cc <- read_csv("data/raw/BLASTp/073CALML6/CALML6-c-carcharias.csv")

# CALML5
# no reliable hits according to microsynteny in c. carcharias or p. marinus

# CALML4
calml4.hs <- read_csv("data/raw/BLASTp/075CALML4/CALML4-h-sapiens.csv")
calml4.cc <- read_csv("data/raw/BLASTp/075CALML4/CALML4-c-carcharias.csv")

# PPP3CA
ppp3ca.hs <- read_csv("data/raw/BLASTp/076PPP3CA/PPP3CA-h-sapiens.csv")
ppp3ca.cc <- read_csv("data/raw/BLASTp/076PPP3CA/PPP3CA-c-carcharias.csv")

# PPP3CB
ppp3cb.hs <- read_csv("data/raw/BLASTp/077PPP3CB/PPP3CB-h-sapiens.csv")
ppp3cb.pm <- read_csv("data/raw/BLASTp/077PPP3CB/PPP3CB-p-marinus.csv")

# PPP3CC
ppp3cc.hs <- read_csv("data/raw/BLASTp/078PPP3CC/PPP3CC-h-sapiens.csv")
ppp3cc.cc <- read_csv("data/raw/BLASTp/078PPP3CC/PPP3CC-c-carcharias.csv")

# PPP3R1
ppp3r1.hs <- read_csv("data/raw/BLASTp/079PPP3R1/PPP3R1-h-sapiens.csv")
ppp3r1.cc <- read_csv("data/raw/BLASTp/079PPP3R1/PPP3R1-c-carcharias.csv")

# PPP3R2
# no reliable hits according to microsynteny in c. carcharias or p. marinus

# NFATC1
nfatc1.hs <- read_csv("data/raw/BLASTp/081NFATC1/NFATC1-h-sapiens.csv")
nfatc1.cc <- read_csv("data/raw/BLASTp/081NFATC1/NFATC1-c-carcharias.csv")

# NFATC2
nfatc2.hs <- read_csv("data/raw/BLASTp/082NFATC2/NFATC2-h-sapiens.csv")
nfatc2.cc <- read_csv("data/raw/BLASTp/082NFATC2/NFATC2-c-carcharias.csv")

# NFATC3
nfatc3.hs <- read_csv("data/raw/BLASTp/083NFATC3/NFATC3-h-sapiens.csv")
nfatc3.cc <- read_csv("data/raw/BLASTp/083NFATC3/NFATC3-c-carcharias.csv")

# NFATC4
# no reliable hits according to microsynteny in c. carcharias or p. marinus

# RGS2
# no reliable hits according to microsynteny in c. carcharias or p. marinus

# RCAN1
# no reliable hits according to microsynteny in c. carcharias or p. marinus

# CAMK2
camk2.hs <- read_csv("data/raw/BLASTp/087CAMK2/CAMK2-h-sapiens.csv")
camk2.cc <- read_csv("data/raw/BLASTp/087CAMK2/CAMK2-c-carcharias.csv")

# PRKAA1
prkaa1.hs <- read_csv("data/raw/BLASTp/088PRKAA1/PRKAA1-h-sapiens.csv")
prkaa1.cc <- read_csv("data/raw/BLASTp/088PRKAA1/PRKAA1-c-carcharias.csv")

# PRKAA2
prkaa2.hs <- read_csv("data/raw/BLASTp/089PRKAA2/PRKAA2-h-sapiens.csv")
prkaa2.cc <- read_csv("data/raw/BLASTp/089PRKAA2/PRKAA2-c-carcharias.csv")

# PRKAB1
prkab1.hs <- read_csv("data/raw/BLASTp/090PRKAB1/PRKAB1-h-sapiens.csv")
prkab1.cc <- read_csv("data/raw/BLASTp/090PRKAB1/PRKAB1-c-carcharias.csv")

# PRKAB2
# no reliable hits according to microsynteny in c. carcharias or p. marinus

# PRKAG1
prkag1.hs <- read_csv("data/raw/BLASTp/092PRKAG1/PRKAG1-h-sapiens.csv")
prkag1.cc <- read_csv("data/raw/BLASTp/092PRKAG1/PRKAG1-c-carcharias.csv")

# PRKAG3
# no reliable hits according to microsynteny in c. carcharias or p. marinus

# PRKAG2
prkag2.hs <- read_csv("data/raw/BLASTp/094PRKAG2/PRKAG2-h-sapiens.csv")
prkag2.cc <- read_csv("data/raw/BLASTp/094PRKAG2/PRKAG2-c-carcharias.csv")

# CAMK1D
camk1d.hs <- read_csv("data/raw/BLASTp/095CAMK1D/CAMK1D-h-sapiens.csv")
camk1d.cc <- read_csv("data/raw/BLASTp/095CAMK1D/CAMK1D-c-carcharias.csv")

# CAMK1G
camk1g.hs <- read_csv("data/raw/BLASTp/096CAMK1G/CAMK1G-h-sapiens.csv")
camk1g.cc <- read_csv("data/raw/BLASTp/096CAMK1G/CAMK1G-c-carcharias.csv")

# CAMK1
camk1.hs <- read_csv("data/raw/BLASTp/097CAMK1/CAMK1-h-sapiens.csv")
camk1.cc <- read_csv("data/raw/BLASTp/097CAMK1/CAMK1-c-carcharias.csv")

# CAMK2A
camk2a.hs <- read_csv("data/raw/BLASTp/098CAMK2A/CAMK2A-h-sapiens.csv")
camk2a.cc <- read_csv("data/raw/BLASTp/098CAMK2A/CAMK2A-c-carcharias.csv")

# CAMK2D
camk2d.hs <- read_csv("data/raw/BLASTp/099CAMK2D/CAMK2D-h-sapiens.csv")
camk2d.cc <- read_csv("data/raw/BLASTp/099CAMK2D/CAMK2D-c-carcharias.csv")

# CAMK2B
# no reliable hits according to microsynteny in c. carcharias or p. marinus

# CAMK2G
camk2g.hs <- read_csv("data/raw/BLASTp/101CAMK2G/CAMK2G-h-sapiens.csv")
camk2g.cc <- read_csv("data/raw/BLASTp/101CAMK2G/CAMK2G-c-carcharias.csv")

# CAMK4
camk4.hs <- read_csv("data/raw/BLASTp/102CAMK4/CAMK4-h-sapiens.csv")
camk4.cc <- read_csv("data/raw/BLASTp/102CAMK4/CAMK4-c-carcharias.csv")

# NOS3
# no reliable hits according to microsynteny in c. carcharias or p. marinus

# GUCY1A2
gucy1a2.hs <- read_csv("data/raw/BLASTp/104GUCY1A2/GUCY1A2-h-sapiens.csv")
gucy1a2.pm <- read_csv("data/raw/BLASTp/104GUCY1A2/GUCY1A2-p-marinus.csv")

# GUCY1A1
# no reliable hits according to microsynteny in c. carcharias or p. marinus

# GUCY1B1
gucy1b1.hs <- read_csv("data/raw/BLASTp/106GUCY1B1/GUCY1B1-h-sapiens.csv")
gucy1b1.cc <- read_csv("data/raw/BLASTp/106GUCY1B1/GUCY1B1-c-carcharias.csv")

# NPPA
# no reliable hits according to microsynteny in c. carcharias or p. marinus

# NPR1
npr1.hs <- read_csv("data/raw/BLASTp/108NPR1/NPR1-h-sapiens.csv")
npr1.cc <- read_csv("data/raw/BLASTp/108NPR1/NPR1-c-carcharias.csv")

# NPR2
# no reliable hits according to microsynteny in c. carcharias or p. marinus

# MYLK
# no reliable hits according to microsynteny in c. carcharias or p. marinus

# MYLK2
# no reliable hits according to microsynteny in c. carcharias or p. marinus

# MYLK3
mylk3.hs <- read_csv("data/raw/BLASTp/112MYLK3/MYLK3-h-sapiens.csv")
mylk3.cc <- read_csv("data/raw/BLASTp/112MYLK3/MYLK3-c-carcharias.csv")

# MYLK4
# no reliable hits according to microsynteny in c. carcharias or p. marinus

# MYL6B
# no reliable hits according to microsynteny in c. carcharias or p. marinus

# MYL6
# no reliable hits according to microsynteny in c. carcharias or p. marinus

# MYL9
myl9.hs <- read_csv("data/raw/BLASTp/116MYL9/MYL9-h-sapiens.csv")
myl9.cc <- read_csv("data/raw/BLASTp/116MYL9/MYL9-c-carcharias.csv")

# ACTG1
# no reliable hits according to microsynteny in c. carcharias or p. marinus

# ACTB
# no reliable hits according to microsynteny in c. carcharias or p. marinus

# RHOA
rhoa.hs <- read_csv("data/raw/BLASTp/119RHOA/RHOA-h-sapiens.csv")
rhoa.cc <- read_csv("data/raw/BLASTp/119RHOA/RHOA-c-carcharias.csv")

# ROCK1
rock1.hs <- read_csv("data/raw/BLASTp/120ROCK1/ROCK1-h-sapiens.csv")
rock1.cc <- read_csv("data/raw/BLASTp/120ROCK1/ROCK1-c-carcharias.csv")

# ROCK2
rock2.hs <- read_csv("data/raw/BLASTp/121ROCK2/ROCK2-h-sapiens.csv")
rock2.cc <- read_csv("data/raw/BLASTp/121ROCK2/ROCK2-c-carcharias.csv")

# PPP1CA
# no reliable hits according to microsynteny in c. carcharias or p. marinus

# PPP1CB
ppp1cb.hs <- read_csv("data/raw/BLASTp/123PPP1CB/PPP1CB-h-sapiens.csv")
ppp1cb.cc <- read_csv("data/raw/BLASTp/123PPP1CB/PPP1CB-c-carcharias.csv")

# PPP1CC
ppp1cc.hs <- read_csv("data/raw/BLASTp/124PPP1CC/PPP1CC-h-sapiens.csv")
ppp1cc.pm <- read_csv("data/raw/BLASTp/124PPP1CC/PPP1CC-p-marinus.csv")

# PPP1R12A
ppp1r12a.hs <- read_csv("data/raw/BLASTp/125PPP1R12A/PPP1R12A-h-sapiens.csv")
ppp1r12a.pm <- read_csv("data/raw/BLASTp/125PPP1R12A/PPP1r12A-p-marinus.csv")

# PPP1R12B
# no reliable hits according to microsynteny in c. carcharias or p. marinus

# PPP1R12C
# no reliable hits according to microsynteny in c. carcharias or p. marinus

# GNAS
# no reliable hits according to microsynteny in c. carcharias or p. marinus

# ADCY1
adcy1.hs <- read_csv("data/raw/BLASTp/129ADCY1/ADCY1-h-sapiens.csv")
adcy1.cc <- read_csv("data/raw/BLASTp/129ADCY1/ADCY1-c-carcharias.csv")

# ADCY2
adcy2.hs <- read_csv("data/raw/BLASTp/130ADCY2/ADCY2-h-sapiens.csv")
adcy2.cc <- read_csv("data/raw/BLASTp/130ADCY2/ADCY2-c-carcharias.csv")

# ADCY3
# no reliable hits according to microsynteny in c. carcharias or p. marinus

# ADCY4
# no reliable hits according to microsynteny in c. carcharias or p. marinus

# ADCY5
adcy5.hs <- read_csv("data/raw/BLASTp/133ADCY5/ADCY5-h-sapiens.csv")
adcy5.cc <- read_csv("data/raw/BLASTp/133ADCY5/ADCY5-c-carcharias.csv")

# ADCY6
# no reliable hits according to microsynteny in c. carcharias or p. marinus

# ADCY7
adcy7.hs <- read_csv("data/raw/BLASTp/135ADCY7/ADCY7-h-sapiens.csv")
adcy7.cc <- read_csv("data/raw/BLASTp/135ADCY7/ADCY7-c-carcharias.csv")

# ADCY8
adcy8.hs <- read_csv("data/raw/BLASTp/136ADCY8/ADCY8-h-sapiens.csv")
adcy8.cc <- read_csv("data/raw/BLASTp/136ADCY8/ADCY8-c-carcharias.csv")

# ADCY9
adcy9.hs <- read_csv("data/raw/BLASTp/137ADCY9/ADCY9-h-sapiens.csv")
adcy9.cc <- read_csv("data/raw/BLASTp/137ADCY9/ADCY9-c-carcharias.csv")

# PRKACA
prkaca.hs <- read_csv("data/raw/BLASTp/138PRKACA/PRKACA-h-sapiens.csv")
prkaca.cc <- read_csv("data/raw/BLASTp/138PRKACA/PRKACA-c-carcharias.csv")

# PRKACB
prkacb.hs <- read_csv("data/raw/BLASTp/139PRKACB/PRKACB-h-sapiens.csv")
prkacb.cc <- read_csv("data/raw/BLASTp/139PRKACB/PRKACB-c-carcharias.csv")

# PRKACB
prkacb.hs <- read_csv("data/raw/BLASTp/139PRKACB/PRKACB-h-sapiens.csv")
prkacb.cc <- read_csv("data/raw/BLASTp/139PRKACB/PRKACB-c-carcharias.csv")

# PRKACG
# no reliable hits according to microsynteny in c. carcharias or p. marinus

# GNAI1
gnai1.hs <- read_csv("data/raw/BLASTp/141GNAI1/GNAI1-h-sapiens.csv")
gnai1.cc <- read_csv("data/raw/BLASTp/141GNAI1/GNAI1-c-carcharias.csv")

# GNAI3
# no reliable hits according to microsynteny in c. carcharias or p. marinus

# GNAI2
# no reliable hits according to microsynteny in c. carcharias or p. marinus

# GNAO1
# no reliable hits according to microsynteny in c. carcharias or p. marinus

# PIK3CG
# no reliable hits according to microsynteny in c. carcharias or p. marinus

# PIK3R5
pik3r5.hs <- read_csv("data/raw/BLASTp/146PIK3R5/PIK3R5-h-sapiens.csv")
pik3r5.pm <- read_csv("data/raw/BLASTp/146PIK3R5/PIK3R5-p-marinus.csv")

# PIK3R6
# no reliable hits according to microsynteny in c. carcharias or p. marinus

# SRC
src.hs <- read_csv("data/raw/BLASTp/148SRC/SRC-h-sapiens.csv")
src.pm <- read_csv("data/raw/BLASTp/148SRC/SRC-c-carcharias.csv")

# KCNJ3
kcnj3.hs <- read_csv("data/raw/BLASTp/149KCNJ3/KCNJ3-h-sapiens.csv")
kcnj3.cc <- read_csv("data/raw/BLASTp/149KCNJ3/KCNJ3-c-carcharias.csv")

# KCNJ6
kcnj6.hs <- read_csv("data/raw/BLASTp/150KCNJ6/KCNJ6-h-sapiens.csv")
kcnj6.cc <- read_csv("data/raw/BLASTp/150KCNJ6/KCNJ6-c-carcharias.csv")

# KCNJ9
# no reliable hits according to microsynteny in c. carcharias or p. marinus

# KCNJ5
kcnj5.hs <- read_csv("data/raw/BLASTp/152KCNJ5/KCNJ5-h-sapiens.csv")
kcnj5.cc <- read_csv("data/raw/BLASTp/152KCNJ5/KCNJ5-c-carcharias.csv")

# EGFR
egfr.hs <- read_csv("data/raw/BLASTp/153EGFR/EGFR-h-sapiens.csv")
egfr.cc <- read_csv("data/raw/BLASTp/153EGFR/EGFR-c-carcharias.csv")

# CDKN1A
# no reliable hits according to microsynteny in c. carcharias or p. marinus















