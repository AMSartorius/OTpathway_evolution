######################################
##### SET UP WORKING ENVIRONMENT #####
######################################

rm(list = ls()) # delete all objects in the workspace
gc(reset = T)   # resest memory (especially useful when working with large data sets)

options(stringsAsFactors = F) # disables automatic conversion of char. strings into factors
Sys.setenv(LANG = "en")


### .......... loading packages
library(data.table) # for loading data
library(tidyr)      # always useful
library(ggplot2)    # plotting
library(stringr)    # replacing strings
library(dplyr)
library(readr)
library(readxl)
library(writexl)
library(knitr)

###################################

## ADJUST THESE PATHS TO YOUR LOCAL DIRECTORY STRUCTURE
path = "C:/Users/daedh/Desktop/OTpathway_evolution-main/OTpathway_evolution-main/"
temp <- list.files(path = paste0(path, "data/raw/BLASTp/all/"), full.names = T)
temp3 <- list.files(path = paste0(path, "data/raw/BLASTp/all/"), full.names = F)


L1 <- list()
for (i in 1:length(temp)) {
  L1[[i]] = data.frame(read_csv(temp[i])) # read in csvs into list
  names(L1)[i] <- paste(temp3[i])
  names(L1[[i]]) <- make.names(names(L1[[i]]), unique = TRUE)
  #L1 <- tools::file_path_sans_ext(L1) 
}

# list2env(L1 ,.GlobalEnv)   # unlists L1 list with dataframes and pate in the environment

pat = data.frame(c("Ciona savignyi", "Ciona intestinalis", "Branchiostoma lanceolatum", "Branchiostoma floridae", "Strongylocentrotus purpuratus", 
                   "Asterias rubens", "Saccoglossus kowalevskii", "Xenoturbella bocki", "Symsagittifera roscoffensis", "Caenorhabditis elegans", 
                   "Priapulus caudatus", "Octopus vulgaris", "Hydra vulgaris", "Acropora millepora", "Clytia hemisphaerica", "Ephydatia muelleri", 
                   "Amphimedon queenslandica", "Salpingoeca rosetta", "Monosiga brevicollis", "Fusarium culmorum", "Aspergillus nidulans", 
                   "Saccharomyces cerevisiae", "Dictyostelium discoideum", "Spironucleus salmonicida", "Methanolobus zinderi", "Escherichia coli"))
pat$c..Ciona.savignyi....Ciona.intestinalis....Branchiostoma.lanceolatum... <- as.character(pat$c..Ciona.savignyi....Ciona.intestinalis....Branchiostoma.lanceolatum...)

names <- data.frame(names(L1))
View(names)


# loops through ONE ref gene csv file and returns list with max scores, accession lengths, 
# protein length of ref gene for each species, and results


### ----------------------------------



# -------------------------------------------- OXT
df.001oxt <- data.frame(1:nrow(pat))
L.001oxt  <- list()

for (k in 1:nrow(pat)) {
  L.001oxt[[k]] <- data.frame(L1[[1]][grep(pat[k,], L1[[1]][["Scientific.Name"]]), ])
  df.001oxt[k,] <- max((data.frame(L1[[1]][grep(pat[k,], L1[[1]][["Scientific.Name"]]), ]))$Max.Score)
  df.001oxt[k, 2] <- first(L.001oxt[[k]][["Acc..Len"]][L.001oxt[[k]][["Max.Score"]] == max(L.001oxt[[k]][["Max.Score"]])])
  df.001oxt[, 3] <- rep(126, 26)
  df.001oxt[, 4] <- df.001oxt[, 1] / ((df.001oxt[, 3]) + df.001oxt[, 2])
  colnames(df.001oxt) <- c("maxscore", "p.length", "ref.p.length", "OXT.maxscore.scaled")
}


# -------------------------------------------- OXTR
df.002oxtr <- data.frame(1:nrow(pat))
L.002oxtr  <- list()

for (k in 1:nrow(pat)) {
  L.002oxtr[[k]] <- data.frame(L1[[2]][grep(pat[k,], L1[[2]][["Scientific.Name"]]), ])
  df.002oxtr[k,] <- max((data.frame(L1[[2]][grep(pat[k,], L1[[2]][["Scientific.Name"]]), ]))$Max.Score)
  df.002oxtr[k, 2] <- first(L.002oxtr[[k]][["Acc..Len"]][L.002oxtr[[k]][["Max.Score"]] == max(L.002oxtr[[k]][["Max.Score"]])])
  df.002oxtr[, 3] <- rep(431, 26)
  df.002oxtr[, 4] <- df.002oxtr[, 1] / ((df.002oxtr[, 3]) + df.002oxtr[, 2])
  colnames(df.002oxtr) <- c("maxscore", "p.length", "ref.p.length", "OXTR.maxscore.scaled")
}


# -------------------------------------------- GNAQ
df.003gnaq <- data.frame(1:nrow(pat))
L.003gnaq  <- list()

for (k in 1:nrow(pat)) {
  L.003gnaq[[k]] <- data.frame(L1[[3]][grep(pat[k,], L1[[3]][["Scientific.Name"]]), ])
  df.003gnaq[k,] <- max((data.frame(L1[[3]][grep(pat[k,], L1[[3]][["Scientific.Name"]]), ]))$Max.Score)
  df.003gnaq[k, 2] <- first(L.003gnaq[[k]][["Acc..Len"]][L.003gnaq[[k]][["Max.Score"]] == max(L.003gnaq[[k]][["Max.Score"]])])
  df.003gnaq[, 3] <- rep(359, 26)
  df.003gnaq[, 4] <- df.003gnaq[, 1] / ((df.003gnaq[, 3]) + df.003gnaq[, 2])
  colnames(df.003gnaq) <- c("maxscore", "p.length", "ref.p.length", "GNAQ.maxscore.scaled")
}


# -------------------------------------------- HRAS
df.004hras <- data.frame(1:nrow(pat))
L.004hras <- list()

for (k in 1:nrow(pat)) {
  L.004hras[[k]] <- data.frame(L1[[4]][grep(pat[k,], L1[[4]][["Scientific.Name"]]), ])
  df.004hras[k,] <- max((data.frame(L1[[4]][grep(pat[k,], L1[[4]][["Scientific.Name"]]), ]))$Max.Score)
  df.004hras[k, 2] <- first(L.004hras[[k]][["Acc..Len"]][L.004hras[[k]][["Max.Score"]] == max(L.004hras[[k]][["Max.Score"]])])
  df.004hras[, 3] <- rep(188, 26)
  df.004hras[, 4] <- df.004hras[, 1] / ((df.004hras[, 3]) + df.004hras[, 2])
  colnames(df.004hras) <- c("maxscore", "p.length", "ref.p.length", "HRAS.maxscore.scaled")
}

# -------------------------------------------- KRAS
df.005kras <- data.frame(1:nrow(pat))
L.005kras <- list()

for (k in 1:nrow(pat)) {
  L.005kras[[k]] <- data.frame(L1[[5]][grep(pat[k,], L1[[5]][["Scientific.Name"]]), ])
  df.005kras[k,] <- max((data.frame(L1[[5]][grep(pat[k,], L1[[5]][["Scientific.Name"]]), ]))$Max.Score)
  df.005kras[k, 2] <- first(L.005kras[[k]][["Acc..Len"]][L.005kras[[k]][["Max.Score"]] == max(L.005kras[[k]][["Max.Score"]])])
  df.005kras[, 3] <- rep(207, 26)
  df.005kras[, 4] <- df.005kras[, 1] / ((df.005kras[, 3]) + df.005kras[, 2])
  colnames(df.005kras) <- c("maxscore", "p.length", "ref.p.length", "KRAS.maxscore.scaled")
}


# -------------------------------------------- NRAS
df.006nras <- data.frame(1:nrow(pat))
L.006nras <- list()

for (k in 1:nrow(pat)) {
  L.006nras[[k]] <- data.frame(L1[[6]][grep(pat[k,], L1[[6]][["Scientific.Name"]]), ])
  df.006nras[k,] <- max((data.frame(L1[[6]][grep(pat[k,], L1[[6]][["Scientific.Name"]]), ]))$Max.Score)
  df.006nras[k, 2] <- first(L.006nras[[k]][["Acc..Len"]][L.006nras[[k]][["Max.Score"]] == max(L.006nras[[k]][["Max.Score"]])])
  df.006nras[, 3] <- rep(189, 26)
  df.006nras[, 4] <- df.006nras[, 1] / ((df.006nras[, 3]) + df.006nras[, 2])
  colnames(df.006nras) <- c("maxscore", "p.length", "ref.p.length", "NRAS.maxscore.scaled")
}


# -------------------------------------------- RAF1
df.007raf1 <- data.frame(1:nrow(pat))
L.007raf1 <- list()

for (k in 1:nrow(pat)) {
  L.007raf1[[k]] <- data.frame(L1[[7]][grep(pat[k,], L1[[7]][["Scientific.Name"]]), ])
  df.007raf1[k,] <- max((data.frame(L1[[7]][grep(pat[k,], L1[[7]][["Scientific.Name"]]), ]))$Max.Score)
  df.007raf1[k, 2] <- first(L.007raf1[[k]][["Acc..Len"]][L.007raf1[[k]][["Max.Score"]] == max(L.007raf1[[k]][["Max.Score"]])])
  df.007raf1[, 3] <- rep(650, 26)
  df.007raf1[, 4] <- df.007raf1[, 1] / ((df.007raf1[, 3]) + df.007raf1[, 2])
  colnames(df.007raf1) <- c("maxscore", "p.length", "ref.p.length", "RAF1.maxscore.scaled")
}


# -------------------------------------------- MAP2K1
df.008map2k1 <- data.frame(1:nrow(pat))
L.008map2k1 <- list()

for (k in 1:nrow(pat)) {
  L.008map2k1[[k]] <- data.frame(L1[[8]][grep(pat[k,], L1[[8]][["Scientific.Name"]]), ])
  df.008map2k1[k,] <- max((data.frame(L1[[8]][grep(pat[k,], L1[[8]][["Scientific.Name"]]), ]))$Max.Score)
  df.008map2k1[k, 2] <- first(L.008map2k1[[k]][["Acc..Len"]][L.008map2k1[[k]][["Max.Score"]] == max(L.008map2k1[[k]][["Max.Score"]])])
  df.008map2k1[, 3] <- rep(393, 26)
  df.008map2k1[, 4] <- df.008map2k1[, 1] / ((df.008map2k1[, 3]) + df.008map2k1[, 2])
  colnames(df.008map2k1) <- c("maxscore", "p.length", "ref.p.length",  "MAP2K1.maxscore.scaled")
}


# -------------------------------------------- MAP2K2
df.009map2k2 <- data.frame(1:nrow(pat))
L.009map2k2 <- list()

for (k in 1:nrow(pat)) {
  L.009map2k2[[k]] <- data.frame(L1[[9]][grep(pat[k,], L1[[9]][["Scientific.Name"]]), ])
  df.009map2k2[k,] <- max((data.frame(L1[[9]][grep(pat[k,], L1[[9]][["Scientific.Name"]]), ]))$Max.Score)
  df.009map2k2[k, 2] <- first(L.009map2k2[[k]][["Acc..Len"]][L.009map2k2[[k]][["Max.Score"]] == max(L.009map2k2[[k]][["Max.Score"]])])
  df.009map2k2[, 3] <- rep(400, 26)
  df.009map2k2[, 4] <- df.009map2k2[, 1] / ((df.009map2k2[, 3]) + df.009map2k2[, 2])
  colnames(df.009map2k2) <- c("maxscore", "p.length", "ref.p.length",  "MAP2K2.maxscore.scaled")
}


# -------------------------------------------- MAPK1
df.010mapk1 <- data.frame(1:nrow(pat))
L.010mapk1 <- list()

for (k in 1:nrow(pat)) {
  L.010mapk1[[k]] <- data.frame(L1[[10]][grep(pat[k,], L1[[10]][["Scientific.Name"]]), ])
  df.010mapk1[k,] <- max((data.frame(L1[[10]][grep(pat[k,], L1[[10]][["Scientific.Name"]]), ]))$Max.Score)
  df.010mapk1[k, 2] <- first(L.010mapk1[[k]][["Acc..Len"]][L.010mapk1[[k]][["Max.Score"]] == max(L.010mapk1[[k]][["Max.Score"]])])
  df.010mapk1[, 3] <- rep(364, 26)
  df.010mapk1[, 4] <- df.010mapk1[, 1] / ((df.010mapk1[, 3]) + df.010mapk1[, 2])
  colnames(df.010mapk1) <- c("maxscore", "p.length", "ref.p.length",  "MAPK1.maxscore.scaled")
}


# -------------------------------------------- PLA2G4A
df.013pla2g4a <- data.frame(1:nrow(pat))
L.013pla2g4a <- list()

for (k in 1:nrow(pat)) {
  L.013pla2g4a[[k]] <- data.frame(L1[[11]][grep(pat[k,], L1[[11]][["Scientific.Name"]]), ])
  df.013pla2g4a[k,] <- max((data.frame(L1[[11]][grep(pat[k,], L1[[11]][["Scientific.Name"]]), ]))$Max.Score)
  df.013pla2g4a[k, 2] <- first(L.013pla2g4a[[k]][["Acc..Len"]][L.013pla2g4a[[k]][["Max.Score"]] == max(L.013pla2g4a[[k]][["Max.Score"]])])
  df.013pla2g4a[, 3] <- rep(780, 26)
  df.013pla2g4a[, 4] <- df.013pla2g4a[, 1] / ((df.013pla2g4a[, 3]) + df.013pla2g4a[, 2])
  colnames(df.013pla2g4a) <- c("maxscore", "p.length", "ref.p.length",  "PLA2G4A.maxscore.scaled")
}


# -------------------------------------------- PLA2G4F
df.018pla2g4f <- data.frame(1:nrow(pat))
L.018pla2g4f <- list()

for (k in 1:nrow(pat)) {
  L.018pla2g4f[[k]] <- data.frame(L1[[12]][grep(pat[k,], L1[[12]][["Scientific.Name"]]), ])
  df.018pla2g4f[k,] <- max((data.frame(L1[[12]][grep(pat[k,], L1[[12]][["Scientific.Name"]]), ]))$Max.Score)
  df.018pla2g4f[k, 2] <- first(L.018pla2g4f[[k]][["Acc..Len"]][L.018pla2g4f[[k]][["Max.Score"]] == max(L.018pla2g4f[[k]][["Max.Score"]])])
  df.018pla2g4f[, 3] <- rep(856, 26)
  df.018pla2g4f[, 4] <- df.018pla2g4f[, 1] / ((df.018pla2g4f[, 3]) + df.018pla2g4f[, 2])
  colnames(df.018pla2g4f) <- c("maxscore", "p.length", "ref.p.length",  "PLA2G4F.maxscore.scaled")
}


# -------------------------------------------- PTGS2
df.019ptgs2 <- data.frame(1:nrow(pat))
L.019ptgs2 <- list()

for (k in 1:nrow(pat)) {
  L.019ptgs2[[k]] <- data.frame(L1[[13]][grep(pat[k,], L1[[13]][["Scientific.Name"]]), ])
  df.019ptgs2[k,] <- max((data.frame(L1[[13]][grep(pat[k,], L1[[13]][["Scientific.Name"]]), ]))$Max.Score)
  df.019ptgs2[k, 2] <- first(L.019ptgs2[[k]][["Acc..Len"]][L.019ptgs2[[k]][["Max.Score"]] == max(L.019ptgs2[[k]][["Max.Score"]])])
  df.019ptgs2[, 3] <- rep(731, 26)
  df.019ptgs2[, 4] <- df.019ptgs2[, 1] / ((df.019ptgs2[, 3]) + df.019ptgs2[, 2])
  colnames(df.019ptgs2) <- c("maxscore", "p.length", "ref.p.length",  "PTGS2.maxscore.scaled")
}


# -------------------------------------------- MAP2K5
df.020map2k5 <- data.frame(1:nrow(pat))
L.020map2k5 <- list()

for (k in 1:nrow(pat)) {
  L.020map2k5[[k]] <- data.frame(L1[[14]][grep(pat[k,], L1[[14]][["Scientific.Name"]]), ])
  df.020map2k5[k,] <- max((data.frame(L1[[14]][grep(pat[k,], L1[[14]][["Scientific.Name"]]), ]))$Max.Score)
  df.020map2k5[k, 2] <- first(L.020map2k5[[k]][["Acc..Len"]][L.020map2k5[[k]][["Max.Score"]] == max(L.020map2k5[[k]][["Max.Score"]])])
  df.020map2k5[, 3] <- rep(444, 26)
  df.020map2k5[, 4] <- df.020map2k5[, 1] / ((df.020map2k5[, 3]) + df.020map2k5[, 2])
  colnames(df.020map2k5) <- c("maxscore", "p.length", "ref.p.length",  "MAP2K5.maxscore.scaled")
}


# -------------------------------------------- JUN
df.022jun <- data.frame(1:nrow(pat))
L.022jun <- list()

for (k in 1:nrow(pat)) {
  L.022jun[[k]] <- data.frame(L1[[15]][grep(pat[k,], L1[[15]][["Scientific.Name"]]), ])
  df.022jun[k,] <- max((data.frame(L1[[15]][grep(pat[k,], L1[[15]][["Scientific.Name"]]), ]))$Max.Score)
  df.022jun[k, 2] <- first(L.022jun[[k]][["Acc..Len"]][L.022jun[[k]][["Max.Score"]] == max(L.022jun[[k]][["Max.Score"]])])
  df.022jun[, 3] <- rep(318, 26)
  df.022jun[, 4] <- df.022jun[, 1] / ((df.022jun[, 3]) + df.022jun[, 2])
  colnames(df.022jun) <- c("maxscore", "p.length", "ref.p.length",  "JUN.maxscore.scaled")
}


# -------------------------------------------- FOS
df.023fos <- data.frame(1:nrow(pat))
L.023fos <- list()

for (k in 1:nrow(pat)) {
  L.023fos[[k]] <- data.frame(L1[[16]][grep(pat[k,], L1[[16]][["Scientific.Name"]]), ])
  df.023fos[k,] <- max((data.frame(L1[[16]][grep(pat[k,], L1[[16]][["Scientific.Name"]]), ]))$Max.Score)
  df.023fos[k, 2] <- first(L.023fos[[k]][["Acc..Len"]][L.023fos[[k]][["Max.Score"]] == max(L.023fos[[k]][["Max.Score"]])])
  df.023fos[, 3] <- rep(301, 26)
  df.023fos[, 4] <- df.023fos[, 1] / ((df.023fos[, 3]) + df.023fos[, 2])
  colnames(df.023fos) <- c("maxscore", "p.length", "ref.p.length",  "FOS.maxscore.scaled")
}


# -------------------------------------------- 024MEF2C
df.024mef2c <- data.frame(1:nrow(pat))
L.024mef2c <- list()

for (k in 1:nrow(pat)) {
  L.024mef2c[[k]] <- data.frame(L1[[17]][grep(pat[k,], L1[[17]][["Scientific.Name"]]), ])
  df.024mef2c[k,] <- max((data.frame(L1[[17]][grep(pat[k,], L1[[17]][["Scientific.Name"]]), ]))$Max.Score)
  df.024mef2c[k, 2] <- first(L.024mef2c[[k]][["Acc..Len"]][L.024mef2c[[k]][["Max.Score"]] == max(L.024mef2c[[k]][["Max.Score"]])])
  df.024mef2c[, 3] <- rep(465, 26)
  df.024mef2c[, 4] <- df.024mef2c[, 1] / ((df.024mef2c[, 3]) + df.024mef2c[, 2])
  colnames(df.024mef2c) <- c("maxscore", "p.length", "ref.p.length",  "MEF2C.maxscore.scaled")
}


# -------------------------------------------- 025CCND1
df.025ccnd1 <- data.frame(1:nrow(pat))
L.025ccnd1 <- list()

for (k in 1:nrow(pat)) {
  L.025ccnd1[[k]] <- data.frame(L1[[18]][grep(pat[k,], L1[[18]][["Scientific.Name"]]), ])
  df.025ccnd1[k,] <- max((data.frame(L1[[18]][grep(pat[k,], L1[[18]][["Scientific.Name"]]), ]))$Max.Score)
  df.025ccnd1[k, 2] <- first(L.025ccnd1[[k]][["Acc..Len"]][L.025ccnd1[[k]][["Max.Score"]] == max(L.025ccnd1[[k]][["Max.Score"]])])
  df.025ccnd1[, 3] <- rep(292, 26)
  df.025ccnd1[, 4] <- df.025ccnd1[, 1] / ((df.025ccnd1[, 3]) + df.025ccnd1[, 2])
  colnames(df.025ccnd1) <- c("maxscore", "p.length", "ref.p.length",  "CCND1.maxscore.scaled")
}


# -------------------------------------------- 026ELK1
df.026elk1 <- data.frame(1:nrow(pat))
L.026elk1 <- list()

for (k in 1:nrow(pat)) {
  L.026elk1[[k]] <- data.frame(L1[[19]][grep(pat[k,], L1[[19]][["Scientific.Name"]]), ])
  df.026elk1[k,] <- max((data.frame(L1[[19]][grep(pat[k,], L1[[19]][["Scientific.Name"]]), ]))$Max.Score)
  df.026elk1[k, 2] <- first(L.026elk1[[k]][["Acc..Len"]][L.026elk1[[k]][["Max.Score"]] == max(L.026elk1[[k]][["Max.Score"]])])
  df.026elk1[, 3] <- rep(372, 26)
  df.026elk1[, 4] <- df.026elk1[, 1] / ((df.026elk1[, 3]) + df.026elk1[, 2])
  colnames(df.026elk1) <- c("maxscore", "p.length", "ref.p.length",  "ELK1.maxscore.scaled")
}


# -------------------------------------------- 028RYR2
df.028ryr2 <- data.frame(1:nrow(pat))
L.028ryr2 <- list()

for (k in 1:nrow(pat)) {
  L.028ryr2[[k]] <- data.frame(L1[[20]][grep(pat[k,], L1[[20]][["Scientific.Name"]]), ])
  df.028ryr2[k,] <- max((data.frame(L1[[20]][grep(pat[k,], L1[[20]][["Scientific.Name"]]), ]))$Max.Score)
  df.028ryr2[k, 2] <- first(L.028ryr2[[k]][["Acc..Len"]][L.028ryr2[[k]][["Max.Score"]] == max(L.028ryr2[[k]][["Max.Score"]])])
  df.028ryr2[, 3] <- rep(4953, 26)
  df.028ryr2[, 4] <- df.028ryr2[, 1] / ((df.028ryr2[, 3]) + df.028ryr2[, 2])
  colnames(df.028ryr2) <- c("maxscore", "p.length", "ref.p.length",  "RYR2.maxscore.scaled")
}


# -------------------------------------------- 029RYR3
df.029RYR3 <- data.frame(1:nrow(pat))
L.029RYR3 <- list()

for (k in 1:nrow(pat)) {
  L.029RYR3[[k]] <- data.frame(L1[[21]][grep(pat[k,], L1[[21]][["Scientific.Name"]]), ])
  df.029RYR3[k,] <- max((data.frame(L1[[21]][grep(pat[k,], L1[[21]][["Scientific.Name"]]), ]))$Max.Score)
  df.029RYR3[k, 2] <- first(L.029RYR3[[k]][["Acc..Len"]][L.029RYR3[[k]][["Max.Score"]] == max(L.029RYR3[[k]][["Max.Score"]])])
  df.029RYR3[, 3] <- rep(4881, 26)
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
  df.030CD38[, 3] <- rep(291, 26)
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
  df.032KCNJ2[, 3] <- rep(426, 26)
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
  df.035KCNJ4[, 3] <- rep(426, 26)
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
  df.037PLCB1[, 3] <- rep(1260, 26)
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
  df.040PLCB4[, 3] <- rep(1207, 26)
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
  df.041PRKCA[, 3] <- rep(664, 26)
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
  df.042PRKCB[, 3] <- rep(671, 26)
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
  df.044EEF2K[, 3] <- rep(715, 26)
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
  df.045EEF2[, 3] <- rep(858, 26)
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
  df.046CACNA1C[, 3] <- rep(2200, 26)
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
  df.047CACNA1D[, 3] <- rep(2130, 26)
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
  df.050CACNB1[, 3] <- rep(800, 26)
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
  df.051CACNB2[, 3] <- rep(613, 26)
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
  df.053CACNB4[, 3] <- rep(476, 26)
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
  df.054CACNA2D1[, 3] <- rep(1086, 26)
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
  df.056CACNA2D3[, 3] <- rep(1085, 26)
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
  df.057CACNA2D4[, 3] <- rep(1195, 26)
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
  df.058CACNG1[, 3] <- rep(227, 26)
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
  df.059CACNG2[, 3] <- rep(322, 26)
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
  df.060CACNG3[, 3] <- rep(317, 26)
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
  df.061CACNG4[, 3] <- rep(328, 26)
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
  df.062CACNG5[, 3] <- rep(279, 26)
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
  df.066ITPR1[, 3] <- rep(2757, 26)
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
  df.067ITPR2[, 3] <- rep(2691, 26)
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
  df.070CALM2[, 3] <- rep(149, 26)
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
  df.072CALM1[, 3] <- rep(149, 26)
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
  df.073CALML6[, 3] <- rep(226, 26)
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
  df.075CALML4[, 3] <- rep(153, 26)
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
  df.076PPP3CA[, 3] <- rep(515, 26)
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
  df.077PPP3CB[, 3] <- rep(515, 26)
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
  df.078PPP3CC[, 3] <- rep(581, 26)
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
  df.079PPP3R1[, 3] <- rep(177, 26)
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
  df.081NFATC1[, 3] <- rep(922, 26)
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
  df.082NFATC2[, 3] <- rep(960, 26)
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
  df.083NFATC3[, 3] <- rep(981, 26)
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
  df.087CAMK2[, 3] <- rep(557, 26)
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
  df.088PRKAA1[, 3] <- rep(558, 26)
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
  df.089PRKAA2[, 3] <- rep(553, 26)
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
  df.090PRKAB1[, 3] <- rep(274, 26)
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
  df.092PRKAG1[, 3] <- rep(304, 26)
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
  df.094PRKAG2[, 3] <- rep(551, 26)
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
  df.095CAMK1D[, 3] <- rep(392, 26)
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
  df.096CAMK1G[, 3] <- rep(468, 26)
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
  df.097CAMK1[, 3] <- rep(386, 26)
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
  df.098CAMK2A[, 3] <- rep(478, 26)
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
  df.099CAMK2D[, 3] <- rep(494, 26)
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
  df.101CAMK2G[, 3] <- rep(546, 26)
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
  df.102CAMK4[, 3] <- rep(428, 26)
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
  df.104GUCY1A2[, 3] <- rep(702, 26)
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
  df.106GUCY1B1[, 3] <- rep(502, 26)
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
  df.108NPR1[, 3] <- rep(1054, 26)
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
  df.112MYLK3[, 3] <- rep(895, 26)
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
  df.116MYL9[, 3] <- rep(172, 26)
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
  df.119RHOA[, 3] <- rep(194, 26)
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
  df.120ROCK1[, 3] <- rep(1364, 26)
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
  df.121ROCK2[, 3] <- rep(1403, 26)
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
  df.123PPP1CB[, 3] <- rep(327, 26)
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
  df.124PPP1CC[, 3] <- rep(325, 26)
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
  df.125PPP1R12A[, 3] <- rep(1117, 26)
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
  df.129ADCY1[, 3] <- rep(1112, 26)
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
  df.130ADCY2[, 3] <- rep(1093, 26)
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
  df.133ADCY5[, 3] <- rep(1189, 26)
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
  df.135ADCY7[, 3] <- rep(1062, 26)
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
  df.136ADCY8[, 3] <- rep(1217, 26)
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
  df.137ADCY9[, 3] <- rep(1552, 26)
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
  df.138PRKACA[, 3] <- rep(403, 26)
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
  df.139PRKACB[, 3] <- rep(403, 26)
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
  df.141GNAI1[, 3] <- rep(354, 26)
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
  df.146PIK3R5[, 3] <- rep(904, 26)
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
  df.148SRC[, 3] <- rep(540, 26)
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
  df.149KCNJ3[, 3] <- rep(493, 26)
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
  df.150KCNJ6[, 3] <- rep(407, 26)
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
  df.152KCNJ5[, 3] <- rep(444, 26)
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
  df.153EGFR[, 3] <- rep(1216, 26)
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

# uncomment to export numeric results
# write_xlsx(res.c.ex, paste0(path, "output/BLASTp/scaled_max_scores.xlsx"))   


# get boolean results
L.hom <- list()

for (i in 1:ncol(res.c)) {
  L.hom[i] <- data.frame(res.c[[i]] > res.c[["rowM"]])
  names(L.hom)[i] <- paste(colnames(res.c[i]))
}

gene.ages <- as.data.frame(L.hom) %>% 
  mutate(pat, .before = "rowM") %>% 
  rename(species.name = c..Ciona.savignyi....Ciona.intestinalis....Branchiostoma.lanceolatum...) 

# uncomment to export boolean results
# write_xlsx(gene.ages, paste0(path, "output/BLASTp/homology_boolean.xlsx"))



### --------------------------

# load results and add columns for gene age categories
OTpthwy.ress <- read_excel(paste0(path, "data/processed/gene_ages.xlsx"), sheet = 1)

# add column containing information on whether a gene is an ortholog or not
for (i in 1:10) {
  OTpthwy.ress$Ortholog[OTpthwy.ress$Phylostratum == i] <- "n"
}
for (i in 11:20) {
  OTpthwy.ress$Ortholog[OTpthwy.ress$Phylostratum == i] <- "y"
}


# add column containing information on whether a gene is an homolog or not
for (i in 1:10) {
  OTpthwy.ress$Homolog[OTpthwy.ress$Phylostratum == i] <- "y"
}
for (i in 11:20) {
  OTpthwy.ress$Homolog[OTpthwy.ress$Phylostratum == i] <- "n"
}


# add column containing information on whether a gene is ancient (PS1-3), modern non-vertebrate (PS4-10), or modern vertebrate (PS11-20)
OTpthwy.ress$age.3.levels.vert <- NA

for (i in 1:3) {
  OTpthwy.ress$age.3.levels.vert[OTpthwy.ress$Phylostratum == i] <- "ancient"
}
for (i in 4:10) {
  OTpthwy.ress$age.3.levels.vert[OTpthwy.ress$Phylostratum == i] <- "mod.non.vertebrate"
}
for (i in 11:20) {
  OTpthwy.ress$age.3.levels.vert[OTpthwy.ress$Phylostratum == i] <- "mod.vertebrate"
}

OTpthwy.ress %>% 
  count(age.3.levels.vert)

# export results
write_xlsx(OTpthwy.ress, paste0(path, "output/BLASTp/OTpthwy_res02_vertebrate_threshold.xlsx"))



######################################
########### END OF SCRIPT ############
######################################



