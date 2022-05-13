

####################################################################
## DATA LOADING, general preparations
####################################################################


## generate a references list with all genes in the genome, containing their phylostratum, ensemble gene ID, and gene symbol
# load phylo map  for phylostratum info
HomoSapiens.PhyloMap <- read_excel(paste0(BASE, "RProject/OT-bd-analyses_loc/MBE_2008_Homo_Sapiens_PhyloMap.xls"), 
                                   sheet = 1, skip = 1)
HomoSapiens.PhyloMap <- HomoSapiens.PhyloMap[ , c("Phylostratum","Gene_ID")]          # 22845 genes


# get gene symbols for all genes in the genome but only those that are available in the phylo map 
genensid <- ensembldb::select(EnsDb.Hsapiens.v79, 
                              keys= HomoSapiens.PhyloMap$Gene_ID, 
                              keytype = "GENEID", 
                              columns = c("SYMBOL","GENEID"))
genensid <- genensid %>% distinct(SYMBOL, .keep_all = TRUE)                          

# join the two lists to annotate list of gene symbols with phylostrata
genensid_ps <- dplyr::inner_join(HomoSapiens.PhyloMap, genensid, by = c("Gene_ID" = "GENEID")) # 18844 genes
## done



# get young OT gene ids and their gene symbols for later subsetting
yPS <- c(4, 5, 6, 7, 10, 12)
M1Cmm <- M1C_tai$mm
yOTids <- toupper(M1Cmm[M1Cmm$Phylostratum %in% yPS, "GeneID"])   # get ensemble gene IDs for young OT pathway genes
yOT <- genensid_ps[genensid_ps$Gene_ID %in% yOTids, ]  # annotate with gene symbols



#####################################################################
## ANALYSIS 1: TITLE
#####################################################################


### >>>>>>>>>>>>>>>>>>> Data loading and prep

## load DK data, subset, calculate means, merge, annotate (I'm sure this can be looped somehow but I'll deal with that later...)

# load DK region annotation
dk_info <- read_csv(paste0(BASE, "data/processed/abagen_DK/atlas-desikankilliany.csv"))



# load Donor 1 DK 
dkD1 <- read_csv(paste0(BASE, "/data/processed/abagen_DK/AHBAdkD1.csv"))
dkD1_y <- dkD1[, names(dkD1) %in% c("label", yOT$SYMBOL)]                # extract subset of modern OT genes, 19 modern genes, 83 regions, lateralized 

ml_dkD1 <- data.frame(rowMeans(dkD1_y[,-1])) 

dkD1_lm <- data.frame(cbind(dkD1_y$label, ml_dkD1$rowMeans.dkD1_y....1..))
colnames(dkD1_lm) <- c("label", "meanD1")


# load Donor 2 DK 
dkD2 <- read_csv(paste0(BASE, "/data/processed/abagen_DK/AHBAdkD2.csv"))
dkD2_y <- dkD2[, names(dkD2) %in% c("label", yOT$SYMBOL)]                # extract subset of modern OT genes, 19 modern genes, 83 regions, lateralized 

ml_dkD2 <- data.frame(rowMeans(dkD2_y[,-1])) 

dkD2_lm <- data.frame(cbind(dkD2_y$label, ml_dkD2$rowMeans.dkD2_y....1..))
colnames(dkD2_lm) <- c("label", "meanD2")


# load Donor 3 DK 
dkD3 <- read_csv(paste0(BASE, "/data/processed/abagen_DK/AHBAdkD3.csv"))
dkD3_y <- dkD3[, names(dkD3) %in% c("label", yOT$SYMBOL)]                # extract subset of modern OT genes, 19 modern genes, 83 regions, lateralized 

ml_dkD3 <- data.frame(rowMeans(dkD3_y[,-1])) 

dkD3_lm <- data.frame(cbind(dkD3_y$label, ml_dkD3$rowMeans.dkD3_y....1..))
colnames(dkD3_lm) <- c("label", "meanD3")


# load Donor 4 DK 
dkD4 <- read_csv(paste0(BASE, "/data/processed/abagen_DK/AHBAdkD4.csv"))
dkD4_y <- dkD4[, names(dkD4) %in% c("label", yOT$SYMBOL)]                # extract subset of modern OT genes, 19 modern genes, 83 regions, lateralized 

ml_dkD4 <- data.frame(rowMeans(dkD4_y[,-1])) 

dkD4_lm <- data.frame(cbind(dkD4_y$label, ml_dkD4$rowMeans.dkD4_y....1..))
colnames(dkD4_lm) <- c("label", "meanD4")


# load Donor 5 DK 
dkD5 <- read_csv(paste0(BASE, "/data/processed/abagen_DK/AHBAdkD5.csv"))
dkD5_y <- dkD5[, names(dkD5) %in% c("label", yOT$SYMBOL)]                # extract subset of modern OT genes, 19 modern genes, 83 regions, lateralized 

ml_dkD5 <- data.frame(rowMeans(dkD5_y[,-1])) 

dkD5_lm <- data.frame(cbind(dkD5_y$label, ml_dkD5$rowMeans.dkD5_y....1..))
colnames(dkD5_lm) <- c("label", "meanD5")


# load Donor 6 DK 
dkD6 <- read_csv(paste0(BASE, "/data/processed/abagen_DK/AHBAdkD6.csv"))
dkD6_y <- dkD6[, names(dkD6) %in% c("label", yOT$SYMBOL)]                # extract subset of modern OT genes, 19 modern genes, 83 regions, lateralized 

ml_dkD6 <- data.frame(rowMeans(dkD6_y[,-1])) 

dkD6_lm <- data.frame(cbind(dkD6_y$label, ml_dkD6$rowMeans.dkD6_y....1..))
colnames(dkD6_lm) <- c("label", "meanD6")



DK_all_d <- dkD1_lm %>%                   # join all donor averages into one df
  left_join(dkD2_lm, by='label') %>%
  left_join(dkD3_lm, by='label') %>%
  left_join(dkD4_lm, by='label') %>%
  left_join(dkD5_lm, by='label') %>%
  left_join(dkD6_lm, by='label') 

DK_ad_i <- dplyr::right_join(dk_info, DK_all_d, by=c("id"="label"))       # annotate with brain region labels
DK_ad_i <- data.frame(DK_ad_i)




### >>>>>>>>>>>>>>>>>>> analysis


## analyses only for the left hemisphere !!


# prepare destrieux data, remove right hemisphere
pat_ana1 <- c("L|B")
DK_ad_il <- subset(DK_ad_i, grepl(pat_ana1, DK_ad_i$hemisphere))  


# convert to long format
DK_ad_ilL <- DK_ad_il %>% pivot_longer(c("meanD1", "meanD2", "meanD3", "meanD4", "meanD5", "meanD6"),        # transform into long format, the way t.test likes it :)
                                         names_to = "donors", values_to = "expression")
DK_ad_ilL <- data.frame(DK_ad_ilL)



mean_dklL <- mean(DK_ad_ilL$expression)             
regions_dklL <- unique(DK_ad_ilL$label)
pvals_dklL <- numeric()
tStatistic_dklL <- numeric()
cohensD_dklL <- numeric()
estimate_dklL <- numeric()
sd_dklL <- numeric()

for (i in 1:length(regions_dklL)){
  Temp_dklL <- t.test(DK_ad_ilL[DK_ad_ilL$label==regions_dklL[i],"expression"], 
                      mu = mean_dklL, alternative = "two.sided")
  CohD_dklL <- ES.t.one(m=Temp_dklL$estimate,
                        sd=sd(DK_ad_ilL[DK_ad_ilL$label==regions_dklL[i],"expression"]),
                        mu=mean_dklL,
                        alternative = "two.sided")
  pvals_dklL <- c(pvals_dklL,Temp_dklL$p.value)
  tStatistic_dklL <- c(tStatistic_dklL,Temp_dklL$statistic)
  cohensD_dklL <- c(cohensD_dklL,CohD_dklL$d)
  estimate_dklL <- c(estimate_dklL, Temp_dklL$estimate)
  sd_dklL <- c(sd_dklL, sd(DK_ad_ilL[DK_ad_ilL$label==regions_dklL[i],"expression"]))
}

SD_mean_dklL <- sd(DK_ad_ilL$expression)

results_dklL <- data.frame(Region=regions_dklL,P_unadjusted=pvals_dklL,
                           t_stat=tStatistic_dklL,CohD=cohensD_dklL)
results_dklL$P_fdr <- p.adjust(results_dklL$P_unadjusted, method="fdr", n = length(results_dklL$P_unadjusted))

results_dklL$stat_sign <- ""
results_dklL$stat_sign[results_dklL$P_fdr <= 0.05]  <- "*"

results_dklL  

#library(writexl)
#write_xlsx(results_dklL, paste0(BASE, "supp_mat02.xlsx"))

# --------------------------------------------------------------------> 3 significant regions


## CONCLUSION: Compared to the average expression of the young OT gene set across all regions, the young OT gene set is significantly 
#              lower expressed in 3 subcortical regions.



# PLOTTING

# plot expression values in dot graph
meanS_dklL <- data.frame(estimate_dklL)
meanS_dklL$mean_allstruc <- Temp_dklL$null.value
meanS_dklL <- cbind(regions_dklL, meanS_dklL)
meanS_dklL$sd <- sd_dklL
meanS_dklL$structure <- "cortical"
meanS_dklL$structure[35:42] <- "subcortical_brainstem"



## ---------------- plotting t-values on cortical map

# plot t values in cortical brain map
library(DataCombine)

dk_dat <- dk$data
DKreg_regions <- data.frame(dk_dat$region)
DKreg_regions$labs <- dk_dat$label
DKreg_regions <- DKreg_regions %>% distinct(dk_dat.region, .keep_all = TRUE)      # remove duplicates
DKreg_regions <- data.frame(DKreg_regions)


# data preparation
DKresdat_l <- results_dklL
DKresdat_l$Region <- Map(paste0, 'lh_', DKresdat_l$Region)           # add "lh_" in front of each region name
DKresdat_l <- DKresdat_l[-c(35:42),]                                 # remove sub-cortical regions
DKresdat_l <- InsertRow(DKresdat_l, c(NA, NA, NA, NA, NA, NA), RowNum = 1)        # add rows to match the data frame structure of "DKreg_regions"
DKresdat_l <- InsertRow(DKresdat_l, c("lh_corpuscallosum", NA, NA, NA, NA, NA), RowNum = 24)    # add rows to match the data frame structure of "DKreg_regions"
DKresdat_l <- data.frame(DKresdat_l[match(DKreg_regions$labs, DKresdat_l$Region),])    # reorder rows according to the order of regions in the atlas provided by ggseg
DKresdat_l <- cbind(DKreg_regions$dk_dat.region, DKresdat_l)          # finally add region column

names(DKresdat_l)[names(DKresdat_l) == "DKreg_regions$dk_dat.region"] <- "region"

DKresdat_l$t_stat <- as.numeric(DKresdat_l$t_stat)

# plotting

# vertical

p127 <-  ggplot() + theme_void() +
  geom_brain(data = DKresdat_l, 
             atlas = dk, 
             colour = "white",
             position = position_brain(side + hemi ~ .),
             aes(fill = t_stat)) +
  viridis::scale_fill_viridis(option="mako", discrete=F) +
  labs(fill='t statistic') 

p127

p127 <- p127 + ylab("\nlateral   |   medial\n") + labs(title="Cortical expression") +
  theme(plot.title = element_text(hjust=.5, size=25, face="bold"),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_text(size=23, angle=90, face="bold"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.title = element_text(size=23),
        legend.text = element_text(size=16),
        legend.key.size = unit(1, 'cm'))

p127



#ggsave(filename = paste0(BASE, "RProject/OT-bd-analyses_loc/output/Ana1_cortical_vert_v2.pdf"), p712,
#       width = 5, height = 10, units = "in", device='pdf')




##  ---------------- plot expression values for subcortical regions in flat violin plot, vertical and horizontal violins

# resdat_left_sc <- resdat_left[-c(1:74),]   ??????

# prep data

# assign young genes dataframe to new objects to not overwrite the original objects, add col with donor id

ydkD1 <- dkD1_y
ydkD1$donor <- "D1"

ydkD2 <- dkD2_y
ydkD2$donor <- "D2"

ydkD3 <- dkD3_y
ydkD3$donor <- "D3"

ydkD4 <- dkD4_y
ydkD4$donor <- "D4"

ydkD5 <- dkD5_y
ydkD5$donor <- "D5"

ydkD6 <- dkD6_y
ydkD6$donor <- "D6"



# average expression for each gene across the six different datasets/donors and collect in a new data frame
# fetch gene colnames
cnames <- colnames(dkD1_y[,c(2:20)])           # any df will work, because they have the same gene colnames
e_list <- list()
for (i in cnames) {
  
  g_means <- (ydkD1[[i]] + ydkD2[[i]] + ydkD3[[i]] + ydkD4[[i]] + ydkD5[[i]] + ydkD6[[i]])/6
  e_list[[i]] <- g_means
  g_means_l <- do.call(rbind, e_list)
}

g_means_lt <- data.frame(t(g_means_l))             # transpose df

g_means_lt <- cbind(dk_info, g_means_lt)
g_means_lt_sc <- subset(g_means_lt, grepl("subcortex/brainstem", g_means_lt$structure))  
g_means_lt_sc <- subset(g_means_lt_sc, grepl(pat_ana1, g_means_lt_sc$hemisphere))


g_means_lt_sc$mean1 <-  meanS_dklL[meanS_dklL$structure=="subcortical_brainstem","estimate_dklL"]        # add means from the t-test
g_means_lt_sc$mean2 <-  rowMeans(g_means_lt_sc[5:23])                                                    # calculate means manually
 
write_xlsx(g_means_lt_sc, paste0(BASE, "supp_mat03.xlsx"))

# --> identical :) Phew.


# transform to long format
g_means_lt_sc_L <- g_means_lt_sc %>% pivot_longer(c("CACNB1", "CACNB2", "CACNB3", "CACNB4", "CACNG1", "CACNG2", "CACNG3", "CACNG6", "CACNG7", 
                                                    "CD38", "CDKN1A", "ELK1", "FOS", "NFATC3", "NFATC4", "NPPA", "OXT", "OXTR", "PIK3R5"),        # transform into long format, the way t.test likes it :)
                                        names_to = "genes", values_to = "expression")
g_means_lt_sc_L <- data.frame(g_means_lt_sc_L)
g_means_lt_sc_L$label <- as.factor(g_means_lt_sc_L$label)



# ------------------ half violin/raincloud plots


# ..... vertical version

# plot specific data prep
# order by mean expression, ASCENDING
g_means_lt_sc_L_asc <- g_means_lt_sc_L %>% 
  arrange(mean2)

# manually set levels to make sure they are plotted in ascending order on the x-axis
order_h <- unique(g_means_lt_sc_L_asc$label)
g_means_lt_sc_L_asc$label <- factor(g_means_lt_sc_L_asc$label, levels = order_h)


# separate values for OXT, OXTR, and CD38
g_means_lt_sc_cd38 <- subset(g_means_lt_sc_L_asc, grepl("^CD38$", g_means_lt_sc_L_asc$genes)) 
g_means_lt_sc_oxt <- subset(g_means_lt_sc_L_asc, grepl("^OXT$", g_means_lt_sc_L_asc$genes)) 
g_means_lt_sc_oxtr <- subset(g_means_lt_sc_L_asc, grepl("^OXTR$", g_means_lt_sc_L_asc$genes)) 
g_means_lt_sc_ot <- rbind(g_means_lt_sc_cd38, g_means_lt_sc_oxt, g_means_lt_sc_oxtr)         # combine
g_means_lt_sc_ot <- g_means_lt_sc_ot %>%                                           # it should still be in the right order, but just to be on the very safe side
  arrange((mean2))

# ... and remove them from base jitter (they are kept for the violin plot to ensure appropriate curving of the plots)
pat_ot <- c("OXT|OXTR|CD38")
g_means_lt_sc_L_short <- subset(g_means_lt_sc_L_asc, !grepl(pat_ot, g_means_lt_sc_L_asc$genes))


# custom y axis labels
ytext <- c("Thalamus proper", "Brain stem",  "Pallidum", "Accumbens area", "Putamen", "Caudate", "Amygdala", "Hippocampus")


# plotting
g222 <- 
  ggplot() +          
  geom_density_ridges(g_means_lt_sc_L_asc, mapping = aes(x = expression, y = factor(label),              # use full dataset with all genes for violin plots
                                                    fill = mean2, colour = mean2), 
                      alpha=.15, scale= .85) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 23, face="bold"),
        axis.title.y = element_text(size = 23, face="bold"),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),
        legend.position="none") +
  ylab("Subcortical regions\n") + xlab("\nmRNA intensity") +
  scale_colour_gradientn(colours = c("#9893DD", "#779ED7", "#5DC3C3", "#3CBC75", "#B8DE29")) +
  scale_fill_gradientn(colours = c("#9893DD", "#779ED7", "#5DC3C3", "#3CBC75", "#B8DE29")) 

g222 <- g222 + geom_point(g_means_lt_sc_L_short, mapping = aes(x = expression, y = factor(label),        # use dataset w/o OT genes for base jitter
                                                        fill = mean2, colour = mean2), 
                          position = position_jitter(height = .15), size = 4, alpha = .65)

g222 <- g222 + geom_point(g_means_lt_sc_ot,    mapping = aes(x = expression, y = factor(label)),         # use dataset with only OT genes for highlight
                          position = position_jitter(height = .15), size = 4, alpha=.6) 

g222 <- g222 + geom_vline(xintercept = mean(g_means_lt_sc_L_asc$mean2), 
                          linetype="dashed", size=.75, colour="gray60")                                  # add mean line of mean expression points

g222 <- g222 + geom_point(g_means_lt_sc_L_asc,      mapping = aes(x = mean2, y = factor(label)),         # add mean expression points
                          position = position_nudge(y = .3), size = 7) 

g222 <- g222 + geom_vline(xintercept = mean(meanS_dklL$mean_allstruc), linetype="dashed", size=.75)      # add mean line for all regions, cortical and subcortical

g222 <- g222 + 
  annotate('text', x = 1.25, y = 1, label = '*', colour= "black", size=11) +
  annotate('text', x = 1.25, y = 2, label = '*', colour= "black", size=11) +
  annotate('text', x = 1.25, y = 3, label = '*', colour= "black", size=11) 

g222 <- g222 + 
  scale_y_discrete(labels= ytext)

g222

#ggsave(filename = paste0(BASE, "RProject/OT-bd-analyses_loc/output/Ana1_subcortical_vert_v3.pdf"), g222,
#       width = 10, height = 15, units = "in", device='pdf')


# arrange in a grid
g222_p127 <- plot_grid(g222, p127,
                       ncol = 2, nrow = 1,
                       rel_widths = c(1.75, 1),
                       labels = c('A', 'B'),
                       label_size = 23) 

g222_p127

#ggsave(filename = paste0(BASE, "RProject/OT-bd-analyses_loc/output/Ana1DK_ScC_combined_a.pdf"), g222_p127,
#       width = 15, height = 10.5, units = "in", device='pdf')




g_means_lt_sc_phylo <- g_means_lt_sc_L_asc

test <- g_means_lt_sc_phylo %>%
  mutate(phylostratum = case_when(
    endsWith(genes, "CACNB2") ~ "PS4",
    endsWith(genes, "CACNB3") ~ "PS4",
    endsWith(genes, "CACNB4") ~ "PS4",
    endsWith(genes, "CACNB1") ~ "PS4",
    endsWith(genes, "ELK1") ~ "PS5",
    endsWith(genes, "OXTR") ~ "PS6",
    endsWith(genes, "CD38") ~ "PS6",
    endsWith(genes, "NFATC4") ~ "PS6",
    endsWith(genes, "NFATC2") ~ "PS6",
    endsWith(genes, "NFATC3") ~ "PS6",
    endsWith(genes, "FOS") ~ "PS6",
    endsWith(genes, "CDKN1A") ~ "PS6",
    endsWith(genes, "OXT") ~ "PS7",
    endsWith(genes, "CACNG7") ~ "PS7",
    endsWith(genes, "CACNG2") ~ "PS7",
    endsWith(genes, "CACNG3") ~ "PS7",
    endsWith(genes, "PIK3R5") ~ "PS10",
    endsWith(genes, "NPPA") ~ "PS12",
    endsWith(genes, "CACNG1") ~ "PS12",
    endsWith(genes, "CACNG6") ~ "PS12"))

patstructures <- c("amygdala|accumbensarea|pallidum|thalamusproper|brainstem")
test <- subset(test, grepl(patstructures, test$label))

###############################################################################
## ANALYSIS 2: one-sample t-test for each region with different population mean
###############################################################################


## compare the expression of young OT genset to the average expression of all protein coding genes for each region separately, no relation to the average expression across the brain.

# .... population mean is the average expression of all genes in region x, sample mean is the average expression of young OT genes subset in region x



### >>>>>>>>>>>>>>>>>>> DESTRIEUX cortical

## subset destrieux data, calculate means, merge, annotate (I'm sure this can be looped somehow but I'll deal with that later...)


## alternative mean calculation procedure (nice to have for double checking):

#des9861_t <- data.frame(t(des9861[,-1]))
#coln9861ad <- colnames(des9861_t)
#ml_9861ad <- numeric()
#for (i in coln9861ad) {
#  means9861ad <- mean(des9861_t[[i]])
#  ml_9861ad <- c(ml_9861ad, means9861ad)
#}

#des9861_ALLs_alt <- data.frame(ml_9861ad) 

# OR

# test <- data.frame(rowSums(des9861_2[,2:14957])/14956)



# 9861 - average expression across all genes for each cortical region 
dkD1_2 <- dkD1
dkD1_ALLs <- data.frame(rowMeans(dkD1_2[,-1])) 
dkD1_ALLs$label <- dkD1_2$label


# 10021 - average expression across all genes for each cortical region
dkD2_2 <- dkD2
dkD2_ALLs <- data.frame(rowMeans(dkD2_2[,-1])) 
dkD2_ALLs$label <- dkD2_2$label


# 12876 - average expression across all genes for each cortical region
dkD3_2 <- dkD3
dkD3_ALLs <- data.frame(rowMeans(dkD3_2[,-1])) 
dkD3_ALLs$label <- dkD3_2$label


# 14380 - average expression across all genes for each cortical region
dkD4_2 <- dkD4
dkD4_ALLs <- data.frame(rowMeans(dkD4_2[,-1])) 
dkD4_ALLs$label <- dkD4_2$label


# 15496 - average expression across all genes for each cortical region
dkD5_2 <- dkD5
dkD5_ALLs <- data.frame(rowMeans(dkD5_2[,-1])) 
dkD5_ALLs$label <- dkD5_2$label


# 15697 - average expression across all genes for each cortical region
dkD6_2 <- dkD6
dkD6_ALLs <- data.frame(rowMeans(dkD6_2[,-1])) 
dkD6_ALLs$label <- dkD6_2$label



DK_all_d_ALL <- dkD1_ALLs %>%                                     # join all donor averages into one df
  left_join(dkD2_ALLs, by='label') %>%
  left_join(dkD3_ALLs, by='label') %>%
  left_join(dkD4_ALLs, by='label') %>%
  left_join(dkD5_ALLs, by='label') %>%
  left_join(dkD6_ALLs, by='label') 


DK_all_d_ALL <- cbind(DK_all_d_ALL$label, DK_all_d_ALL)
DK_all_d_ALL <- DK_all_d_ALL[, -3]

colnames(DK_all_d_ALL) <- c("label", "meanALLD1", "meanALLD2", "meanALLD3", "meanALLD4", "meanALLD5", "meanALLD6")


# average across donors. This is just for info and not needed in the analysis. Originally, it was implemented but is now deprecated.
#des_all_donors_ALL_t <- data.frame(t(des_all_donors_ALL[,-1]))
#colndada <- colnames(des_all_donors_ALL_t)
#ml_dada <- numeric()
#for (i in colndada) {
#  meansdada <- mean(des_all_donors_ALL_t[[i]])
#  ml_dada <- c(ml_dada, meansdada)
#}
#######


DK_all_d_ALLi <- dplyr::right_join(dk_info, DK_all_d_ALL, by=c("id"="label"))          # annotate with region names

DK_all_d_ALLiL <- subset(DK_all_d_ALLi, grepl(pat_ana1, DK_all_d_ALLi$hemisphere))     # select only left hemisphere and brain stem
DK_all_d_ALLiL <- data.frame(DK_all_d_ALLiL)

# transform into long format, the way t.test likes it :)
DK_all_d_ALLiL_l <- DK_all_d_ALLiL %>% pivot_longer(c("meanALLD1", "meanALLD2", "meanALLD3", "meanALLD4", "meanALLD5", "meanALLD6"),        
                                                      names_to = "donors", values_to = "expressionALL")
DK_all_d_ALLiL_l <- data.frame(DK_all_d_ALLiL_l)


# OT dataframe: DK_ad_ilL


# analysis

regions_dk2 <- unique(DK_ad_ilL$label)
pvals_dk2 <- numeric()
tStatistic_dk2 <- numeric()
cohensD_dk2 <- numeric()

estimateOT_dk2 <- numeric()
estimate_dk2 <- numeric()
sdOT_dk2 <- numeric()
sd_dk2 <- numeric()


for (i in 1:length(regions_dk2)){
  Temp_dk2 <- t.test(DK_ad_ilL[DK_ad_ilL$label==regions_dk2[i],"expression"], 
                      mu = mean(DK_all_d_ALLiL_l[DK_all_d_ALLiL_l$label==regions_dk2[i],"expressionALL"]), 
                      alternative = "two.sided")
  CohD_dk2 <- ES.t.one(m=Temp_dk2$estimate,
                        sd=sd(DK_ad_ilL[DK_ad_ilL$label==regions_dk2[i],"expression"]),
                        mu=mean(DK_all_d_ALLiL_l[DK_all_d_ALLiL_l$label==regions_dk2[i],"expressionALL"]),
                        alternative = "two.sided")
  pvals_dk2 <- c(pvals_dk2,Temp_dk2$p.value)
  tStatistic_dk2 <- c(tStatistic_dk2,Temp_dk2$statistic)
  cohensD_dk2 <- c(cohensD_dk2,CohD_dk2$d)
  
  estimateOT_dk2 <- c(estimateOT_dk2, Temp_dk2$estimate)
  estimate_dk2 <- c(estimate_dk2, Temp_dk2$null.value)
  sdOT_dk2 <- c(sdOT_dk2, sd(DK_ad_ilL[DK_ad_ilL$label==regions_dk2[i],"expression"]))
  sd_dk2 <- c(sd_dk2, sd(DK_all_d_ALLiL_l[DK_all_d_ALLiL_l$label==regions_dk2[i],"expressionALL"]))
}

results_dk2 <- data.frame(Region=regions_dk2,P_unadjusted=pvals_dk2,
                           t_stat=tStatistic_dk2,CohD=cohensD_dk2)

results_dk2$P_fdr <- p.adjust(results_dk2$P_unadjusted, method="fdr", n = length(results_dk2$P_unadjusted))

results_dk2$stat_sign <- ""
results_dk2$stat_sign[results_dk2$P_fdr <= 0.05]  <- "*"

results_dk2

# combine test statistics and parameters and regions

test_stats_dk2 <- data.frame(OTmean=estimateOT_dk2, mean=estimate_dk2, OTSD=sdOT_dk2, SD=sd_dk2)
test_stats_dk2 <- cbind(results_dk2$Region, test_stats_dk2)

# --------------------------------------------------------------------> significant results for some regions


# plotting

# ----- plot expression values in lolli pop graph

names(test_stats_dk2)[names(test_stats_dk2) == "results_dk2$Region"] <- "region"
test_stats_dk2$lobe <- "Frontal"
test_stats_dk2$lobe[c(1, 5, 6, 8, 14, 15, 29, 32, 33)] <- "Temporal"
test_stats_dk2$lobe[c(7, 21, 22, 24, 28, 30)] <- "Parietal"
test_stats_dk2$lobe[c(4, 10, 20)] <- "Occipital"
test_stats_dk2$lobe[34] <- "Insula"
test_stats_dk2$lobe[c(2, 9, 22, 25)] <- "Limbic"
test_stats_dk2$lobe[c(35:42)] <- "XSub-cortical"
test_stats_dk2$LongNames <- c("Banks of superior temporal S", "Caudal anterior cingulate C", "Caudal middle frontal G", "Cuneus C", 
                              "Entorhinal C",  "Fusiform G", "Inferior parietal C", "Inferior temporal G", "Isthmus cingulate C", 
                              "Lateral occipital C", "Lateral orbitofrontal C",  "Lingual G", "Medial orbitofrontal C", 
                              "Middle temporal G", "Parahippocampal G", "Paracentral L", "Pars opercularis", "Pars orbitalis", 
                              "Pars triangularis", "Pericalcarine C", "Postcentral G", "Posterior cingulate C", "Precentral G", 
                              "Precuneus C", "Rostral anterior cingulate C", "Rostral middle frontal G", "Superior frontal G", 
                              "Superior parietal C", "Superior temporal G", "Supramarginal G", "Frontal pole", "Temporal pole",
                              "Transverse temporal C", "Insula", "Thalamus proper", "Caudate", "Putamen", "Pallidum", "Accumbens area", 
                              "Hippocampus", "Amygdala", "Brain stem")


hlinOT <- mean(test_stats_dk2$OTmean) 
hlinAll <- mean(test_stats_dk2$mean)



test_stats_dk2_o <- test_stats_dk2 %>%           # order descending by mean
  arrange(desc(mean))

test_stats_dk2_o <- test_stats_dk2_o %>%       # order ascending by lobe label
  arrange((lobe))


# manually set levels to make sure they are plotted in descending order on the x-axis
order_dk2 <- unique(test_stats_dk2_o$LongNames)
test_stats_dk2_o$LongNames <- factor(test_stats_dk2_o$LongNames, levels = order_dk2)


# Plot
lpp <- ggplot(test_stats_dk2_o) + theme_classic() +
  geom_hline(yintercept=hlinAll, linetype="dashed", color="grey65", size=1) +
  geom_hline(yintercept=hlinOT, linetype="dashed", color="black", size=1) +
  
  geom_segment( aes(x=LongNames, xend=LongNames, y=OTmean, yend=mean), color="grey") +
  geom_point( aes(x=LongNames, y=OTmean, color = lobe), size=9, alpha=.85) +
  geom_point( aes(x=LongNames, y=mean, color = lobe), size=9, alpha=.25) +
  viridis::scale_colour_viridis(option="turbo", discrete=T, 
                                labels = c("Frontal", "Insula", "Limbic", "Occipital",
                                           "Parietal", "Temporal", "Sub-cortical"))

lpp <- lpp + ylab("mRNA intensity\n") + xlab("\nBrain regions") +
  theme(axis.title = element_text(size=25, face="bold"),
        axis.text.x = element_text(angle=50, size=16, vjust =1, hjust = 1),
        # axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=15),
        legend.position = c(0.25, 0.2),
        legend.direction = "horizontal",
        legend.text = element_text(size=22),
        legend.title = element_blank()) 

lpp <- lpp + 
  annotate('text', x = 17, y = 0.5417715, label = '*', colour= "black", size=14) +
  annotate('text', x = 26, y = 0.5500894, label = '*', colour= "black", size=14) +
  annotate('text', x = 28, y = 0.5595212, label = '*', colour= "black", size=14) +
  annotate('text', x = 35, y = 0.4395165, label = '*', colour= "black", size=14) +
  annotate('text', x = 36, y = 0.4275628, label = '*', colour= "black", size=14) 

lpp

#ggsave(filename = paste0(BASE, "RProject/OT-bd-analyses_loc/output/Ana2DK_cSc3.pdf"), lpp,
#       width = 20, height = 10, units = "in", device='pdf')












##################################################################
# combine plots from analysis 1 and plot from analysis 2 in a grid
##################################################################

#g22 <- g2 + theme(plot.margin = unit(c(0, 2, 0, 0), "cm"))
#g22

#g2_p098x <- plot_grid(g22, p098,
#                       ncol = 2, nrow = 1,
#                       rel_widths = c(1.4, 1),
#                       labels = c('A', 'B'),
#                       label_size = 25) 

#g2_p098x

#pa1pa2 <- plot_grid(g2_p098x, pa2, nrow =2, 
#                    rel_heights = c(1.45, 1),
#                    labels = c('', 'C'),
#                    label_size = 25) 

#pa1pa2

#ggsave(filename = paste0(BASE, "RProject/OT-bd-analyses_loc/output/pa1pa2_2.pdf"), pa1pa2,
#       width = 25, height = 25, units = "in", device='pdf')





