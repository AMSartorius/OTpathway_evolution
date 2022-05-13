

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
## ANALYSIS 1: Replicate analysis in Dans Nat Com paper, 2019, Fig. 1
#####################################################################



# the first two following subsections are for cortical and subcortical expression separately. Refer to subsection three for combined analysis. 
# I think this was important to compare to a population mean across ALL regions - cortical AND subcortical



### >>>>>>>>>>>>>>>>>>> 1 DESTRIEUX cortical

## load destrieux data, subset, calculate means, merge, annotate (I'm sure this can be looped somehow but I'll deal with that later...)

# load destrieux region annotation
des_info <- read_csv("C:/Users/alina/nilearn_data/destrieux_2009/destrieux2009_rois_labels_lateralized.csv")
des_info_s <- des_info[-1,]



# load Donor 1 destrieux 
desD1 <- read_csv(paste0(BASE, "/data/processed/abagen_DesPauli_rm/AHBAdesD0.csv"))
desD1_y <- desD1[, names(desD1) %in% c("label", yOT$SYMBOL)]                # extract subset of modern OT genes, 19 modern genes, 148 cortical regions, lateralized 

ml_D1d <- data.frame(rowMeans(desD1_y[,-1])) 

desD1_lm <- data.frame(cbind(desD1_y$label, ml_D1d$rowMeans.desD1_y....1..))
colnames(desD1_lm) <- c("label", "meanD1")


# load Donor 2 destrieux 
desD2 <- read_csv(paste0(BASE, "/data/processed/abagen_DesPauli_rm/AHBAdesD1.csv"))
desD2_y <- desD2[, names(desD2) %in% c("label", yOT$SYMBOL)]                # extract subset of modern OT genes, 19 modern genes, 148 cortical regions, lateralized

ml_D2d <- data.frame(rowMeans(desD2_y[,-1]))

desD2_lm <- data.frame(cbind(desD2_y$label, ml_D2d$rowMeans.desD2_y....1..))
colnames(desD2_lm) <- c("label", "meanD2")


# load Donor 3 destrieux 
desD3 <- read_csv(paste0(BASE, "/data/processed/abagen_DesPauli_rm/AHBAdesD2.csv"))
desD3_y <- desD3[, names(desD3) %in% c("label", yOT$SYMBOL)]                # extract subset of modern OT genes, 19 modern genes, 148 cortical regions, lateralized

ml_D3d <- data.frame(rowMeans(desD3_y[,-1]))

desD3_lm <- data.frame(cbind(desD3_y$label, ml_D3d$rowMeans.desD3_y....1..))
colnames(desD3_lm) <- c("label", "meanD3")


# load Donor 4 destrieux 
desD4 <- read_csv(paste0(BASE, "/data/processed/abagen_DesPauli_rm/AHBAdesD3.csv"))
desD4_y <- desD4[, names(desD4) %in% c("label", yOT$SYMBOL)]                # extract subset of modern OT genes, 19 modern genes, 148 cortical regions, lateralized

ml_D4d <- data.frame(rowMeans(desD4_y[,-1]))

desD4_lm <- data.frame(cbind(desD4_y$label, ml_D4d$rowMeans.desD4_y....1..))
colnames(desD4_lm) <- c("label", "meanD4")


# load Donor 5 destrieux 
desD5 <- read_csv(paste0(BASE, "/data/processed/abagen_DesPauli_rm/AHBAdesD4.csv"))
desD5_y <- desD5[, names(desD5) %in% c("label", yOT$SYMBOL)]                # extract subset of modern OT genes, 19 modern genes, 148 cortical regions, lateralized

ml_D5d <- data.frame(rowMeans(desD5_y[,-1]))

desD5_lm <- data.frame(cbind(desD5_y$label, ml_D5d$rowMeans.desD5_y....1..))
colnames(desD5_lm) <- c("label", "meanD5")


# load Donor 6 destrieux 
desD6 <- read_csv(paste0(BASE, "/data/processed/abagen_DesPauli_rm/AHBAdesD5.csv"))
desD6_y <- desD6[, names(desD6) %in% c("label", yOT$SYMBOL)]                # extract subset of modern OT genes, 19 modern genes, 148 cortical regions, lateralized

ml_D6d <- data.frame(rowMeans(desD6_y[,-1]))

desD6_lm <- data.frame(cbind(desD6_y$label, ml_D6d$rowMeans.desD6_y....1..))
colnames(desD6_lm) <- c("label", "meanD6")



## INFO: "Medial wall" is missing in both hemispheres



des_all_d <- desD1_lm %>%                   # join all donor averages into one df
  left_join(desD2_lm, by='label') %>%
  left_join(desD3_lm, by='label') %>%
  left_join(desD4_lm, by='label') %>%
  left_join(desD5_lm, by='label') %>%
  left_join(desD6_lm, by='label') 

des_ad_i <- dplyr::right_join(des_info_s, des_all_d, by=c("index"="label"))       # annotate with brain region labels
des_ad_i <- data.frame(des_ad_i)

des_ad_i$hemisphere <- "L"                                         # add hemisphere column with "L" and "R"
des_ad_i$hemisphere[75:148] <- "R"

des_ad_i$name <- gsub('L ', '', des_ad_i$name)                  # remove "L " in front of each region
des_ad_i$name <- gsub('R ', '', des_ad_i$name)                  # remove "R " in front of each region





### >>>>>>>>>>>>>>>>>>> 2 PAULI subcortical

## load pauli data, subset, calculate means, merge, annotate (I'm sure this can be looped somehow but I'll deal with that later...)

# load pauli region annotation
pau_info <- read.table("C:/Users/alina/nilearn_data/pauli_2017/labels.txt")
pau_info$V1 <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16)
pau_info$node_info <-  c("Putamen", "Caudate nucleus", "Nucleus accumbens", "Extended amygdala", "Globus pallidus, external", "Globus pallidus, internal", "Substantia nigra, pars compacta", 
                         "Red nucleus", "Substantia nigra, pars reticulata", "Parabrachial pigmented nucleus", "Ventral tegmental area", "Ventral pallidum", "Habenular nuclei", "Hypothalamus", 
                         "Mammillary nucleus", "Subthalamic nucleus")


# load Donor 1 pauli 
pauD1 <- read_csv(paste0(BASE, "/data/processed/abagen_DesPauli_rm/AHBApauD0.csv"))
pauD1_y <- pauD1[, names(pauD1) %in% c("label", yOT$SYMBOL)]                         # extract subset of modern OT genes, 18 modern genes, 16 subcortical regions !! ONE GENES LESS !!

ml_D1p <- data.frame(rowMeans(pauD1_y[,-1])) 

pauD1_lm <- data.frame(cbind(pauD1_y$label, ml_D1p$rowMeans.pauD1_y....1..))
colnames(pauD1_lm) <- c("label", "meanD1")


# load Donor 2 pauli 
pauD2 <- read_csv(paste0(BASE, "/data/processed/abagen_DesPauli_rm/AHBApauD1.csv"))
pauD2_y <- pauD2[, names(pauD2) %in% c("label", yOT$SYMBOL)]                # extract subset of modern OT genes, 19 modern genes, 148 cortical regions, lateralized

ml_D2p <- data.frame(rowMeans(pauD2_y[,-1])) 

pauD2_lm <- data.frame(cbind(pauD2_y$label, ml_D2p$rowMeans.pauD2_y....1..))
colnames(pauD2_lm) <- c("label", "meanD2")


# load Donor 3 pauli 
pauD3 <- read_csv(paste0(BASE, "/data/processed/abagen_DesPauli_rm/AHBApauD2.csv"))
pauD3_y <- pauD3[, names(pauD3) %in% c("label", yOT$SYMBOL)]                # extract subset of modern OT genes, 19 modern genes, 148 cortical regions, lateralized

ml_D3p <- data.frame(rowMeans(pauD3_y[,-1])) 

pauD3_lm <- data.frame(cbind(pauD3_y$label, ml_D3p$rowMeans.pauD3_y....1..))
colnames(pauD3_lm) <- c("label", "meanD3")


# load Donor 4 pauli 
pauD4 <- read_csv(paste0(BASE, "/data/processed/abagen_DesPauli_rm/AHBApauD3.csv"))
pauD4_y <- pauD4[, names(pauD4) %in% c("label", yOT$SYMBOL)]                # extract subset of modern OT genes, 19 modern genes, 148 cortical regions, lateralized

ml_D4p <- data.frame(rowMeans(pauD4_y[,-1])) 

pauD4_lm <- data.frame(cbind(pauD4_y$label, ml_D4p$rowMeans.pauD4_y....1..))
colnames(pauD4_lm) <- c("label", "meanD4")


# load Donor 5 pauli 
pauD5 <- read_csv(paste0(BASE, "/data/processed/abagen_DesPauli_rm/AHBApauD4.csv"))
pauD5_y <- pauD5[, names(pauD5) %in% c("label", yOT$SYMBOL)]                # extract subset of modern OT genes, 19 modern genes, 148 cortical regions, lateralized

ml_D5p <- data.frame(rowMeans(pauD5_y[,-1])) 

pauD5_lm <- data.frame(cbind(pauD5_y$label, ml_D5p$rowMeans.pauD5_y....1..))
colnames(pauD5_lm) <- c("label", "meanD5")


# load Donor 6 pauli 
pauD6 <- read_csv(paste0(BASE, "/data/processed/abagen_DesPauli_rm/AHBApauD5.csv"))
pauD6_y <- pauD6[, names(pauD6) %in% c("label", yOT$SYMBOL)]                # extract subset of modern OT genes, 19 modern genes, 148 cortical regions, lateralized

ml_D6p <- data.frame(rowMeans(pauD6_y[,-1])) 

pauD6_lm <- data.frame(cbind(pauD6_y$label, ml_D6p$rowMeans.pauD6_y....1..))
colnames(pauD6_lm) <- c("label", "meanD6")



pau_all_d <- pauD1_lm %>%                   # join all donor averages into one df
  left_join(pauD2_lm, by='label') %>%
  left_join(pauD3_lm, by='label') %>%
  left_join(pauD4_lm, by='label') %>%
  left_join(pauD5_lm, by='label') %>%
  left_join(pauD6_lm, by='label') 

pau_ad_i <- dplyr::right_join(pau_info, pau_all_d, by=c("V1"="label"))       # annotate with brain region labels

pau_ad_i$hemisphere <- "B"                                                        # add hemisphere column with "B"




### >>>>>>>>>>>>>>>>>>> 3 Combined cortical and subcortical analysis


## analyses only for the left hemisphere !!


# prepare destrieux data, remove right hemisphere
des_ad_ileft <- subset(des_ad_i, grepl("L", des_ad_i$hemisphere))  


# prepare pauli data
pau_ad_i_noV2 <- pau_ad_i[, -2]                      # remove second column with abbreviations of regions
names(pau_ad_i_noV2)[names(pau_ad_i_noV2) == "V1"] <- "index" 
names(pau_ad_i_noV2)[names(pau_ad_i_noV2) == "node_info"] <- "name" 



des_pau_l <- rbind(des_ad_ileft, pau_ad_i_noV2)   # combine into one df

# convert to long format
des_pau_lL <- des_pau_l %>% pivot_longer(c("meanD1", "meanD2", "meanD3", "meanD4", "meanD5", "meanD6"),        # transform into long format, the way t.test likes it :)
                                            names_to = "donors", values_to = "expression")
des_pau_lL <- data.frame(des_pau_lL)



mean_dplL <- mean(des_pau_lL$expression)             
regions_dplL <- unique(des_pau_lL$name)
pvals_dplL <- numeric()
tStatistic_dplL <- numeric()
cohensD_dplL <- numeric()
estimate_dplL <- numeric()
sd_dplL <- numeric()

for (i in 1:length(regions_dplL)){
  Temp_dplL <- t.test(des_pau_lL[des_pau_lL$name==regions_dplL[i],"expression"], 
                     mu = mean_dplL, alternative = "two.sided")
  CohD_dplL <- ES.t.one(m=Temp_dplL$estimate,
                       sd=sd(des_pau_lL[des_pau_lL$name==regions_dplL[i],"expression"]),
                       mu=mean_dplL,
                       alternative = "two.sided")
  pvals_dplL <- c(pvals_dplL,Temp_dplL$p.value)
  tStatistic_dplL <- c(tStatistic_dplL,Temp_dplL$statistic)
  cohensD_dplL <- c(cohensD_dplL,CohD_dplL$d)
  estimate_dplL <- c(estimate_dplL, Temp_dplL$estimate)
  sd_dplL <- c(sd_dplL, sd(des_pau_lL[des_pau_lL$name==regions_dplL[i],"expression"]))
}

SD_mean_dplL <- sd(des_pau_lL$expression)

results_dplL <- data.frame(Region=regions_dplL,P_unadjusted=pvals_dplL,
                          t_stat=tStatistic_dplL,CohD=cohensD_dplL)
results_dplL$P_fdr <- p.adjust(results_dplL$P_unadjusted, method="fdr", n = length(results_dplL$P_unadjusted))

results_dplL$stat_sign <- ""
results_dplL$stat_sign[results_dplL$P_fdr <= 0.05]  <- "*"

results_dplL  


# --------------------------------------------------------------------> 5 significant regions


## CONCLUSION: Compared to the average expression of the young OT gene set across all regions, the young OT gene set is not significantly higher or
#              lower expressed in any cortical or subcortical region.



# PLOTTING

# plot expression values in dot graph
meanS_dplL <- data.frame(estimate_dplL)
meanS_dplL$mean_allstruc <- Temp_dplL$null.value
meanS_dplL <- cbind(regions_dplL, meanS_dplL)
meanS_dplL$sd <- sd_dplL
meanS_dplL$structure <- "cortical"
meanS_dplL$structure[75:90] <- "subcortical"




## ---------------- plotting t-values on cortical map

# plot t values in cortical brain map

###### fetch destrieux atlas map

ggseg_atlas_repos("Desterieux")
des_repo <- ggseg_atlas_repos("desterieux", ignore.case = TRUE)
if (!requireNamespace("Desterieux", quietly = TRUE)) {
  install_ggseg_atlas(repo = des_repo$repo, source = des_repo$source)
}


library(ggsegDesterieux)

# inspect atlas
ggseg(atlas = ggsegDesterieux::desterieux, 
      mapping=aes(fill = region)) +
  scale_fill_brain("desterieux", package="ggsegDesterieux")

#####





library(ggsegDesterieux)
des_dat <- desterieux$data
reg_regions <- unique(des_dat$region)


# data preparation
resdat_l <- results_dplL
resdat_l$Region <- gsub('_', ' ', resdat_l$Region)    # remove "_" in region names
resdat_l_c <- resdat_l[-c(75:90),]                    # remove subcortical regions because ggseg destrieux only supports cortical regions

resdat_l_ro <- data.frame(resdat_l_c[match(reg_regions, resdat_l_c$Region),])    # reorder rows according to the order of regions in the atlas provided by ggseg

names(resdat_l_ro)[names(resdat_l_ro) == "Region"] <- "region"

# plotting

# vertical

p712 <- ggplot() +
  geom_brain(data = resdat_l_ro,
             atlas = desterieux,
             mapping = aes(fill = t_stat),
             colour = "black",
             position = position_brain(side + hemi ~ .))
p712 <- p712 + scale_fill_viridis_c(option = "mako", direction = -1) +
  theme_void() +
  labs(fill='t statistic') 

p712 <- p712 + ylab("\nlateral   |   medial\n") + labs(title="Cortical expression") +
  theme(plot.title = element_text(hjust=.5, size=25, face="bold"),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_text(size=23, angle=90, face="bold"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.title = element_text(size=23),
        legend.text = element_text(size=16),
        legend.key.size = unit(1, 'cm'))
p712

#ggsave(filename = paste0(BASE, "RProject/OT-bd-analyses_loc/output/Ana1_cortical_vert_v2.pdf"), p712,
#       width = 5, height = 10, units = "in", device='pdf')




##  ---------------- plot expression values for subcortical regions in flat violin plot, vertical and horizontal violins

# resdat_left_sc <- resdat_left[-c(1:74),]   ??????

# prep data

# assign young genes dataframe to new objects to not overwrite the original objects, add col with donor id

pD1 <- pauD1_y
pD1$donor <- "D1"

pD2 <- pauD2_y
pD2$donor <- "D2"

pD3 <- pauD3_y
pD3$donor <- "D3"

pD4 <- pauD4_y 
pD4$donor <- "D4"

pD5 <- pauD5_y
pD5$donor <- "D5"

pD6 <- pauD6_y
pD6$donor <- "D6"



# average expression for each gene across the six different datasets/donors and collect in a new data frame
# fetch gene colnames
cnames <- colnames(pauD1_y[,c(2:20)])           # any df will work, because they have the same gene colnames
e_list <- list()
for (i in cnames) {
  
  g_means <- (pD1[[i]] + pD2[[i]] + pD3[[i]] + pD4[[i]] + pD5[[i]] + pD6[[i]])/6
  e_list[[i]] <- g_means
  g_means_l <- do.call(rbind, e_list)
}

g_means_lt <- data.frame(t(g_means_l))             # transpose df

index <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16)      # add index column with region label indices
g_means_lt <- cbind(index, g_means_lt)

pauli_aR <- g_means_lt                                       # assign into a new object to again not overwrite old objects

pauli_aR <- cbind(pau_info$node_info, pauli_aR)                             # add long region names
names(pauli_aR)[names(pauli_aR) == "pau_info$node_info"] <- "region"

pauli_aR$mean1 <-  meanS_dplL[meanS_dplL$structure=="subcortical","estimate_dplL"]        # add means from the t-test
pauli_aR$mean2 <-  rowMeans(pauli_aR[3:21])                                               # calculate means manually

# --> identical :) Phew.


# transform to long format
pauli_aR_L <- pauli_aR %>% pivot_longer(c("CACNB1", "CACNB2", "CACNB3", "CACNB4", "CACNG1", "CACNG2", "CACNG3", "CACNG6", "CACNG7", "CD38", "CDKN1A", "ELK1",
                                              "FOS", "NFATC3", "NFATC4", "NPPA", "OXT", "OXTR", "PIK3R5"),        # transform into long format, the way t.test likes it :)
                                            names_to = "genes", values_to = "expression")
pauli_aR_L <- data.frame(pauli_aR_L)
pauli_aR_L$region <- as.factor(pauli_aR_L$region)



# ------------------ half violin/raincloud plots


# ..... vertical version

# plot specific data prep
# order by mean expression, ASCENDING
pauli_aR_L_asc <- pauli_aR_L %>% 
  arrange(mean2)

# manually set levels to make sure they are plotted in ascending order on the x-axis
order_h <- unique(pauli_aR_L_asc$region)
pauli_aR_L_asc$region <- factor(pauli_aR_L_asc$region, levels = order_h)


# separate values for OXT, OXTR, and CD38
pauli_aR_L_cd38 <- subset(pauli_aR_L_asc, grepl("^CD38$", pauli_aR_L_asc$genes)) 
pauli_aR_L_oxt <- subset(pauli_aR_L_asc, grepl("^OXT$", pauli_aR_L_asc$genes)) 
pauli_aR_L_oxtr <- subset(pauli_aR_L_asc, grepl("^OXTR$", pauli_aR_L_asc$genes)) 
pauli_aR_L_ot <- rbind(pauli_aR_L_cd38, pauli_aR_L_oxt, pauli_aR_L_oxtr)         # combine
pauli_aR_L_ot <- pauli_aR_L_ot %>%                                           # it should still be in the right order, but just to be on the very safe side
  arrange((mean2))

# ... and remove them from base jitter (they are kept for the violin plot to ensure appropriate curving of the plots)
pat_ot <- c("OXT|OXTR|CD38")
pauli_aR_L_short <- subset(pauli_aR_L_asc, !grepl(pat_ot, pauli_aR_L_asc$genes))



# order OT subset ascending
pauli_aR_L_asc_ot <- pauli_aR_L_ot %>%                           
  arrange(mean2)

# order subset w/o OT ascending
pauli_aR_L_asc_short <- pauli_aR_L_short %>%                           
  arrange(mean2)



# plotting
g22 <- 
  ggplot() +          
  geom_density_ridges(pauli_aR_L_asc, mapping = aes(x = expression, y = factor(region),                           # use full dataset with all genes for violin plots
                                                       fill = mean2, colour = mean2), 
                      alpha=.1, scale= .9) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 23, face="bold"),
        axis.title.y = element_text(size = 23, face="bold"),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),
        legend.position="none") +
  ylab("Subcortical regions\n") + xlab("\nmRNA intensity") +
  viridis::scale_colour_viridis(option="viridis", discrete = F, begin = .35, end = .9) +
  viridis::scale_fill_viridis(option="viridis", discrete = F, begin = .35, end = .9) 

g22 <- g22 + geom_point(pauli_aR_L_short, mapping = aes(x = expression, y = region,                               # use dataset w/o OT genes for base jitter
                                                            fill = mean2, colour = mean2), 
                      position = position_jitter(height = .15), size = 4, alpha = .6)

g22 <- g22 + geom_point(pauli_aR_L_ot,    mapping = aes(x = expression, y = factor(region)),                      # use dataset with only OT genes for highlight
                      position = position_jitter(height = .15), size = 4, alpha=.6) 

g22 <- g22 + geom_vline(xintercept = mean(pauli_aR_L_asc$mean2), linetype="dashed", size=.75, colour="gray60")       # add mean line of mean expression points

g22 <- g22 + geom_point(pauli_aR_L_asc,      mapping = aes(x = mean2, y = factor(region)),                           # add mean expression points
                      position = position_nudge(y = .3), size = 7) 

g22 <- g22 + geom_vline(xintercept = mean(meanS_dplL$mean_allstruc), linetype="dashed", size=.75)           # add mean line for all regions, cortical and subcortical

g22 <- g22 + 
  annotate('text', x = 1.35, y = 2, label = '*', colour= "black", size=11) +
  annotate('text', x = 1.35, y = 5, label = '*', colour= "black", size=11) +
  annotate('text', x = 1.35, y = 6, label = '*', colour= "black", size=11) +
  annotate('text', x = 1.35, y = 10, label = '*', colour= "black", size=11) +
  annotate('text', x = 1.35, y = 12, label = '*', colour= "black", size=11)

g22

#ggsave(filename = paste0(BASE, "RProject/OT-bd-analyses_loc/output/Ana1_subcortical_vert_v2.pdf"), g22,
#       width = 10, height = 15, units = "in", device='pdf')


# arrange in a grid
g22_p712 <- plot_grid(g22, p712,
                       ncol = 2, nrow = 1,
                       rel_widths = c(1.75, 1),
                       labels = c('A', 'B'),
                       label_size = 23) 

g22_p712

ggsave(filename = paste0(BASE, "RProject/OT-bd-analyses_loc/output/Ana1_ScC_combined2.pdf"), g22_p712,
       width = 15, height = 11, units = "in", device='pdf')







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
desD1_2 <- desD1
desD1_ALLs <- data.frame(rowMeans(desD1_2[,-1])) 
desD1_ALLs$label <- desD1_2$label


# 10021 - average expression across all genes for each cortical region
desD2_2 <- desD2
desD2_ALLs <- data.frame(rowMeans(desD2_2[,-1])) 
desD2_ALLs$label <- desD2_2$label


# 12876 - average expression across all genes for each cortical region
desD3_2 <- desD3
desD3_ALLs <- data.frame(rowMeans(desD3_2[,-1])) 
desD3_ALLs$label <- desD3_2$label


# 14380 - average expression across all genes for each cortical region
desD4_2 <- desD4
desD4_ALLs <- data.frame(rowMeans(desD4_2[,-1])) 
desD4_ALLs$label <- desD4_2$label


# 15496 - average expression across all genes for each cortical region
desD5_2 <- desD5
desD5_ALLs <- data.frame(rowMeans(desD5_2[,-1])) 
desD5_ALLs$label <- desD5_2$label


# 15697 - average expression across all genes for each cortical region
desD6_2 <- desD6
desD6_ALLs <- data.frame(rowMeans(desD6_2[,-1])) 
desD6_ALLs$label <- desD6_2$label



des_all_d_ALL <- desD1_ALLs %>%                                     # join all donor averages into one df
  left_join(desD2_ALLs, by='label') %>%
  left_join(desD3_ALLs, by='label') %>%
  left_join(desD4_ALLs, by='label') %>%
  left_join(desD5_ALLs, by='label') %>%
  left_join(desD6_ALLs, by='label') 


des_all_d_ALL <- cbind(des_all_d_ALL$label, des_all_d_ALL)
des_all_d_ALL <- des_all_d_ALL[, -3]

colnames(des_all_d_ALL) <- c("label", "meanALLD1", "meanALLD2", "meanALLD3", "meanALLD4", "meanALLD5", "meanALLD6")


# average across donors. This is just for info and not needed in the analysis. Originally, it was implemented but is now deprecated.
#des_all_donors_ALL_t <- data.frame(t(des_all_donors_ALL[,-1]))
#colndada <- colnames(des_all_donors_ALL_t)
#ml_dada <- numeric()
#for (i in colndada) {
#  meansdada <- mean(des_all_donors_ALL_t[[i]])
#  ml_dada <- c(ml_dada, meansdada)
#}
#######


des_all_d_ALLi <- dplyr::right_join(des_info_s, des_all_d_ALL, by=c("index"="label"))   # annotate with region names

des_all_d_ALLiL <- subset(des_all_d_ALLi, grepl("^L ", des_all_d_ALLi$name))       # select only left hemisphere
des_all_d_ALLiL$name <- gsub('L ', '', des_all_d_ALLiL$name)                       # remove left hemisphere indicator in the region names column
des_all_d_ALLiL <- data.frame(des_all_d_ALLiL)

# transform into long format, the way t.test likes it :)
des_all_d_ALLiL_l <- des_all_d_ALLiL %>% pivot_longer(c("meanALLD1", "meanALLD2", "meanALLD3", "meanALLD4", "meanALLD5", "meanALLD6"),        
                                                                names_to = "donors", values_to = "expressionALL")
des_all_d_ALLiL_l <- data.frame(des_all_d_ALLiL_l)



### OT gene set subset into long format
des_ad_ileft_L <- des_ad_ileft %>% pivot_longer(c("meanD1", "meanD2", "meanD3", "meanD4", "meanD5", "meanD6"),        # transform into long format, the way t.test likes it :)
                                                      names_to = "donors", values_to = "expression")
des_ad_ileft_L <- data.frame(des_ad_ileft_L)
## ----------


# analysis

regions_des2 <- unique(des_ad_ileft_L$name)
pvals_des2 <- numeric()
tStatistic_des2 <- numeric()
cohensD_des2 <- numeric()

estimateOT_des2 <- numeric()
estimate_des2 <- numeric()
sdOT_des2 <- numeric()
sd_des2 <- numeric()


for (i in 1:length(regions_des2)){
  Temp_des2 <- t.test(des_ad_ileft_L[des_ad_ileft_L$name==regions_des2[i],"expression"], 
                       mu = mean(des_all_d_ALLiL_l[des_all_d_ALLiL_l$name==regions_des2[i],"expressionALL"]), 
                       alternative = "two.sided")
  CohD_des2 <- ES.t.one(m=Temp_des2$estimate,
                         sd=sd(des_ad_ileft_L[des_ad_ileft_L$name==regions_des2[i],"expression"]),
                         mu=mean(des_all_d_ALLiL_l[des_all_d_ALLiL_l$name==regions_des2[i],"expressionALL"]),
                         alternative = "two.sided")
  pvals_des2 <- c(pvals_des2,Temp_des2$p.value)
  tStatistic_des2 <- c(tStatistic_des2,Temp_des2$statistic)
  cohensD_des2 <- c(cohensD_des2,CohD_des2$d)
  
  estimateOT_des2 <- c(estimateOT_des2, Temp_des2$estimate)
  estimate_des2 <- c(estimate_des2, Temp_des2$null.value)
  sdOT_des2 <- c(sdOT_des2, sd(des_ad_ileft_L[des_ad_ileft_L$name==regions_des2[i],"expression"]))
  sd_des2 <- c(sd_des2, sd(des_all_d_ALLiL_l[des_all_d_ALLiL_l$name==regions_des2[i],"expressionALL"]))
}

results_des2 <- data.frame(Region=regions_des2,P_unadjusted=pvals_des2,
                            t_stat=tStatistic_des2,CohD=cohensD_des2)

test_stats_des2 <- data.frame(OTmean=estimateOT_des2, mean=estimate_des2, OTSD=sdOT_des2, SD=sd_des2)


# combine with subcortical and then FDR it


# --------------------------------------------------------------------> nothing significant




### >>>>>>>>>>>>>>>>>>> PAULI subcortical

## prepare data for average expression of all coding genes for each subcortical region

## alternative mean calculation procedure (nice to have for double checking):
#pau9861_t <- data.frame(t(pau9861[,-1]))     
#coln9861ap <- colnames(pau9861_t)
#ml_9861ap <- numeric()
#for (i in coln9861ap) {
#  means9861ap <- mean(pau9861_t[[i]])
#  ml_9861ap <- c(ml_9861ap, means9861ap)
#}

#pau9861_ALLs_alt <- data.frame(ml_9861ap) 



# Donor 1 - average expression across all genes for each sub cortical region 
pauD1_2 <- pauD1
pauD1_ALLs <- data.frame(rowMeans(pauD1_2[,-1]))
pauD1_ALLs$label <- pauD1_2$label


# Donor 2 - average expression across all genes for each sub cortical region 
pauD2_2 <- pauD2
pauD2_ALLs <- data.frame(rowMeans(pauD2_2[,-1]))
pauD2_ALLs$label <- pauD2_2$label


# Donor 3 - average expression across all genes for each sub cortical region 
pauD3_2 <- pauD3
pauD3_ALLs <- data.frame(rowMeans(pauD3_2[,-1]))
pauD3_ALLs$label <- pauD3_2$label


# Donor 4 - average expression across all genes for each sub cortical region 
pauD4_2 <- pauD4
pauD4_ALLs <- data.frame(rowMeans(pauD4_2[,-1]))
pauD4_ALLs$label <- pauD4_2$label


# Donor 5 - average expression across all genes for each sub cortical region 
pauD5_2 <- pauD5
pauD5_ALLs <- data.frame(rowMeans(pauD5_2[,-1]))
pauD5_ALLs$label <- pauD5_2$label


# Donor 6 - average expression across all genes for each sub cortical region  
pauD6_2 <- pauD6
pauD6_ALLs <- data.frame(rowMeans(pauD6_2[,-1]))
pauD6_ALLs$label <- pauD6_2$label



pau_all_d_ALL <- pauD1_ALLs %>%                                     # join all donor averages into one df
  left_join(pauD2_ALLs, by='label') %>%
  left_join(pauD3_ALLs, by='label') %>%
  left_join(pauD4_ALLs, by='label') %>%
  left_join(pauD5_ALLs, by='label') %>%
  left_join(pauD6_ALLs, by='label') 


pau_all_d_ALL <- cbind(pau_all_d_ALL$label, pau_all_d_ALL)
pau_all_d_ALL <- pau_all_d_ALL[, -3]

colnames(pau_all_d_ALL) <- c("label", "meanALLD1", "meanALLD2", "meanALLD3", "meanALLD4", "meanALLD5", "meanALLD6")


# average across donors. This is just for info and not needed in the analysis. Originally, it was implemented but is now deprecated.
#pau_all_donors_ALL_t <- data.frame(t(pau_all_donors_ALL[,-1]))
#colnpada <- colnames(pau_all_donors_ALL_t)
#ml_pada <- numeric()
#for (i in colnpada) {
#  meanspada <- mean(pau_all_donors_ALL_t[[i]])
#  ml_pada <- c(ml_pada, meanspada)
#}

########


pau_all_d_ALL_info <- dplyr::right_join(pau_info, pau_all_d_ALL, by=c("V1"="label"))     # annotate with region names
pau_all_d_ALLi_L <- pau_all_d_ALL_info %>% pivot_longer(c("meanALLD1", "meanALLD2", "meanALLD3", "meanALLD4", "meanALLD5", "meanALLD6"),        # transform into long format, the way t.test likes it :)
                                                                  names_to = "donors", values_to = "expressionALL")
pau_all_d_ALLi_L <- data.frame(pau_all_d_ALLi_L)


# transform young OT geneset data to long format
pau_ad_i_l <- pau_ad_i_noV2 %>% pivot_longer(c("meanD1", "meanD2", "meanD3", "meanD4", "meanD5", "meanD6"),        # transform into long format, the way t.test likes it :)
                                              names_to = "donors", values_to = "expression")
pau_ad_i_l <- data.frame(pau_ad_i_l)



# analysis

regions_pau2 <- unique(pau_ad_i_l$name)
pvals_pau2 <- numeric()
tStatistic_pau2 <- numeric()
cohensD_pau2 <- numeric()

estimateOT_pau2 <- numeric()
estimate_pau2 <- numeric()
sdOT_pau2 <- numeric()
sd_pau2 <- numeric()



for (i in 1:length(regions_pau2)){
  Temp_pau2 <- t.test(pau_ad_i_l[pau_ad_i_l$name==regions_pau2[i],"expression"], 
                       mu = mean(pau_all_d_ALLi_L[pau_all_d_ALLi_L$node_info==regions_pau2[i],"expressionALL"]), 
                       alternative = "two.sided")
  CohD_pau2 <- ES.t.one(m=Temp_pau2$estimate,
                         sd=sd(pau_ad_i_l[pau_ad_i_l$name==regions_pau2[i],"expression"]),
                         mu=mean(pau_all_d_ALLi_L[pau_all_d_ALLi_L$node_info==regions_pau2[i],"expressionALL"]),
                         alternative = "two.sided")
  pvals_pau2 <- c(pvals_pau2,Temp_pau2$p.value)
  tStatistic_pau2 <- c(tStatistic_pau2,Temp_pau2$statistic)
  cohensD_pau2 <- c(cohensD_pau2,CohD_pau2$d)
  
  estimateOT_pau2 <- c(estimateOT_pau2, Temp_pau2$estimate)
  estimate_pau2 <- c(estimate_pau2, Temp_pau2$null.value)
  sdOT_pau2 <- c(sdOT_pau2, sd(pau_ad_i_l[pau_ad_i_l$name==regions_pau2[i],"expression"]))
  sd_pau2 <- c(sd_pau2, sd(pau_all_d_ALLi_L[pau_all_d_ALLi_L$node_info==regions_pau2[i],"expressionALL"]))
  
}

results_pau2 <- data.frame(Region=regions_pau2,P_unadjusted=pvals_pau2,
                            t_stat=tStatistic_pau2,CohD=cohensD_pau2)

test_stats_pau2 <- data.frame(OTmean=estimateOT_pau2, mean=estimate_pau2, OTSD=sdOT_pau2, SD=sd_pau2)


# combine with cortical and FDR it

## INFO ##
# It doesn't make a difference whether these are separately analyzed and then stacked or if all data is in one df. The t-tests work with multiple population means
# of ROWS (=brain regions), not one COLUMN (=whole brain), each t-test has its own population mean (as opposed to analysis #1). Therefore the population means do not
# depend on the mean of one column and the values in that column. If they would rely on the mean of the same column and would not change with each iteration based on 
# rows, then indeed it would make a difference because the values in the column the mean is based on would change (but again, that is not the case here). Anyways, see 
# end of script for analysis with all the data in one df.


res_dp2 <- rbind(results_des2, results_pau2)

res_dp2$P_fdr <- p.adjust(res_dp2$P_unadjusted, method="fdr", n = length(res_dp2$P_unadjusted))

res_dp2$stat_sign <- ""
res_dp2$stat_sign[res_dp2$P_fdr <= 0.05]  <- "*"

res_dp2

# combine test statistics and parameters

test_stats_dp2 <- rbind(test_stats_des2, test_stats_pau2)
test_stats_dp2 <- cbind(res_dp2$Region, test_stats_dp2)

# --------------------------------------------------------------------> significant results for some regions


# plotting

# ----- plot expression values in lolli pop graph

names(test_stats_dp2)[names(test_stats_dp2) == "res_dp2$Region"] <- "region"

# annotate with full region names
des_ln <- read_excel(paste0(BASE, "RProject/OT-bd-analyses_loc/Data/destrieux_fullnames.xlsx"))    # load long names table
df_dp2 <- dplyr::full_join(des_ln, test_stats_dp2, by=c("ShortName"="region"))
df_dp2$LongName[75:90] <- df_dp2$ShortName[75:90]                                                # add names for subcortical regions
df_dp2$Lobe[75:90] <- "x_subcortical"
df_dp2 <- data.frame(df_dp2[, -c(1,2)])                                                            # remove index and short names
names(df_dp2)[names(df_dp2) == "LongName"] <- "region"
df_dp2$region <- factor(df_dp2$region, levels=unique(df_dp2$region))

df_dp2$structure_OTmean <- "cortical"
df_dp2$structure_OTmean[75:90] <- "subcortical"


hlinOT <- c(mean(df_dp2$OTmean) - mean(df_dp2$OTSD), mean(df_dp2$OTmean), mean(df_dp2$OTmean) + mean(df_dp2$OTSD))
hlinAll <- c(mean(df_dp2$mean) - mean(df_dp2$SD), mean(df_dp2$mean), mean(df_dp2$mean) + mean(df_dp2$SD))



df_dp2_ord <- df_dp2 %>%           # order descending by mean
  arrange(desc(mean))

df_dp2_ord <- df_dp2_ord %>%       # order ascending by lobe label
  arrange((Lobe))


# manually set levels to make sure they are plotted in descending order on the x-axis
order_2 <- unique(df_dp2_ord$region)
df_dp2_ord$region <- factor(df_dp2_ord$region, levels = order_2)

# custom x axis labels
xtext <- c("", "", "", "", "", "", "", "", "", "",  "", "", "", "", "", "", "", "", "", "", 
           "", "", "", "", "", "", "", "", "", "",  "", "", "", "", "", "", "", "", "", "",
           "", "", "", "", "", "", "", "", "", "",  "", "Post-central sulcus", "", "", "", "", "", "", "", "",
           "", "", "", "", "", "", "", "", "", "",  "", "", "", "", "External globus pallidus", "", "", "", "", "",
           "Extended amygdala", "", "", "", "", "", "", "", "", "")


# Plot
lpp <- ggplot(df_dp2_ord) + theme_classic() +
       geom_segment( aes(x=region, xend=region, y=OTmean, yend=mean), color="grey") +
       geom_point( aes(x=region, y=OTmean, color = Lobe), size=6, alpha=.7) +
       geom_point( aes(x=region, y=mean, color = Lobe), size=6, alpha=.3) +
       viridis::scale_colour_viridis(option="turbo", discrete=T, 
                                     labels = c("Frontal", "Fronto-parietal", "Insula", "Lateral sulcus", "Limbic", "Occipital",
                                                "Parietal", "Parieto-limbic", "Parieto-occipital", "Temporal", "Temporo-occipital",
                                                "Sub-cortical"))

lpp <- lpp + ylab("mRNA intensity\n") + xlab("\nBrain regions") +
  theme(axis.title = element_text(size=25, face="bold"),
        axis.text.x = element_text(angle=60, size=14, vjust =1, hjust = 1),
       # axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=14),
        legend.position = c(0.37, 0.2),
        legend.direction = "horizontal",
        legend.text = element_text(size=20),
        legend.title = element_blank()) 
 
lpp <- lpp + 
  annotate('text', x = 52, y = .557, label = '*', colour= "black", size=11) +
  annotate('text', x = 75, y = .427, label = '*', colour= "black", size=11) +
  annotate('text', x = 81, y = .452, label = '*', colour= "black", size=11) 

lpp <- lpp + 
  scale_x_discrete(labels= xtext)

lpp

ggsave(filename = paste0(BASE, "RProject/OT-bd-analyses_loc/output/Ana2_cSc2.pdf"), lpp,
       width = 24, height = 9.5, units = "in", device='pdf')





#pa2 <- ggplot(df_dpa2, aes(x=region)) + theme_classic() +
#  geom_hline(yintercept=hlinesOT, linetype=c("dashed", "solid", "dashed"), color=c("#EA7E98", "#DE2056", "#EA7E98"), size=c(1, 2, 1)) +  
#  #geom_hline(yintercept=hlinesOT, linetype=c("dashed", "solid", "dashed"), color=c("#97E1CE", "#16CC8F", "#97E1CE"), size=c(.5, 1, .5))  # for plasma
#  geom_hline(yintercept=hlinesAll, linetype=c("dashed", "solid", "dashed"), color=c("gray50", "gray35", "gray50"), size=c(.5, 1, .5)) +
#  geom_point(size=6, aes(y=mean))  +
#  geom_errorbar(aes(ymin = mean - SD, ymax = mean + SD, width = .5)) 

#pa2
#pa2 <- pa2 +
#  geom_point(size=6, aes(y=OTmean, colour = region)) +
#  geom_errorbar(aes(ymin = OTmean - OTSD, ymax = OTmean + OTSD, color=factor(region), width = .5)) 
##  viridis::scale_colour_viridis(option="viridis", discrete = T, direction = 1, begin=.3, end=.9)  
## viridis::scale_colour_viridis(option="plasma", discrete = T, direction = 1, begin=.2, end=.875)         # for plasma 

#pa2

#pa2 <- pa2 + ylab("mRNA intensity\n") + xlab("\nBrain regions")
#pa2 <- pa2 + theme(axis.text.x = element_text(angle=50, vjust=1, hjust=1, size=16),
#                   axis.title.x = element_text(size=33, face="bold"),
#                   axis.title.y = element_text(size=30, face="bold"),
#                   axis.text.y = element_text(size=12),
#                   legend.position ="none")  

#pa2 <- pa2 + annotate("text", x = 66, y = .59, label = "*", size=8) +
#  annotate("text", x = 67, y = .58, label = "*", size=8) +
#  annotate("text", x = 78, y = .433, label = "*", size=8) +
#  annotate("text", x = 79, y = .392, label = "*", size=8) 
#pa2 

#ggsave(filename = paste0(BASE, "RProject/OT-bd-analyses_loc/output/test2005.pdf"), pa2,
#       width = 24, height = 9.5, units = "in", device='pdf')










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










############################################################################

## ----------------------- sanity check analysis #2 - all regions in one df

pau_v <- pau_all_donors_ALL_info
pau_v <- pau_v[,-2]
names(pau_v)[names(pau_v) == "V1"] <- "index"
names(pau_v)[names(pau_v) == "node_info"] <- "name"


des_v <- des_all_donors_ALLiL


despau <- rbind(des_v, pau_v)
despau_l <- despau %>% pivot_longer(c("meanALL9861", "meanALL10021", "meanALL12876", "meanALL14380", "meanALL15496", "meanALL15697"),        # transform into long format, the way t.test likes it :)
                                    names_to = "donors", values_to = "expressionALL")
despau_l <- data.frame(despau_l)







regions_x <- unique(des_pau_leftL$name)
pvals_x <- numeric()
tStatistic_x <- numeric()
cohensD_x <- numeric()

for (i in 1:length(regions_x)){
  Temp_x <- t.test(des_pau_leftL[des_pau_leftL$name==regions_x[i],"expression"], 
                   mu = mean(despau_l[despau_l$name==regions_x[i],"expressionALL"]), 
                   alternative = "two.sided")
  CohD_x <- ES.t.one(m=Temp_x$estimate,
                     sd=sd(des_pau_leftL[des_pau_leftL$name==regions_x[i],"expression"]),
                     mu=mean(despau_l[despau_l$name==regions_x[i],"expressionALL"]),
                     alternative = "two.sided")
  pvals_x <- c(pvals_x,Temp_x$p.value)
  tStatistic_x <- c(tStatistic_x,Temp_x$statistic)
  cohensD_x <- c(cohensD_x,CohD_x$d)
}

res_x <- data.frame(Region=regions_x,P_unadjusted=pvals_x,
                    t_stat=tStatistic_x,CohD=cohensD_x)

res_x

res_x$P_fdr <- p.adjust(res_x$P_unadjusted, method="fdr", n = length(res_x$P_unadjusted))

res_x$stat_sign <- ""
res_x$stat_sign[res_x$P_fdr <= 0.05]  <- "*"

res_x





## ----------------------- quality control: check whether means that i calculate when fetching data for the 6 donors separately are the same or correlate with 
#                                           the means that one gets when fetching the data average across all donors


# load destrieux dataset where all donors have been fetched together
AHBAdesALL <- read_csv(paste0(BASE, "data/processed/AHBAdesALL.csv"))
AHBAdesALL_y <- AHBAdesALL[, names(AHBAdesALL) %in% c("label", yOT$SYMBOL)]   # extract young OT genes

# calculate mean across genes
AHBAdesALL_y_t <- data.frame(t(AHBAdesALL_y[,-1]))
colnadall  <- colnames(AHBAdesALL_y_t)
ml_adall<- numeric()
for (i in colnadall) {
  meansadall <- mean(AHBAdesALL_y_t[[i]])
  ml_adall <- c(ml_adall, meansadall)
}

mean1 <- data.frame(ml_adall) 
mean1L <- data.frame(mean1[1:74,]) # keep only left hemisphere


ju <- c(1:74)          # create column with label ids
ju <- data.frame(ju)

mean1L <- cbind(mean1L, ju)    
names(mean1L)[names(mean1L) == "mean1.1.74..."] <- "MeanM"
names(mean1L)[names(mean1L) == "ju"] <- "lab"





# load pauli dataset where all donors have been fetched together
AHBApauALL <- read_csv(paste0(BASE, "data/processed/AHBApauliALL.csv"))
AHBApauALL_y <- AHBApauALL[, names(AHBApauALL) %in% c("label", yOT$SYMBOL)]    # extract young OT genes

AHBApauALL_y_t <- data.frame(t(AHBApauALL_y[,-1]))       # calculate mean across genes
colnadall2  <- colnames(AHBApauALL_y_t)
ml_adall2 <- numeric()
for (i in colnadall2) {
  meansadall2 <- mean(AHBApauALL_y_t[[i]])
  ml_adall2 <- c(ml_adall2, meansadall2)
}

mean2 <- data.frame(ml_adall2) 


ji <- c(1:16)             # generate column with region indices
ji <- data.frame(ji)

mean2 <- cbind(mean2, ji)
names(mean2)[names(mean2) == "ml_adall2"] <- "MeanM"
names(mean2)[names(mean2) == "ji"] <- "lab"





mean12 <- rbind(mean1L, mean2)         # combine cortical (destrieux) and subcortical (pauli)
mean12_a <- cbind(mean12, meanS_dpL$regions_dpL)                            # add long region names
names(mean12_a)[names(mean12_a) == "meanS_dpL$regions_dpL"] <- "region"

mean12_aa <- cbind(mean12_a ,meanS_dpL$estimate_dpL)                              # add means from my analyses
names(mean12_aa)[names(mean12_aa) == "meanS_dpL$estimate_dpL"] <- "mean2"
mean12_aa$structure <- "cortical"
mean12_aa$structure[75:90] <- "subcortical"

gp_pp <- ggplot(mean12_aa, aes(x=region)) + theme_minimal() +                                # plot both and visually compare
  geom_point(size=6, aes(y=MeanM, colour = factor(structure))) 

gp_pp <- gp_pp +
  geom_point(size=6, shape=1, aes(y=mean2, colour = factor(structure))) +
  theme(axis.text.x = element_text(angle=50, vjust=1, hjust=1, size=12))


gp_pp

ggsave(filename = paste0(BASE, "RProject/OT-bd-analyses_loc/output/dotplotAllvsSix.pdf"), gp_pp,
       width = 20, height = 10, units = "in", device='pdf')

# plot as correlation plot

corr_p <- ggplot(mean12_aa, aes(x=MeanM, y=mean2)) +
  geom_point()
corr_p
ggsave(filename = paste0(BASE, "RProject/OT-bd-analyses_loc/output/corrplotAllvsSix.pdf"), corr_p,
       width = 10, height = 10, units = "in", device='pdf')

# calculate linear model
model <- lm(MeanM ~ mean2, mean12_aa)
summary(model)



## ----> they correlate quite strongly (r2 = .96, p < 2.2e-16). But they are not completely identical. Why do they differ, even if its just a bit?

# check expression patterns in desikan atlas - maybe its due to the atlas??

# dk labeling info
DK_info <- read_csv(paste0(BASE, "RProject/OT-bd-analyses_loc/Data/atlas-desikankilliany_space.csv"))

AHBAdk <- read_csv(paste0(BASE, "data/processed/expression_main_mis.csv"))
AHBAdk_y <- AHBAdk[, names(AHBAdk) %in% c("label", yOT$SYMBOL)]   # extract young OT genes

AHBAdk_y_a <- cbind(DK_info$label, DK_info$structure, DK_info$hemisphere, AHBAdk_y)

# calculate mean across genes
# version #1
ml_DK1 <- data.frame(rowSums(AHBAdk_y[2:20])/19)

# version #2
AHBAdk_y_t <- data.frame(t(AHBAdk_y[,-1]))
colnDK  <- colnames(AHBAdk_y_t)
ml_DK <- numeric()
for (i in colnDK) {
  meansDK <- mean(AHBAdk_y_t[[i]])
  ml_DK2 <- c(ml_DK, meansDK)
}

meanDK = ml_DK1 
meanDK_labs <- cbind(meanDK, DK_info$label, DK_info$hemisphere, DK_info$structure)
names(meanDK_labs)[names(meanDK_labs) == "DK_info$label"] <- "region"
names(meanDK_labs)[names(meanDK_labs) == "DK_info$hemisphere"] <- "hemi"
names(meanDK_labs)[names(meanDK_labs) == "DK_info$structure"] <- "structure"
names(meanDK_labs)[names(meanDK_labs) == "rowSums.AHBAdk_y.2.20...19"] <- "exp_mean"

pattern <- c("L|B")
meanDK_labsL <- subset(meanDK_labs, grepl(pattern, meanDK_labs$hemi))



pp <- ggplot(meanDK_labsL, aes(x=region, y=exp_mean)) + theme_minimal() +                                # plot 
  geom_point(size=6, aes(colour = factor(structure))) +
  theme(axis.text.x = element_text(angle=50, vjust=1, hjust=1, size=12))


pp

ggsave(filename = paste0(BASE, "RProject/OT-bd-analyses_loc/output/brainexpression_DK.pdf"), pp,
       width = 20, height = 10, units = "in", device='pdf')


# compare plot next to each other

pp.grid <- pp + theme(legend.position = "none") + xlab("\n\nregion")
gp_pp.grid <- gp_pp +  theme(legend.position = "none")


arr <- pp.grid + gp_pp.grid

ggsave(filename = paste0(BASE, "RProject/OT-bd-analyses_loc/output/brainexpression_DKvsDesPauli.pdf"), pp_gppp,
       width = 20, height = 10, units = "in", device='pdf')











