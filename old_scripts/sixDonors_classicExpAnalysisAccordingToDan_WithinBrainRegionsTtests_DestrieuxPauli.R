

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



# load 9861 destrieux 
des9861 <- read_csv(paste0(BASE, "/data/processed/abagen_destrieux_pauli/AHBAdes9861.csv"))
des9861_y <- des9861[, names(des9861) %in% c("label", yOT$SYMBOL)]                   # extract subset of modern OT genes, 18 modern genes, 148 cortical regions, lateralized !! ONE GENE LESS !!

ml_9861d <- data.frame(rowMeans(des9861_y[,-1])) 

des9861_lm <- data.frame(cbind(des9861_y$label, ml_9861d$rowMeans.des9861_y....1..))
colnames(des9861_lm) <- c("label", "mean9861")


# load 10021 destrieux 
des10021 <- read_csv(paste0(BASE, "/data/processed/abagen_destrieux_pauli/AHBAdes10021.csv"))
des10021_y <- des10021[, names(des10021) %in% c("label", yOT$SYMBOL)]                # extract subset of modern OT genes, 19 modern genes, 148 cortical regions, lateralized

ml_10021d <- data.frame(rowMeans(des10021_y[,-1]))

des10021_lm <- data.frame(cbind(des10021_y$label, ml_10021d$rowMeans.des10021_y....1..))
colnames(des10021_lm) <- c("label", "mean10021")


# load 12876 destrieux 
des12876 <- read_csv(paste0(BASE, "/data/processed/abagen_destrieux_pauli/AHBAdes12876.csv"))
des12876_y <- des12876[, names(des12876) %in% c("label", yOT$SYMBOL)]                # extract subset of modern OT genes, 19 modern genes, 148 cortical regions, lateralized

ml_12876d <- data.frame(rowMeans(des12876_y[,-1]))

des12876_lm <- data.frame(cbind(des12876_y$label, ml_12876d$rowMeans.des12876_y....1..))
colnames(des12876_lm) <- c("label", "mean12876")


# load 14380 destrieux 
des14380 <- read_csv(paste0(BASE, "/data/processed/abagen_destrieux_pauli/AHBAdes14380.csv"))
des14380_y <- des14380[, names(des14380) %in% c("label", yOT$SYMBOL)]                # extract subset of modern OT genes, 19 modern genes, 148 cortical regions, lateralized

ml_14380d <- data.frame(rowMeans(des14380_y[,-1]))

des14380_lm <- data.frame(cbind(des14380_y$label, ml_14380d$rowMeans.des14380_y....1..))
colnames(des14380_lm) <- c("label", "mean14380")


# load 15496 destrieux 
des15496 <- read_csv(paste0(BASE, "/data/processed/abagen_destrieux_pauli/AHBAdes15496.csv"))
des15496_y <- des15496[, names(des15496) %in% c("label", yOT$SYMBOL)]                # extract subset of modern OT genes, 19 modern genes, 148 cortical regions, lateralized

ml_15496d <- data.frame(rowMeans(des15496_y[,-1]))

des15496_lm <- data.frame(cbind(des15496_y$label, ml_15496d$rowMeans.des15496_y....1..))
colnames(des15496_lm) <- c("label", "mean15496")


# load 15697 destrieux 
des15697 <- read_csv(paste0(BASE, "/data/processed/abagen_destrieux_pauli/AHBAdes15697.csv"))
des15697_y <- des15697[, names(des15697) %in% c("label", yOT$SYMBOL)]                # extract subset of modern OT genes, 19 modern genes, 148 cortical regions, lateralized

ml_15697d <- data.frame(rowMeans(des15697_y[,-1]))

des15697_lm <- data.frame(cbind(des15697_y$label, ml_15697d$rowMeans.des15697_y....1..))
colnames(des15697_lm) <- c("label", "mean15697")



## INFO: "Medial wall" is missing in both hemispheres



des_all_donors <- des9861_lm %>%                   # join all donor averages into one df
  left_join(des10021_lm, by='label') %>%
  left_join(des12876_lm, by='label') %>%
  left_join(des14380_lm, by='label') %>%
  left_join(des15496_lm, by='label') %>%
  left_join(des15697_lm, by='label') 

des_ad_info <- dplyr::right_join(des_info_s, des_all_donors, by=c("index"="label"))       # annotate with brain region labels
des_ad_info <- data.frame(des_ad_info)

#des_ad_infoL <- subset(des_ad_info, grepl("^L ", des_ad_info$name))                       # keep only left hemisphere and remove right

des_ad_info$hemisphere <- "L"                                         # add hemisphere column with "L" and "R"
des_ad_info$hemisphere[75:148] <- "R"

des_ad_info$name <- gsub('L ', '', des_ad_info$name)                  # remove "L " in front of each region
des_ad_info$name <- gsub('R ', '', des_ad_info$name)                  # remove "R " in front of each region





### >>>>>>>>>>>>>>>>>>> 2 PAULI subcortical

## load pauli data, subset, calculate means, merge, annotate (I'm sure this can be looped somehow but I'll deal with that later...)

# load pauli region annotation
pau_info <- read.table("C:/Users/alina/nilearn_data/pauli_2017/labels.txt")
pau_info$V1 <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16)
pau_info$node_info <-  c("Putamen", "Caudate nucleus", "Nucleus accumbens", "Extended amygdala", "Globus pallidus, external", "Globus pallidus, internal", "Substantia nigra, pars compacta", 
                             "Red nucleus", "Substantia nigra, pars reticulata", "Parabrachial pigmented nucleus", "Ventral tegmental area", "Ventral pallidum", "Habenular nuclei", "Hypothalamus", 
                             "Mammillary nucleus", "Subthalamic nucleus")


# load 9861 pauli 
pau9861 <- read_csv(paste0(BASE, "/data/processed/abagen_destrieux_pauli/AHBApauli9861.csv"))
pau9861_y <- pau9861[, names(pau9861) %in% c("label", yOT$SYMBOL)] # extract subset of modern OT genes, 18 modern genes, 16 subcortical regions !! ONE GENES LESS !!

ml_9861p <- data.frame(rowMeans(pau9861_y[,-1])) 

pau9861_lm <- data.frame(cbind(pau9861_y$label, ml_9861p$rowMeans.pau9861_y....1..))
colnames(pau9861_lm) <- c("label", "mean9861")


# load 10021 pauli 
pau10021 <- read_csv(paste0(BASE, "/data/processed/abagen_destrieux_pauli/AHBApauli10021.csv"))
pau10021_y <- pau10021[, names(pau10021) %in% c("label", yOT$SYMBOL)]                # extract subset of modern OT genes, 19 modern genes, 148 cortical regions, lateralized

ml_10021p <- data.frame(rowMeans(pau10021_y[,-1])) 

pau10021_lm <- data.frame(cbind(pau10021_y$label, ml_10021p$rowMeans.pau10021_y....1..))
colnames(pau10021_lm) <- c("label", "mean10021")


# load 12876 pauli 
pau12876 <- read_csv(paste0(BASE, "/data/processed/abagen_destrieux_pauli/AHBApauli12876.csv"))
pau12876_y <- pau12876[, names(pau12876) %in% c("label", yOT$SYMBOL)]                # extract subset of modern OT genes, 19 modern genes, 148 cortical regions, lateralized

ml_12876p <- data.frame(rowMeans(pau12876_y[,-1])) 

pau12876_lm <- data.frame(cbind(pau12876_y$label, ml_12876p$rowMeans.pau12876_y....1..))
colnames(pau12876_lm) <- c("label", "mean12876")


# load 14380 pauli 
pau14380 <- read_csv(paste0(BASE, "/data/processed/abagen_destrieux_pauli/AHBApauli14380.csv"))
pau14380_y <- pau14380[, names(pau14380) %in% c("label", yOT$SYMBOL)]                # extract subset of modern OT genes, 19 modern genes, 148 cortical regions, lateralized

ml_14380p <- data.frame(rowMeans(pau14380_y[,-1])) 

pau14380_lm <- data.frame(cbind(pau14380_y$label, ml_14380p$rowMeans.pau14380_y....1..))
colnames(pau14380_lm) <- c("label", "mean14380")


# load 15496 pauli 
pau15496 <- read_csv(paste0(BASE, "/data/processed/abagen_destrieux_pauli/AHBApauli15496.csv"))
pau15496_y <- pau15496[, names(pau15496) %in% c("label", yOT$SYMBOL)]                # extract subset of modern OT genes, 19 modern genes, 148 cortical regions, lateralized

ml_15496p <- data.frame(rowMeans(pau15496_y[,-1])) 

pau15496_lm <- data.frame(cbind(pau15496_y$label, ml_15496p$rowMeans.pau15496_y....1..))
colnames(pau15496_lm) <- c("label", "mean15496")


# load 15697 pauli 
pau15697 <- read_csv(paste0(BASE, "/data/processed/abagen_destrieux_pauli/AHBApauli15697.csv"))
pau15697_y <- pau15697[, names(pau15697) %in% c("label", yOT$SYMBOL)]                # extract subset of modern OT genes, 19 modern genes, 148 cortical regions, lateralized

ml_15697p <- data.frame(rowMeans(pau15697_y[,-1])) 

pau15697_lm <- data.frame(cbind(pau15697_y$label, ml_15697p$rowMeans.pau15697_y....1..))
colnames(pau15697_lm) <- c("label", "mean15697")



pau_all_donors <- pau9861_lm %>%                   # join all donor averages into one df
  left_join(pau10021_lm, by='label') %>%
  left_join(pau12876_lm, by='label') %>%
  left_join(pau14380_lm, by='label') %>%
  left_join(pau15496_lm, by='label') %>%
  left_join(pau15697_lm, by='label') 

pau_ad_info <- dplyr::right_join(pau_info, pau_all_donors, by=c("V1"="label"))       # annotate with brain region labels

pau_ad_info$hemisphere <- "B"                                                        # add hemisphere column with "B"




### >>>>>>>>>>>>>>>>>>> 3 Combined cortical and subcortical analysis


## only left hemisphere !!

# prepare destrieux data, remove right hemisphere
des_ad_infoleft <- subset(des_ad_info, grepl("L", des_ad_info$hemisphere))  


# prepare pauli data
pau_ad_info_noV2 <- pau_ad_info[, -2]
names(pau_ad_info_noV2)[names(pau_ad_info_noV2) == "V1"] <- "index" 
names(pau_ad_info_noV2)[names(pau_ad_info_noV2) == "node_info"] <- "name" 



des_pau_left <- rbind(des_ad_infoleft, pau_ad_info_noV2)   # combine into one df

# convert to long format
des_pau_leftL <- des_pau_left %>% pivot_longer(c("mean9861", "mean10021", "mean12876", "mean14380", "mean15496", "mean15697"),        # transform into long format, the way t.test likes it :)
                                               names_to = "donors", values_to = "expression")
des_pau_leftL <- data.frame(des_pau_leftL)



mean_dpL <- mean(des_pau_leftL$expression)             
regions_dpL <- unique(des_pau_leftL$name)
pvals_dpL <- numeric()
tStatistic_dpL <- numeric()
cohensD_dpL <- numeric()
estimate_dpL <- numeric()
sd_dpL <- numeric()

for (i in 1:length(regions_dpL)){
  Temp_dpL <- t.test(des_pau_leftL[des_pau_leftL$name==regions_dpL[i],"expression"], 
                    mu = mean_dpL, alternative = "two.sided")
  CohD_dpL <- ES.t.one(m=Temp_dpL$estimate,
                      sd=sd(des_pau_leftL[des_pau_leftL$name==regions_dpL[i],"expression"]),
                      mu=mean_dpL,
                      alternative = "two.sided")
  pvals_dpL <- c(pvals_dpL,Temp_dpL$p.value)
  tStatistic_dpL <- c(tStatistic_dpL,Temp_dpL$statistic)
  cohensD_dpL <- c(cohensD_dpL,CohD_dpL$d)
  estimate_dpL <- c(estimate_dpL, Temp_dpL$estimate)
  sd_dpL <- c(sd_dpL, sd(des_pau_leftL[des_pau_leftL$name==regions_dpL[i],"expression"]))
}

SD_mean_dpL <- sd(des_pau_leftL$expression)

results_dpL <- data.frame(Region=regions_dpL,P_unadjusted=pvals_dpL,
                          t_stat=tStatistic_dpL,CohD=cohensD_dpL)
results_dpL$P_fdr <- p.adjust(results_dpL$P_unadjusted, method="fdr", n = length(results_dpL$P_unadjusted))

results_dpL$stat_sign <- ""
results_dpL$stat_sign[results_dpL$P_fdr <= 0.05]  <- "*"

results_dpL  


# --------------------------------------------------------------------> 5 subcortical regions are significant


## RESULTS: Compared to the average expression of the young OT gene set across all regions, the young OT gene set is significantly 
#           lower expressed in the extended amygdala, the external globus pallidus, the substantia nigra, pars reticulata, the 
#           parabrachial pigmented nucleus, and the hypothalamus.



# PLOTTING

# plot expression values in dot graph
meanS_dpL <- data.frame(estimate_dpL)
meanS_dpL$mean_allstruc <- Temp_dpL$null.value
meanS_dpL <- cbind(regions_dpL, meanS_dpL)
meanS_dpL$sd <- sd_dpL
meanS_dpL$structure <- "cortical"
meanS_dpL$structure[75:90] <- "subcortical"


#######################

#hlinesL <- c(unique(meanS_dpL$mean_allstruc - SD_mean_dpL), unique(meanS_dpL$mean_allstruc), unique(meanS_dpL$mean_allstruc + SD_mean_dpL))

#p2 <- ggplot(meanS_dpL, aes(x=regions_dpL, y=estimate_dpL)) + theme_classic() +
#  geom_hline(yintercept=hlinesL,  linetype=c("dashed", "solid", "dashed"), color=c("gray50", "black", "gray50"), size=c(1, 2, 1))  
#p2 <- p2 +
#  geom_point(size=6, aes(colour = regions_dpL)) +
#  geom_errorbar(aes(ymin = estimate_dpL - sd, ymax = estimate_dpL + sd, color=factor(regions_dpL), width = .5)) +
#  viridis::scale_colour_viridis(option="viridis", discrete = T)
#p2 <- p2 + ylab("mRNA intensity\n") + xlab("\nBrain regions")
#p2 <- p2 + theme(axis.text.x = element_text(angle=50, vjust=1, hjust=1, size=12),
#                axis.title.x = element_text(size=25, face="bold"),
#                axis.title.y = element_text(size=22, face="bold"),
#                axis.text.y = element_text(size=12),
#                legend.position = "none")  
#p2 

# important: this way of annotation only work if you keep the width and height of the pdf file EXACTLY the same

#pa <- grid.arrange(p2, right = textGrob("+1SD\n\n\n\nMean\n\n\n\n-1SD", vjust =-.23, gp=gpar(fontsize=15,font=3)))

#ggsave(filename = paste0(BASE, "RProject/OT-bd-analyses_loc/output/test1004.pdf"), pa,
#       width = 24, height = 9.5, units = "in", device='pdf')




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
resdat_left <- results_dpL
resdat_left$Region <- gsub('_', ' ', resdat_left$Region)
resdat_left_c <- resdat_left[-c(75:90),]                    # remove subcortical regions because ggseg destrieux only supports cortical regions

resdat_left_ro <- data.frame(resdat_left_c[match(reg_regions, resdat_left_c$Region),])    # reorder rows according to the order of regions in the atlas

names(resdat_left_ro)[names(resdat_left_ro) == "Region"] <- "region"

# plotting

# vertical

p098 <- ggplot() +
  geom_brain(data = resdat_left_ro,
             atlas = desterieux,
             mapping = aes(fill = t_stat),
             colour = "black",
             position = position_brain(side + hemi ~ .))
p098 <- p098 + scale_fill_viridis_c(option = "mako", direction = -1) +
  theme_void() +
  labs(fill='t statistic') 

p098 <- p098 + ylab("\nlateral   |   medial\n") + labs(title="Cortical expression") +
  theme(plot.title = element_text(hjust=.5, size=25, face="bold"),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_text(size=23, angle=90, face="bold"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.title = element_text(size=23),
        legend.text = element_text(size=16),
        legend.key.size = unit(1, 'cm'))
p098

ggsave(filename = paste0(BASE, "RProject/OT-bd-analyses_loc/output/t_statAna1_cortical_v.pdf"), p098,
       width = 5, height = 10, units = "in", device='pdf')



# horizontal

#p987 <- ggseg(.data = resdat_left_ro,
#              atlas = desterieux,
#              hemisphere = "left", 
#              mapping = aes(fill = t_stat),
#              colour = "black",
#              size = .1, 
#              position = "dispersed")
#p987 <- p987 + scale_fill_viridis_c(option = "mako", direction = -1) +
#  theme_void() +
#  labs(fill='t statistic') 
#p987 <- p987 + labs(title= "Left hemisphere", subtitle="lateral | medial\n") +
#  theme(plot.title = element_text(hjust=.5, size=25, face="bold"),
#        plot.subtitle = element_text(hjust=.5, size=18))

#p987

#ggsave(filename = paste0(BASE, "RProject/OT-bd-analyses_loc/output/t_statAna1_cortical.pdf"), p987,
#       width = 10, height = 4, units = "in", device='pdf')




##  ---------------- plot expression values for subcortical regions in flat violin plot, vertical and horizontal violins

# resdat_left_sc <- resdat_left[-c(1:74),]   ??????

# prep data

# assign young genes dataframe to new objects to not overwrite the original objects, add col with donor id

p9861 <- pau9861_y
p9861$donor <- "d9861"
p9861 <- add_column(p9861, ELK1 = NA, .after = "CDKN1A")    # ELK1 is not present in 9861, therefore include it as a col with NAs
#p9861 <- subset(p9861, select=-c(mean9861))

p10021 <- pau10021_y
p10021$donor <- "d10021"

p12876 <- pau12876_y
p12876$donor <- "d12876"

p14380 <- pau14380_y 
p14380$donor <- "d14380"

p15496 <- pau15496_y
p15496$donor <- "d15496"

p15697 <- pau15697_y
p15697$donor <- "d15697"



# average expression for each gene across the six different datasets/donors and collect in a new data frame
colns <- colnames(pau10021_y[,c(2:12,14:20)])
plk <- list()
for (i in colns) {
  
  sumer <- (p9861[[i]] + p10021[[i]] + p12876[[i]] + p14380[[i]] + p15496[[i]] + p15697[[i]])/6
  plk[[i]] <- sumer
  jop <- do.call(rbind, plk)
}

jopt <- data.frame(t(jop))    # transpose df

ELK1m <- data.frame((p10021$ELK1 + p12876$ELK1 + p14380$ELK1 + p15496$ELK1 + p15697$ELK1)/5) # add ELK1 manually because it is not available
                                                                                             # for all donors (id 9861)

jopt <- add_column(jopt, ELK1 = ELK1m$X.p10021.ELK1...p12876.ELK1...p14380.ELK1...p15496.ELK1...p15697.ELK1..5, .after = "CDKN1A")

# add index column with region label indices
index <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16)
jopt <- cbind(index, jopt)

pauli_allR <- jopt # assign into a new object to again not overwrite old objects

pauli_allR <- cbind(pau_info$node_info, pauli_allR)                                      # add long region names
names(pauli_allR)[names(pauli_allR) == "pau_info$node_info"] <- "region"

pauli_allR$mean1 <-  meanS_dpL[meanS_dpL$structure=="subcortical","estimate_dpL"]        # add means from the t-test
pauli_allR$mean2 <-  rowSums(pauli_allR[, c(3:21)])/19                                   # calculate means manually
# pauli_allR$sd    <-  meanS_dpL[meanS_dpL$structure=="subcortical","sd"]                # add column with SD, EDIT: these are the wrong SDs

# the means differ slightly, the values in mean2 are about .002-.003 smaller than in mean1 - why?
# for mean1, I averaged across genes per donor, and then average across donors
# for mean2, I averaged across donors per gene, and then average across genes



# transform to long format
pauli_allR_L <- pauli_allR %>% pivot_longer(c("CACNB1", "CACNB2", "CACNB3", "CACNB4", "CACNG1", "CACNG2", "CACNG3", "CACNG6", "CACNG7", "CD38", "CDKN1A", "ELK1",
                                             "FOS", "NFATC3", "NFATC4", "NPPA", "OXT", "OXTR", "PIK3R5"),        # transform into long format, the way t.test likes it :)
                                               names_to = "genes", values_to = "expression")
pauli_allR_L <- data.frame(pauli_allR_L)
pauli_allR_L$region <- as.factor(pauli_allR_L$region)



# ------------------ half violin/raincloud plots


# ..... vertical version

# plot specific data prep
# order by mean expression, DESCENDING
pauli_all_R_L_ord <- pauli_allR_L %>% 
  arrange(desc(mean2))

# manually set levels to make sure they are plotted in descending order on the x-axis
order <- unique(pauli_all_R_L_ord$region)
pauli_all_R_L_ord$region <- factor(pauli_all_R_L_ord$region, levels = order)


# separate values for OXT, OXTR, and CD38
pauli_allR_L_cd38 <- subset(pauli_all_R_L_ord, grepl("^CD38$", pauli_all_R_L_ord$genes)) 
pauli_allR_L_oxt <- subset(pauli_all_R_L_ord, grepl("^OXT$", pauli_all_R_L_ord$genes)) 
pauli_allR_L_oxtr <- subset(pauli_all_R_L_ord, grepl("^OXTR$", pauli_all_R_L_ord$genes)) 
pauli_allR_L_ot <- rbind(pauli_allR_L_cd38, pauli_allR_L_oxt, pauli_allR_L_oxtr)
pauli_allR_L_ot <- pauli_allR_L_ot %>%                            # it should still be in the right order, but just to be on the very safe side
  arrange(desc(mean2))


# ... and remove them from base jitter (they are kept for the violin plot to ensure appropriate curving of the plots)
pat_ot <- c("OXT|OXTR|CD38")
pauli_allR_L_short <- subset(pauli_all_R_L_ord, !grepl(pat_ot, pauli_all_R_L_ord$genes))


# plotting
#g1 <- 
#  ggplot() +
#  geom_flat_violin(pauli_all_R_L_ord, mapping = aes(x = factor(region), y = expression,                         # use full dataset with all genes for violin plots
#                                                    fill = factor(region), color = factor(region)), 
#                   position = position_nudge(x = .2, y = 0), trim = F, alpha = .1, scale = "width") +
#  theme_classic() +
#  theme(axis.title.x = element_text(size = 22, face="bold"),
#        axis.title.y = element_text(size = 22, face="bold"),
#        axis.text.x = element_text(size=15, angle=50, vjust=1, hjust=1),
#        axis.text.y = element_text(size=15),
#        legend.position="none") +
#  ylab("mRNA intensity") + xlab("\nSubcortical brain regions")

#g1 <- g1 + geom_point(pauli_allR_L_short, mapping = aes(x= factor(region), y = expression,                     # use dataset w/o OT genes for base jitter
#                                                        fill = factor(region), color = factor(region)), 
#                      position = position_jitter(width = .15), size = 4, alpha = .6)

#g1 <- g1 + geom_point(pauli_allR_L_ot,    mapping = aes(x = factor(region), y = expression),                   # use dataset with only OT genes for highlight
#                      position = position_jitter(width = .15), size = 4, alpha = .5) 

#g1 <- g1 + geom_point(pauli_all_R_L_ord, mapping = aes(x = factor(region), y = mean2),                         # add mean expression points
#                    position = position_nudge(x = 0.3), size = 7) 

#g1 <- g1 + geom_hline(yintercept=mean(pauli_all_R_L_ord$mean2), linetype="dashed", size=.75)                   # add mean line of mean expression points

#g1

#ggsave(filename = paste0(BASE, "RProject/OT-bd-analyses_loc/output/subcortical_expression_analysis1_7.pdf"), g1,
#       width = 30, height = 10, units = "in", device='pdf')



# ..... horizontal version

# plot specific data prep
# order by mean expression, ASCENDING
pauli_all_R_L_asc <- pauli_allR_L %>% 
  arrange(mean2)

# manually set levels to make sure they are plotted in ascending order on the x-axis
order_2 <- unique(pauli_all_R_L_asc$region)
pauli_all_R_L_asc$region <- factor(pauli_all_R_L_asc$region, levels = order_2)

# order OT subset ascending
pauli_allR_L_asc_ot <- pauli_allR_L_ot %>%                           
  arrange(mean2)

# order subset w/o OT ascending
pauli_allR_L_asc_short <- pauli_allR_L_short %>%                           
  arrange(mean2)



# plotting
g2 <- 
  ggplot() +
  geom_density_ridges(pauli_all_R_L_asc, mapping = aes(x = expression, y = factor(region),                           # use full dataset with all genes for violin plots
                                                       fill = factor(region), colour = factor(region)), 
                      alpha=.1, scale= .9) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 23, face="bold"),
        axis.title.y = element_text(size = 23, face="bold"),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),
        legend.position="none") +
  ylab("Subcortical regions\n") + xlab("\nmRNA intensity") +
  viridis::scale_colour_viridis(option="viridis", discrete = T, begin = .2, end = .9) +
  viridis::scale_fill_viridis(option="viridis", discrete = T, begin = .2, end = .9) 

g2 <- g2 + geom_point(pauli_allR_L_asc_short, mapping = aes(x = expression, y = region,                               # use dataset w/o OT genes for base jitter
                                                            fill = factor(region), colour = factor(region)), 
                      position = position_jitter(height = .15), size = 4, alpha = .6)

g2 <- g2 + geom_point(pauli_allR_L_asc_ot,    mapping = aes(x = expression, y = factor(region)),                      # use dataset with only OT genes for highlight
                      position = position_jitter(height = .15), size = 4, alpha=.6) 

g2 <- g2 + geom_point(pauli_all_R_L_asc,      mapping = aes(x = mean2, y = factor(region)),                           # add mean expression points
                      position = position_nudge(y = .3), size = 7) 

g2 <- g2 + geom_vline(xintercept = mean(pauli_all_R_L_asc$mean2), linetype="dashed", size=.75)                        # add mean line of mean expression points

g2



#ggsave(filename = paste0(BASE, "RProject/OT-bd-analyses_loc/output/subcortical_expression_analysis1_hori3.pdf"), g2,
#       width = 15, height = 10, units = "in", device='pdf')


# arrange in a grid
g2_p098_2 <- plot_grid(g2, p098,
                       ncol = 2, nrow = 1,
                       rel_widths = c(1.75, 1),
                       labels = c('A', 'B'),
                       label_size = 23) 

g2_p098_2

ggsave(filename = paste0(BASE, "RProject/OT-bd-analyses_loc/output/ScC3_expression_analysis1.pdf"), g2_p098_2,
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
des9861_2 <- des9861
des9861_ALLs <- data.frame(rowMeans(des9861_2[,-1])) 
des9861_ALLs$label <- des9861_2$label


# 10021 - average expression across all genes for each cortical region
des10021_2 <- des10021
des10021_ALLs <- data.frame(rowMeans(des10021_2[,-1])) 
des10021_ALLs$label <- des10021_2$label


# 12876 - average expression across all genes for each cortical region
des12876_2 <- des12876
des12876_ALLs <- data.frame(rowMeans(des12876_2[,-1])) 
des12876_ALLs$label <- des12876_2$label


# 14380 - average expression across all genes for each cortical region
des14380_2 <- des14380
des14380_ALLs <- data.frame(rowMeans(des14380_2[,-1])) 
des14380_ALLs$label <- des14380_2$label


# 15496 - average expression across all genes for each cortical region
des15496_2 <- des15496
des15496_ALLs <- data.frame(rowMeans(des15496_2[,-1])) 
des15496_ALLs$label <- des15496_2$label


# 15697 - average expression across all genes for each cortical region
des15697_2 <- des15697
des15697_ALLs <- data.frame(rowMeans(des15697_2[,-1])) 
des15697_ALLs$label <- des15697_2$label



des_all_donors_ALL <- des9861_ALLs %>%                                     # join all donor averages into one df
  left_join(des10021_ALLs, by='label') %>%
  left_join(des12876_ALLs, by='label') %>%
  left_join(des14380_ALLs, by='label') %>%
  left_join(des15496_ALLs, by='label') %>%
  left_join(des15697_ALLs, by='label') 


des_all_donors_ALL <- cbind(des_all_donors_ALL$label, des_all_donors_ALL)
des_all_donors_ALL <- des_all_donors_ALL[, -3]

colnames(des_all_donors_ALL) <- c("label", "meanALL9861", "meanALL10021", "meanALL12876", "meanALL14380", "meanALL15496", "meanALL15697")


# average across donors. This is just for info and not needed in the analysis. Originally, it was implemented but is now deprecated.
#des_all_donors_ALL_t <- data.frame(t(des_all_donors_ALL[,-1]))
#colndada <- colnames(des_all_donors_ALL_t)
#ml_dada <- numeric()
#for (i in colndada) {
#  meansdada <- mean(des_all_donors_ALL_t[[i]])
#  ml_dada <- c(ml_dada, meansdada)
#}
#######


des_all_donors_ALLi <- dplyr::right_join(des_info_s, des_all_donors_ALL, by=c("index"="label"))   # annotate with region names

des_all_donors_ALLiL <- subset(des_all_donors_ALLi, grepl("^L ", des_all_donors_ALLi$name))  # select only left hemisphere
des_all_donors_ALLiL$name <- gsub('L ', '', des_all_donors_ALLiL$name)                       # remove left hemisphere indicator in the region names column
des_all_donors_ALLiL <- data.frame(des_all_donors_ALLiL)

# transform into long format, the way t.test likes it :)
des_all_donors_ALLiL_l <- des_all_donors_ALLiL %>% pivot_longer(c("meanALL9861", "meanALL10021", "meanALL12876", "meanALL14380", "meanALL15496", "meanALL15697"),        
                                               names_to = "donors", values_to = "expressionALL")
des_all_donors_ALLiL_l <- data.frame(des_all_donors_ALLiL_l)



### OT gene set subset into long format
des_ad_infoleft_L <- des_ad_infoleft %>% pivot_longer(c("mean9861", "mean10021", "mean12876", "mean14380", "mean15496", "mean15697"),        # transform into long format, the way t.test likes it :)
                                                           names_to = "donors", values_to = "expression")
des_ad_infoleft_L <- data.frame(des_ad_infoleft_L)
## ----------


# analysis

regions_desa2 <- unique(des_ad_infoleft_L$name)
pvals_desa2 <- numeric()
tStatistic_desa2 <- numeric()
cohensD_desa2 <- numeric()

estimateOT_desa2 <- numeric()
estimate_desa2 <- numeric()
sdOT_desa2 <- numeric()
sd_desa2 <- numeric()


for (i in 1:length(regions_desa2)){
  Temp_desa2 <- t.test(des_ad_infoleft_L[des_ad_infoleft_L$name==regions_desa2[i],"expression"], 
                 mu = mean(des_all_donors_ALLiL_l[des_all_donors_ALLiL_l$name==regions_desa2[i],"expressionALL"]), 
                 alternative = "two.sided")
  CohD_desa2 <- ES.t.one(m=Temp_desa2$estimate,
                   sd=sd(des_ad_infoleft_L[des_ad_infoleft_L$name==regions_desa2[i],"expression"]),
                   mu=mean(des_all_donors_ALLiL_l[des_all_donors_ALLiL_l$name==regions_desa2[i],"expressionALL"]),
                   alternative = "two.sided")
  pvals_desa2 <- c(pvals_desa2,Temp_desa2$p.value)
  tStatistic_desa2 <- c(tStatistic_desa2,Temp_desa2$statistic)
  cohensD_desa2 <- c(cohensD_desa2,CohD_desa2$d)
  
  estimateOT_desa2 <- c(estimateOT_desa2, Temp_desa2$estimate)
  estimate_desa2 <- c(estimate_desa2, Temp_desa2$null.value)
  sdOT_desa2 <- c(sdOT_desa2, sd(des_ad_infoleft_L[des_ad_infoleft_L$name==regions_desa2[i],"expression"]))
  sd_desa2 <- c(sd_desa2, sd(des_all_donors_ALLiL_l[des_all_donors_ALLiL_l$name==regions_desa2[i],"expressionALL"]))
}

results_desa2 <- data.frame(Region=regions_desa2,P_unadjusted=pvals_desa2,
                      t_stat=tStatistic_desa2,CohD=cohensD_desa2)

test_stats_desa2 <- data.frame(OTmean=estimateOT_desa2, mean=estimate_desa2, OTSD=sdOT_desa2, SD=sd_desa2)


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

 

# 9861 - average expression across all genes for each sub cortical region 
pau9861_2 <- pau9861
pau9861_ALLs <- data.frame(rowMeans(pau9861_2[,-1]))
pau9861_ALLs$label <- pau9861$label


# 10021 - average expression across all genes for each sub cortical region 
pau10021_2 <- pau10021
pau10021_ALLs <- data.frame(rowMeans(pau10021_2[,-1]))
pau10021_ALLs$label <- pau10021$label


# 12876 - average expression across all genes for each sub cortical region 
pau12876_2 <- pau12876
pau12876_ALLs <- data.frame(rowMeans(pau12876_2[,-1]))
pau12876_ALLs$label <- pau12876$label


# 14380 - average expression across all genes for each sub cortical region 
pau14380_2 <- pau14380
pau14380_ALLs <- data.frame(rowMeans(pau14380_2[,-1]))
pau14380_ALLs$label <- pau14380$label


# 15496 - average expression across all genes for each sub cortical region 
pau15496_2 <- pau15496
pau15496_ALLs <- data.frame(rowMeans(pau15496_2[,-1]))
pau15496_ALLs$label <- pau15496$label


# 15697 - average expression across all genes for each sub cortical region  
pau15697_2 <- pau15697
pau15697_ALLs <- data.frame(rowMeans(pau15697_2[,-1]))
pau15697_ALLs$label <- pau15697$label



pau_all_donors_ALL <- pau9861_ALLs %>%                                     # join all donor averages into one df
  left_join(pau10021_ALLs, by='label') %>%
  left_join(pau12876_ALLs, by='label') %>%
  left_join(pau14380_ALLs, by='label') %>%
  left_join(pau15496_ALLs, by='label') %>%
  left_join(pau15697_ALLs, by='label') 


pau_all_donors_ALL <- cbind(pau_all_donors_ALL$label, pau_all_donors_ALL)
pau_all_donors_ALL <- pau_all_donors_ALL[, -3]

colnames(pau_all_donors_ALL) <- c("label", "meanALL9861", "meanALL10021", "meanALL12876", "meanALL14380", "meanALL15496", "meanALL15697")


# average across donors. This is just for info and not needed in the analysis. Originally, it was implemented but is now deprecated.
#pau_all_donors_ALL_t <- data.frame(t(pau_all_donors_ALL[,-1]))
#colnpada <- colnames(pau_all_donors_ALL_t)
#ml_pada <- numeric()
#for (i in colnpada) {
#  meanspada <- mean(pau_all_donors_ALL_t[[i]])
#  ml_pada <- c(ml_pada, meanspada)
#}

########


pau_all_donors_ALL_info <- dplyr::right_join(pau_info, pau_all_donors_ALL, by=c("V1"="label"))     # annotate with region names
pau_all_donors_ALLi_L <- pau_all_donors_ALL_info %>% pivot_longer(c("meanALL9861", "meanALL10021", "meanALL12876", "meanALL14380", "meanALL15496", "meanALL15697"),        # transform into long format, the way t.test likes it :)
                                                                              names_to = "donors", values_to = "expressionALL")
pau_all_donors_ALLi_L <- data.frame(pau_all_donors_ALLi_L)


# transform young OT geneset data to long format
pau_ad_info_l <- pau_ad_info %>% pivot_longer(c("mean9861", "mean10021", "mean12876", "mean14380", "mean15496", "mean15697"),        # transform into long format, the way t.test likes it :)
                                                                     names_to = "donors", values_to = "expression")
pau_ad_info_l <- data.frame(pau_ad_info_l)



# analysis

regions_paua2 <- unique(pau_ad_info_l$node_info)
pvals_paua2 <- numeric()
tStatistic_paua2 <- numeric()
cohensD_paua2 <- numeric()

estimateOT_paua2 <- numeric()
estimate_paua2 <- numeric()
sdOT_paua2 <- numeric()
sd_paua2 <- numeric()



for (i in 1:length(regions_paua2)){
  Temp_paua2 <- t.test(pau_ad_info_l[pau_ad_info_l$node_info==regions_paua2[i],"expression"], 
                 mu = mean(pau_all_donors_ALLi_L[pau_all_donors_ALLi_L$node_info==regions_paua2[i],"expressionALL"]), 
                 alternative = "two.sided")
  CohD_paua2 <- ES.t.one(m=Temp_paua2$estimate,
                   sd=sd(pau_ad_info_l[pau_ad_info_l$node_info==regions_paua2[i],"expression"]),
                   mu=mean(pau_all_donors_ALLi_L[pau_all_donors_ALLi_L$node_info==regions_paua2[i],"expressionALL"]),
                   alternative = "two.sided")
  pvals_paua2 <- c(pvals_paua2,Temp_paua2$p.value)
  tStatistic_paua2 <- c(tStatistic_paua2,Temp_paua2$statistic)
  cohensD_paua2 <- c(cohensD_paua2,CohD_paua2$d)
  
  estimateOT_paua2 <- c(estimateOT_paua2, Temp_paua2$estimate)
  estimate_paua2 <- c(estimate_paua2, Temp_paua2$null.value)
  sdOT_paua2 <- c(sdOT_paua2, sd(pau_ad_info_l[pau_ad_info_l$node_info==regions_paua2[i],"expression"]))
  sd_paua2 <- c(sd_paua2, sd(pau_all_donors_ALLi_L[pau_all_donors_ALLi_L$node_info==regions_paua2[i],"expressionALL"]))
  
}

results_paua2 <- data.frame(Region=regions_paua2,P_unadjusted=pvals_paua2,
                         t_stat=tStatistic_paua2,CohD=cohensD_paua2)

test_stats_paua2 <- data.frame(OTmean=estimateOT_paua2, mean=estimate_paua2, OTSD=sdOT_paua2, SD=sd_paua2)


# combine with cortical and FDR it

## INFO ##
# It doesn't make a difference whether these are separately analyzed and then stacked or if all data is in one df. The t-tests work with multiple population means
# of ROWS (=brain regions), not one COLUMN (=whole brain), each t-test has its own population mean (as opposed to analysis #1). Therefore the population means do not
# depend on the mean of one column and the values in that column. If they would rely on the mean of the same column and would not change with each iteration based on 
# rows, then indeed it would make a difference because the values in the column the mean is based on would change (but again, that is not the case here). Anyways, see 
# end of script for analysis with all the data in one df.


res_dpa2 <- rbind(results_desa2, results_paua2)

res_dpa2$P_fdr <- p.adjust(res_dpa2$P_unadjusted, method="fdr", n = length(res_dpa2$P_unadjusted))

res_dpa2$stat_sign <- ""
res_dpa2$stat_sign[res_dpa2$P_fdr <= 0.05]  <- "*"

res_dpa2

# combine test statistics and parameters

test_stats_dpa2 <- rbind(test_stats_desa2, test_stats_paua2)
test_stats_dpa2 <- cbind( res_dpa2$Region, test_stats_dpa2)

# --------------------------------------------------------------------> significant results for some regions


# plotting

# plot expression values in dot graph
names(test_stats_dpa2)[names(test_stats_dpa2) == "res_dpa2$Region"] <- "region"

# annotate with full region names
des_ln <- read_excel(paste0(BASE, "RProject/OT-bd-analyses_loc/Data/destrieux_longnames.xlsx"))
df_dpa2 <- dplyr::full_join(des_ln, test_stats_dpa2, by=c("Short Name"="region"))
df_dpa2$`Long name (TA nomenclature is bold typed)`[75:90] <- df_dpa2$`Short Name`[75:90]
df_dpa2 <- data.frame(df_dpa2[, -c(1,2)])
names(df_dpa2)[names(df_dpa2) == "Long.name..TA.nomenclature.is.bold.typed."] <- "region"
df_dpa2$region <- factor(df_dpa2$region, levels=unique(df_dpa2$region))


hlinesOT <- c(mean(df_dpa2$OTmean) - mean(df_dpa2$OTSD), mean(df_dpa2$OTmean), mean(df_dpa2$OTmean) + mean(df_dpa2$OTSD))
hlinesAll <- c(mean(df_dpa2$mean) - mean(df_dpa2$SD), mean(df_dpa2$mean), mean(df_dpa2$mean) + mean(df_dpa2$SD))

pa2 <- ggplot(df_dpa2, aes(x=region)) + theme_classic() +
  geom_hline(yintercept=hlinesOT, linetype=c("dashed", "solid", "dashed"), color=c("#EA7E98", "#DE2056", "#EA7E98"), size=c(1, 2, 1)) +  
  #geom_hline(yintercept=hlinesOT, linetype=c("dashed", "solid", "dashed"), color=c("#97E1CE", "#16CC8F", "#97E1CE"), size=c(.5, 1, .5))  # for plasma
  geom_hline(yintercept=hlinesAll, linetype=c("dashed", "solid", "dashed"), color=c("gray50", "gray35", "gray50"), size=c(.5, 1, .5)) +
  geom_point(size=6, aes(y=mean))  +
  geom_errorbar(aes(ymin = mean - SD, ymax = mean + SD, width = .5)) 
  

  
pa2
pa2 <- pa2 +
  geom_point(size=6, aes(y=OTmean, colour = region)) +
  geom_errorbar(aes(ymin = OTmean - OTSD, ymax = OTmean + OTSD, color=factor(region), width = .5)) 
#  viridis::scale_colour_viridis(option="viridis", discrete = T, direction = 1, begin=.3, end=.9)  
# viridis::scale_colour_viridis(option="plasma", discrete = T, direction = 1, begin=.2, end=.875)         # for plasma 

pa2

pa2 <- pa2 + ylab("mRNA intensity\n") + xlab("\nBrain regions")
pa2 <- pa2 + theme(axis.text.x = element_text(angle=50, vjust=1, hjust=1, size=16),
                 axis.title.x = element_text(size=33, face="bold"),
                 axis.title.y = element_text(size=30, face="bold"),
                 axis.text.y = element_text(size=12),
                 legend.position ="none")  

pa2 <- pa2 + annotate("text", x = 66, y = .59, label = "*", size=8) +
             annotate("text", x = 67, y = .58, label = "*", size=8) +
             annotate("text", x = 78, y = .433, label = "*", size=8) +
             annotate("text", x = 79, y = .392, label = "*", size=8) 
pa2 

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











