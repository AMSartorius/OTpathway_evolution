
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



### >>>>>>>>>>>>>>>>>>> 1 HO cortical

## load HO data, subset, calculate means, merge, annotate (I'm sure this can be looped somehow but I'll deal with that later...)

# load HO region annotation
hoc_info <- read_csv(paste0(BASE, "data/processed/abagen_harvard-oxford/labels_ho_c.csv"))
hoc_info_s <- data.frame(hoc_info[-1,])                   # remove "Background" row



# load Donor 1 HO cortical 
hocD1 <- read_csv(paste0(BASE, "/data/processed/abagen_harvard-oxford/AHBAhocD1.csv"))
hocD1_y <- hocD1[, names(hocD1) %in% c("label", yOT$SYMBOL)]                # extract subset of modern OT genes, 19 modern genes, 148 cortical regions, lateralized 

ml_D1hoc <- data.frame(rowMeans(hocD1_y[,-1])) 

hocD1_lm <- data.frame(cbind(hocD1_y$label, ml_D1hoc$rowMeans.hocD1_y....1..))
colnames(hocD1_lm) <- c("label", "meanD1")


# load Donor 2 HO cortical 
hocD2 <- read_csv(paste0(BASE, "/data/processed/abagen_harvard-oxford/AHBAhocD2.csv"))
hocD2_y <- hocD2[, names(hocD2) %in% c("label", yOT$SYMBOL)]                # extract subset of modern OT genes, 19 modern genes, 148 cortical regions, lateralized 

ml_D2hoc <- data.frame(rowMeans(hocD2_y[,-1])) 

hocD2_lm <- data.frame(cbind(hocD2_y$label, ml_D2hoc$rowMeans.hocD2_y....1..))
colnames(hocD2_lm) <- c("label", "meanD2")


# load Donor 3 HO cortical  
hocD3 <- read_csv(paste0(BASE, "/data/processed/abagen_harvard-oxford/AHBAhocD3.csv"))
hocD3_y <- hocD3[, names(hocD3) %in% c("label", yOT$SYMBOL)]                # extract subset of modern OT genes, 19 modern genes, 148 cortical regions, lateralized 

ml_D3hoc <- data.frame(rowMeans(hocD3_y[,-1])) 

hocD3_lm <- data.frame(cbind(hocD3_y$label, ml_D3hoc$rowMeans.hocD3_y....1..))
colnames(hocD3_lm) <- c("label", "meanD3")


# load Donor 4 HO cortical  
hocD4 <- read_csv(paste0(BASE, "/data/processed/abagen_harvard-oxford/AHBAhocD4.csv"))
hocD4_y <- hocD4[, names(hocD4) %in% c("label", yOT$SYMBOL)]                # extract subset of modern OT genes, 19 modern genes, 148 cortical regions, lateralized 

ml_D4hoc <- data.frame(rowMeans(hocD4_y[,-1])) 

hocD4_lm <- data.frame(cbind(hocD4_y$label, ml_D4hoc$rowMeans.hocD4_y....1..))
colnames(hocD4_lm) <- c("label", "meanD4")



# load Donor 5 HO cortical  
hocD5 <- read_csv(paste0(BASE, "/data/processed/abagen_harvard-oxford/AHBAhocD5.csv"))
hocD5_y <- hocD5[, names(hocD5) %in% c("label", yOT$SYMBOL)]                # extract subset of modern OT genes, 19 modern genes, 148 cortical regions, lateralized 

ml_D5hoc <- data.frame(rowMeans(hocD5_y[,-1])) 

hocD5_lm <- data.frame(cbind(hocD5_y$label, ml_D5hoc$rowMeans.hocD5_y....1..))
colnames(hocD5_lm) <- c("label", "meanD5")


# load Donor 6 HO cortical  
hocD6 <- read_csv(paste0(BASE, "/data/processed/abagen_harvard-oxford/AHBAhocD6.csv"))
hocD6_y <- hocD6[, names(hocD6) %in% c("label", yOT$SYMBOL)]                # extract subset of modern OT genes, 19 modern genes, 148 cortical regions, lateralized 

ml_D6hoc <- data.frame(rowMeans(hocD6_y[,-1])) 

hocD6_lm <- data.frame(cbind(hocD6_y$label, ml_D6hoc$rowMeans.hocD6_y....1..))
colnames(hocD6_lm) <- c("label", "meanD6")



## INFO: "Medial wall" is missing in both hemispheres



hoc_all_d <- hocD1_lm %>%                   # join all donor averages into one df
  left_join(hocD2_lm, by='label') %>%
  left_join(hocD3_lm, by='label') %>%
  left_join(hocD4_lm, by='label') %>%
  left_join(hocD5_lm, by='label') %>%
  left_join(hocD6_lm, by='label') 

hoc_ad_i <- dplyr::right_join(hoc_info_s, hoc_all_d, by=c("X1"="label"))       # annotate with brain region labels
hoc_ad_i <- data.frame(hoc_ad_i)

hoc_ad_i <- hoc_ad_i %>%           # order ascending by region name to group into right and left
  arrange((X0))

hoc_ad_i$hemisphere <- "L"         # add hemisphere column with "L" and "R"
hoc_ad_i$hemisphere[49:96] <- "R"

hoc_ad_i <- hoc_ad_i %>%           # re-order descending by index
  arrange((X1))

hoc_ad_i$X0 <- gsub('Left ', 'lh_', hoc_ad_i$X0)                  # replace "L " in front of each region with "lh_"
hoc_ad_i$X0 <- gsub('Right ', 'rh_', hoc_ad_i$X0)                  # replace "R " in front of each region with "rh_"





### >>>>>>>>>>>>>>>>>>> 2 HO subcortical

## load HO data, subset, calculate means, merge, annotate (I'm sure this can be looped somehow but I'll deal with that later...)

# load pauli region annotation
hoSc_info <- read_csv(paste0(BASE, "data/processed/abagen_harvard-oxford/labels_ho_sc.csv"))
hoSc_info_s <- hoSc_info[-1,]   # remove "Background" row
hoSc_info_s <- data.frame(hoSc_info_s)


# load Donor 1 HO subcortical
hoScD1 <- read_csv(paste0(BASE, "/data/processed/abagen_harvard-oxford/AHBAhoScD1.csv"))
hoScD1_y <- hoScD1[, names(hoScD1) %in% c("label", yOT$SYMBOL)]                         # extract subset of modern OT genes, 18 modern genes, 16 subcortical regions !! ONE GENES LESS !!

ml_D1hoSc <- data.frame(rowMeans(hoScD1_y[,-1])) 

hoScD1_lm <- data.frame(cbind(hoScD1_y$label, ml_D1hoSc$rowMeans.hoScD1_y....1..))
colnames(hoScD1_lm) <- c("label", "meanD1")


# load Donor 2 HO subcortical
hoScD2 <- read_csv(paste0(BASE, "/data/processed/abagen_harvard-oxford/AHBAhoScD2.csv"))
hoScD2_y <- hoScD2[, names(hoScD2) %in% c("label", yOT$SYMBOL)]                         # extract subset of modern OT genes, 18 modern genes, 16 subcortical regions !! ONE GENES LESS !!

ml_D2hoSc <- data.frame(rowMeans(hoScD2_y[,-1])) 

hoScD2_lm <- data.frame(cbind(hoScD2_y$label, ml_D2hoSc$rowMeans.hoScD2_y....1..))
colnames(hoScD2_lm) <- c("label", "meanD2")


# load Donor 3 HO subcortical 
hoScD3 <- read_csv(paste0(BASE, "/data/processed/abagen_harvard-oxford/AHBAhoScD3.csv"))
hoScD3_y <- hoScD3[, names(hoScD3) %in% c("label", yOT$SYMBOL)]                         # extract subset of modern OT genes, 18 modern genes, 16 subcortical regions !! ONE GENES LESS !!

ml_D3hoSc <- data.frame(rowMeans(hoScD3_y[,-1])) 

hoScD3_lm <- data.frame(cbind(hoScD3_y$label, ml_D3hoSc$rowMeans.hoScD3_y....1..))
colnames(hoScD3_lm) <- c("label", "meanD3")


# load Donor 4 HO subcortical 
hoScD4 <- read_csv(paste0(BASE, "/data/processed/abagen_harvard-oxford/AHBAhoScD4.csv"))
hoScD4_y <- hoScD4[, names(hoScD4) %in% c("label", yOT$SYMBOL)]                         # extract subset of modern OT genes, 18 modern genes, 16 subcortical regions !! ONE GENES LESS !!

ml_D4hoSc <- data.frame(rowMeans(hoScD4_y[,-1])) 

hoScD4_lm <- data.frame(cbind(hoScD4_y$label, ml_D4hoSc$rowMeans.hoScD4_y....1..))
colnames(hoScD4_lm) <- c("label", "meanD4")


# load Donor 5 HO subcortical 
hoScD5 <- read_csv(paste0(BASE, "/data/processed/abagen_harvard-oxford/AHBAhoScD5.csv"))
hoScD5_y <- hoScD5[, names(hoScD5) %in% c("label", yOT$SYMBOL)]                         # extract subset of modern OT genes, 18 modern genes, 16 subcortical regions !! ONE GENES LESS !!

ml_D5hoSc <- data.frame(rowMeans(hoScD5_y[,-1])) 

hoScD5_lm <- data.frame(cbind(hoScD5_y$label, ml_D5hoSc$rowMeans.hoScD5_y....1..))
colnames(hoScD5_lm) <- c("label", "meanD5")


# load Donor 6 HO subcortical 
hoScD6 <- read_csv(paste0(BASE, "/data/processed/abagen_harvard-oxford/AHBAhoScD6.csv"))
hoScD6_y <- hoScD6[, names(hoScD6) %in% c("label", yOT$SYMBOL)]                         # extract subset of modern OT genes, 18 modern genes, 16 subcortical regions !! ONE GENES LESS !!

ml_D6hoSc <- data.frame(rowMeans(hoScD6_y[,-1])) 

hoScD6_lm <- data.frame(cbind(hoScD6_y$label, ml_D6hoSc$rowMeans.hoScD6_y....1..))
colnames(hoScD6_lm) <- c("label", "meanD6")



hoSc_all_d <- hoScD1_lm %>%                   # join all donor averages into one df
  left_join(hoScD2_lm, by='label') %>%
  left_join(hoScD3_lm, by='label') %>%
  left_join(hoScD4_lm, by='label') %>%
  left_join(hoScD5_lm, by='label') %>%
  left_join(hoScD6_lm, by='label') 

hoSc_ad_i <- dplyr::right_join(hoSc_info_s, hoSc_all_d, by=c("X1"="label"))       # annotate with brain region labels

hoSc_ad_i <- data.frame(hoSc_ad_i)                                                     

hoSc_ad_i$hemisphere <- "L"                                  # add hemisphere column with "L", "R" and "B"
hoSc_ad_i$hemisphere[8] <- "B"
hoSc_ad_i$hemisphere[12:21] <- "R"

hoSc_ad_i$X0 <- gsub('Left ', '', hoSc_ad_i$X0)                  # remove "L " in front of each region
hoSc_ad_i$X0 <- gsub('Right ', '', hoSc_ad_i$X0)                  # remove "R " in front of each region




### >>>>>>>>>>>>>>>>>>> 3 Combined cortical and subcortical analysis


## analyses only for the left hemisphere !!


# prepare HO cortical data, remove right hemisphere
hoc_ad_ileft <- subset(hoc_ad_i, grepl("L", hoc_ad_i$hemisphere))  


# prepare HO sub-cortical data, remove right hemisphere
hosc_pat <- c("L|B")
hoSc_ad_ileft <- subset(hoSc_ad_i, grepl(hosc_pat, hoSc_ad_i$hemisphere)) 


HO_cSc_l <- rbind(hoc_ad_ileft, hoSc_ad_ileft)   # combine into one df


names(HO_cSc_l)[names(HO_cSc_l) == "X1"] <- "index" 
names(HO_cSc_l)[names(HO_cSc_l) == "X0"] <- "region" 



# convert to long format
HO_cSc_lL <- HO_cSc_l %>% pivot_longer(c("meanD1", "meanD2", "meanD3", "meanD4", "meanD5", "meanD6"),        # transform into long format, the way t.test likes it :)
                                         names_to = "donors", values_to = "expression")
HO_cSc_lL <- data.frame(HO_cSc_lL)



mean_holL <- mean(HO_cSc_lL$expression)             
regions_holL <- unique(HO_cSc_lL$region)
pvals_holL <- numeric()
tStatistic_holL <- numeric()
cohensD_holL <- numeric()
estimate_holL <- numeric()
sd_holL <- numeric()

for (i in 1:length(regions_holL)){
  Temp_holL <- t.test(HO_cSc_lL[HO_cSc_lL$region==regions_holL[i],"expression"], 
                      mu = mean_holL, alternative = "two.sided")
  CohD_holL <- ES.t.one(m=Temp_holL$estimate,
                        sd=sd(HO_cSc_lL[HO_cSc_lL$region==regions_holL[i],"expression"]),
                        mu=mean_holL,
                        alternative = "two.sided")
  pvals_holL <- c(pvals_holL,Temp_holL$p.value)
  tStatistic_holL <- c(tStatistic_holL,Temp_holL$statistic)
  cohensD_holL <- c(cohensD_holL,CohD_holL$d)
  estimate_holL <- c(estimate_holL, Temp_holL$estimate)
  sd_holL <- c(sd_holL, sd(HO_cSc_lL[HO_cSc_lL$region==regions_holL[i],"expression"]))
}

SD_mean_holL <- sd(HO_cSc_lL$expression)

results_holL <- data.frame(Region=regions_holL,P_unadjusted=pvals_holL,
                           t_stat=tStatistic_holL,CohD=cohensD_holL)
results_holL$P_fdr <- p.adjust(results_holL$P_unadjusted, method="fdr", n = length(results_holL$P_unadjusted))

results_holL$stat_sign <- ""
results_holL$stat_sign[results_holL$P_fdr <= 0.05]  <- "*"

results_holL  


# --------------------------------------------------------------------> 1 significant region (brainstem)

meanS_holL <- data.frame(estimate_holL)
meanS_holL$mean_allstruc <- Temp_holL$null.value
meanS_holL <- cbind(regions_holL, meanS_holL)
meanS_holL$sd <- sd_holL
meanS_holL$structure <- "cortical"
meanS_holL$structure[49:59] <- "subcortical"

## CONCLUSION: Compared to the average expression of the young OT gene set across all regions, the young OT gene set is significantly lower 
#              expressed in 1 subcortical region.



# PLOTTING

## ---------------- plotting t-values on cortical map

# plot t values in cortical brain map


###### fetch HO cortical atlas map

ggseg_atlas_repos("HO")
ggsegHO <- ggseg_atlas_repos("HO", ignore.case = TRUE)
if (!requireNamespace("HO", quietly = TRUE)) {
  install_ggseg_atlas(repo = ggsegHO$repo, source = ggsegHO$source, force=T)
}


library(ggsegHO)

# inspect atlas
ggseg(atlas = ggsegHO::hoCort, 
      mapping=aes(fill = region)) +
  scale_fill_brain("hoCort", package="ggsegHO")

#####

HO_dat <- hoCort$data
HO_dat <- HO_dat %>% distinct(label, .keep_all = TRUE)    # remove duplicates
HO_ggseg_regions <- unique(HO_dat$label)
HO_ggseg_regions <- data.frame(HO_ggseg_regions)


# data preparation
res_holL <- results_holL
res_holL$Region[1:48] <- gsub(',', '.', res_holL$Region[1:48])                   # remove special characters and replace with "." in region names
res_holL$Region[1:48] <- gsub(' ', '.', res_holL$Region[1:48]) 
res_holL$Region[1:48] <- gsub("'", '.', res_holL$Region[1:48]) 
res_holL$Region[1:48] <- gsub("[(]", '.', res_holL$Region[1:48]) 
res_holL$Region[1:48] <- gsub("[)]", '.', res_holL$Region[1:48]) 
res_holL_c <- res_holL[-c(49:59),]             # remove subcortical regions because ggseg HO only supports cortical regions


# extract region names from the ggseg atlas because data is merged by the region column, not the label column
# later add the region names to the results df so ggseg has a column to merge by
HO_dat_s <- data.frame(HO_dat$region)
HO_dat_s <- cbind(HO_dat_s, HO_dat$label)


res_holL_j <- dplyr::full_join(HO_dat_s, res_holL_c, by=c("HO_dat$label"="Region"))            # merge region names, labels and results data
res_holL_j <- res_holL_j[-c(50:98),]                                                           # remove right hemisphere

names(res_holL_j)[names(res_holL_j) == "HO_dat$label"] <- "label"
names(res_holL_j)[names(res_holL_j) == "HO_dat.region"] <- "region"

# plotting

# vertical

p712 <- ggplot() +
  geom_brain(data = res_holL_j,
             atlas = hoCort,
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

ggsave(filename = paste0(BASE, "RProject/OT-bd-analyses_loc/output/apple.pdf"), p712,
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
cnames <- colnames(pauD1_y[,c(2:20)])       # any df will work, because they have the same gene colnames
e_list <- list()
for (i in cnames) {
  
  g_means <- (pD1[[i]] + pD2[[i]] + pD3[[i]] + pD4[[i]] + pD5[[i]] + pD6[[i]])/6
  e_list[[i]] <- g_means
  g_means_l <- do.call(rbind, e_list)
}

g_means_lt <- data.frame(t(g_means_l))    # transpose df

# add index column with region label indices
index <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16)
g_means_lt <- cbind(index, g_means_lt)

pauli_aR <- g_means_lt # assign into a new object to again not overwrite old objects

pauli_aR <- cbind(pau_info$node_info, pauli_aR)                                      # add long region names
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


# ..... horizontal version

# plot specific data prep
# order by mean expression, DESCENDING
pauli_aR_L_ord <- pauli_aR_L %>% 
  arrange(desc(mean2))

# manually set levels to make sure they are plotted in descending order on the x-axis
order <- unique(pauli_aR_L_ord$region)
pauli_aR_L_ord$region <- factor(pauli_aR_L_ord$region, levels = order)


# separate values for OXT, OXTR, and CD38
pauli_aR_L_cd38 <- subset(pauli_aR_L_ord, grepl("^CD38$", pauli_aR_L_ord$genes)) 
pauli_aR_L_oxt <- subset(pauli_aR_L_ord, grepl("^OXT$", pauli_aR_L_ord$genes)) 
pauli_aR_L_oxtr <- subset(pauli_aR_L_ord, grepl("^OXTR$", pauli_aR_L_ord$genes)) 
pauli_aR_L_ot <- rbind(pauli_aR_L_cd38, pauli_aR_L_oxt, pauli_aR_L_oxtr)     # combine
pauli_aR_L_ot <- pauli_aR_L_ot %>%                                               # it should still be in the right order, but just to be on the very safe side
  arrange(desc(mean2))


# ... and remove them from base jitter (they are kept for the violin plot to ensure appropriate curving of the plots)
pat_ot <- c("OXT|OXTR|CD38")
pauli_aR_L_short <- subset(pauli_aR_L_ord, !grepl(pat_ot, pauli_aR_L_ord$genes))




# ..... vertical version

# plot specific data prep
# order by mean expression, ASCENDING
pauli_aR_L_asc <- pauli_aR_L %>% 
  arrange(mean2)

# manually set levels to make sure they are plotted in ascending order on the x-axis
order_2 <- unique(pauli_aR_L_asc$region)
pauli_aR_L_asc$region <- factor(pauli_aR_L_asc$region, levels = order_2)

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

g22 <- g22 + geom_point(pauli_aR_L_asc_short, mapping = aes(x = expression, y = region,                               # use dataset w/o OT genes for base jitter
                                                            fill = factor(region), colour = factor(region)), 
                        position = position_jitter(height = .15), size = 4, alpha = .6)

g22 <- g22 + geom_point(pauli_aR_L_asc_ot,    mapping = aes(x = expression, y = factor(region)),                      # use dataset with only OT genes for highlight
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

ggsave(filename = paste0(BASE, "RProject/OT-bd-analyses_loc/output/Ana1_subcortical_vert_v2.pdf"), g22,
       width = 10, height = 15, units = "in", device='pdf')


# arrange in a grid
g22_p712 <- plot_grid(g22, p712,
                      ncol = 2, nrow = 1,
                      rel_widths = c(1.75, 1),
                      labels = c('A', 'B'),
                      label_size = 23) 

g22_p712

ggsave(filename = paste0(BASE, "RProject/OT-bd-analyses_loc/output/Ana1_ScC_combined.pdf"), g22_p712,
       width = 15, height = 11, units = "in", device='pdf')







###############################################################################
## ANALYSIS 2: one-sample t-test for each region with different population mean
###############################################################################


## compare the expression of young OT genset to the average expression of all protein coding genes for each region separately, no relation to the average expression across the brain.

# .... population mean is the average expression of all genes in region x, sample mean is the average expression of young OT genes subset in region x



### >>>>>>>>>>>>>>>>>>> DESTRIEUX cortical

## subset HO cortical data, calculate means, merge, annotate (I'm sure this can be looped somehow but I'll deal with that later...)


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



# Donor 1 - average expression across all genes for each cortical region 
hocD1_2 <- hocD1
hocD1_2_ALLs <- data.frame(rowMeans(hocD1_2[,-1])) 
hocD1_2_ALLs$label <- hocD1_2$label


# Donor 2 - average expression across all genes for each cortical region
hocD2_2 <- hocD2
hocD2_2_ALLs <- data.frame(rowMeans(hocD2_2[,-1])) 
hocD2_2_ALLs$label <- hocD2_2$label


# Donor 3 - average expression across all genes for each cortical region
hocD3_2 <- hocD3
hocD3_2_ALLs <- data.frame(rowMeans(hocD3_2[,-1])) 
hocD3_2_ALLs$label <- hocD3_2$label


# Donor 4 - average expression across all genes for each cortical region
hocD4_2 <- hocD4
hocD4_2_ALLs <- data.frame(rowMeans(hocD4_2[,-1])) 
hocD4_2_ALLs$label <- hocD4_2$label


# Donor 5 - average expression across all genes for each cortical region
hocD5_2 <- hocD5
hocD5_2_ALLs <- data.frame(rowMeans(hocD5_2[,-1])) 
hocD5_2_ALLs$label <- hocD5_2$label


# Donor 6 - average expression across all genes for each cortical region
hocD6_2 <- hocD6
hocD6_2_ALLs <- data.frame(rowMeans(hocD6_2[,-1])) 
hocD6_2_ALLs$label <- hocD6_2$label




hoc_all_d_ALL <- hocD1_2_ALLs %>%                                     # join all donor averages into one df
  left_join(hocD2_2_ALLs, by='label') %>%
  left_join(hocD3_2_ALLs, by='label') %>%
  left_join(hocD4_2_ALLs, by='label') %>%
  left_join(hocD5_2_ALLs, by='label') %>%
  left_join(hocD6_2_ALLs, by='label') 


hoc_all_d_ALL <- cbind(hoc_all_d_ALL$label, hoc_all_d_ALL)
hoc_all_d_ALL <- hoc_all_d_ALL[, -3]

colnames(hoc_all_d_ALL) <- c("label", "meanALLD1", "meanALLD2", "meanALLD3", "meanALLD4", "meanALLD5", "meanALLD6")


# average across donors. This is just for info and not needed in the analysis. Originally, it was implemented but is now deprecated.
#des_all_donors_ALL_t <- data.frame(t(des_all_donors_ALL[,-1]))
#colndada <- colnames(des_all_donors_ALL_t)
#ml_dada <- numeric()
#for (i in colndada) {
#  meansdada <- mean(des_all_donors_ALL_t[[i]])
#  ml_dada <- c(ml_dada, meansdada)
#}
#######


hoc_all_d_ALLi <- dplyr::right_join(hoc_info_s, hoc_all_d_ALL, by=c("X1"="label"))   # annotate with region names

hoc_all_d_ALLiL <- subset(hoc_all_d_ALLi, grepl("^Left ", hoc_all_d_ALLi$X0))       # select only left hemisphere
hoc_all_d_ALLiL$X0 <- gsub('Left ', '', hoc_all_d_ALLiL$X0)                        # remove left hemisphere indicator in the region names column
hoc_all_d_ALLiL <- data.frame(hoc_all_d_ALLiL)

# transform into long format, the way t.test likes it :)
hoc_all_d_ALLiL_l <- hoc_all_d_ALLiL %>% pivot_longer(c("meanALLD1", "meanALLD2", "meanALLD3", "meanALLD4", "meanALLD5", "meanALLD6"),        
                                                      names_to = "donors", values_to = "expressionALL")
hoc_all_d_ALLiL_l <- data.frame(hoc_all_d_ALLiL_l)



### OT gene set subset into long format
hoc_ad_ileft_2 <- hoc_ad_ileft
hoc_ad_ileft_2$X0 <- gsub('lh_', '', hoc_ad_ileft_2$X0)                  

hoc_ad_ileft_L <- hoc_ad_ileft_2 %>% pivot_longer(c("meanD1", "meanD2", "meanD3", "meanD4", "meanD5", "meanD6"),        # transform into long format, the way t.test likes it :)
                                                names_to = "donors", values_to = "expression")
hoc_ad_ileft_L <- data.frame(hoc_ad_ileft_L)
## ----------


# analysis

regions_hoc2 <- unique(hoc_ad_ileft_L$X0)
pvals_hoc2 <- numeric()
tStatistic_hoc2 <- numeric()
cohensD_hoc2 <- numeric()

estimateOT_hoc2 <- numeric()
estimate_hoc2 <- numeric()
sdOT_hoc2 <- numeric()
sd_hoc2 <- numeric()


for (i in 1:length(regions_hoc2)){
  Temp_hoc2 <- t.test(hoc_ad_ileft_L[hoc_ad_ileft_L$X0==regions_hoc2[i],"expression"], 
                      mu = mean(hoc_all_d_ALLiL_l[hoc_all_d_ALLiL_l$X0==regions_hoc2[i],"expressionALL"]), 
                      alternative = "two.sided")
  CohD_hoc2 <- ES.t.one(m=Temp_hoc2$estimate,
                        sd=sd(hoc_ad_ileft_L[hoc_ad_ileft_L$X0==regions_hoc2[i],"expression"]),
                        mu=mean(hoc_all_d_ALLiL_l[hoc_all_d_ALLiL_l$X0==regions_hoc2[i],"expressionALL"]),
                        alternative = "two.sided")
  pvals_hoc2 <- c(pvals_hoc2,Temp_hoc2$p.value)
  tStatistic_hoc2 <- c(tStatistic_hoc2,Temp_hoc2$statistic)
  cohensD_hoc2 <- c(cohensD_hoc2,CohD_hoc2$d)
  
  estimateOT_hoc2 <- c(estimateOT_hoc2, Temp_hoc2$estimate)
  estimate_hoc2 <- c(estimate_hoc2, Temp_hoc2$null.value)
  sdOT_hoc2 <- c(sdOT_hoc2, sd(hoc_ad_ileft_L[hoc_ad_ileft_L$X0==regions_hoc2[i],"expression"]))
  sd_hoc2 <- c(sd_hoc2, sd(hoc_all_d_ALLiL_l[hoc_all_d_ALLiL_l$X0==regions_hoc2[i],"expressionALL"]))
}

results_hoc2 <- data.frame(Region=regions_hoc2,P_unadjusted=pvals_hoc2,
                           t_stat=tStatistic_hoc2,CohD=cohensD_hoc2)

test_stats_hoc2 <- data.frame(OTmean=estimateOT_hoc2, mean=estimate_hoc2, OTSD=sdOT_hoc2, SD=sd_hoc2)


# combine with subcortical and then FDR it





### >>>>>>>>>>>>>>>>>>> HO subcortical

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
hoScD1_2 <- hoScD1
hoScD1_ALLs <- data.frame(rowMeans(hoScD1_2[,-1]))
hoScD1_ALLs$label <- hoScD1_2$label


# Donor 2 - average expression across all genes for each sub cortical region 
hoScD2_2 <- hoScD2
hoScD2_ALLs <- data.frame(rowMeans(hoScD2_2[,-1]))
hoScD2_ALLs$label <- hoScD2_2$label


# Donor 3 - average expression across all genes for each sub cortical region 
hoScD3_2 <- hoScD3
hoScD3_ALLs <- data.frame(rowMeans(hoScD3_2[,-1]))
hoScD3_ALLs$label <- hoScD3_2$label


# Donor 4 - average expression across all genes for each sub cortical region 
hoScD4_2 <- hoScD4
hoScD4_ALLs <- data.frame(rowMeans(hoScD4_2[,-1]))
hoScD4_ALLs$label <- hoScD4_2$label


# Donor 5 - average expression across all genes for each sub cortical region 
hoScD5_2 <- hoScD5
hoScD5_ALLs <- data.frame(rowMeans(hoScD5_2[,-1]))
hoScD5_ALLs$label <- hoScD5_2$label


# Donor 6 - average expression across all genes for each sub cortical region  
hoScD6_2 <- hoScD6
hoScD6_ALLs <- data.frame(rowMeans(hoScD6_2[,-1]))
hoScD6_ALLs$label <- hoScD6_2$label



hoSc_all_d_ALL <- hoScD1_ALLs %>%                                     # join all donor averages into one df
  left_join(hoScD2_ALLs, by='label') %>%
  left_join(hoScD3_ALLs, by='label') %>%
  left_join(hoScD4_ALLs, by='label') %>%
  left_join(hoScD5_ALLs, by='label') %>%
  left_join(hoScD6_ALLs, by='label') 


hoSc_all_d_ALL <- cbind(hoSc_all_d_ALL$label, hoSc_all_d_ALL)
hoSc_all_d_ALL <- hoSc_all_d_ALL[, -3]

colnames(hoSc_all_d_ALL) <- c("label", "meanALLD1", "meanALLD2", "meanALLD3", "meanALLD4", "meanALLD5", "meanALLD6")


# average across donors. This is just for info and not needed in the analysis. Originally, it was implemented but is now deprecated.
#pau_all_donors_ALL_t <- data.frame(t(pau_all_donors_ALL[,-1]))
#colnpada <- colnames(pau_all_donors_ALL_t)
#ml_pada <- numeric()
#for (i in colnpada) {
#  meanspada <- mean(pau_all_donors_ALL_t[[i]])
#  ml_pada <- c(ml_pada, meanspada)
#}

########

hoSc_all_d_ALLi <- dplyr::right_join(hoSc_info_s, hoSc_all_d_ALL, by=c("X1"="label"))   # annotate with region names

hosc_pat2 <- c("Left|Brain")
hoSc_all_d_ALLiL <- subset(hoSc_all_d_ALLi, grepl(hosc_pat2, hoSc_all_d_ALLi$X0))       # select only left hemisphere
hoSc_all_d_ALLiL$X0 <- gsub('Left ', '', hoSc_all_d_ALLiL$X0)                           # remove left hemisphere indicator in the region names column
hoSc_all_d_ALLiL <- data.frame(hoSc_all_d_ALLiL)


# transform into long format, the way t.test likes it :)
hoSc_all_d_ALLiL_l <- hoSc_all_d_ALLiL %>% pivot_longer(c("meanALLD1", "meanALLD2", "meanALLD3", "meanALLD4", "meanALLD5", "meanALLD6"),        
                                                        names_to = "donors", values_to = "expressionALL")
hoSc_all_d_ALLiL_l <- data.frame(hoSc_all_d_ALLiL_l)


# transform young OT geneset data to long format
hoSc_ad_ileft_l <- hoSc_ad_ileft %>% pivot_longer(c("meanD1", "meanD2", "meanD3", "meanD4", "meanD5", "meanD6"),        # transform into long format, the way t.test likes it :)
                                             names_to = "donors", values_to = "expression")
hoSc_ad_ileft_l <- data.frame(hoSc_ad_ileft_l)



# analysis

regions_hoSc2 <- unique(hoSc_ad_ileft_l$X0)
pvals_hoSc2 <- numeric()
tStatistic_hoSc2 <- numeric()
cohensD_hoSc2 <- numeric()

estimateOT_hoSc2 <- numeric()
estimate_hoSc2 <- numeric()
sdOT_hoSc2 <- numeric()
sd_hoSc2 <- numeric()



for (i in 1:length(regions_hoSc2)){
  Temp_hoSc2 <- t.test(hoSc_ad_ileft_l[hoSc_ad_ileft_l$X0==regions_hoSc2[i],"expression"], 
                       
                      mu = mean(hoSc_all_d_ALLiL_l[hoSc_all_d_ALLiL_l$X0==regions_hoSc2[i],"expressionALL"]), 
                      
                      alternative = "two.sided")
  
  CohD_hoSc2 <- ES.t.one(m=Temp_hoSc2$estimate,
                         
                        sd=sd(hoSc_ad_ileft_l[hoSc_ad_ileft_l$X0==regions_hoSc2[i],"expression"]),
                        
                        mu=mean(hoSc_all_d_ALLiL_l[hoSc_all_d_ALLiL_l$X0==regions_hoSc2[i],"expressionALL"]),
                        
                        alternative = "two.sided")
  
  pvals_hoSc2 <- c(pvals_hoSc2,Temp_hoSc2$p.value)
  
  tStatistic_hoSc2 <- c(tStatistic_hoSc2,Temp_hoSc2$statistic)
  
  cohensD_hoSc2 <- c(cohensD_hoSc2,CohD_hoSc2$d)
  
  
  
  estimateOT_hoSc2 <- c(estimateOT_hoSc2, Temp_hoSc2$estimate)
  
  estimate_hoSc2 <- c(estimate_hoSc2, Temp_hoSc2$null.value)
  
  sdOT_hoSc2 <- c(sdOT_hoSc2, sd(hoSc_ad_ileft_l[hoSc_ad_ileft_l$X0==regions_hoSc2[i],"expression"]))
  
  sd_hoSc2 <- c(sd_hoSc2, sd(hoSc_all_d_ALLiL_l[hoSc_all_d_ALLiL_l$X0==regions_hoSc2[i],"expressionALL"]))
  
}

results_hoSc2 <- data.frame(Region=regions_hoSc2,P_unadjusted=pvals_hoSc2,
                           t_stat=tStatistic_hoSc2,CohD=cohensD_hoSc2)

test_stats_hoSc2 <- data.frame(OTmean=estimateOT_hoSc2, mean=estimate_hoSc2, OTSD=sdOT_hoSc2, SD=sd_hoSc2)


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
  viridis::scale_colour_viridis(option="turbo", discrete=T)

lpp <- lpp + ylab("mRNA intensity\n") + xlab("\nBrain regions") +
  theme(axis.title = element_text(size=25, face="bold"),
        axis.text.x = element_text(angle=60, size=14, vjust =1, hjust = 1),
        # axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=14),
        legend.position = c(0.37, 0.2),
        legend.direction = "horizontal") 

lpp <- lpp + scale_fill_discrete(labels = c("Frontal", "Fronto-parietal", "Insula", "Lateral sulcus", "Limbic", "Occipital",
                                            "Parietal", "Parieto-limbic", "Parieto-occipital", "Temporal", "Temporo-occipital",
                                            "Sub-cortical"))
lpp <- lpp + 
  annotate('text', x = 52, y = .557, label = '*', colour= "black", size=11) +
  annotate('text', x = 75, y = .427, label = '*', colour= "black", size=11) +
  annotate('text', x = 81, y = .452, label = '*', colour= "black", size=11) 

lpp <- lpp + 
  scale_x_discrete(labels= xtext)

lpp








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











