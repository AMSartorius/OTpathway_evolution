######################################
##### SET UP WORKING ENVIRONMENT #####
######################################

rm(list = ls()) # delete all objects in the workspace
gc(reset = T)   # resest memory (especially useful when working with large data sets)

options(stringsAsFactors = F) # disables automatic conversion of char. strings into factors
Sys.setenv(LANG = "en")

### .......... loading packages
# library(powerAnalysis)   # for "ES.t.one" function
library(writexl)         # export .xlsx files 
library(readxl)          # load/read in .xlsx files 
library(readr)           # load/read in .csv files
library(tidyverse)       # for everything and pipes
library(DataCombine)     # for "InsertRow" function

# for data visualization
library(ggseg)
# library(ggsegDefaultExtra)
library(ggridges)   
library(cowplot)



# get modern vertebrate OT gene ids and their gene symbols for later subsetting
OTpthwy.res <- data.frame(read_excel("output/BLASTp/OTpthwy_res02.xlsx", sheet = 1))
mv.PS <- c(11, 12, 13, 14, 15, 16, 17, 19, 20)
mvOTids <- data.frame(toupper(OTpthwy.res[OTpthwy.res$Phylostratum %in% mv.PS, "Gene"]))   # get gene labels for modern vertebrate OT pathway genes
colnames(mvOTids) <- "mvGenes"


#####################################################################
## ANALYSIS 1
#####################################################################


### >>>>>>>>>>>>>>>>>>> Data loading and prep

## load DK data, subset, calculate means, merge, annotate (I'm sure this can be looped somehow but I'll deal with that later...)

# 17 of the 98 genes are not available in DK abagen output, hence there are 81 genes


# load DK region annotation
dk.info <- read_csv("data/raw/abagen_dk/atlas-desikankilliany.csv")

pat.dk <- c("9861", "10021", "12876", "14380", "15496", "15697")

list.dk1 <- list()
list.dk2 <- list()
list.dk3 <- list()
list.dk4 <- list()

for (i in pat.dk) {
  list.dk1[[i]] <- read_csv(paste0("output/abagen_dk/AHBAdk_220513_", i, ".csv"))          # fetch original files with all genes
}

for (i in 1:length(pat.dk)) {
  list.dk2[[i]] <- list.dk1[[i]][, names(list.dk1[[i]]) %in% c("label", mvOTids$mvGenes)]  # extract subset of modern vertebrate OT genes, 98 modern genes, 83 regions, lateralized
  list.dk3[[i]] <- data.frame(rowMeans(list.dk2[[i]][, -1]))                               # calculate row means
  list.dk4[[i]] <- data.frame(cbind(list.dk2[[i]][["label"]], list.dk3[[i]][[1]]))         # merge again with label column
  colnames(list.dk4[[i]]) <- c("label", paste0("meanDK", i))                               # change column names for readability
  names(list.dk4)[i] <- paste0("D", pat.dk[i]) 
  names(list.dk2)[i] <- paste0("D.org", pat.dk[i])
  
 }    


list2env(list.dk4, envir = .GlobalEnv)   # unlists list with dataframes and pastes in the environment
list2env(list.dk2, envir = .GlobalEnv)


DK.all <- D9861 %>%                      # join all donor averages into one df
  left_join(D10021, by = 'label') %>%
  left_join(D12876, by = 'label') %>%
  left_join(D14380, by = 'label') %>%
  left_join(D15496, by = 'label') %>%
  left_join(D15697, by = 'label') 


DK.all.i <- data.frame(dplyr::right_join(dk.info, DK.all, by = c("id" = "label")))       # annotate with brain region labels




### >>>>>>>>>>>>>>>>>>> analysis


## analyses only for the left hemisphere !!


# prepare desikan killiany data, remove right hemisphere
pat.ana1 <- c("L|B")                                                  # pattern for hemispheres to keep, L = left, B = subcortex/brainstem
DK.all.il <- subset(DK.all.i, grepl(pat.ana1, DK.all.i$hemisphere))   # DK_ad_il

# convert to long format
DK.all.ilL <- data.frame(DK.all.il %>% pivot_longer(c("meanDK1", "meanDK2", "meanDK3", "meanDK4", "meanDK5", "meanDK6"),        # transform into long format, the way t.test likes it :)
                                       names_to = "donors", values_to = "expression"))


mean.dklL       <- mean(DK.all.ilL$expression)        # mean across all donors and all brain regions      
regions.dklL    <- unique(DK.all.ilL$label)           # string with unique region names to avoid duplicates
pvals.dklL      <- numeric() 
tStatistic.dklL <- numeric()
cohensD.dklL    <- numeric()
estimate.dklL   <- numeric()
sd.dklL         <- numeric()

for (i in 1:length(regions.dklL)){
  Temp.dklL <- t.test(DK.all.ilL[DK.all.ilL$label == regions.dklL[i], "expression"], 
                      mu = mean.dklL, alternative = "two.sided")
# CohD.dklL <- ES.t.one(m           = Temp.dklL$estimate,
#                       sd          = sd(DK.all.ilL[DK.all.ilL$label == regions.dklL[i], "expression"]),
#                       mu          = mean.dklL,
#                       alternative = "two.sided")
  pvals.dklL      <- c(pvals.dklL, Temp.dklL$p.value)
  tStatistic.dklL <- c(tStatistic.dklL, Temp.dklL$statistic)
#  cohensD.dklL    <- c(cohensD.dklL, CohD.dklL$d)
  estimate.dklL   <- c(estimate.dklL, Temp.dklL$estimate)
  sd.dklL         <- c(sd.dklL, sd(DK.all.ilL[DK.all.ilL$label == regions.dklL[i], "expression"]))
}

SD.mean.dklL <- sd(DK.all.ilL$expression)            # overall standard deviation

results.dklL <- data.frame(Region = regions.dklL, P.unadjusted = pvals.dklL,
                           t.stat = tStatistic.dklL) #, CohD = cohensD.dklL)

results.dklL$P.fdr <- p.adjust(results.dklL$P.unadjusted, method = "fdr", n = length(results.dklL$P.unadjusted))

results.dklL$stat.sign <- ""
results.dklL$stat.sign[results.dklL$P.fdr <= 0.05]  <- "*"

results.dklL  

# ....... for supplementary materials
# write_xlsx(results.dklL, "output/supp_mat_new/sup_mat_07.xlsx")

# --------------------------------------------------------------------> 9 significant regions


# PLOTTING

# plot expression values in dot graph
meanS.dklL <- data.frame(estimate.dklL)
meanS.dklL$mean.all.struc <- Temp.dklL$null.value
meanS.dklL <- cbind(regions.dklL, meanS.dklL)
meanS.dklL$sd <- sd.dklL
meanS.dklL$structure <- "cortical"
meanS.dklL$structure[35:42] <- "subcortical.brainstem"



## ---------------- plotting t-values on cortical map

# plot t values in cortical brain map

dk.dat <- dk$data                                 # "dk" is a dataframe with coordinates from the ggseg package universe
DKreg.regions <- data.frame(dk.dat$region)
DKreg.regions$labs <- dk.dat$label
DKreg.regions <- DKreg.regions %>% distinct(dk.dat.region, .keep_all = TRUE)      # remove duplicates
DKreg.regions <- data.frame(DKreg.regions)


# data preparation
DKresdat.l <- results.dklL
DKresdat.l$Region <- Map(paste0, 'lh_', DKresdat.l$Region)                      # add "lh_" in front of each region name
DKresdat.l <- DKresdat.l[-c(35:42),]                                            # remove sub-cortical regions
DKresdat.l <- InsertRow(DKresdat.l, c(NA, NA, NA, NA, NA, NA), RowNum = 1)        # add row in the beginning of the df to match the df structure of "DKreg.regions"
DKresdat.l <- InsertRow(DKresdat.l, c("lh_corpuscallosum", NA, NA, NA, NA, NA), RowNum = 24)    # add row at position 24 to match the df structure of "DKreg.regions"
DKresdat.l <- data.frame(DKresdat.l[match(DKreg.regions$labs, DKresdat.l$Region),])    # reorder rows according to the order of regions in the atlas provided by ggseg
DKresdat.l <- cbind(DKreg.regions$dk.dat.region, DKresdat.l)                           # finally add region column

names(DKresdat.l)[names(DKresdat.l) == "DKreg.regions$dk.dat.region"] <- "region"      # rename to match "DKreg.regions"

DKresdat.l$t.stat <- as.numeric(DKresdat.l$t.stat)



# plotting

# vertical

p127 <-  ggplot() + theme_void() +
  geom_brain(data = DKresdat.l, 
             atlas = dk, 
             colour = "white",
             position = position_brain(side + hemi ~ .),
             aes(fill = t.stat)) +
  viridis::scale_fill_viridis(option = "mako", discrete = F) +
  labs(fill = 't statistic') 

p127

p127 <- p127 + ylab("\nlateral   |   medial\n") + labs(title = "Cortical expression") +
  theme(plot.title   = element_text(hjust = .5, size = 25, face = "bold"),
        axis.ticks.y = element_blank(),
        axis.text.y  = element_blank(),
        axis.title.y = element_text(size = 23, angle = 90, face = "bold"),
        axis.ticks.x = element_blank(),
        axis.text.x  = element_blank(),
        legend.title = element_text(size = 23),
        legend.text  = element_text(size = 16),
        legend.key.size = unit(1, 'cm'))

p127

# ....... for supplementary materials
# write_xlsx(DKresdat.l, "output/supp_mat_new/sup_mat_10a.xlsx")

# ggsave(filename = "output/abagen_analysis01.pdf", p127,
#        width = 5, height = 10, units = "in", device = 'pdf')




##  ---------------- plot expression values for subcortical regions in flat violin plot, vertical and horizontal violins

# resdat_left_sc <- resdat_left[-c(1:74),]   ??????

# prep data

# assign young genes df to new objects to not overwrite the original objects, add col with donor id

mvdk9861 <- D.org9861
mvdk9861$donor <- "9861"

mvdk10021 <- D.org10021
mvdk10021$donor <- "10021"

mvdk12876 <- D.org12876
mvdk12876$donor <- "12876"

mvdk14380 <- D.org14380
mvdk14380$donor <- "14380"

mvdk15496 <- D.org15496
mvdk15496$donor <- "15496"

mvdk15697 <- D.org15697
mvdk15697$donor <- "15697"



# average expression for each gene across the six different datasets/donors and collect in a new data frame
# fetch gene colnames
cnames <- colnames(mvdk9861[,c(2:82)])           # any df will work, because they have the same gene colnames
e.list <- list()
for (i in cnames) {
  g.means <- (mvdk9861[[i]] + mvdk10021[[i]] + mvdk12876[[i]] + mvdk14380[[i]] + mvdk15496[[i]] + mvdk15697[[i]])/6
  e.list[[i]] <- g.means
  g.means.l <- do.call(rbind, e.list)
}

g.means.lt <- data.frame(t(g.means.l))             # transpose df

g.means.lt <- cbind(dk.info, g.means.lt)
g.means.lt.sc <- subset(g.means.lt, grepl("subcortex/brainstem", g.means.lt$structure))  # subset into subcortical structures only
g.means.lt.sc <- subset(g.means.lt.sc, grepl(pat.ana1, g.means.lt.sc$hemisphere))        # keep only "L" and "B" structures (remove right hemisphere)


g.means.lt.sc$mean1 <-  meanS.dklL[meanS.dklL$structure == "subcortical.brainstem", "estimate.dklL"]        # add means from the t-test
g.means.lt.sc$mean2 <-  rowMeans(g.means.lt.sc[5:85])                                                       # calculate means manually (sanity check)
# --> identical :) Phew.


# ....... for supplementary materials
# write_xlsx(g.means.lt.sc, "output/supp_mat_new/sup_mat_10b.xlsx")



# transform to long format
g.means.lt.scl <- data.frame(g.means.lt.sc %>% pivot_longer(names(g.means.lt.sc[,5:85]), 
                                                 names_to = "genes", values_to = "expression"))    # transform into long format, the way t.test likes it :)

g.means.lt.scl$label <- as.factor(g.means.lt.scl$label)



# ------------------ half violin/raincloud plots


# ..... vertical version

# plot specific data prep
# order by mean expression, ASCENDING
g.means.lt.scl.asc <- g.means.lt.scl %>% 
  arrange(mean2)

# manually set levels to make sure they are plotted in ascending order on the x-axis
order_h <- unique(g.means.lt.scl.asc$label)
g.means.lt.scl.asc$label <- factor(g.means.lt.scl.asc$label, levels = order_h)


# separate values for OXT, OXTR, and CD38
g.means.lt.sc.cd38 <- subset(g.means.lt.scl.asc, grepl("^CD38$", g.means.lt.scl.asc$genes)) 
g.means.lt.sc.oxt  <- subset(g.means.lt.scl.asc, grepl("^OXT$", g.means.lt.scl.asc$genes)) 
g.means.lt.sc.oxtr <- subset(g.means.lt.scl.asc, grepl("^OXTR$", g.means.lt.scl.asc$genes)) 
g.means.lt.sc.ot <- rbind(g.means.lt.sc.cd38, g.means.lt.sc.oxt, g.means.lt.sc.oxtr)         # combine
g.means.lt.sc.ot <- g.means.lt.sc.ot %>%                                                     # it should still be in the right order, but just to be on the very safe side
  arrange((mean2))

# ... and remove them from base jitter (they are kept for the violin plot to ensure appropriate curving of the plots)
pat.ot <- c("OXT|OXTR|CD38")
g.means.lt.scl.short <- subset(g.means.lt.scl.asc, !grepl(pat.ot, g.means.lt.scl.asc$genes))


# custom y axis labels
ytext <- c("Brain stem", "Accumbens area", "Pallidum",  "Thalamus proper",  "Amygdala", "Hippocampus", "Putamen", "Caudate")  


# plotting
g222 <- 
  ggplot() +          
  geom_density_ridges(g.means.lt.scl.asc, mapping = aes(x = expression, y = factor(label),              # use full dataset with all genes for violin plots
                                                        fill = mean2, colour = mean2), 
                      alpha = .15, scale = .8) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 23, face = "bold"),
        axis.title.y = element_text(size = 23, face = "bold"),
        axis.text.x  = element_text(size = 16),
        axis.text.y  = element_text(size = 16),
        legend.position = "none") +
  ylab("Subcortical regions\n") + xlab("\nmRNA intensity") +
  scale_colour_gradientn(colours = c("#9893DD", "#779ED7", "#5DC3C3", "#3CBC75", "#B8DE29")) +
  scale_fill_gradientn(colours = c("#9893DD", "#779ED7", "#5DC3C3", "#3CBC75", "#B8DE29")) 

g222 <- g222 + geom_point(g.means.lt.scl.short, mapping = aes(x = expression, y = factor(label),    # use dataset w/o OT genes for base jitter
                                                              fill = mean2, colour = mean2), 
                          position = position_jitter(height = .15), size = 4, alpha = .65)

g222 <- g222 + geom_point(g.means.lt.sc.ot,    mapping = aes(x = expression, y = factor(label)),    # use dataset with only OT genes for highlight
                          position = position_jitter(height = .15), size = 4, alpha=.6) 

g222 <- g222 + geom_vline(xintercept = mean(g.means.lt.scl.asc$mean2), 
                          linetype = "dashed", size = .75, colour = "gray60")                       # add LINE for mean subcortical expression points

g222 <- g222 + geom_point(g.means.lt.scl.asc, mapping = aes(x = mean2, y = factor(label)),          # add mean subcortical expression POINTS
                          position = position_nudge(y = .3), size = 7) 

g222 <- g222 + geom_vline(xintercept = mean(meanS.dklL$mean.all.struc),                             # add mean line for all regions, cortical and sub-cortical
                          linetype = "dashed", size = .75)      

g222 <- g222 + 
  annotate('text', x = 1.25, y = 1, label = '*', colour = "black", size = 11) +
  annotate('text', x = 1.25, y = 2, label = '*', colour = "black", size = 11) +
  annotate('text', x = 1.25, y = 3, label = '*', colour = "black", size = 11) +
  annotate('text', x = 1.25, y = 4, label = '*', colour = "black", size = 11) +
  annotate('text', x = 1.25, y = 6, label = '*', colour = "black", size = 11) 

g222 <- g222 + 
  scale_y_discrete(labels = ytext)

g222

# ggsave(filename = "output/abagen_analysis01_subcortical_vert.pdf", g222,
#        width = 10, height = 15, units = "in", device = 'pdf')


# arrange in a grid
g222_p127 <- plot_grid(g222, p127,
                       ncol       = 2, 
                       nrow       = 1,
                       rel_widths = c(1.75, 1),
                       labels     = c('A', 'B'),
                       label_size = 23) 

g222_p127

# ggsave(filename = "output/abagen_analysis01_ScC_combined02.pdf", g222_p127,
#        width = 15, height = 10.5, units = "in", device = 'pdf')





###############################################################################
## ANALYSIS 2
###############################################################################


pat.dk <- c("9861", "10021", "12876", "14380", "15496", "15697")

list.dk5 <- list()


for (i in pat.dk) {
  list.dk1[[i]] <- read_csv(paste0("output/abagen_dk/AHBAdk_220513_", i, ".csv"))
}

for (i in 1:length(pat.dk)) {
  list.dk5[[i]] <- data.frame(rowMeans(list.dk1[[i]][, -1]))                       # calculate row means
  list.dk5[[i]][["label"]] <- list.dk1[[i]][["label"]]                             # merge again with label column                              # change column names for readability
  names(list.dk5)[i] <- paste0("D.all.", pat.dk[i])                                                  
}    


list2env(list.dk5, envir = .GlobalEnv)   # unlists list with dataframes and pastes in the environment



DK.all.genes <- D.all.9861 %>%                                     # join all donor averages into one df
  left_join(D.all.10021, by = 'label') %>%
  left_join(D.all.12876, by = 'label') %>%
  left_join(D.all.14380, by = 'label') %>%
  left_join(D.all.15496, by = 'label') %>%
  left_join(D.all.15697, by = 'label') 


DK.all.genes <- cbind(DK.all.genes$label, DK.all.genes)
DK.all.genes <- DK.all.genes[, -3]

colnames(DK.all.genes) <- c("label", "meanALLDK1", "meanALLDK2", "meanALLDK3", "meanALLDK4", "meanALLDK5", "meanALLDK6")



DK.all.genes.i <- dplyr::right_join(dk.info, DK.all.genes, by = c("id" = "label"))                       # annotate with region names

DK.all.genes.iL <- data.frame(subset(DK.all.genes.i, grepl(pat.ana1, DK.all.genes.i$hemisphere)))        # select only left hemisphere and brain stem

# transform into long format, the way t.test likes it :)
DK.all.genes.iLl <- data.frame(DK.all.genes.iL %>% pivot_longer(c("meanALLDK1", "meanALLDK2", "meanALLDK3", "meanALLDK4", "meanALLDK5", "meanALLDK6"),        
                                                    names_to = "donors", values_to = "expression.all.genes"))


# OT dataframe: DK.all.ilL


# analysis

regions.dk2    <- unique(DK.all.genes.iLl$label)
pvals.dk2      <- numeric()
tStatistic.dk2 <- numeric()
cohensD.dk2    <- numeric()

estimateOT.dk2 <- numeric()
estimate.dk2   <- numeric()
sdOT.dk2       <- numeric()
sd.dk2         <- numeric()


for (i in 1:length(regions.dk2)){
  Temp.dk2 <- t.test(DK.all.ilL[DK.all.ilL$label == regions.dk2[i], "expression"], 
                     mu = mean(DK.all.genes.iLl[DK.all.genes.iLl$label == regions.dk2[i], "expression.all.genes"]), 
                     alternative = "two.sided")
#  CohD.dk2 <- ES.t.one(m  = Temp.dk2$estimate,
#                       sd = sd(DK.all.ilL[DK.all.ilL$label == regions.dk2[i], "expression"]),
#                       mu = mean(DK.all.genes.iLl[DK.all.genes.iLl$label == regions.dk2[i], "expression.all.genes"]),
#                       alternative = "two.sided")
  pvals.dk2      <- c(pvals.dk2, Temp.dk2$p.value)
  tStatistic.dk2 <- c(tStatistic.dk2,Temp.dk2$statistic)
#  cohensD.dk2    <- c(cohensD.dk2,CohD.dk2$d)
  
  estimateOT.dk2 <- c(estimateOT.dk2, Temp.dk2$estimate)
  estimate.dk2   <- c(estimate.dk2, Temp.dk2$null.value)
  sdOT.dk2       <- c(sdOT.dk2, sd(DK.all.ilL[DK.all.ilL$label == regions.dk2[i], "expression"]))
  sd.dk2         <- c(sd.dk2,   sd(DK.all.genes.iLl[DK.all.genes.iLl$label == regions.dk2[i], "expression.all.genes"]))
}

results.dk2 <- data.frame(Region = regions.dk2, P.unadjusted = pvals.dk2,
                          t.stat = tStatistic.dk2) #, CohD = cohensD.dk2)

results.dk2$P.fdr <- p.adjust(results.dk2$P.unadjusted, method = "fdr", n = length(results.dk2$P.unadjusted))

results.dk2$stat.sign <- ""
results.dk2$stat.sign[results.dk2$P.fdr <= 0.05]  <- "*"

results.dk2

# ....... for supplementary materials
# write_xlsx(g.means.lt.sc, "output/supp_mat_new/within_brain_test_stats.xlsx")



# combine test statistics and parameters and regions

test.stats.dk2 <- data.frame(OTmean = estimateOT.dk2, mean = estimate.dk2, OTSD = sdOT.dk2, SD = sd.dk2)
test.stats.dk2 <- cbind(results.dk2$Region, test.stats.dk2)

# --------------------------------------------------------------------> significant results for some regions

# sanity checks for means
# OT genes banksstats
check01 <- data.frame(c(0.5247907, 0.5346995 , 0.5320750 , 0.5285898 , .5457009 , 0.5185269))
mean(check01$c.0.5247907..0.5346995..0.532075..0.5285898..0.5457009..0.5185269.)
# [1] 0.5307305

# OT genes cuneus
check02 <- data.frame(c(0.5055762, 0.5315147 , 0.5269086 , 0.5050939 , 0.5077779 , 0.5350095))
mean(check02$c.0.5055762..0.5315147..0.5269086..0.5050939..0.5077779..0.5350095.)
# [1] 0.5186468
# everything ok :) 

# ALL genes banksstats
check03 <- data.frame(c(0.5109677, 0.5064424, 0.5000834, 0.5036631, 0.5081212, 0.5052236))
mean(check03$c.0.5109677..0.5064424..0.5000834..0.5036631..0.5081212..0.5052236.)

# ALL genes cuneus
check04 <- data.frame(c(0.5074279, 0.5114647, 0.4849419, 0.4962632, 0.4928239, 0.5195357))
mean(check04$c.0.5074279..0.5114647..0.4849419..0.4962632..0.4928239..0.5195357.)
# everything ok :)


# plotting

# ----- plot expression values in lolli pop graph

names(test.stats.dk2)[names(test.stats.dk2) == "results.dk2$Region"] <- "region"
test.stats.dk2$lobe <- "Frontal"
test.stats.dk2$lobe[c(1, 5, 6, 8, 14, 15, 29, 32, 33)] <- "Temporal"
test.stats.dk2$lobe[c(7, 21, 22, 24, 28, 30)] <- "Parietal"
test.stats.dk2$lobe[c(4, 10, 20)] <- "Occipital"
test.stats.dk2$lobe[34] <- "Insula"
test.stats.dk2$lobe[c(2, 9, 22, 25)] <- "Limbic"
test.stats.dk2$lobe[c(35:42)] <- "XSub-cortical"
test.stats.dk2$LongNames <- c("Banks of superior temporal S", "Caudal anterior cingulate C", "Caudal middle frontal G", "Cuneus C", 
                              "Entorhinal C",  "Fusiform G", "Inferior parietal C", "Inferior temporal G", "Isthmus cingulate C", 
                              "Lateral occipital C", "Lateral orbitofrontal C",  "Lingual G", "Medial orbitofrontal C", 
                              "Middle temporal G", "Parahippocampal G", "Paracentral L", "Pars opercularis", "Pars orbitalis", 
                              "Pars triangularis", "Pericalcarine C", "Postcentral G", "Posterior cingulate C", "Precentral G", 
                              "Precuneus C", "Rostral anterior cingulate C", "Rostral middle frontal G", "Superior frontal G", 
                              "Superior parietal C", "Superior temporal G", "Supramarginal G", "Frontal pole", "Temporal pole",
                              "Transverse temporal C", "Insula", "Thalamus proper", "Caudate", "Putamen", "Pallidum", "Accumbens area", 
                              "Hippocampus", "Amygdala", "Brain stem")


hlinOT <- mean(test.stats.dk2$OTmean) 
hlinAll <- mean(test.stats.dk2$mean)



test.stats.dk2.desc <- test.stats.dk2 %>%           # order descending by mean
  arrange(desc(mean))

test.stats.dk2.desc <- test.stats.dk2.desc %>%      # order ascending by lobe label
  arrange((lobe))


# manually set levels to make sure they are plotted in descending order on the x-axis
order.dk2 <- unique(test.stats.dk2.desc$LongNames)
test.stats.dk2.desc$LongNames <- factor(test.stats.dk2.desc$LongNames, levels = order.dk2)


# ...... for supplementary materials
write_xlsx(test.stats.dk2.desc, "output/supp_mat_new/within_brain_plot_stats.xlsx")


# Plot
lpp <- ggplot(test.stats.dk2.desc) + theme_classic() +
  geom_hline(yintercept = hlinAll, linetype = "dashed", color = "grey65", size = 1) +
  geom_hline(yintercept = hlinOT,  linetype = "dashed", color = "black",  size = 1) +
  
  geom_segment(aes(x = LongNames, xend = LongNames, y = OTmean, yend = mean), color = "grey") +
  geom_point(aes(x = LongNames, y = OTmean, color = lobe), size = 9, alpha = .85) +
  geom_point(aes(x = LongNames, y = mean, color = lobe), size = 9, alpha = .25) +
  viridis::scale_colour_viridis(option = "turbo", discrete = T, 
                                labels = c("Frontal", "Insula", "Limbic", "Occipital",
                                           "Parietal", "Temporal", "Sub-cortical"))

lpp <- lpp + ylab("mRNA intensity\n") + xlab("\nBrain regions") +
  theme(axis.title       = element_text(size = 25, face = "bold"),
        axis.text.x      = element_text(angle = 50, size = 20, vjust = 1, hjust = 1),
        #axis.text.x  = element_blank(),
        axis.ticks.x     = element_blank(),
        axis.text.y      = element_text(size = 15),
        legend.position  = c(0.25, 0.2),
        legend.direction = "horizontal",
        legend.text      = element_text(size = 22),
        legend.title     = element_blank()) 

lpp <- lpp + 
  annotate('text', x =  1, y = 0.5574277, label = '*', colour = "black", size = 14) +
  annotate('text', x =  2, y = 0.5557939, label = '*', colour = "black", size = 14) +
  annotate('text', x =  3, y = 0.5596114, label = '*', colour = "black", size = 14) +
  annotate('text', x =  4, y = 0.5572229, label = '*', colour = "black", size = 14) +
  annotate('text', x =  5, y = 0.5570432, label = '*', colour = "black", size = 14) +
  annotate('text', x =  6, y = 0.5503407, label = '*', colour = "black", size = 14) +
  annotate('text', x =  7, y = 0.5525331, label = '*', colour = "black", size = 14) +
  annotate('text', x =  8, y = 0.5508206, label = '*', colour = "black", size = 14) +
  annotate('text', x =  9, y = 0.5518614, label = '*', colour = "black", size = 14) +
  annotate('text', x = 10, y = 0.5569483, label = '*', colour = "black", size = 14) +
  annotate('text', x = 11, y = 0.5532204, label = '*', colour = "black", size = 14) +
  annotate('text', x = 12, y = 0.5315804, label = '*', colour = "black", size = 14) +
  annotate('text', x = 13, y = 0.5521486, label = '*', colour = "black", size = 14) +
  annotate('text', x = 14, y = 0.5553737, label = '*', colour = "black", size = 14) +
  annotate('text', x = 15, y = 0.5511881, label = '*', colour = "black", size = 14) +
  annotate('text', x = 16, y = 0.5407926, label = '*', colour = "black", size = 14) +
  annotate('text', x = 19, y = 0.539877,  label = '*', colour = "black", size = 14) +
  annotate('text', x = 20, y = 0.5286468, label = '*', colour = "black", size = 14) +
  annotate('text', x = 21, y = 0.5431154, label = '*', colour = "black", size = 14) +
  annotate('text', x = 22, y = 0.5443513, label = '*', colour = "black", size = 14) +
  annotate('text', x = 23, y = 0.543591,  label = '*', colour = "black", size = 14) +
  annotate('text', x = 24, y = 0.5507908, label = '*', colour = "black", size = 14) +
  annotate('text', x = 25, y = 0.5416664, label = '*', colour = "black", size = 14) +
  annotate('text', x = 26, y = 0.5517053, label = '*', colour = "black", size = 14) +
  annotate('text', x = 27, y = 0.5584234, label = '*', colour = "black", size = 14) +
  annotate('text', x = 28, y = 0.5545811, label = '*', colour = "black", size = 14) +
  annotate('text', x = 29, y = 0.5501405, label = '*', colour = "black", size = 14) +
  annotate('text', x = 30, y = 0.5496329, label = '*', colour = "black", size = 14) +
  annotate('text', x = 31, y = 0.5407305, label = '*', colour = "black", size = 14) +
  annotate('text', x = 32, y = 0.5458628, label = '*', colour = "black", size = 14) +
  annotate('text', x = 33, y = 0.5469003, label = '*', colour = "black", size = 14) +
  annotate('text', x = 34, y = 0.541717,  label = '*', colour = "black", size = 14) +
  annotate('text', x = 35, y = 0.458435,  label = '*', colour = "black", size = 14) +   
  annotate('text', x = 36, y = 0.4658219, label = '*', colour = "black", size = 14) +
  annotate('text', x = 37, y = 0.4570534, label = '*', colour = "black", size = 14) +
  annotate('text', x = 38, y = 0.5433302, label = '*', colour = "black", size = 14) +
  annotate('text', x = 42, y = 0.4453109, label = '*', colour = "black", size = 14) +
     scale_y_continuous(breaks = c(.45, .475, .5, .525, .55))

lpp


# ggsave(filename = "output/abagen_analysis02_cSc_03.pdf", lpp,
#        width = 20, height = 10, units = "in", device = 'pdf')


#### END OF SCRIPT ####


