#########################
# setting up working environment 
#########################

rm(list=ls()) # delete all objects in the workspace
gc(reset=T) # resest memory (especially useful when working with large data sets)

options(stringsAsFactors=F) # disables automatic conversion of char. strings into factors

Sys.setenv(LANG = "en")

#install.packages("installr")

# check for R version updates
library(installr)
updateR()

# PRO TIP (jk, it's a rookie tip): Ctr + left click on an object to view it in a new R tab

BASE="C:/01Dokumente/07Erwerbstaetigkeit/01PhD_UiO/projects/first_author/oxytocinHealth/RQ4/"      ### IMPORTANT: adjust and update wd

# Base dir for laptop
BASE="D:/eigene_dateien/02Dokumente/07Erwerbstaetigkeit/01PhD_UiO/projects/first_author/OTpathway_evolution/RProject/OT-bd-analyses_loc/"

#########################
# setting up working environment 
#########################


# Manuscript analysis 

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("ABAEnrichment")
#BiocManager::install("EnsDb.Hsapiens.v79")
#BiocManager::install("ggtree")


# general
library(tidyverse)
library(tibble)
library(readxl)
library(data.table)
library(readr)
library(tidyr)
library(powerAnalysis)
library(gage)          # not avail, check up
library(KEGGREST)      # not avail, check up
library(DataCombine)

# data visualization
library(ggplot2)
library(ggplotify)
library(ggpubr)
library(ggrepel)
library(cowplot)
library(gridExtra)
library(grid)
library(PupillometryR)
library(ggridges)
library(hrbrthemes)
library(ggforce)

# myTAI, phylostratigraphy, phylogenetics
library(ABAEnrichment) 
#library(myTAI)
require(ABAData)
library(EnsDb.Hsapiens.v79)
library(ggtree)
library(ape)

# brain expression enrichment
library(ggseg)
library(ggsegDefaultExtra)
#library(ggsegExtra)
# INFO: Install w/o vignettes, see documention: https://github.com/LCBC-UiO/ggseg/issues/36, https://github.com/LCBC-UiO/ggseg3d/issues/12


# INFO: If "permission denied" when installing packages: Remove affected package manually from the folder it is
# stored in. Open R as admin and start installation again.


#####

# https://cran.r-project.org/web/packages/myTAI/vignettes/Introduction.html

# https://github.com/drostlab/myTAI#tutorials

# https://drostlab.github.io/myTAI/articles/Introduction.html#scientific-introduction-performing-evolutionary-transcriptomics-with-r

# https://github.com/HajkD/published_phylomaps

# https://academic-oup-com.ezproxy.uio.no/gbe/article/8/6/1812/2574026

# https://bioconductor.org/packages/3.11/bioc/vignettes/ABAEnrichment/inst/doc/ABAEnrichment.html

# https://guangchuangyu.github.io/ggtree-book/short-introduction-to-r.html 

# https://github.com/drostlab/orthologr 

# https://academic.oup.com/bioinformatics/article/20/2/289/204981 

#####


## download structure ontology at: http://help.brain-map.org/display/api/Atlas+Drawings+and+Ontologies

all_struct=c('Allen:10163','Allen:10173',
             'Allen:10185','Allen:10194',
             'Allen:10209','Allen:10225',
             'Allen:10236','Allen:10243',
             'Allen:10252','Allen:10269',
             'Allen:10278','Allen:10294',
             'Allen:10333','Allen:10361',
             'Allen:10398','Allen:10657')
#########

# retrieve gene list from KEGG
# NOTE: genes belonging to a pathway are updated constantly, new genes are added on a regular basis, thus the
# list you retrieve might vary slightly from the list we initially retrieved

# first fetch information from KEGG to fetch label of the pathway we're looking for
kg.hsa = kegg.gsets("hsa") # retrieve gene set data from KEGG for Homo SApiens
kg.hsa.sets <- kg.hsa$kg.sets # fetch gene sets and search for "oxytocin signaling pathway"
# ----> pathway label: hsa04921

# get gene names
OTpthwy <- keggGet("hsa04921")[[1]]$GENE    
OTpthwy_nms_base <- OTpthwy[seq(2,length(OTpthwy),2)]
OTpthwy_nms <- gsub("\\;.*","",OTpthwy_nms_base)
OTpthwy_df <- data.frame(OTpthwy_nms)

#########

# this is our list
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

OTpthwyDF <- data.frame(OTpw.genes)

# Create OT pathway genes list




## =================================  Microsynteny analyses

# see R script "BLASTp_evaluation.R"





## ================================= Phylostratigraphy and myTAI analyses


###### load phylo map





#####

# left hemisphere hypothalamic structures

left_hyp=c('Allen:4598',  # anterior hypothalamic area, left
           'Allen:4596',  # lateral hypothalamic area, anterior region, left
           'Allen:4597',  # pallidohypothalamic nucleus, left
           'Allen:4573',  # paraventricular nucleus of the hypothalamus, left
           'Allen:4593',  # supraoptic nucleus, left
           
           'Allen:13005', # lateral hypothalamic area, mammillary region, left
           'Allen:4671',  # mammillary body, left
           'Allen:4675',  # lateral mammillary nucleus, left
           'Allen:4672',  # medial mammillary nucleus, left
           'Allen:4679',  # posterior hypothalamic area, left
           
           'Allen:4669',  # supramammillary nucleus, left
           'Allen:4667',  # tuberomammillary nucleus, left
           'Allen:4542',  # preoptic region, left
           'Allen:4644',  # arcuate nucleus of the hypothalamus, left
           'Allen:4641',  # dorsomedial hypothalamic nucleus, left
           
           'Allen:13008', # lateral hypothalamic area, tuberal region, left
           'Allen:10151', # lateral tuberal nucleus, left
           'Allen:4648',  # perifornical nucleus, left
           'Allen:4637'   # ventromedial hypothalamic nucleus, left
)


OTpthwy.res <- data.frame(read_excel(paste0(path, "output/BLASTp/OTpthwy_res.xlsx"), sheet = 1))

counts <- table(OTpthwy.res$Phylostratum)

barplot(counts, main = "Phylostratum levels",
        xlab = "Phylostratum")


# calculate relative frequencies for OT genes
rel_fOT <- table(OTpthwy.res$Phylostratum)/length(OTpthwy.res$Phylostratum)
rel_fOT_r <- round((rel_fOT*100), digits = 1)
abs_fOT <- table(OTpthwy.res$Phylostratum)


c1 <- c(12, 15, 1, 9, 0, 5, 8, 3, 1, 2, 5, 34, 11, 19, 14, 10, 2, 0, 1, 2 )                                           # absolute freqs  
c2 <- c(0.077922078, 0.097402597, 0.006493506, 0.058441558, 0, 0.032467532, 0.051948052,
        0.019480519, 0.006493506, 0.012987013, 0.032467532, 0.220779221, 0.071428571,
        0.123376623, 0.090909091, 0.064935065, 0.012987013, 0, 0.006493506, 0.012987013)   # relative freqs
c3 <- c(7.8, 9.7, 0.6, 5.8, 0, 3.2, 5.2, 1.9, 0.6, 1.3, 3.2, 22.1, 7.1, 12.3, 9.1, 6.5, 1.3, 0, 0.6, 1.3)                         # relative freqs, rounded



### ------- comparison with all genes not possible anymore since we ditched the old inaccurate method and
### ------- repeating the new method for all protein-coding genes is not possible....

#exp_all <- get_expression(structure_ids = c("Allen:10389"),      # AHBA structure ontology ID for the diencephalon. The other IDs for the hypothalamus did not work. 
#                          gene_ids      = HomoSapiens.PhyloMap2013$EnsemblGeneID, 
#                          dataset       = '5_stages')

#exp_all_df <- as.data.frame(do.call(rbind, exp_all)) 
#rownames(exp_all_df) <- NULL
#rownames(exp_all_df) <- c("Prenatal", "Infant", "Child", "Adolescent", "Adult")
#exp_all_t <- as.data.frame(t(exp_all_df))
#exp_all_t <- rownames_to_column(exp_all_t, var = "gene")

#names(exp_all_t)[names(exp_all_t) == "gene"] <- "Gene_ID"   
#exp_all_t <- exp_all_t %>% distinct(Gene_ID, .keep_all = TRUE) 

#mmALLE <- MatchMap(HomoSapiens.PhyloMap2013, exp_all_t)         # merge expression and phylogenetic data
# final available genes: 17027

# calculate relative frequencies for ALL genes
#rel_fALL <- data.frame(table(mmALLE$Phylostratum)/length(mmALLE$Phylostratum))
#rel_fALL_r <- data.frame(round((rel_fALL$Freq*100), digits = 1))
#abs_fALL <- data.frame(table(mmALLE$Phylostratum))

#d1 <- abs_fALL$Freq[1:12]               # absolute freqs
#d2 <- rel_fALL$Freq[1:12]               # relative freqs
#d3 <- rel_fALL_r$round..rel_fALL.Freq...100...digits...1.[1:12]  # relative freqs, rounded

# mention that the remaining percentages are distributed among PS 13-19

# combine ALL genes and OT genes freqs
labs <-c("1. Cellular organisms",
         "2. Eukaryota",
         "3. Opisthokonta",
         "4. Holozoa",
         "5. Metazoa",
         "6. Eumetazoa",
         "7. Bilateria",
         "8. Deuterostomia",
         "9. Chordata",
         "10. Craniata",
         "11. Vertebrata",
         "12. Gnathostomata",
         "13. Osteichthyes/Sarcopterygii",
         "14. Tetrapoda", 
         "15. Amniota",
         "16. Mammalia/Theria",
         "17. Eutheria/Boreoeutheria",
         "18. Euarchontaglires",
         "19. Primates",
         "20. Homo")

mm_absDF <- data.frame(labs = labs, counts_abs_OT = c1)

mm_relDF <- data.frame(labs = labs, counts_rel_OT = c2)

mm_relDF_round <- data.frame(labs = labs, counts_rel_OT_round = c3)

mm_relDF$labs <- factor(mm_relDF$labs, levels = mm_relDF$labs) 


# merge
mm_relDF2 <- cbind(mm_absDF, mm_relDF, mm_relDF_round)
mm_relDF2$labs_short <- c("PS1", "PS2", "PS3", "PS4", "PS5", "PS6", "PS7", "PS8", "PS9", "PS10", "PS11", 
                          "PS12", "PS13", "PS14", "PS15", "PS16", "PS17", "PS18", "PS19", "PS20")

mm_relDF2$labs_short <- factor(mm_relDF2$labs_short , levels = c("PS1", "PS2", "PS3", "PS4", "PS5", "PS6", "PS7", 
                                                                 "PS8", "PS9", "PS10", "PS11", "PS12", "PS13", 
                                                                 "PS14", "PS15", "PS16", "PS17", "PS18", "PS19", 
                                                                 "PS20"))

mm_relDF2 <- mm_relDF2[, -c(3,5)]




## ================================= FUMA analyses

# refer to "FUMA_results" script!






## ================================= dN/dS values from Dumas et al., 2021 --- GenEvo tool

# load file generated by https://genevo.pasteur.fr/, entered the OT pathway geneset as gene list. 
# Settings: Exact match, only 1-to-1 orthologs

dnds_primates_OT <- fread(paste0(BASE, "RProject/OT-bd-analyses_loc/Data/dNdS/genevo_2021-04-30_1to1.tsv"))
#View(dnds_primates_OT)

ids_OT <- data.frame(dnds_primates_OT$`Entrez ID`)

# generally, the majority of genes from the OT pathway geneset in HUMANS compared to a common primate ancestor have a dN/dS ratio < 1, indicating evolutionary constraint and high conservation.
# only 107 out of the 153 initial genes could be included in the analysis.

# For direct comparison against the neanderthals, we can use the NSS score. 


### ------ Box plot

dnds_primates_OT$Gene <- as.factor(dnds_primates_OT$Gene)

head(dnds_primates_OT)


## functions for whisker length adjustment

o <- function(x) {
  subset(x, x == max(x) | x == min(x))
}

f <- function(x) {
  r <- quantile(x, probs = c(0.00, 0.25, 0.5, 0.75, 1))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}


## phylogenetic tree generation

# >>>>>>>>>>>>>> INFO: check out https://www.r-graph-gallery.com/dendrogram.html and maybe prettify the script with this

# Guangchuang Yu. Using ggtree to visualize data on tree-like structures. Current Protocols in Bioinformatics, 2020, 69:e96. doi:10.1002/cpbi.96
# Guangchuang Yu, Tommy Tsan-Yuk Lam, Huachen Zhu, Yi Guan. Two methods for mapping and visualizing associated data on phylogeny using ggtree. Molecular Biology and Evolution 2018, 35(12):3041-3043. doi:10.1093/molbev/msy194
# Guangchuang Yu, David Smith, Huachen Zhu, Yi Guan, Tommy Tsan-Yuk Lam. ggtree: an R package for visualization and annotation of phylogenetic trees with their covariates and other associated data. Methods in Ecology and Evolution 2017, 8(1):28-36. doi:10.1111/2041-210X.12628


# old tree
#ptree_tx <- "(((C. jacchus, M. mulatta), P. abelii),(G. gorilla,(P. troglodytes,(H. denisova,(H. neanderthalensis, H. sapiens)))));"
#ptree <-read.tree(text=ptree_tx)


# generate tree w/o branch labels
tree_test_tx <- "(((,),),(,(,(,(,)))));"
ptree2 <-read.tree(text=tree_test_tx)

# plot tree
gg_ptr2 <- ggtree(ptree2, ladderize=FALSE) + geom_tiplab(align=TRUE) +
  geom_text(aes(label=c("C. jacchus", "M. mulatta", "P. abelii", "G. gorilla", "P. troglodytes",
                        "H. denisova", "H. neanderthalensis",  "H. sapiens",
                        "","", "", "","","", "")), 
            size=5, color="black", hjust = 0) + 
  scale_x_continuous(expand=expansion(1)) +
  xlab("") +
  theme(axis.title.x=element_text(margin = margin(t = 40, r = 0, b = 0, l = 0))) 

gg_ptr2


# Troubleshooting: If error message "Error in DataMask$new(.data, caller_env) : argument "caller_env" is missing, with no default", follow this thread: 
#                  https://github.com/YuLab-SMU/ggtree/issues/399. tl;dr: Update tidytree and set dplyr back to version 1.0.5



# ..... high quality
subset_hq <- dnds_primates_OT[, c(2, 6, 9, 12, 15, 18, 21, 24, 27)]

colnames(subset_hq) <- c("Gene",                                                                # rename col names
                         "H. sapiens", 
                         "H. neanderthalensis", 
                         "H. denisova", 
                         "P. troglodytes", 
                         "G. gorilla",
                         "P. abelii",
                         "M. mulatta",
                         "C. jacchus")

subset_hq_l <- subset_hq %>% pivot_longer(c("H. sapiens",                                       # transform into long format
                                            "H. neanderthalensis", 
                                            "H. denisova", 
                                            "P. troglodytes", 
                                            "G. gorilla",
                                            "P. abelii",
                                            "M. mulatta",
                                            "C. jacchus"),
                                          names_to = "Species", values_to = "value")

subset_hq_l$Species <- factor(subset_hq_l$Species , levels = c("H. sapiens",                     # make sure the order of the species is kept
                                                               "H. neanderthalensis", 
                                                               "H. denisova", 
                                                               "P. troglodytes", 
                                                               "G. gorilla",
                                                               "P. abelii",
                                                               "M. mulatta",
                                                               "C. jacchus"))

subset_hq_l$Gene <- as.factor(subset_hq_l$Gene)

p_hq <- ggplot(subset_hq_l, aes(x = value, y = Species, fill = Species, colour = Species)) + theme_minimal()

p_hq <- p_hq +  
  stat_summary(fun.data=f, geom="boxplot", alpha=0.45, show.legend = FALSE, width=.75) +
  stat_summary(fun = o, geom="point") 
p_hq <- p_hq + 
  scale_y_discrete(limits = rev(levels(subset_hq_l$Species)))
p_hq <- p_hq + 
  geom_point(aes(fill=Species), colour = "black", size = 4, shape = 21, alpha = .6, position = position_jitterdodge(jitter.width = 1))
p_hq <- p_hq + 
  xlab("dN/dS ratio") + xlim(0, 1.5)
p_hq <- p_hq + 
  theme(axis.title.x = element_text(size=25, face="bold", margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.text.x  = element_text(size = 15),
        axis.ticks.y = element_blank(),       
        axis.text.y  = element_blank(),         
        axis.title.y = element_blank()) 
p_hq <- p_hq + 
  viridis::scale_fill_viridis(option = 'viridis', discrete = T, begin = 0, end = .9) +
  viridis::scale_colour_viridis(option = 'viridis', discrete = T, begin = 0, end = .9)
p_hq <- p_hq + 
  geom_vline(xintercept = 1, linetype = "dashed", size = .75, color="gray50")
p_hq <- p_hq + geom_label_repel(data=(subset_hq_l %>% dplyr::filter(value > 1)), 
                                aes(label     = Gene), 
                                nudge_x       = 0,
                                nudge_y       = .4,
                                label.padding =.25,
                                colour        = "gray50",
                                fill          = alpha(c("white"),0.5),
                                size          = 4) 
p_hq <- p_hq + geom_label_repel(data=(subset_hq_l %>% dplyr::filter(Gene %in% c("OXTR"))), 
                                aes(label     = Gene), 
                                nudge_x       = .05, 
                                nudge_y       = .4,
                                label.padding = .25,
                                colour        = "steelblue4",
                                fill          = alpha(c("steelblue2"),.4),
                                size          = 4) 
p_hq <- p_hq + geom_label_repel(data=(subset_hq_l %>% dplyr::filter(Gene %in% c("OXT"))), 
                                aes(label     = Gene), 
                                nudge_x       = -.05, 
                                nudge_y       = .4,
                                label.padding = .25,
                                colour        = "steelblue4",
                                fill          = alpha(c("steelblue2"),.4),
                                size          = 4)
p_hq <- p_hq + geom_label_repel(data=(subset_hq_l %>% dplyr::filter(Gene %in% c("CD38"))), 
                                aes(label      = Gene), 
                                nudge_x        = 0, 
                                nudge_y        = .4,
                                label.padding  = .25,
                                colour         = "steelblue4",
                                fill           = alpha(c("steelblue2"),.4),
                                size           = 4)
p_hq <- p_hq + theme(legend.position = "none")

p_hq

#ggsave(filename = paste0(BASE, "RProject/OT-bd-analyses_loc/output/dnds_OTgeneset_hq.pdf"), p_hq,
#       width = 18, height = 8, units = "in", device='pdf')

# combine high quality barplot and phylogenetic tree
dnds_OT_tree_hq <- plot_grid(gg_ptr2, p_hq, rel_widths = c(1.1, 3)) 
dnds_OT_tree_hq


# save plot
ggsave(filename = paste0(BASE, "RProject/OT-bd-analyses_loc/output/dndsOT_tree_hq_phylo2.pdf"), dnds_OT_tree_hq,
       width = 20, height = 8, units = "in", device='pdf')



# ..... medium quality
subset_mq <- dnds_primates_OT[, c(2, 7, 10, 13, 16, 19, 22, 25, 28)]                            # subset with only medium quality columns, order of three last species changed
# to match phylo tree

colnames(subset_mq) <- c("Gene",                                                                # rename col names
                         "H. sapiens", 
                         "H. neanderthalensis", 
                         "H. denisova", 
                         "P. troglodytes", 
                         "G. gorilla",
                         "P. abelii",
                         "M. mulatta",
                         "C. jacchus")

subset_mq_l <- subset_mq %>% pivot_longer(c("H. sapiens",                                       # transform into long format
                                            "H. neanderthalensis", 
                                            "H. denisova", 
                                            "P. troglodytes", 
                                            "G. gorilla",
                                            "P. abelii",
                                            "M. mulatta",
                                            "C. jacchus"),
                                          names_to = "Species", values_to = "value")


subset_mq_l$Species <- factor(subset_mq_l$Species , levels = c("H. sapiens",                     # make sure the order of the species is kept
                                                               "H. neanderthalensis", 
                                                               "H. denisova", 
                                                               "P. troglodytes", 
                                                               "G. gorilla",
                                                               "P. abelii",
                                                               "M. mulatta",
                                                               "C. jacchus"))

subset_mq_l$Gene <- as.factor(subset_mq_l$Gene)

p_mq <- ggplot(subset_mq_l, aes(x = value, y = Species, fill = Species, colour = Species)) + theme_minimal()
p_mq <- p_mq +  
  stat_summary(fun.data=f, geom="boxplot", alpha=0.45, show.legend = FALSE, width=.75) +
  stat_summary(fun = o, geom="point") 
p_mq <- p_mq + 
  scale_y_discrete(limits = rev(levels(subset_mq_l$Species)))
p_mq <- p_mq + 
  geom_point(aes(fill=Species), colour = "black", size = 4, shape = 21, alpha = .6, position = position_jitterdodge(jitter.width = 1)) 
p_mq <- p_mq + 
  xlab("dN/dS ratio") + xlim(0, 1.5)
p_mq <- p_mq + 
  theme(axis.title.x = element_text(size=25, face="bold", margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.text.x  = element_text(size = 15),
        axis.ticks.y = element_blank(),       
        axis.text.y  = element_blank(),         
        axis.title.y = element_blank()) 
p_mq <- p_mq + 
  viridis::scale_fill_viridis(option = 'mako', discrete = T, begin = .3, end = .95) +
  viridis::scale_colour_viridis(option = 'mako', discrete = T, begin = .3, end = .95)
p_mq <- p_mq + 
  geom_vline(xintercept = 1, linetype = "dashed", size = .75, color="gray50") 
p_mq <- p_mq + geom_label_repel(data=(subset_mq_l %>% dplyr::filter(value > 1)), 
                                aes(label     = Gene), 
                                nudge_x       = 0,
                                nudge_y       = .4,
                                label.padding = .25,
                                colour        = "gray50",
                                fill          = alpha(c("white"),0.5),
                                size          = 4) 
p_mq <- p_mq + geom_label_repel(data=(subset_mq_l %>% dplyr::filter(Gene %in% c("OXTR"))), 
                                aes(label     = Gene), 
                                nudge_x       = .05, 
                                nudge_y       = .4,
                                label.padding = .25,
                                colour        = "#FF0564",
                                fill          = alpha(c("#FF0564"),.25),
                                size          = 4) 
p_mq <- p_mq + geom_label_repel(data=(subset_mq_l %>% dplyr::filter(Gene %in% c("OXT"))), 
                                aes(label     = Gene), 
                                nudge_x       = -.05, 
                                nudge_y       = .4,
                                label.padding = .25,
                                colour        = "#FF0564",
                                fill          = alpha(c("#FF0564"),.25),
                                size          = 4)
p_mq <- p_mq + geom_label_repel(data=(subset_mq_l %>% dplyr::filter(Gene %in% c("CD38"))), 
                                aes(label      = Gene), 
                                nudge_x        = 0, 
                                nudge_y        = .4,
                                label.padding  = .25,
                                colour         = "#FF0564",
                                fill           = alpha(c("#FF0564"),.25),
                                size           = 4)
p_mq <- p_mq + theme(legend.position = "none")

p_mq

#ggsave(filename = paste0(BASE, "RProject/OT-bd-analyses_loc/output/testy.pdf"), p_mq,
#       width = 18, height = 8, units = "in", device='pdf')

# combine med quality barplot and phylogenetic tree
dnds_OT_tree_mq <- plot_grid(gg_ptr2, p_mq, rel_widths = c(1.1, 3)) 
dnds_OT_tree_mq

# save plot
ggsave(filename = paste0(BASE, "RProject/OT-bd-analyses_loc/output/dndsOT_tree_mq_phylo6.pdf"), dnds_OT_tree_mq,
       width = 20, height = 8, units = "in", device='pdf')



###########################
###########################

# plot for ISPNE conference
p_mq <- ggplot(subset_mq_l, aes(x = value, y = Species, fill = Species, colour = Species)) + theme_minimal()
p_mq <- p_mq +  
  stat_summary(fun.data=f, geom="boxplot", alpha=0.45, show.legend = FALSE, width=.75) +
  stat_summary(fun = o, geom="point") 
p_mq <- p_mq + 
  scale_y_discrete(limits = rev(levels(subset_mq_l$Species)))
p_mq <- p_mq + 
  geom_point(aes(fill=Species), colour = "black", size = 5, shape = 21, alpha = .6, position = position_jitterdodge(jitter.width = 1)) 
p_mq <- p_mq + 
  xlab("dN/dS ratio") + xlim(0, 1.5)
p_mq <- p_mq + 
  theme(axis.title.x = element_text(size=25, face="bold", margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.text.x  = element_text(size = 15),
        axis.ticks.y = element_blank(),       
        axis.text.y  = element_blank(),         
        axis.title.y = element_blank()) 
p_mq <- p_mq + 
  viridis::scale_fill_viridis(option = 'magma', discrete = T, begin = .2, end = .85) +
  viridis::scale_colour_viridis(option = 'magma', discrete = T, begin = .2, end = .85)
p_mq <- p_mq + 
  geom_vline(xintercept = 1, linetype = "dashed", size = .75, color="gray50")

p_mq <- p_mq + theme(legend.position = "none")

p_mq


# save plot
ggsave(filename = paste0(BASE, "RProject/OT-bd-analyses_loc/output/dndsOT_tree_mq_phyloISPNE.pdf"), dnds_OT_tree_mq,
       width = 20, height = 8, units = "in", device='pdf')

###########################
###########################


# ..... low quality
subset_lq <- dnds_primates_OT[, c(2, 8, 11, 14, 17, 20, 23, 26, 29)]

colnames(subset_lq) <- c("Gene",                                                                # rename col names
                         "H. sapiens", 
                         "H. neanderthalensis", 
                         "H. denisova", 
                         "P. troglodytes", 
                         "G. gorilla",
                         "P. abelii",
                         "M. mulatta",
                         "C. jacchus")

subset_lq_l <- subset_lq %>% pivot_longer(c("H. sapiens",                                       # transform into long format
                                            "H. neanderthalensis", 
                                            "H. denisova", 
                                            "P. troglodytes", 
                                            "G. gorilla",
                                            "P. abelii",
                                            "M. mulatta",
                                            "C. jacchus"),
                                          names_to = "Species", values_to = "value")

subset_lq_l$Species <- factor(subset_lq_l$Species , levels = c("H. sapiens",                     # make sure the order of the species is kept
                                                               "H. neanderthalensis", 
                                                               "H. denisova", 
                                                               "P. troglodytes", 
                                                               "G. gorilla",
                                                               "P. abelii",
                                                               "M. mulatta",
                                                               "C. jacchus"))

subset_lq_l$Gene <- as.factor(subset_lq_l$Gene)

p_lq <- ggplot(subset_lq_l, aes(x = value, y = Species, fill = Species, colour = Species)) + theme_minimal()
p_lq <- p_lq +  
  stat_summary(fun.data=f, geom="boxplot", alpha=0.45, show.legend = FALSE, width=.75) +
  stat_summary(fun = o, geom="point") 
p_lq <- p_lq + 
  scale_y_discrete(limits = rev(levels(subset_lq_l$Species))) 
p_lq <- p_lq + 
  geom_point(aes(fill=Species), colour = "black", size = 4, shape = 21, alpha = .6, position = position_jitterdodge(jitter.width = 1)) 
p_lq <- p_lq + 
  xlab("dN/dS ratio") + xlim(0, 1.5) 
p_lq <- p_lq + 
  theme(axis.title.x = element_text(size=25, face="bold", margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.text.x  = element_text(size = 15),
        axis.ticks.y = element_blank(),       
        axis.text.y  = element_blank(),         
        axis.title.y = element_blank()) 
p_lq <- p_lq +
  viridis::scale_fill_viridis(option = 'cividis', discrete = T, end = .8) +
  viridis::scale_colour_viridis(option = 'cividis', discrete = T, end = .8)
p_lq <- p_lq + 
  geom_vline(xintercept = 1, linetype = "dashed", size = .75, color="gray50")
p_lq <- p_lq + geom_label_repel(data=(subset_lq_l %>% dplyr::filter(value > 1)), 
                                aes(label     = Gene), 
                                nudge_x       = 0,
                                nudge_y       = .4,
                                label.padding = .25,
                                colour        = "gray50",
                                fill          = alpha(c("white"),0.5),
                                size          = 4) 
p_lq <- p_lq + geom_label_repel(data=(subset_lq_l %>% dplyr::filter(Gene %in% c("OXTR"))), 
                                aes(label     = Gene), 
                                nudge_x       = .05, 
                                nudge_y       = .4,
                                label.padding = .25,
                                colour        = "steelblue4",
                                fill          = alpha(c("steelblue2"),.4),
                                size          = 4)  
p_lq <- p_lq + geom_label_repel(data=(subset_lq_l %>% dplyr::filter(Gene %in% c("OXT"))), 
                                aes(label     = Gene), 
                                nudge_x       = -.05, 
                                nudge_y       = .4,
                                label.padding = .25,
                                colour        = "steelblue4",
                                fill          = alpha(c("steelblue2"),.4),
                                size          = 4)
p_lq <- p_lq + geom_label_repel(data=(subset_lq_l %>% dplyr::filter(Gene %in% c("CD38"))), 
                                aes(label      = Gene), 
                                nudge_x        = 0, 
                                nudge_y        = .4,
                                label.padding  = .25,
                                colour         = "steelblue4",
                                fill           = alpha(c("steelblue2"),.4),
                                size           = 4)
p_lq <- p_lq + theme(legend.position = "none")

p_lq

#ggsave(filename = paste0(BASE, "RProject/OT-bd-analyses_loc/output/dnds_OTgeneset_lq.pdf"), p_lq,
#       width = 18, height = 8, units = "in", device='pdf')


# combine low quality barplot and phylogenetic tree
dnds_OT_tree_lq <- plot_grid(gg_ptr2, p_lq, rel_widths = c(1.1, 3)) 
dnds_OT_tree_lq


# save plot
ggsave(filename = paste0(BASE, "RProject/OT-bd-analyses_loc/output/dndsOT_tree_lq_phylo2.pdf"), dnds_OT_tree_lq,
       width = 20, height = 8, units = "in", device='pdf')


# ..... positioning OXT, OXTR, CD38 in background set (supplements)


data = c("genevo_2021-09-01_1to1000","genevo_2021-09-01_1001to2000", "genevo_2021-09-01_2001to3000", "genevo_2021-09-01_3001to4000", "genevo_2021-09-01_4001to5000", 
         "genevo_2021-09-01_5001to6000", "genevo_2021-09-01_6001to7000", "genevo_2021-09-01_7001to8000", "genevo_2021-09-02_8001to9000", "genevo_2021-09-02_9001to10000",
         "genevo_2021-09-02_10001to11000", "genevo_2021-09-02_11001to12000", "genevo_2021-09-02_12001to13000", "genevo_2021-09-02_13001to14000", "genevo_2021-09-02_14001to15000",
         "genevo_2021-09-02_15001to16000", "genevo_2021-09-02_16001to17000", "genevo_2021-09-02_17001to18000", "genevo_2021-09-02_18001to18445")

for (i in data) {
  df = fread(paste0(BASE, "RProject/OT-bd-analyses_loc/Data/dNdS/", i, ".tsv")) # read in dfs
  df = as.data.frame(df) #change colname
  assign(paste0("df_",i), df) #make new df
}


df_list1 = list(`df_genevo_2021-09-01_1to1000`, `df_genevo_2021-09-01_1001to2000`, `df_genevo_2021-09-01_2001to3000`, `df_genevo_2021-09-01_3001to4000`, `df_genevo_2021-09-01_4001to5000`,
                `df_genevo_2021-09-01_5001to6000`, `df_genevo_2021-09-01_6001to7000`, `df_genevo_2021-09-01_7001to8000`, `df_genevo_2021-09-02_8001to9000`, `df_genevo_2021-09-02_9001to10000`, 
                `df_genevo_2021-09-02_10001to11000`, `df_genevo_2021-09-02_11001to12000`, `df_genevo_2021-09-02_12001to13000`, `df_genevo_2021-09-02_13001to14000`, `df_genevo_2021-09-02_14001to15000`,
                `df_genevo_2021-09-02_15001to16000`, `df_genevo_2021-09-02_16001to17000`, `df_genevo_2021-09-02_17001to18000`, `df_genevo_2021-09-02_18001to18445`)
#empty_li <- list()
newdf <- NULL

for (i in df_list1) {
  df2 <- data.frame(i[, c(2, 7, 10, 13, 16, 19, 22, 25, 28)])
  #empty_li[[length(empty_li)+1]] <- df2                                  # if you want to create a list and from that list
  #assign(paste0("set", length(empty_li)), df2) # saves to environment    # 'export' data frames into the environment
  newdf <- rbind(df2, newdf)
}

colnames(newdf) <- c("Gene",                                                                # rename col names
                     "H. sapiens", 
                     "H. neanderthalensis", 
                     "H. denisova", 
                     "P. troglodytes", 
                     "G. gorilla",
                     "P. abelii",
                     "M. mulatta",
                     "C. jacchus")
newdf_l <- newdf %>% pivot_longer(c("H. sapiens",                                       # transform into long format
                                    "H. neanderthalensis", 
                                    "H. denisova", 
                                    "P. troglodytes", 
                                    "G. gorilla",
                                    "P. abelii",
                                    "M. mulatta",
                                    "C. jacchus"),
                                  names_to = "Species", values_to = "value")


newdf_l$Species <- factor(newdf_l$Species , levels = c("H. sapiens",                     # make sure the order of the species is kept
                                                       "H. neanderthalensis", 
                                                       "H. denisova", 
                                                       "P. troglodytes", 
                                                       "G. gorilla",
                                                       "P. abelii",
                                                       "M. mulatta",
                                                       "C. jacchus"))

newdf_l$Gene <- as.factor(newdf_l$Gene)

pALL_mq <- ggplot(newdf_l, aes(x = value, y = Species, fill = Species, colour = Species)) + theme_minimal()
pALL_mq <- pALL_mq +  
  stat_summary(fun.data=f, geom="boxplot", alpha=0.45, show.legend = FALSE, width=.75) +
  stat_summary(fun = o, geom="point") 
pALL_mq <- pALL_mq + 
  scale_y_discrete(limits = rev(levels(newdf_l$Species)))
pALL_mq <- pALL_mq + 
  geom_point(aes(fill=Species), colour = "black", size = 4, shape = 21, alpha = .6, position = position_jitterdodge(jitter.width = 1)) 
pALL_mq <- pALL_mq + 
  xlab("dN/dS ratio") #+ xlim(0, 1.5)
pALL_mq <- pALL_mq + 
  theme(axis.title.x = element_text(size=25, face="bold", margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.text.x  = element_text(size = 15),
        axis.ticks.y = element_blank(),       
        axis.text.y  = element_blank(),         
        axis.title.y = element_blank()) 
pALL_mq <- pALL_mq + 
  viridis::scale_fill_viridis(option = 'magma', discrete = T, begin = .3, end = .85) +
  viridis::scale_colour_viridis(option = 'magma', discrete = T, begin = .3, end = .85)
pALL_mq <- pALL_mq + 
  geom_vline(xintercept = 1, linetype = "dashed", size = .75, color="gray50") +
  scale_x_continuous(breaks = seq(0, 14, by = 1))
pALL_mq <- pALL_mq + geom_label_repel(data=(newdf_l %>% dplyr::filter(Gene %in% c("OXTR"))), 
                                      aes(label     = Gene), 
                                      nudge_x       = .3, 
                                      nudge_y       = .4,
                                      label.padding = .25,
                                      colour        = "steelblue4",
                                      fill          = alpha(c("steelblue2"),.4),
                                      size          = 4) 
pALL_mq <- pALL_mq + geom_label_repel(data=(newdf_l %>% dplyr::filter(Gene %in% c("OXT"))), 
                                      aes(label     = Gene), 
                                      nudge_x       = -.3, 
                                      nudge_y       = .4,
                                      label.padding = .25,
                                      colour        = "steelblue4",
                                      fill          = alpha(c("steelblue2"),.4),
                                      size          = 4)
pALL_mq <- pALL_mq + geom_label_repel(data=(newdf_l %>% dplyr::filter(Gene %in% c("CD38"))), 
                                      aes(label      = Gene), 
                                      nudge_x        = 0.6, 
                                      nudge_y        = .4,
                                      label.padding  = .25,
                                      colour         = "steelblue4",
                                      fill           = alpha(c("steelblue2"),.4),
                                      size           = 4)
pALL_mq <- pALL_mq + theme(legend.position = "none")

pALL_mq

#combine med quality barplot and phylogenetic tree
dnds_ALL_tree_mq <- plot_grid(gg_ptr2, pALL_mq, rel_widths = c(1.1, 3.4)) 
dnds_ALL_tree_mq


# save plot
ggsave(filename = paste0(BASE, "RProject/OT-bd-analyses_loc/output/dndsALL_tree_mq4.pdf"), dnds_ALL_tree_mq,
       width = 22, height = 8, units = "in", device='pdf')







## ================================= Abagen toolbox expression data pre-processing

# >>>>>>>>>>>>>>> python <<<<<<<<<<<<<<<<<<

## IMPORTANT INFO ##
# This part of the script is using the toolbox 'abagen' (https://github.com/rmarkello/abagen), which is written in PYTHON.
# Do not attempt to run in R, it will not work. Run in python shell in the terminal or in Jupyter Lab/Notebook.
# Use 'abagen_analysis.ipynb' for detailed information, pre-requisites, dependencies, etc.
# Most importantly, read https://abagen.readthedocs.io/en/stable/ carefully for optimal use. 


# NOT RUN {
# PREPARE
# download microarray data
import abagen    # initiate abagen
files = abagen.fetch_microarray(donors='all', verbose=0)

# ANALYSIS
# fetch atlas/parcellation
import abagen   # initiate abagen
atlasDK = abagen.fetch_desikan_killiany()

# get expression data and do the parcellation 
# --> IF data is already downloaded to local machine, specify the directory it is stored in with " data_dir='your/local/path/microarray/' " parameter
# --> otherwise, the function will download the data per default
AHBAdk, reportAHBAdk = abagen.get_expression_data(atlasDK['image'], atlasDK['info'],
                                                  norm_matched  = False,           # necessary for the missing parameter to work
                                                  missing       = 'interpolate',   # fill missing values
                                                  return_donors = True,            # returns all 6 donors separately instead of one matrix 1
                                                  donors        = 'all',           # returns all 6 donors separately instead of one matrix 2
                                                  return_report = True,            # returns report about what was done
                                                  data_dir      = 'C:\\Users\\daedh\\abagen-data\\microarray\\')  # directory where the data is stored

stabilityAMS = abagen.correct.keep_stable_genes(list(AHBAdk.values()), return_stability=True)[1]   # get differential stability values

# check if everything worked
print(reportAHBAdk)

print(AHBAdk)

print(stabilityAMS)


# export donor matrices
expression_main.to_csv('C:/01Dokumente/07Erwerbstaetigkeit/01PhD_UiO/projects/first_author/oxytocinHealth/RQ4/data/processed/expression_main_mis.csv', 
                       index = False, header=True)


# tbc....


# export DS values
import pandas as pd 
pd.DataFrame(stabilityAMS).to_csv('D:\\eigene_dateien\\02Dokumente\\07Erwerbstaetigkeit\\01PhD_UiO\\projects\\first_author\\oxytocinHealth\\RQ4\\analyses\\differential_stability\\DS_values_abagen.csv',
                                  index = True)
# }

#############################
##### end python script #####
#############################


## INFO: must run script "4sixDonors_classicExpAnalysisAccordingToDan_WithinBrainRegionsTtests_DK2013.R" first


# Differential stability analysis

# load data generated with abagen
DSabagen <- read_csv(paste0(BASE, "/analyses/differential_stability/DS_values_abagen.csv"))   # each row corresponds to a gene column from the expression matrix, i.e., rows are genes here
gene_names <- data.frame(colnames(dkD1))     # get gene names
gene_namess <- data.frame(gene_names[-1,])
DSabagen_anno <- cbind(DSabagen, gene_namess)   # annotate w/ gene names
15633/2                                          # estimate top 50% cut-off value
# cut off: 7817, top 50% are N = 7816

# order by DS value
DSabagen_anno_desc <- DSabagen_anno %>% 
  arrange(desc(`0`))

DStop50 <- DSabagen_anno_desc[1:7816,]             # extract top 50%
DSlower50 <- DSabagen_anno_desc[7817:15633,]        # extract lowest 50%


DSOTtop50 <- dplyr::inner_join(OTpthwyDF, DStop50, by=c("OTpw.genes"="gene_names..1..."))       # match with *all* OT genes
DSOTlower50 <- dplyr::inner_join(OTpthwyDF, DSlower50, by=c("OTpw.genes"="gene_names..1..."))    # match with OT *all* genes

OTgenes.out <- dplyr::anti_join(OTpthwyDF, DSabagen_anno_desc, by=c("OTpw.genes"="gene_names..1..."))    # check whether some genes were not available

# no data available for 18 genes
# out of 136 genes, 105 were among the top 50% genes with highest DS values, 31 were among the genes with lower DS values
# out of the genes for which data was available (n = 136), 77.21% have high DS values and are thus consistently and reliabley expressed across brain regions independet of donor characteristics. 


# prepare supp mat excel files
DSabagen_anno_desc$OTpathwaygene <- "n"
DSabagen_anno_desc[DSabagen_anno_desc$gene_names..1... %in% OTpw.genes, "OTpathwaygene"] <- "*YES"
DSabagen_anno_desc$division <- "bottom50"
DSabagen_anno_desc$division[1:7816] <- "top50"

names(DSabagen_anno_desc)[names(DSabagen_anno_desc) == "X1"] <- "index"
names(DSabagen_anno_desc)[names(DSabagen_anno_desc) == "0"] <- "DS"
names(DSabagen_anno_desc)[names(DSabagen_anno_desc) == "gene_names..1..."] <- "GeneLabel"


library("writexl")
write_xlsx(DSabagen_anno_desc, paste0(BASE, "supp_mat07.xlsx"))




## ================================= Visualization of AHBA cortical data with ggseg











## ================================= Trial TDI and dN/dS analysis with human VS macaque

## IMPORTANT: this code is not working yet ##

# UPDATE1: The code is either not working or it is taking incredibly long to run.
# UPDATE2: Inquire Constantina about this?

## helpful links and websites ##
# https://github.com/drostlab/orthologr 
# https://drostlab.github.io/orthologr/articles/Install.html 
# https://drostlab.github.io/orthologr/reference/index.html 
# https://github.com/ropensci/biomartr 
# https://drostlab.github.io/orthologr/articles/dNdS_estimation.html 

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install()



### -- probably these only work for Linux/Unix/OS -- ###

# install orthologr from GitHub
#devtools::install_github("HajkD/metablastr")

# install orthologr from GitHub
#devtools::install_github("HajkD/orthologr")

### ------------------------------------------------ ###



# Install package dependencies
BiocManager::install(c("GenomicRanges", "GenomicFeatures", "Rsamtools", "rtracklayer"))

BiocManager::install("Biostrings")

# issues with downloading "rtracklayer" and "Biostrings": remove the package from the library directory manually and run R as administrator, then install



### ------------- for WINDOWS users ------------ ###



install.packages("rtools")   # it says this doesn't work or is not available for this version of R
# which is honestly a straight up lie. I installed R tools manually 
# via: https://cran.r-project.org/bin/windows/Rtools/. Why is R lying to me?

install.packages("devtools")

devtools::install_github("HajkD/orthologr", build_vignettes = TRUE, dependencies = TRUE, force=T)

library("orthologr", lib.loc = "C:/Program Files/R/R-3.1.1/library")

### -------------------------------------------- ###


BiocManager::install("biomaRt")
install.packages("biomartr", dependencies = TRUE)

install.packages("readr")

# R tools is strictly necessary here apparently
# also install (important)....
# BLAST: ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ !!OPEN IN FIREFOX!! (info:https://www.ncbi.nlm.nih.gov/books/NBK52637/)
# C++: https://visualstudio.microsoft.com/de/vs/features/cplusplus/
# Strawberry Perl: https://www.perl.org/get.html#win32 
# KaKs calculator: https://sourceforge.net/projects/kakscalculator2/ 

# installing older version of "readr" because the download function for the genome data seems to have a problem with paths and readr as of version 1.4.0??
require(devtools)
install_version("readr", version = "1.3.1", repos = "http://cran.us.r-project.org")

library(biomaRt)
library(orthologr)



## analyses ##

# download reference genomes 

# get h. sapiens genome
hSapiens <- biomartr::getCDS(organism = "Homo sapiens", path = getwd())

# get m. mulatta genome 
mMulatta <- biomartr::getCDS(organism = "Macaca mulatta", path = getwd())


# Detect orthologous genes between a query species and a subject species
# and compute the synonymous versus non-synonymous substitution rates (dN/dS)
# following this paradigm:
# 1) reciprocal best hit for orthology inference (RBH)
# 2) Needleman-Wunsch for pairwise amino acid alignments
# 3) pal2nal for codon alignments
# 4) Comeron for dNdS estimation
# 5) multi-core processing 'comp_cores = 1'
sap_vs_mac <- 
  dNdS(query_file      = hSapiens,
       subject_file    = mMulatta,
       delete_corrupt_cds = TRUE,    # coding sequences that cannot be divided by 3 (triplets) will be removed
       ortho_detection = "RBH",      # perform BLAST best reciprocal hit orthology inference
       aa_aln_type     = "pairwise", # perform pairwise global alignments of AA seqs 
       aa_aln_tool     = "NW",       # using Needleman-Wunsch
       codon_aln_tool  = "pal2nal",  # perform codon alignments using the tool Pal2Nal
       dnds_est.method = "Comeron",  # use Comeron's method for dN/dS inference
       comp_cores      = 1,
       kaks_calc_path  = "C:/Program Files/KaKs_Calculator1.2_Windows.tar/Windows/"    # specifying path to KaKs/dNdS calculator
       # blast_path      = "C:/Program Files/NCBI/blast-2.11.0+/bin/"                  # specifying path to BLAST program    
  )

# Fehler in value[[3L]](cond) : 
# File blastresult_Homo_sapiens_cds_from_genomic_refseq.fna.gz.csv could not be read correctly. 
# Please check the correct path to blastresult_Homo_sapiens_cds_from_genomic_refseq.fna.gz.csv or whether BLAST did write the resulting hit table correctly.

# Re-installed all necessary packages as admin, installed older version of readr() as required, specified paths to KaKs calculator and BLAST

# UPDATE: still not working...? WHY?

# UPDATE: removed the specified path to BLAST and now it's running for hours again....

View(sap_vs_mac)

readr::write_excel_csv(sap_vs_mac, "hs_vs_mm_dNdS.csv")

