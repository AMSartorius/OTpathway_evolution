######################################
##### SET UP WORKING ENVIRONMENT #####
######################################

# rm(list = ls()) # delete all objects in the workspace, USE WITH CAUTION
# gc(reset = T)   # resest memory (especially useful when working with large data sets), USE WITH CAUTION

options(stringsAsFactors = F) # disables automatic conversion of char. strings into factors
Sys.setenv(LANG = "en")


# check for R version updates
# library(installr)  # install with: install.packages("installr")
# updateR()


### .......... loading packages
# general
library(tidyverse)     # for everything and pipes
library(tibble)
library(readxl)        # read in/load .xlsx files
library(data.table)
library(readr)         # read in/load .csv files
# library(powerAnalysis) # ISSUE: not avail for R version 4.2.1.
library(gage)          # install with BiocManager (Bioconductor)
library(KEGGREST)      # install with BiocManager (Bioconductor)
library(DataCombine)   # for "InsertRow" function
library(writexl)       # export as .xlsx files 

# data visualization
library(ggplot2)
library(ggplotify)
library(viridis)       # for colour palette 
library(ggpubr)
library(ggrepel)       # for manually labeling data points in plots
library(cowplot)
library(PupillometryR) # for violin plots??
library(ggridges)
library(hrbrthemes)
library(ggforce)

# phylogenetics, dn/ds ratio
library(ggtree)        # install with BiostcManager (Bioconductor)
library(ape)

# brain expression enrichment
library(ggseg)
# INFO: Install w/o vignettes, see documention: https://github.com/LCBC-UiO/ggseg/issues/36, https://github.com/LCBC-UiO/ggseg3d/issues/12


# INFO: If "permission denied" when installing packages: Remove affected package manually from the folder it is
# stored in. Open R as admin and start installation again.


###################################

# GAGE: Luo, Weijun, Friedman, Michael, Shedden, Kerby, Hankenson, Kurt, Woolf, Peter (2009). 
#       “GAGE: generally applicable gene set enrichment for pathway analysis.” BMC Bioinformatics, 10, 161.

# https://academic-oup-com.ezproxy.uio.no/gbe/article/8/6/1812/2574026

# https://bioconductor.org/packages/3.11/bioc/vignettes/ABAEnrichment/inst/doc/ABAEnrichment.html

# https://guangchuangyu.github.io/ggtree-book/short-introduction-to-r.html 

# https://github.com/drostlab/orthologr 

# https://academic.oup.com/bioinformatics/article/20/2/289/204981 

#####

# Retrieve OT pathway gene list from KEGG
# NOTE: genes belonging to a pathway are updated constantly, new genes are added on a regular basis, thus the
# list you retrieve might vary slightly from the list we initially retrieved

# first fetch information from KEGG to fetch label of the pathway we're looking for
kg.hsa = kegg.gsets("hsa")    # retrieve gene set data from KEGG for Homo SApiens
kg.hsa.sets <- kg.hsa$kg.sets # fetch gene sets and search for "oxytocin signaling pathway"
# ----> pathway label: hsa04921

# get gene names
OTpthwy <- keggGet("hsa04921")[[1]]$GENE    
OTpthwy_nms_base <- OTpthwy[seq(2,length(OTpthwy),2)]
OTpthwy_nms <- gsub("\\;.*","",OTpthwy_nms_base)
OTpthwy_df <- data.frame(OTpthwy_nms)

#####

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

OTpthwyDF <- data.frame(OTpw.genes) # Create OT pathway genes list


## ================================= BLASTp/microsynteny results visualisation prep

OTpthwy.res <- data.frame(read_excel("output/BLASTp/OTpthwy_res02_vertebrate_threshold.xlsx", sheet = 1))
   

# inspect results, gene count per phylostratum
counts <- OTpthwy.res %>% 
  count(Phylostratum)                                                 

# template for adobe illustrator figure
bp <- ggplot(counts, aes(Phylostratum, n)) + 
  geom_bar(aes(fill = n), stat = "identity", colour = "black") + 
  scale_fill_gradient() +
  theme_classic()
bp

# save for adobe illustrator
# ggsave(filename = "output/phylotree_template.pdf", bp,
#        width = 18, height = 6, units = "in", device = 'pdf')
  
OTpthwy.res %>% 
  count(Phylostratum)

OTpthwy.res %>% 
  count(age.3.levels.vert)

# continue with manual visualization in Adobe Illustrator, project: phylogeny_circular_02_v2, 
# current path (fixed): D:\eigene_dateien\02Dokumente\07Erwerbstaetigkeit\01PhD_UiO\projects\first_author\OTpathway_evolution\images_figures\circular_phylogeny\phylogeny_circular_02_v2.ai



## ================================= dN/dS values from Dumas et al., 2021 --- GenEvo tool

# load file generated by https://genevo.pasteur.fr/ 
# If you want to manually replicate the output from the online tool, enter the entire OT pathway geneset as gene list and then export it
# from the tool. Settings: Exact match, only 1-to-1 orthologs

dnds_primates_OT <- fread("data/raw/dNdS/genevo_2021-04-30_1to1.tsv")
#View(dnds_primates_OT)

ids_OT <- data.frame(dnds_primates_OT$`Entrez ID`)


### ------ Visualization: Box plot

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

# >>>>>>>>>>>>>> INFO: check out https://www.r-graph-gallery.com/dendrogram.html if you want to prettify the script

# Guangchuang Yu. Using ggtree to visualize data on tree-like structures. Current Protocols in Bioinformatics, 2020, 69:e96. doi:10.1002/cpbi.96
# Guangchuang Yu, Tommy Tsan-Yuk Lam, Huachen Zhu, Yi Guan. Two methods for mapping and visualizing associated data on phylogeny using ggtree. Molecular Biology and Evolution 2018, 35(12):3041-3043. doi:10.1093/molbev/msy194
# Guangchuang Yu, David Smith, Huachen Zhu, Yi Guan, Tommy Tsan-Yuk Lam. ggtree: an R package for visualization and annotation of phylogenetic trees with their covariates and other associated data. Methods in Ecology and Evolution 2017, 8(1):28-36. doi:10.1111/2041-210X.12628



# generate tree w/o branch labels
tree_test_tx <- "(((,),),(,(,(,(,)))));"
ptree2 <- read.tree(text=tree_test_tx)

# plot tree
gg_ptr2 <- ggtree(ptree2, ladderize = FALSE) + geom_tiplab(align = TRUE) +
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



# ..... Supplementary materials: high quality
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

# ggsave(filename = "output/supp_mat_new/dnds_OTgeneset_hq.pdf", p_hq,
#       width = 18, height = 8, units = "in", device='pdf')

# combine high quality barplot and phylogenetic tree
dnds_OT_tree_hq <- plot_grid(gg_ptr2, p_hq, rel_widths = c(1.1, 3)) 
dnds_OT_tree_hq


# save plot
ggsave(filename = "output/supp_mat_new/dnds_OTgeneset_hq.pdf", dnds_OT_tree_hq,
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

#ggsave(filename = "output/dnds_OTgeneset_hq.pdf", p_mq,
#       width = 18, height = 8, units = "in", device='pdf')

# combine med quality barplot and phylogenetic tree
dnds_OT_tree_mq <- plot_grid(gg_ptr2, p_mq, rel_widths = c(1.1, 3)) 
dnds_OT_tree_mq


# create directory
res_dir <- "output/figures/"
dir.create(res_dir)

# save plot
ggsave(filename = "output/figures/dndsOT_tree_mq_phylo6.pdf", dnds_OT_tree_mq,
       width = 20, height = 8, units = "in", device = 'pdf')




# ..... Supplementary materials: low quality
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

#ggsave(filename = output/dnds_OTgeneset_lq.pdf", p_lq,
#       width = 18, height = 8, units = "in", device='pdf')


# combine low quality barplot and phylogenetic tree
dnds_OT_tree_lq <- plot_grid(gg_ptr2, p_lq, rel_widths = c(1.1, 3)) 
dnds_OT_tree_lq


# save plot
ggsave(filename = "output/supp_mat_new/dndsOT_tree_lq_phylo2.pdf", dnds_OT_tree_lq,
      width = 20, height = 8, units = "in", device = 'pdf')




# ..... Supplementary materials: positioning OXT, OXTR, CD38 in background set (supplements)


data = c("genevo_2021-09-01_1to1000","genevo_2021-09-01_1001to2000", "genevo_2021-09-01_2001to3000", "genevo_2021-09-01_3001to4000", "genevo_2021-09-01_4001to5000", 
         "genevo_2021-09-01_5001to6000", "genevo_2021-09-01_6001to7000", "genevo_2021-09-01_7001to8000", "genevo_2021-09-02_8001to9000", "genevo_2021-09-02_9001to10000",
         "genevo_2021-09-02_10001to11000", "genevo_2021-09-02_11001to12000", "genevo_2021-09-02_12001to13000", "genevo_2021-09-02_13001to14000", "genevo_2021-09-02_14001to15000",
         "genevo_2021-09-02_15001to16000", "genevo_2021-09-02_16001to17000", "genevo_2021-09-02_17001to18000", "genevo_2021-09-02_18001to18445")

for (i in data) {
  df = fread(paste0("data/raw/dNdS/", i, ".tsv")) # read in dfs
  df = as.data.frame(df)                      # change colname
  assign(paste0("df_", i), df)                # make new df
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
ggsave(filename = "output/supp_mat_new/dndsALL_tree_mq4.pdf", dnds_ALL_tree_mq,
       width = 22, height = 8, units = "in", device='pdf')




######################################
########### END OF SCRIPT ############
######################################

