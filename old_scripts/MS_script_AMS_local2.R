#########################
# setting up working environment 
#########################

rm(list=ls()) # delete all objects in the workspace
gc(reset=T) # resest memory (especially useful when working with large data sets)

options(stringsAsFactors=F) # disables automatic conversion of char. strings into factors

Sys.setenv(LANG = "en")

#install.packages("installr")

# check for R version updates
# library(installr)
# updateR()

# PRO TIP (jk, it's a rookie tip): Ctr + left click on an object to view it in a new R tab

BASE="C:/01Dokumente/07Erwerbstaetigkeit/01PhD_UiO/projects/first_author/oxytocinHealth/RQ4/"      ### IMPORTANT: adjust and update wd

# Base dir for laptop
BASE="D:/eigene_dateien/02Dokumente/07Erwerbstaetigkeit/01PhD_UiO/projects/first_author/oxytocinHealth/RQ4/"

#########################
# setting up working environment 
#########################


# Manuscript analysis 

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("ABAEnrichment")

#BiocManager::install("EnsDb.Hsapiens.v79")


# general
library(tidyverse)
library(tibble)
library(readxl)
library(data.table)
library(readr)
library(tidyr)
library(powerAnalysis)
library(gage)
library(KEGGREST)

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
library(myTAI)
require(ABAData)
library(EnsDb.Hsapiens.v79)
library(ggtree)
library(ape)

# brain expression enrichment
library(ggseg)
library(ggsegExtra)
# INFO: Install w/o vignettes, see documention: https://github.com/LCBC-UiO/ggseg/issues/36, https://github.com/LCBC-UiO/ggseg3d/issues/12


# INFO: If "permission denied" when installing packages: Remove affected package manually from the folder it is
# stored in. Open R as admin and start install again.


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
#######

# retrieve gene list from KEGG
# NOTE: genes belonging to a pathway are updated constantly, new genes are added on a regular basis, thus the
# list you retrieve might vary slightly from the list we initially retrieved

# first fetch information from KEGG to fetch label of the pathway we're looking for
kg.hsa=kegg.gsets("hsa") # retrieve gene set data from KEGG for Homo SApiens
kg.hsa.sets <- kg.hsa$kg.sets # fetch gene sets and search for "oxytocin signaling pathway"
# ----> pathway label: hsa04921

# get gene names
OTpthwy <- keggGet("hsa04921")[[1]]$GENE    
OTpthwy_nms_base <- OTpthwy[seq(2,length(OTpthwy),2)]
OTpthwy_nms <- gsub("\\;.*","",OTpthwy_nms_base)
OTpthwy_df <- data.frame(OTpthwy_nms)

# Create OT pathway genes list

kop = c('GUCY1A2', 'GUCY1A1', 'CACNA2D2', 'CACNA2D1', 'CACNA2D4',  'CACNA2D3','PRKAB1', 'PRKAB2', 'PPP1CA', 'PPP1CB',
        'MYL6', 'MYL9', 'CALM2','CALM3','RAF1', 'CALM1','PRKAA2', 'ROCK1','ADCY1','ROCK2',
        'CALML6', 'ADCY5', 'CALML3', 'ADCY4','CALML4', 'ADCY3','ADCY2','CALML5', 'ADCY9','EGFR',   
        'PRKAA1', 'ADCY8','ADCY7','ADCY6','PPP1CC', 'NRAS', 'ACTB',  'CCND1', 'MAP2K5', 'PRKCG',
        'MAP2K2', 'PLA2G4E','PLA2G4F','MAP2K1', 'PLA2G4C','PLA2G4D','PLA2G4A','PLA2G4B','PRKCB','EEF2',
        'PRKCA',  'CACNB2', 'RCAN1',  'CACNB3', 'CACNB4', 'GNAQ', 'GNAS','CACNB1','OXTR', 'GUCY1B1',
        'CACNA1C','ACTG1','CACNA1D','ELK1',  'CACNA1F', 'MYL6B','MYLK', 'TRPM2','PPP3R1', 'PRKACG',
        'EEF2K','PPP3R2', 'RGS2', 'NPPA', 'CACNA1S','PRKACB', 'PRKACA', 'KCNJ5','RHOA', 'KCNJ6',
        'PPP1R12A', 'PPP1R12B', 'JUN', 'KCNJ12', 'KCNJ9', 'KCNJ14', 'NFATC4', 'NFATC3', 'NFATC2', 'NFATC1',
        'GNAO1', 'CAMK4', 'CAMK1', 'PPP1R12C', 'NPR1', 'NPR2', 'PIK3R5', 'ITPR2', 'ITPR3', 'OXT',
        'ITPR1', 'PIK3R6',  'PPP3CA', 'PPP3CB', 'PPP3CC', 'CD38', 'CAMK1D', 'MEF2C', 'NOS3', 'FOS',
        'PLCB4', 'JMJD7-PLA2G4B', 'KRAS', 'PLCB2', 'CAMK1G', 'PLCB3', 'PLCB1', 'MYLK2', 'RYR2', 'RYR3',
        'CAMK2D', 'MYLK3', 'SRC', 'RYR1', 'CDKN1A', 'CAMK2A', 'CAMK2B', 'PRKAG2', 'PRKAG3', 'PRKAG1',
        'GNAI1', 'GNAI2', 'CAMKK2', 'MYLK4', 'PTGS2', 'GNAI3', 'PIK3CG', 'CACNG7', 'CACNG8', 'MAPK3',
        'CACNG1', 'MAPK1', 'HRAS', 'CACNG2', 'CACNG3', 'MAPK7', 'CAMK2G', 'CACNG4', 'KCNJ2', 'KCNJ3',
        'CACNG5', 'CACNG6', 'KCNJ4')

######

# download the Phylostratigraphic Map of Homo sapiens from Tomislav Domazet-Loso and Diethard Tautz, 2008

# 1. download supplementary data zip file from: https://academic.oup.com/mbe/article/25/12/2699/1114618?searchresult=1#supplementary-data
# 2. extract "mbe-08-0522-File008_msn214.xls" file and save as "MBE_2008_Homo_Sapiens_PhyloMap.xls" in your BASE directory



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

hyp_ex <- get_expression(structure_ids=left_hyp,      # fetch expression data
                         gene_ids=OTpthwy_df$OTpthwy_nms, 
                         dataset='adult')

hyp_ex <- as.data.frame(hyp_ex)                       # prepare data: transpose, create column with gene labels

hyp_ex_t <- as.data.frame(t(hyp_ex))

hyp_ex_t <- rownames_to_column(hyp_ex_t, var = "gene")

hyp_ex_t$gene

ens_gene_hyp <- ensembldb::select(EnsDb.Hsapiens.v79,                   # fetch ensemble gene names
                                  keys= hyp_ex_t$gene, 
                                  keytype = "SYMBOL", 
                                  columns = c("SYMBOL","GENEID"))

ens_gene_hyp <- ens_gene_hyp %>% distinct(SYMBOL, .keep_all = TRUE)     # remove duplicates (?)
names(ens_gene_hyp)[names(ens_gene_hyp) == "SYMBOL"] <- "gene"
hyp_ex_t_f <- dplyr::full_join(ens_gene_hyp, hyp_ex_t, by = "gene")     # join with expression data
hyp_ex_t_f <- dplyr::select(hyp_ex_t_f, -gene)

names(hyp_ex_t_f)[names(hyp_ex_t_f) == "GENEID"] <- "Gene_ID"

HomoSapiens.PhyloMap <- read_excel(paste0(BASE, "RProject/OT-bd-analyses_loc/MBE_2008_Homo_Sapiens_PhyloMap.xls"),  # load phylo map
                                   sheet = 1, skip = 1)

HomoSapiens.PhyloMap <- HomoSapiens.PhyloMap[ , c("Phylostratum","Gene_ID")]

counts <- table(HomoSapiens.PhyloMap$Phylostratum)

barplot(counts, main="Phylostratum levels",
        xlab="Phylostratum")

mm_hyp <- MatchMap(HomoSapiens.PhyloMap, hyp_ex_t_f)
#write.csv(mm_hyp, paste0(BASE, "mm_hyp.csv"))


# ---- get supp. table with genes and phylostrata
mm_hyp_reduced <- mm_hyp[,1:2] 
mm_hyp_reduced$GeneID <- toupper(mm_hyp_reduced$GeneID)

supp_mat01 <- dplyr::left_join(mm_hyp_reduced, ens_gene_hyp, by=c("GeneID"="GENEID"))
names(supp_mat01)[names(supp_mat01) == "gene"] <- "GeneLabel"
supp_mat01 <- supp_mat01 %>%       # order ascending by lobe label
  arrange((Phylostratum))

col1_sm <- NA
col2_sm <- NA
supp_mat01.2 <- dplyr::anti_join(OTpthwy_df, supp_mat01, by= c("OTpthwy_nms"="GeneLabel"))
supp_mat01.22 <- cbind(col1_sm, col2_sm, supp_mat01.2)
names(supp_mat01.22)[names(supp_mat01.22) == "col1_sm"] <- "Phylostratum"
names(supp_mat01.22)[names(supp_mat01.22) == "col2_sm"] <- "GeneID"
names(supp_mat01.22)[names(supp_mat01.22) == "OTpthwy_nms"] <- "GeneLabel"

supp_mat01_fin <- rbind(supp_mat01, supp_mat01.22)

library("writexl")
write_xlsx(supp_mat01_fin, paste0(BASE, "supp_mat01_c.xlsx"))
# ---- ---- ----


# calculate relative frequencies for OT genes
rel_fOT <- table(mm_hyp$Phylostratum)/length(mm_hyp$Phylostratum)
rel_fOT_r <- round((rel_fOT*100), digits = 1)
abs_fOT <- table(mm_hyp$Phylostratum)


c1 <- c(94, 19, 5, 4, 1, 7, 4, 0, 0, 1, 0, 3)                                           # absolute freqs  
c2 <- c(0.681159420, 0.137681159, 0.036231884, 0.028985507, 0.007246377, 0.050724638,   # relative freqs
        0.028985507, 0, 0, 0.007246377, 0, 0.021739130) 
c3 <- c(68.1, 13.8, 3.6, 2.9, 0.7, 5.1, 2.9, 0, 0, 0.7, 0, 2.2)                         # relative freqs, rounded
            







#### ALL genes
# get expression data from all ontogenetic stages for each brain region for as many genes as possible (ens gene ids in the phylo map as reference)
exp_all <- get_expression(structure_ids=left_hyp,      # fetch expression data
                          gene_ids=HomoSapiens.PhyloMap$Gene_ID, # I'm using the gene ids for genes that are avail in the phylo map because these are the ones we need and will use in the end. No need to call genes for which there is no phylogenetic data avail anyways.
                          dataset='adult')

exp_all <- as.data.frame(exp_all)                       # prepare data: transpose, create column with gene labels

exp_all_t <- as.data.frame(t(exp_all))

exp_all_t <- rownames_to_column(exp_all_t, var = "gene")

exp_all_t <- exp_all_t %>% distinct(gene, .keep_all = T)

# merge expression and phylogenetic data
mmALLE <- MatchMap(HomoSapiens.PhyloMap, exp_all_t)
# final available genes: 14957




# calculate relative frequencies for ALL genes
rel_fALL <- data.frame(table(mmALLE$Phylostratum)/length(mmALLE$Phylostratum))
rel_fALL_r <- data.frame(round((rel_fALL$Freq*100), digits = 1))
abs_fALL <- data.frame(table(mmALLE$Phylostratum))


d1 <- abs_fALL$Freq[1:12]               # absolute freqs
d2 <- rel_fALL$Freq[1:12]               # relative freqs
d3 <- rel_fALL_r$round..rel_fALL.Freq...100...digits...1.[1:12]  # relative freqs, rounded


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
         "9. Olfactores",
         "10. Chordata",
         "11. Craniata",
         "12. Euteleostomi")

mm_absDF <- data.frame(labs=labs, counts_rel_ALL=d1, counts_rel_OT=c1)

mm_relDF <- data.frame(labs=labs, counts_rel_ALL=d2, counts_rel_OT=c2)

mm_relDF$labs <- factor(mm_relDF$labs, levels = mm_relDF$labs) 


# transform into long format
mm_relDF_L <- mm_relDF %>% pivot_longer(c("counts_rel_ALL", "counts_rel_OT"),        
                                        names_to = "counts", values_to = "values")
mm_relDF_L <- data.frame(mm_relDF_L)



# df for within plot bar labs with rounded ("short") values
rel_barlabs <- data.frame(labs=labs, barALL=d3, barOT=c3)

# transform into long format
rel_barlabs_L <- rel_barlabs %>% pivot_longer(c("barALL", "barOT"),        
                                              names_to = "barCounts", values_to = "vals_short")
rel_barlabs_L <- data.frame(rel_barlabs_L)



# merge
mm_relDF2 <- cbind(mm_relDF_L, rel_barlabs_L)
mm_relDF2 <- mm_relDF2[, -c(4,5)]
mm_relDF2$labs_short <- c("PS1", "PS1", "PS2", "PS2", "PS3", "PS3" ,"PS4", "PS4", "PS5", "PS5", "PS6", "PS6", 
                          "PS7", "PS7", "PS8", "PS8", "PS9", "PS9", "PS10", "PS10", "PS11", "PS11", "PS12", "PS12")

mm_relDF2$labs_short <- factor(mm_relDF2$labs_short , levels = c("PS1", "PS2", "PS3", "PS4", "PS5", "PS6",
                                                                 "PS7", "PS8", "PS9", "PS10", "PS11", "PS12"))

# OLD VERSION -- plotting Overleaf version

#mmoverleaf <- ggplot(mm_relDF2, aes(values, labs_short, fill = counts, width=.75)) +
#  theme_classic() +
#  geom_col(aes(fill=counts), position="dodge") +
#  geom_text(aes(values, labs_short, colour=counts, label = vals_short), 
#            position = position_dodge(width = .8), hjust = -.2, size=4) 

#mmoverleaf <- mmoverleaf +
#  theme(axis.text = element_text(size=15),
#        axis.title.y = element_text(size=17.5, face="bold"),
#        axis.title.x = element_text(size=15, face="bold"),
#        axis.ticks.y = element_blank(),
#        legend.position = "none") +
#  xlim(0,.75) + 
#  ylab("Phylostratum \n (phylogenetic stage)\n") + xlab("Relative gene frequencies \n per phylostratum ") +
#  scale_fill_manual(values=c("#1675AA", "#93CFF1")) +
#  scale_colour_manual(values=c("black", "gray50"))

#mmoverleaf

#ggsave(filename = paste0(BASE, "RProject/OT-bd-analyses_loc/output/genes_per_PS_relV2.1.pdf"), mmoverleaf,
#       width = 4.2, height = 4.7, units = "in", device='pdf')





# NEW VERSION --


tt_tx <- "(((,(,(,(,(,(,(,(,(,(,(,(,(())))))))))))))));"
tt <-read.tree(text=tt_tx)
#ggtree(tt) + geom_tiplab(align=TRUE)

# plot tree
gg_tt <- ggtree(tt, ladderize=FALSE) + geom_tiplab(align=TRUE) +
  geom_text(aes(label=c("  PS1", "  PS2", "  PS3", "  PS4", "  PS5", "  PS6", "  PS7", "  PS8", "  PS9", "  PS10", "  PS11", "  PS12", "", "", "",
                        "","", "", "","", "", "", "", "", "", "", "", "", "")), 
            size=14, color="black", hjust = 0) +
  xlab("") +
  theme(axis.title.x=element_text(margin = margin(t = 115, r = 0, b = 0, l = 0))) 

#gg_tt

mmoverleaf_t <- ggplot(mm_relDF2, aes(values, labs_short, fill = counts, width=.75)) +
  theme_bw() +
  geom_col(aes(fill=counts), position="dodge") +
  geom_text(aes(values, labs_short, colour=counts, label = vals_short), 
            position = position_dodge(width = .8), hjust = -.2, size=13) 

mmoverleaf_t <- mmoverleaf_t +
  theme(axis.text.x = element_text(size=37),
        axis.title.x = element_text(size=42, face="bold"),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(color = 'black'))+
  xlim(0,.75) + 
  xlab("Relative gene frequencies \n per phylostratum ") +
  scale_fill_manual(values=c("#1675AA", "#93CFF1")) +
  scale_colour_manual(values=c("black", "gray50"))


#mmoverleaf_t

blank <- ggplot() + theme_void()

right_side <- plot_grid(blank, mmoverleaf_t, 
                        ncol=1, rel_heights  = c(.14, 2)) 

right_side

testxb <- plot_grid(gg_tt, right_side, 
                    ncol=2, rel_widths = c(2.5, 1)) 
testxb

#ggsave(filename = paste0(BASE, "RProject/OT-bd-analyses_loc/output/fig1template.jpg"), testxb,
#       width = 46, height = 20, units = "in", device='jpg')


# ---------> proceed with bioRender





#########################
# Function to run relevant myTAI functions
#########################


tai_dev <-
  function(region) 
  {
    exp <- get_expression(structure_ids=c(region),
                          gene_ids=kop, 
                          dataset='5_stages') # Get expression values
    
    
    exp_df <- as.data.frame(do.call(rbind, exp)) 
    rownames(exp_df) <- NULL
    rownames(exp_df) <- c("Prenatal", "Infant", "Child", "Adolescent", "Adult")
    exp_t <- as.data.frame(t(exp_df))
    exp_t <- rownames_to_column(exp_t, var = "gene")                         
    
    ens_gene_dev <- ensembldb::select(EnsDb.Hsapiens.v79, 
                                      keys= exp_t$gene, 
                                      keytype = "SYMBOL", 
                                      columns = c("SYMBOL","GENEID"))
    
    ens_gene_dev <- ens_gene_dev %>% distinct(SYMBOL, .keep_all = TRUE)
    
    names(ens_gene_dev)[names(ens_gene_dev) == "SYMBOL"] <- "gene"
    
    ex_t_f <- dplyr::full_join(ens_gene_dev, exp_t, by = "gene")
    
    ex_t_f <- dplyr::select(ex_t_f, -gene)
    
    names(ex_t_f)[names(ex_t_f) == "GENEID"] <- "Gene_ID"
    
    HomoSapiens.PhyloMap <- read_excel(paste0(BASE, "RProject/OT-bd-analyses_loc/MBE_2008_Homo_Sapiens_PhyloMap.xls"), 
                                       sheet = 1, skip = 1)
    
    HomoSapiens.PhyloMap <- HomoSapiens.PhyloMap[ , c("Phylostratum","Gene_ID")]
    
    mm <- MatchMap(HomoSapiens.PhyloMap, ex_t_f)
    
    
    ps <- PlotSignature( ExpressionSet = mm,
                         measure       = "TAI", 
                         TestStatistic = "FlatLineTest",
                         p.value       = FALSE,
                         xlab          = "Developmental stage", 
                         ylab          = "TAI" )
    
    ps_dat <- ps$data                                               # Plotted
    
    pc <- PlotContribution(mm,                                      # Plotted
                           legendName = "PS")
    
    re <-PlotRE(mm,                                                 # Weirdly and funnily, if the "Groups" parameter is NOT defined correctly, as is the case here, one
                Groups = list(c(1:7)),                              # gets access to the raw data. If, however, two groups are defined as its supposed to be, the
                legendName = "PS")                                  # raw data is not provided in the object "re".
    
    pm <- PlotMeans(mm,                                             # Plotted as average. Weirdly, even though the group has been limited to PS 1-7, it plots all PS?
                    Groups = list(c(1:7)),          
                    legendName = "PS",
                    adjust.range = TRUE)
    
    #pb <- PlotBarRE(ExpressionSet = mm,                             
    #                Groups        = list(c(1:3), c(4:7, 10, 12)),  # PS 1-3 = "old", PS 4-12 = "young" (myTAI documentation p. 32/33)?
    #                p.adjust.method = "fdr")
    
    #pCE <- PlotCategoryExpr(mm,                                    # not sure whether to keep this
    #                        legendName = "PS",
    #                        test.stat = TRUE,
    #                        type = "category-centered",
    #                        distr.type = "dotplot",
    #                        log.expr = F)
    
    
    ################# I am not 100% sure anymore why I wrote this part of the function anymore but I believe it deals with the pMatric function of
    # myTAI that Computes Partial TAI or TDI Values and then it transforms it into long format and 
    # adds a column marking phylostrata as either "old" or "young"... I guess it prepares for data visualization.
    
    pmm <- data.frame(pMatrix(mm))
    pmm <- rownames_to_column(pmm)
    pmm <- pmm %>% 
      mutate(index = 1:138)
    pmm_t <- dplyr::right_join(mm, pmm, by=c("GeneID"="rowname"))
    pmm_t <- pmm_t[,-c(3:7)]
    colnames(pmm_t) <- c("Phylostratum", "GeneID", "Prenatal", "Infant", "Child", "Adolescent", "Adult", "index")
    pmm_t$`Phylostratum class` <- pmm_t$Phylostratum
    
    # add column with phylostratum classes, 1:3 = old, 4:12 = young
    pattern_young = c("4|5|6|7|10|12")
    pmm_t$`Phylostratum class` <- str_replace_all(pmm_t$`Phylostratum class`, pattern_young, "young")
    
    pattern_old = c("1|2|3")
    pmm_t$`Phylostratum class` <- str_replace_all(pmm_t$`Phylostratum class`, pattern_old, "old")
    
        # transform into long format
    pmm_l <- pmm_t %>% pivot_longer(c("Prenatal", "Infant", "Child", "Adolescent", "Adult"),        # transform into long format, the way t.test likes it :)
                                    names_to = "Stages", values_to = "PartialTAI")
    
    #################
    
    
    #gs_oxtr <- SelectGeneSet(ExpressionSet = mm,                    # plotted
    #                     gene.set = "ENSG00000180914",
    #                     use.only.map = F)
    
    # gs_oxt <- SelectGeneSet(ExpressionSet = mm,                     # plotted
    #                         gene.set = "ENSG00000101405",
    #                         use.only.map = F)
    
    # gs_cd38 <- SelectGeneSet(ExpressionSet = mm,                    # plotted
    #                         gene.set = "ENSG00000004468",
    #                         use.only.map = F)
    
    value <- list(
      mm = mm,
      ps = ps,
      ps_dat = ps_dat,
      pc = pc,
      re = re,
      pm = pm,
      pmm_l = pmm_l
      #pCE = pCE,
      #pb = pb,
      #gs_oxtr = gs_oxtr,
      #gs_oxt = gs_oxt,
      #gs_cd38 = gs_cd38
    ) # Create a list of output objects
    attr(value, "class") <- "tai_dev"
    value
  }



## get main data

get_name(10163)
M1C_tai <- tai_dev("Allen:10163") 

get_name(10173)
DFC_tai <- tai_dev("Allen:10173")

get_name(10185)
VFC_tai <- tai_dev("Allen:10185")

get_name(10194)
OFC_tai <- tai_dev("Allen:10194")

get_name(10209)
S1C_tai <-tai_dev("Allen:10209")

get_name(10225)
IPC_tai <- tai_dev("Allen:10225")

get_name(10236)
A1C_tai <- tai_dev("Allen:10236")

get_name(10243)
STC_tai <- tai_dev("Allen:10243")

get_name(10252)
ITC_tai <- tai_dev("Allen:10252")

get_name(10269)
V1C_tai <- tai_dev("Allen:10269")

get_name(10278)
MFC_tai <- tai_dev("Allen:10278")

get_name(10294)
HIP_tai <- tai_dev("Allen:10294")

get_name(10333)
STR_tai <-tai_dev("Allen:10333")

get_name(10361)
AMY_tai <- tai_dev("Allen:10361")

get_name(10398)
MD_tai <- tai_dev("Allen:10398")

get_name(10657)
CBC_tai <- tai_dev("Allen:10657")




## -------------- PLOTTING

## ================================= Combine individual TAI plots 

# Prepare variables

M1Cdat <- M1C_tai$ps_dat
colnames(M1Cdat)[2] <- "M1C"

DFCdat <- DFC_tai$ps_dat
colnames(DFCdat)[2] <- "DFC"

VFCdat <- VFC_tai$ps_dat
colnames(VFCdat)[2] <- "VFC"

OFCdat <- OFC_tai$ps_dat
colnames(OFCdat)[2] <- "OFC"

S1Cdat <- S1C_tai$ps_dat
colnames(S1Cdat)[2] <- "S1C"

IPCdat <- IPC_tai$ps_dat
colnames(IPCdat)[2] <- "IPC"

A1Cdat <- A1C_tai$ps_dat
colnames(A1Cdat)[2] <- "A1C"

STCdat <- STC_tai$ps_dat
colnames(STCdat)[2] <- "STC"

ITCdat <- ITC_tai$ps_dat
colnames(ITCdat)[2] <- "ITC"

V1Cdat <- V1C_tai$ps_dat
colnames(V1Cdat)[2] <- "V1C"

MFCdat <- MFC_tai$ps_dat
colnames(MFCdat)[2] <- "MFC"

HIPdat <- HIP_tai$ps_dat
colnames(HIPdat)[2] <- "HIP"

STRdat <- STR_tai$ps_dat
colnames(STRdat)[2] <- "STR"

AMYdat <- AMY_tai$ps_dat
colnames(AMYdat)[2] <- "AMY"

MDdat <- MD_tai$ps_dat
colnames(MDdat)[2] <- "MD"

CBCdat <- CBC_tai$ps_dat
colnames(CBCdat)[2] <- "CBC"

jt <- M1Cdat %>%
  left_join(DFCdat, by='Stage') %>%
  left_join(VFCdat, by='Stage') %>%
  left_join(OFCdat, by='Stage') %>%
  left_join(S1Cdat, by='Stage') %>%
  left_join(IPCdat, by='Stage') %>%
  left_join(A1Cdat, by='Stage') %>%
  left_join(STCdat, by='Stage') %>%
  left_join(ITCdat, by='Stage') %>%
  left_join(V1Cdat, by='Stage') %>%
  left_join(MFCdat, by='Stage') %>%
  left_join(HIPdat, by='Stage') %>%
  left_join(STRdat, by='Stage') %>%
  left_join(AMYdat, by='Stage') %>%
  left_join(MDdat,  by='Stage') %>%
  left_join(CBCdat, by='Stage')


################
min(jt[,2:17])
max(jt[,2:17])
################


jt$Stage <- as.factor(jt$Stage)

level_order <- c("Prenatal", "Infant",
                 "Child", "Adolescent", 
                 "Adult")

mean_TAIs_stages <- rowMeans(jt[, -1])


jtl <- gather(jt, 
              brain_region, 
              TAI, 
              M1C:CBC, 
              factor_key=TRUE)


# !!!!!!!!!!!!!!!!!!!!!
# !!!!!!!!!!!!!!!!!!!!!


# IMPORTANT:
# execute script "TAIcom-PlotSignature_facetZOOM_allProteinCoding_vs_OT" now and then come back to this script for 
# further analyses. Skip commented part below (do not pass Go. Do not collect $200) directly proceed with "Plot Phylostratum Enrichment" section of *this* script.


# !!!!!!!!!!!!!!!!!!!!!
# !!!!!!!!!!!!!!!!!!!!!



##################################

#jtl$cat <- "Majority pattern"

#jtl[jtl$brain_region %in% c('V1C', 'MD'),]$cat <- "Childhood peak"

#jtl[jtl$brain_region %in% c('M1C', 'IPC', 'AMY'),]$cat <- "Adolescence drop"

#jtl[jtl$brain_region %in% c('STR'),]$cat <- "Combined"

#jtl$Stage <- factor(jtl$Stage , levels = c("Prenatal", "Infant",
#                                           "Child", "Adolescent", 
#                                           "Adult"))

# single plot
#p <- ggplot(jtl, aes(x= Stage, y=TAI, group=brain_region)) +
#  geom_line(aes(color=cat), size = 1) +
#  theme_classic() 
#p <- p + 
#  annotate('text', x = 0.75, y = 1.692, label = 'STR', colour= "#FF0564") +   # #05FFBE
#  annotate('text', x = 0.75, y = 1.628, label = 'V1C', colour= "#81C800") +
#  annotate('text', x = 0.75, y = 1.644, label = 'MD', colour= "#81C800") +
#  annotate('text', x = 0.75, y = 1.72, label = 'AMY', colour= "#05FFBE") +
#  annotate('text', x = 0.75, y = 1.67, label = 'M1C', colour= "#05FFBE") +
#  annotate('text', x = 0.75, y = 1.636, label = 'IPC', colour= "#05FFBE") +
#  scale_color_manual(values = c("#05FFBE", "#81C800", "#FF0564", "#3D477B")) 
#p <- p + labs(color = "TAI trajectories") + 
#        theme(legend.position = c(.77, .75),
#              axis.title.x = element_text(face="bold", size=16),
#              axis.text.x  = element_text(angle=40, hjust=1, vjust = 1, size=12),
#              axis.title.y = element_text(face="bold", size=16),
#              axis.text.y = element_text(size=12))

#p <- p + xlab("Ontogenetic stage") + ylab("Transcriptome Age Index\n")

#p

#ggsave(filename = paste0(BASE, "RProject/OT-bd-analyses_loc/output/TAI_combined.pdf"), p,
#       width = 6, height = 6, units = "in", device='pdf')


# prepare for combined plot (different sizes but generally the same plot) and overleaf
#pcom <- ggplot(jtl, aes(x= Stage, y=TAI, group=brain_region)) +
#  geom_line(aes(color=cat), size = 1.5) +
#  theme_classic() 
#pcom <- pcom + 
#  annotate('text', x = 0.75, y = 1.692, label = 'STR', colour= "#FF0564", size=7) +
#  annotate('text', x = 0.75, y = 1.627, label = 'V1C', colour= "#81C800", size=7) +
#  annotate('text', x = 0.75, y = 1.645, label = 'MD', colour= "#81C800", size=7) +
#  annotate('text', x = 0.75, y = 1.72, label = 'AMY', colour= "#FAA100", size=7) +
#  annotate('text', x = 0.75, y = 1.67, label = 'M1C', colour= "#FAA100", size=7) +
#  annotate('text', x = 0.75, y = 1.636, label = 'IPC', colour= "#FAA100", size=7) +
#  scale_color_manual(values = c("#FAA100", "#81C800", "#FF0564", "#104C7E"))           
#pcom <- pcom + labs(color = "TAI trajectories") + 
#  theme(legend.position = c(.77, .75),
#       legend.title = element_text(size=20, face="bold"),
#        legend.text = element_text(size=20),
#        axis.title.x = element_text(face="bold", size=23),
#        axis.text.x  = element_text(angle=20, vjust=.5, hjust=0.5, size=21),
#        axis.title.y = element_text(face="bold", size=22),
#        axis.text.y = element_text(size=21))

#pcom <- pcom + xlab("\nOntogenetic stage") + ylab("Transcriptome Age Index\n")

#pcom

#ggsave(filename = paste0(BASE, "RProject/OT-bd-analyses_loc/output/TAIcombined_overleaf.pdf"), pcom,
#       width = 10, height = 7, units = "in", device='pdf')

#mmFreqP_2 <- mmFreqP +
#  theme(axis.title.x = element_text(margin = margin(t=0, r=0, b=12, l=0)),
#        legend.position = "none") 

#mmFreqP_2
# combine pcom with pplot from above
#TAIcom_genesperPS <- plot_grid(mmFreqP_2, pcom,
#                               ncol = 2, rel_widths = c(.7, 1),
#                               labels = c('B', 'C'),
#                               label_size = 28) # 12 * 20

#TAIcom_genesperPS

#ggsave(filename = paste0(BASE, "RProject/OT-bd-analyses_loc/output/TAIcom_genesperPS_rel8.pdf"), TAIcom_genesperPS,
#       width = 20, height = 7.5, units = "in", device='pdf')

##################################





## ================================= RELATIVE expression levels of PS age categories per ontogenetic stage per brain region


M1Cre <- M1C_tai$re
M1Cre_dat <- M1Cre$data

DFCre <- DFC_tai$re
DFCre_dat <- DFCre$data

VFCre <- VFC_tai$re
VFCre_dat <- VFCre$data

OFCre <- OFC_tai$re
OFCre_dat <- OFCre$data

S1Cre <- S1C_tai$re
S1Cre_dat <- S1Cre$data

IPCre <- IPC_tai$re
IPCre_dat <- IPCre$data

A1Cre <- A1C_tai$re
A1Cre_dat <- A1Cre$data

STCre <- STC_tai$re
STCre_dat <- STCre$data

ITCre <- ITC_tai$re
ITCre_dat <- ITCre$data

V1Cre <- V1C_tai$re
V1Cre_dat <- V1Cre$data

MFCre <- MFC_tai$re
MFCre_dat <- MFCre$data

HIPre <- HIP_tai$re
HIPre_dat <- HIPre$data

STRre <- STR_tai$re
STRre_dat <- STRre$data

AMYre <- AMY_tai$re
AMYre_dat <- AMYre$data

MDre <- MD_tai$re
MDre_dat <- MDre$data

CBCre <- CBC_tai$re
CBCre_dat <- CBCre$data


## ---------- individual bar plots per region in a grid

pattern_old = c("1", "2", "3")
pattern_you = c("4", "5", "6", "7", "10", "12")


#  function for individual bar plots

plotRE_ams <- function(regionre_dat) {
  m_pre_o <- with(regionre_dat, mean(expr[age %in% pattern_old & stage == "Prenatal"])) # takes the expression values that are in PS 1-3 ("pattern_old") and in a specific stage ("Prenatal") and calculates the mean
  m_pre_y <- with(regionre_dat, mean(expr[age %in% pattern_you & stage == "Prenatal"])) # takes the expression values that are in PS 4-12 ("pattern_young") and in a specific stage ("Prenatal") and calculates the mean
  m_inf_o <- with(regionre_dat, mean(expr[age %in% pattern_old & stage == "Infant"]))   # takes the expression values that are in PS 1-3 ("pattern_old") and in a specific stage ("Infant") and calculates the mean
  m_inf_y <- with(regionre_dat, mean(expr[age %in% pattern_you & stage == "Infant"]))   # takes the expression values that are in PS 4-12 ("pattern_young") and in a specific stage ("Infant") and calculates the mean
  m_chi_o <- with(regionre_dat, mean(expr[age %in% pattern_old & stage == "Child"]))    # ""
  m_chi_y <- with(regionre_dat, mean(expr[age %in% pattern_you & stage == "Child"]))    # ""
  m_ado_o <- with(regionre_dat, mean(expr[age %in% pattern_old & stage == "Adolescent"])) # ""
  m_ado_y <- with(regionre_dat, mean(expr[age %in% pattern_you & stage == "Adolescent"])) # ""
  m_adu_o <- with(regionre_dat, mean(expr[age %in% pattern_old & stage == "Adult"]))    # ""
  m_adu_y <- with(regionre_dat, mean(expr[age %in% pattern_you & stage == "Adult"]))    # ""
  
  exp_re <- m_pre_o %>%
    rbind(m_pre_y) %>%
    rbind(m_inf_o) %>%
    rbind(m_inf_y) %>%
    rbind(m_chi_o) %>%
    rbind(m_chi_y) %>%
    rbind(m_ado_o) %>%
    rbind(m_ado_y) %>%
    rbind(m_adu_o) %>%
    rbind(m_adu_y)
  
  exp_re <- data.frame(exp_re)
  exp_re$PSclass <- rep(c("Old", "Young"), 5)
  exp_re$stage <- c(rep("1Prenatal", 2), rep("2Infant", 2), rep("3Child",2), rep("4Adolescent",2), rep("5Adult",2))
  
  re_p <- 
    exp_re %>% 
    group_by(stage) %>% 
    mutate(sd = sd(exp_re)) %>% 
    ggplot(aes(x = stage, y = exp_re, fill=PSclass)) + 
    geom_bar(stat="identity", position=position_dodge()) 
  re_p <- re_p + ylim(0, 1.25) + theme_classic() + scale_x_discrete(labels=c("Prenatal", "Infant", "Child", "Adolescent", "Adult"))
#  re_p <- re_p + ylab("Relative expression intensity\n") + xlab("Ontogenetic stage") +
#    labs(fill = "Phylostratum \nage category") 
  re_p <- re_p + theme(axis.title.x = element_blank(),
                       axis.text.x = element_text(color="white", size=9, angle=40, hjust=0.8, vjust=.8, face="bold"),
                       axis.title.y = element_blank(),
                       axis.text.y = element_text(color="white", size=11))
  re_p <- re_p + geom_errorbar(aes(ymin=exp_re, ymax=exp_re+sd), width=.2, colour="orange", position=position_dodge(.9))
  #re_p <- re_p + scale_x_discrete(labels=c("Prenatal", "Infant", "Child", "Adolescent", "Adult"))
  re_p <- re_p + theme(legend.position = "none")
  
  value <- list(
    exp_re = exp_re,
    re_p = re_p
  ) # Create a list of output objects
  attr(value, "class") <- "plotRE_ams"
  value
  
}

# 1st row
# M1C
outm1c <- plotRE_ams(M1Cre_dat) 
M1Cre_p <- outm1c$re_p + ggtitle("M1C") + theme(axis.text.y = element_text(face="bold", color="black", size=11, margin = margin(t = 0, r = 10, b = 0, l = 0)),
                                                plot.title = element_text(size=18, face="bold")) 
M1Cre_p <- M1Cre_p + scale_fill_manual(values=c("#104C7E", "#90CEF0"))   #6E0288, #F3C6FE
M1Cre_p

# DFC
outdfc <- plotRE_ams(DFCre_dat) 
DFCre_p <- outdfc$re_p + ggtitle("DFC") + theme(plot.title = element_text(size=18, face="bold"))
DFCre_p <- DFCre_p + scale_fill_manual(values=c("#104C7E", "#90CEF0"))   #5B2187, #D5B5ED
DFCre_p

# VFC
outvfc <- plotRE_ams(VFCre_dat) 
VFCre_p <- outvfc$re_p + ggtitle("VFC") + theme(plot.title = element_text(size=18, face="bold"))
VFCre_p <- VFCre_p + scale_fill_manual(values=c("#104C7E", "#90CEF0"))   #414487, #C4C5E2
VFCre_p

# OFC
outofc <- plotRE_ams(OFCre_dat) 
OFCre_p <- outofc$re_p + ggtitle("OFC") + theme(plot.title = element_text(size=18, face="bold"))
OFCre_p <- OFCre_p + scale_fill_manual(values=c("#104C7E", "#90CEF0"))   #2A788E, #AFDBE7
OFCre_p

# 2nd row
# S1C        
outs1c <- plotRE_ams(S1Cre_dat) 
S1Cre_p <- outs1c$re_p + ggtitle("S1C") + theme(axis.text.y = element_text(face="bold", color="black", size=11, margin = margin(t = 0, r = 10, b = 0, l = 0)),
                                                plot.title = element_text(size=18, face="bold"))
S1Cre_p <- S1Cre_p + scale_fill_manual(values=c("#104C7E", "#90CEF0"))   #51358F, #CEC2E8
S1Cre_p

# IPC
outipc <- plotRE_ams(IPCre_dat) 
IPCre_p <- outipc$re_p + ggtitle("IPC") + theme(plot.title = element_text(size=18, face="bold"))
IPCre_p <- IPCre_p + scale_fill_manual(values=c("#104C7E", "#90CEF0"))   #39568C, #C9D4E9
IPCre_p

# A1C
outa1c <- plotRE_ams(A1Cre_dat) 
A1Cre_p <- outa1c$re_p + ggtitle("A1C") + theme(plot.title = element_text(size=18, face="bold"))
A1Cre_p <- A1Cre_p + scale_fill_manual(values=c("#104C7E", "#90CEF0"))   #23888E, #A2E5E8
A1Cre_p

# STC
outstc <- plotRE_ams(STCre_dat) 
STCre_p <- outstc$re_p + ggtitle("STC") + theme(plot.title = element_text(size=18, face="bold"))
STCre_p <- STCre_p + scale_fill_manual(values=c("#104C7E", "#90CEF0"))   #3CB250, #BAE8C2
STCre_p

# 3rd row
# ITC
outitc <- plotRE_ams(ITCre_dat) 
ITCre_p <- outitc$re_p + ggtitle("ITC") + theme(axis.text.y = element_text(face="bold", color="black", size=11, margin = margin(t = 0, r = 10, b = 0, l = 0)),
                                                plot.title = element_text(size=18, face="bold"))
ITCre_p <- ITCre_p + scale_fill_manual(values=c("#104C7E", "#90CEF0"))   #2B5C7D, #A9CAE1
ITCre_p

# V1C
outv1c <- plotRE_ams(V1Cre_dat) 
V1Cre_p <- outv1c$re_p + ggtitle("V1C") + theme(plot.title = element_text(size=18, face="bold"))
V1Cre_p <- V1Cre_p + scale_fill_manual(values=c("#104C7E", "#90CEF0"))   #1F988B, #9DEBE2
V1Cre_p

# MFC
outmfc <- plotRE_ams(MFCre_dat) 
MFCre_p <- outmfc$re_p + ggtitle("MFC") + theme(plot.title = element_text(size=18, face="bold"))
MFCre_p <- MFCre_p + scale_fill_manual(values=c("#104C7E", "#90CEF0"))   #35B779, #AEE8CC
MFCre_p

# HIP
outhip <- plotRE_ams(HIPre_dat) 
HIPre_p <- outhip$re_p + ggtitle("HIP") + theme(plot.title = element_text(size=18, face="bold"))
HIPre_p <- HIPre_p + scale_fill_manual(values=c("#104C7E", "#90CEF0"))   #95CD25, #D2ED9D
HIPre_p

# 4th row
# STR
outstr <- plotRE_ams(STRre_dat) 
STRre_p <- outstr$re_p + ggtitle("STR") + 
  theme(axis.text.x = element_text(size=9, color="black", angle=40, hjust=0.8, vjust=.8, face="bold"),
        axis.text.y = element_text(face="bold", color="black", size=11, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        plot.title = element_text(size=18, face="bold"))
STRre_p <- STRre_p + scale_fill_manual(values=c("#104C7E", "#90CEF0"))   #22A884, #AFEFDE"
STRre_p

# AMY
outamy <- plotRE_ams(AMYre_dat) 
AMYre_p <- outamy$re_p + ggtitle("AMY") + 
  theme(axis.text.x = element_text(size=9, color="black", angle=40, hjust=0.8, vjust=.8, face="bold"),
        plot.title = element_text(size=18, face="bold"))
AMYre_p <- AMYre_p + scale_fill_manual(values=c("#104C7E", "#90CEF0"))   #60C234, #BBE7A7
AMYre_p

# MD
outmd <- plotRE_ams(MDre_dat) 
MDre_p <- outmd$re_p + ggtitle("MD") + 
  theme(axis.text.x = element_text(size=9, color="black", angle=40, hjust=0.8, vjust=.8, face="bold"),
        plot.title = element_text(size=18, face="bold"))
MDre_p <- MDre_p + scale_fill_manual(values=c("#104C7E", "#90CEF0"))     #B7C317, #EAF18B
MDre_p

# CBC
outcbc <- plotRE_ams(CBCre_dat) 
CBCre_p <- outcbc$re_p + ggtitle("CBC") + 
  theme(axis.text.x = element_text(size=9, color="black", angle=40, hjust=0.8, vjust=.8, face="bold"),
        plot.title = element_text(size=18, face="bold"))
CBCre_p <- CBCre_p + scale_fill_manual(values=c("#104C7E", "#90CEF0"))   #D4C002, #FEF38A
CBCre_p


# ----------------- combine plots in one large grid

grid_re <- plot_grid(M1Cre_p, DFCre_p, VFCre_p, OFCre_p, 
                     S1Cre_p, IPCre_p, A1Cre_p, STCre_p, 
                     ITCre_p, V1Cre_p, MFCre_p, HIPre_p,
                     STRre_p, AMYre_p, MDre_p, CBCre_p,
                     nrow = 4)
grid_re

ggsave(filename = paste0(BASE, "supp_mat06.pdf"), grid_re,    # ehemals "PlotRE_grid.pdf"
       width = 14, height = 12, units = "in", device='pdf')





## ================================= Plot Phylostratum Enrichment 
# myTAI documentation p44 and p 18
# This Phylostratum or Divergence Stratum enrichment analysis is motivated by Sestak and Domazet-
# Loso (2015) who perform Phylostratum or Divergence Stratum enrichment analyses to correlate
# organ evolution with the origin of organ specific genes.

# Gene IDs for OT gene set as reference (should not matter from which mm file these are retrieved...?)
M1C_mm <- M1C_tai$mm
OT_subset <-M1C_mm$GeneID

# !!IMPORTANT!! check whether the genes are *really* the same 

# INFO: This analysis only needs to be performed once, it is not necessary to calculate everything for each brain regions individuals. The function only compares the
# phylostrata and these are the same in each expression set, it does NOT call the expression values. 


# get expression data from all ontogenetic stages for each brain region for *ALL* genes (ens gene ids in the phylo map as reference)
exp_allENRICH <- get_expression(structure_ids=c("Allen:10163"),                       # any brain region, doesn't matter as PS distributions are the same/remain unchanged
                            gene_ids=HomoSapiens.PhyloMap$Gene_ID, 
                            dataset='5_stages') 
  
xallENRICH_df <- as.data.frame(do.call(rbind, exp_allENRICH)) 
  
rownames(xallENRICH_df) <- c("Prenatal", "Infant", "Child", "Adolescent", "Adult")
  
xallENRICH_t <- as.data.frame(t(xallENRICH_df))
xallENRICH_t <- rownames_to_column(xallENRICH_t, var = "Gene_ID")
xallENRICH_t <- xallENRICH_t %>% distinct(Gene_ID, .keep_all = TRUE)                    # remove duplicates?
# available genes: 15858
  
  
# merge expression and phylogenetic data
mm_allENRICH <- MatchMap(HomoSapiens.PhyloMap, xallENRICH_t)
# available genes: 15858
  
enrich_p_dat <- EnrichmentTest(ExpressionSet = mm_allENRICH, 
                               test.set = OT_subset, 
                               measure = "log-foldchange",
                               use.only.map = F, 
                               complete.bg = T)   
  
res_m <- data.frame(enrich_p_dat$enrichment.matrix)
  
res_m$pvals <- enrich_p_dat$p.values
  
res_m$pFDR <-  p.adjust(res_m$pvals, method = "fdr", n = length(res_m$pvals)) 

# Usually, PlotEnrichment() adjust for multiple comparisons with p.adjust(), but EnrichmentTest() just... does not include that parameter
# even though it does essential the same except for the plotting?

res_m <- tibble::rownames_to_column(res_m, "Phylostratum") #library(dplyr)
res_m$Phylostratum <- Map(paste0, 'PS', res_m$Phylostratum)

# thresholds copied from PlotEnrichment() source code
res_m$star <- ""
res_m$star[res_m$pFDR <= .05]  <- "*"
res_m$star[res_m$pFDR <= .005]  <- "**"
res_m$star[res_m$pFDR <= 5e-04] <- "***"


bp <- ggplot(res_m) + 
  geom_col(aes(y=factor(Phylostratum, levels=Phylostratum), x = Test_Set, fill=Test_Set)) +
  theme_classic() 
bp <- bp +  
  scale_fill_gradient(breaks=c((min(res_m$Test_Set) + .3), (max(res_m$Test_Set) - .3)), 
                      labels = c("Under-represented", "Over-represented"), 
                      name="") 
bp <- bp + 
  geom_vline(aes(xintercept = 0), size=2) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5))
bp <- bp + 
  xlab("\n\nRelative enrichment") +  ylab("Phylostratum \n (phylogenetic stage)\n") 
bp <- bp + 
  theme(axis.title  = element_text(size = 23, face="bold"),
        axis.text.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        legend.position = c(.75, .9),
        legend.title = element_blank(),
        legend.text = element_text(size=16))#,
        #panel.grid.major.y = element_blank(),
        #panel.grid.major.x = element_line( size=.1, color="gray90"))
bp2 <- bp + 
  annotate("text", y = 1, x = 0.89897619, label = "*", size=10) +
  annotate("text", y = 2, x = -1.19685124, label = "*", size=10)

bp2


#ggsave(filename = paste0(BASE, "RProject/OT-bd-analyses_loc/output/enrichmentPS_new.pdf"), bp,
#       width = 10, height = 12, units = "in", device='pdf')



# combine pcom with bp (TAI with enrichment)
enrich_TAIcom <- plot_grid(bp, NULL, pg,
                           ncol = 3, rel_widths = c(.4, 0.05, 1),
                           labels = c('A', 'B'),
                           label_size = 28) 

enrich_TAIcom


ggsave(filename = paste0(BASE, "RProject/OT-bd-analyses_loc/output/enrich_TAIcom_v2_6.pdf"), enrich_TAIcom,
       width = 22, height = 10, units = "in", device='pdf')

# alternative red: EA0075

# for solid bar colours w/o gradient
# M1C_res_m <- mutate(M1C_res_m, Enrichment=if_else(Test_Set>0, "Over-represented", "Under-represented"))
# + scale_fill_manual(values=c("#FF0564", "#3D477B")) 




## ================================= Early Conservation & hourglass models

# UPDATE: Possibly run the analyses with non-human data from GEO?

# Probably we cannot compute these because the theories/models are about different stages during embryogenesis. We do not have
# detailed information about specific prenatal stages, we just have one big "prenatal" stage... at least in the data set that is 
# retrieved here. But in the downloadable version, the prenatal stage is available much more detailed...




# >>>>>>>>>>>>>>>>>>>>>> #
##########################
### end myTAI analyses ###
##########################
# >>>>>>>>>>>>>>>>>>>>>> #







## ================================= dN/dS values from Dumas et al., 2021 --- GenEvo tool

# load file generated by https://genevo.pasteur.fr/, entered the OT pathway geneset as gene list. 
# Settings: Exact match, only 1-to-1 orthologs

dnds_primates_OT <- fread(paste0(BASE, "RProject/OT-bd-analyses_loc/Data/genevo_2021-04-30_1to1.tsv"))
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
 
# combine low quality barplot and phylogenetic tree
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
  viridis::scale_fill_viridis(option = 'magma', discrete = T, begin = .2, end = .85) +
  viridis::scale_colour_viridis(option = 'magma', discrete = T, begin = .2, end = .85)
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
                                colour        = "steelblue4",
                                fill          = alpha(c("steelblue2"),.4),
                                size          = 4) 
p_mq <- p_mq + geom_label_repel(data=(subset_mq_l %>% dplyr::filter(Gene %in% c("OXT"))), 
                                aes(label     = Gene), 
                                nudge_x       = -.05, 
                                nudge_y       = .4,
                                label.padding = .25,
                                colour        = "steelblue4",
                                fill          = alpha(c("steelblue2"),.4),
                                size          = 4)
p_mq <- p_mq + geom_label_repel(data=(subset_mq_l %>% dplyr::filter(Gene %in% c("CD38"))), 
                               aes(label      = Gene), 
                               nudge_x        = 0, 
                               nudge_y        = .4,
                               label.padding  = .25,
                               colour         = "steelblue4",
                               fill           = alpha(c("steelblue2"),.4),
                               size           = 4)
p_mq <- p_mq + theme(legend.position = "none")

p_mq

#ggsave(filename = paste0(BASE, "RProject/OT-bd-analyses_loc/output/testy.pdf"), p_mq,
#       width = 18, height = 8, units = "in", device='pdf')

# combine med quality barplot and phylogenetic tree
dnds_OT_tree_mq <- plot_grid(gg_ptr2, p_mq, rel_widths = c(1.1, 3)) 
dnds_OT_tree_mq

# save plot
ggsave(filename = paste0(BASE, "RProject/OT-bd-analyses_loc/output/dndsOT_tree_mq_phylo4.pdf"), dnds_OT_tree_mq,
       width = 20, height = 8, units = "in", device='pdf')



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


DSOTtop50 <- dplyr::full_join(OTpthwy_df, DStop50, by=c("OTpthwy_nms"="gene_names..1..."))       # match with *all* OT genes
DSOTlower50 <- dplyr::full_join(OTpthwy_df, DSlower50, by=c("OTpthwy_nms"="gene_names..1..."))    # match with OT *all* genes

OTgenes.out <- dplyr::anti_join(OTpthwy_df, DSabagen_anno_desc, by=c("OTpthwy_nms"="gene_names..1..."))    # check whether some genes were not available

# no data available for 18 genes
# out of 136 genes, 105 were among the top 50% genes with highest DS values, 31 were among the genes with lower DS values
# out of the genes for which data was available (n = 136), 77.21% have high DS values and are thus consistently and reliabley expressed across brain regions independet of donor characteristics. 



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



