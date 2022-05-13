#########################
# setting up working environment 
#########################

rm(list=ls()) # delete all objects in the workspace
gc(reset=T) # resest memory (especially useful when working with large data sets)

options(stringsAsFactors=F) # disables automatic conversion of char. strings into factors

Sys.setenv(LANG = "en")

#install.packages("installr")

# check for R version updates
#library(installr)
#updateR()

# PRO TIP (jk, it's a rookie tip): Ctr + left click on an object to view it in a new R tab

BASE="C:/01Dokumente/07Erwerbstaetigkeit/01PhD_UiO/projects/first_author/oxytocinHealth/RQ4/"      ### IMPORTANT: adjust and update wd

# Base dir for laptop
BASE="D:/eigene_dateien/02Dokumente/07Erwerbstaetigkeit/01PhD_UiO/projects/first_author/oxytocinHealth/RQ4/"

#########################
# setting up working environment 
#########################


# Manuscript analysis 

# general
library(tidyverse)
library(tibble)
library(readxl)
library(data.table)
library(readr)
library(tidyr)
library(powerAnalysis)
library(gage)           # BiocManager
library(KEGGREST)       # BiocManager
library(DataCombine)
library(writexl)

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

# phylostratigraphy, phylogenetics
library(ABAEnrichment) 
require(ABAData)
library(EnsDb.Hsapiens.v79)
library(ggtree)          # BiocManager
library(ape)

# brain expression enrichment
library(ggseg)
library(ggsegDefaultExtra)
#library(ggsegExtra)
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
#########

# retrieve gene list from KEGG
# NOTE: genes belonging to a pathway are updated constantly, new genes are added on a regular basis, thus the
# list you retrieve might vary slightly from the list we initially retrieved


kg.hsa = kegg.gsets("hsa") # retrieve gene set data from KEGG for Homo SApiens
kg.hsa.sets <- kg.hsa$kg.sets # fetch gene sets and search for "oxytocin signaling pathway"

# get gene names
OTpthwy <- keggGet("hsa04921")[[1]]$GENE    
OTpthwy_nms_base <- OTpthwy[seq(2,length(OTpthwy),2)]
OTpthwy_nms <- gsub("\\;.*","",OTpthwy_nms_base)
OTpthwy_df <- data.frame(OTpthwy_nms)

# export as excel for supps

write_xlsx(OTpthwy_df, paste0(BASE, "OTpthwy_list.xlsx"))




## =================================  Microsynteny analyses

# ................. 



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

#hyp_ex <- get_expression(structure_ids=left_hyp,      # fetch expression data
#                         gene_ids=OTpw.genes, 
#                         dataset='adult')

#hyp_ex <- as.data.frame(hyp_ex)                       # prepare data: transpose, create column with gene labels

#hyp_ex_t <- as.data.frame(t(hyp_ex))

#hyp_ex_t <- rownames_to_column(hyp_ex_t, var = "gene")

#hyp_ex_t$gene

#ens_gene_hyp <- ensembldb::select(EnsDb.Hsapiens.v79,                   # fetch ensemble gene names
#                                  keys= hyp_ex_t$gene, 
#                                  keytype = "SYMBOL", 
#                                  columns = c("SYMBOL","GENEID"))

#ens_gene_hyp <- ens_gene_hyp %>% distinct(SYMBOL, .keep_all = TRUE)     # remove duplicates (?)
#names(ens_gene_hyp)[names(ens_gene_hyp) == "SYMBOL"] <- "gene"
#hyp_ex_t_f <- dplyr::full_join(ens_gene_hyp, hyp_ex_t, by = "gene")     # join with expression data
#hyp_ex_t_f <- dplyr::select(hyp_ex_t_f, -gene)

#names(hyp_ex_t_f)[names(hyp_ex_t_f) == "GENEID"] <- "Gene_ID"

counts <- table(HomoSapiens.PhyloMap2013$Phylostratum)

barplot(counts, main="Phylostratum levels",
        xlab="Phylostratum")

#mm_hyp2013 <- MatchMap(HomoSapiens.PhyloMap2013, hyp_ex_t_f)  # 141 genes, maybe because it is aligned with a different expression set
#write.csv(mm_hyp2013, paste0(BASE, "mm_hyp2013.csv"))


expDie <- get_expression(structure_ids=c("Allen:10389"),       
                         gene_ids=OTpw.genes, 
                         dataset='5_stages')                  # Get expression values, dataset: BrainSpan



expDie_df <- as.data.frame(do.call(rbind, expDie)) 
rownames(expDie_df) <- NULL
rownames(expDie_df) <- c("Prenatal", "Infant", "Child", "Adolescent", "Adult")
expDie_t <- as.data.frame(t(expDie_df))
expDie_t <- rownames_to_column(expDie_t, var = "gene")                         

ens_gene_devDie <- ensembldb::select(EnsDb.Hsapiens.v79, 
                                     keys= expDie_t$gene, 
                                     keytype = "SYMBOL", 
                                     columns = c("SYMBOL","GENEID"))

ens_gene_devDie <- ens_gene_devDie %>% distinct(SYMBOL, .keep_all = TRUE)

names(ens_gene_devDie)[names(ens_gene_devDie) == "SYMBOL"] <- "gene"

expDie_t_f <- dplyr::full_join(ens_gene_devDie, expDie_t, by = "gene")

expDie_t_f <- dplyr::select(expDie_t_f, -gene)

names(expDie_t_f)[names(expDie_t_f) == "GENEID"] <- "Gene_ID"     

mmDie5stages <- MatchMap(HomoSapiens.PhyloMap2013, expDie_t_f)



# ---- get supp. table with genes and phylostrata
mmDie5stages_reduced <- mmDie5stages[,1:2] 
mmDie5stages_reduced$GeneID <- toupper(mmDie5stages_reduced$GeneID)

supp_mat01 <- dplyr::left_join(mmDie5stages_reduced, ens_gene_devDie, by=c("GeneID"="GENEID"))
names(supp_mat01)[names(supp_mat01) == "gene"] <- "GeneLabel"
supp_mat01 <- supp_mat01 %>%       # order ascending by lobe label
  arrange((Phylostratum))


col1_sm <- NA
col2_sm <- NA
supp_mat01.2 <- dplyr::anti_join(data.frame(OTpw.genes), supp_mat01, by= c("OTpw.genes"="GeneLabel"))
supp_mat01.22 <- cbind(col1_sm, col2_sm, supp_mat01.2)
names(supp_mat01.22)[names(supp_mat01.22) == "col1_sm"] <- "Phylostratum"
names(supp_mat01.22)[names(supp_mat01.22) == "col2_sm"] <- "GeneID"
names(supp_mat01.22)[names(supp_mat01.22) == "OTpw.genes"] <- "GeneLabel"

supp_mat01_fin <- rbind(supp_mat01, supp_mat01.22)

library("writexl")
#write_xlsx(supp_mat01_fin, paste0(BASE, "supp_mat01_e2013.xlsx"))
# ---- ---- ----


# calculate relative frequencies for OT genes
rel_fOT <- table(mmDie5stages$Phylostratum)/length(mmDie5stages$Phylostratum)
rel_fOT_r <- round((rel_fOT*100), digits = 1)
abs_fOT <- table(mmDie5stages$Phylostratum)


c1 <- c(94, 24, 3, 6, 1, 7, 4, 0, 0, 0, 0, 3)                                           # absolute freqs  
c2 <- c(0.661971831, 0.169014085, 0.021126761, 0.042253521, 0.007042254, 0.049295775, 0.028169014, 0, 0, 0, 0, 0.021126761)   # relative freqs
c3 <- c(66.2, 16.9, 2.1, 4.2, 0.7, 4.9, 2.8, 0, 0, 0, 0, 2.1)                         # relative freqs, rounded






#### ALL genes
# get expression data from all ontogenetic stages for each brain region for as many genes as possible (ens gene ids in the phylo map as reference)

#ens_gene_devALL <- ensembldb::select(EnsDb.Hsapiens.v79,                                  
#                                     keys    = HomoSapiens.PhyloMap2013$EnsemblGeneID, 
#                                     keytype = "GENEID", 
#                                    columns = c("SYMBOL","GENEID"))

#ens_gene_devALL <- ens_gene_devALL %>% distinct(SYMBOL, .keep_all = TRUE)


exp_all <- get_expression(structure_ids = c("Allen:10389"),      # AHBA structure ontology ID for the diencephalon. The other IDs for the hypothalamus did not work. 
                          gene_ids      = HomoSapiens.PhyloMap2013$EnsemblGeneID, 
                          dataset       = '5_stages')

exp_all_df <- as.data.frame(do.call(rbind, exp_all)) 
rownames(exp_all_df) <- NULL
rownames(exp_all_df) <- c("Prenatal", "Infant", "Child", "Adolescent", "Adult")
exp_all_t <- as.data.frame(t(exp_all_df))
exp_all_t <- rownames_to_column(exp_all_t, var = "gene")

names(exp_all_t)[names(exp_all_t) == "gene"] <- "Gene_ID"   
exp_all_t <- exp_all_t %>% distinct(Gene_ID, .keep_all = TRUE) 

mmALLE <- MatchMap(HomoSapiens.PhyloMap2013, exp_all_t)         # merge expression and phylogenetic data
# final available genes: 17027


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

#ggsave(filename = paste0(BASE, "RProject/OT-bd-analyses_loc/output/fig1template2013c.jpg"), testxb,
#       width = 46, height = 20, units = "in", device='jpg')


# ---------> proceed with bioRender




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


# BACKGROUND SET: mmALLE
# available genes: 17027

enrich_p_dat <- EnrichmentTest(ExpressionSet = mmALLE, 
                               test.set = OT_subset, 
                               measure = "log-foldchange",
                               use.only.map = F, 
                               complete.bg = T)   

res_m <- data.frame(enrich_p_dat$enrichment.matrix)

res_m$pvals <- enrich_p_dat$p.values

res_m$pFDR <-  p.adjust(res_m$pvals, method = "fdr")#, n = length(res_m$pvals)) 

# Usually, PlotEnrichment() adjust for multiple comparisons with p.adjust(), but EnrichmentTest() just... does not include that parameter
# even though it does essential the same except for the plotting?

res_m <- tibble::rownames_to_column(res_m, "Phylostratum") #library(dplyr)
res_m$Phylostratum <- Map(paste0, 'PS', res_m$Phylostratum)

# thresholds copied from PlotEnrichment() source code
res_m$star <- ""
res_m$star[res_m$pFDR <= .05]  <- "*"
res_m$star[res_m$pFDR <= .005]  <- "**"
res_m$star[res_m$pFDR <= 5e-04] <- "***"

res_m_short <- res_m[-c(13:20),] # remove phylostrata that do not contain any OT genes



# plotting

bp <- ggplot(res_m_short) + 
  geom_col(aes(y=factor(Phylostratum, levels=Phylostratum), x = Test_Set, fill=Test_Set)) +
  theme_classic() 
bp <- bp +  
  scale_fill_gradient(breaks=c((min(res_m_short$Test_Set) + .3), (max(res_m_short$Test_Set) - .3)), 
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
        legend.position = c(.3, .9),
        legend.title = element_blank(),
        legend.text = element_text(size=16))#,
#panel.grid.major.y = element_blank(),
#panel.grid.major.x = element_line( size=.1, color="gray90"))
bp <- bp + 
  annotate("text", y = 0.9, x = 1.03629910, label = "*", size=10) +
  annotate("text", y = 1.9, x = -0.93520941, label = "*", size=10) +
  annotate("text", y = 4.9, x = -3.03374016, label = "*", size=10) +
  annotate("text", y = 7.9, x = .25, label = "NA", size=6, color = "grey", fontface = "italic") +
  annotate("text", y = 8.9, x = .25, label = "NA", size=6, color = "grey", fontface = "italic") +
  annotate("text", y = 9.9, x = .25, label = "NA", size=6, color = "grey", fontface = "italic") +
  annotate("text", y = 10.9, x = .25, label = "NA", size=6, color = "grey", fontface = "italic")

bp


#ggsave(filename = paste0(BASE, "RProject/OT-bd-analyses_loc/output/enrichmentPS_new_ps1-12.pdf"), bp,
#       width = 10, height = 12, units = "in", device='pdf')



# combine pcom with bp (TAI with enrichment)
enrich_TAIcom <- plot_grid(bp, NULL, pg,
                           ncol = 3, rel_widths = c(.4, 0.05, 1),
                           labels = c('A', 'B'),
                           label_size = 28) 

enrich_TAIcom


ggsave(filename = paste0(BASE, "RProject/OT-bd-analyses_loc/output/enrich_TAIcom_v2_2013d.pdf"), enrich_TAIcom,
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

