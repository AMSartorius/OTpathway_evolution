######################################
##### SET UP WORKING ENVIRONMENT #####
######################################

# rm(list = ls()) # delete all objects in the workspace, USE WITH CAUTION
# gc(reset = T)   # resest memory (especially useful when working with large data sets), USE WITH CAUTION

options(stringsAsFactors = F) # disables automatic conversion of char. strings into factors
Sys.setenv(LANG = "en")

### .......... loading packages
library(readxl)          # load/read in .xlsx files
library(data.table)      # fread() function
library(stringr)
library(tidyverse)
library(writexl)         # export .xlsx files 


# preprare gene age categories for submission to FUMA

OTpthwy.res <- data.frame(read_excel("output/BLASTp/OTpthwy_res02_vertebrate_threshold.xlsx", sheet = 1))

# modern vertebates
modern <- data.frame(OTpthwy.res[OTpthwy.res$age.3.levels.vert == "mod.vertebrate", "Gene"])
colnames(modern) <- "modern.genes"
# write_xlsx(modern, "data/processed/modern_genes.xlsx")

# modern non-vertebates
med.age <- data.frame(OTpthwy.res[OTpthwy.res$age.3.levels.vert == "mod.non.vertebrate", "Gene"])
colnames(med.age) <- "med.age.genes"
# write_xlsx(med.age, "data/processed/medage_genes.xlsx")

# ancient
ancient <- data.frame(OTpthwy.res[OTpthwy.res$age.3.levels.vert == "ancient", "Gene"])
colnames(ancient) <- "ancient.genes"
# write_xlsx(ancient, "data/processed/ancient_genes.xlsx")




### !!! Use the files in the respective folder for subsequent analyses (or replicate your own FUMA analyses and return to this script, but the files provided work just fine) !!!




###############################
# FUMA Main OT analysis results 
###############################

# ----------------------- ancient gene set 
# settings: 28 genes, exclude MHC, use FDR correction, date: 22/05/18

sum.anc <- fread("data/raw/FUMA/ancient_FUMA_gene2func82286/summary.txt", h = T)

geneTab.anc <- read.table("data/raw/FUMA/ancient_FUMA_gene2func82286/geneTable.txt", h = T)
geneIDs.anc <- read.table("data/raw/FUMA/ancient_FUMA_gene2func82286/geneIDs.txt", h = T)

GS.anc <- fread("data/raw/FUMA/ancient_FUMA_gene2func82286/GS.txt", h = T)

# for supplementary materials
# write_xlsx(GS.anc, "output/supp_mat_new/sup_mat_06a.xlsx")

GS.anc.GWAS <- GS.anc[GS.anc$Category == "GWAScatalog",]

# none              


gtex.30.log.anc  <- read.table("data/raw/FUMA/ancient_FUMA_gene2func82286/gtex_v8_ts_general_avg_log2TPM_exp.txt", h = T)
gtex.30.norm.anc <- read.table("data/raw/FUMA/ancient_FUMA_gene2func82286/gtex_v8_ts_general_avg_normTPM_exp.txt", h = T)
gtex.30.DEG.anc  <-      fread("data/raw/FUMA/ancient_FUMA_gene2func82286/gtex_v8_ts_general_DEG.txt", h = T)


tissueEnrich30anc <- gtex.30.DEG.anc[gtex.30.DEG.anc$Category == "DEG.up" & gtex.30.DEG.anc$adjP < 0.05, ]

# none

# ..... continue working with "gtex.30.DEG.anc"

gtex30.anc.up <- subset(gtex.30.DEG.anc, grepl("DEG.up", gtex.30.DEG.anc$Category)) 
gtex30.anc.up$significant <- "n"
gtex30.anc.up$significant[gtex30.anc.up$adjP <= 0.05]  <- "y"      # no rows with "y", good double check

gtex30.anc.up$log10pval <- -log10(gtex30.anc.up$p)


df.anc <- data.frame(c(gtex30.anc.up[, c(2, 8, 9)]))    # select relevant columns

df.anc$GeneSet <- str_replace_all(df.anc$GeneSet, "_T", " t")
df.anc$GeneSet <- str_replace_all(df.anc$GeneSet, "_V", " v")
df.anc$GeneSet <- str_replace_all(df.anc$GeneSet, "_G", " g")
df.anc$GeneSet <- str_replace_all(df.anc$GeneSet, "_I", " i")
df.anc$GeneSet <- str_replace_all(df.anc$GeneSet, "_Uteri", "")

df.anc <- df.anc %>%                                    # order descending by mean
  arrange(desc(log10pval))

p.df.anc <- ggplot(df.anc, aes(x = GeneSet, y = log10pval, fill = significant)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#999999", "#8D6EFF"),   
                    breaks = c("n", "y")) +
  theme_classic() 

p.df.anc <- p.df.anc + 
  aes(x = fct_inorder(GeneSet)) + xlab("Tissue") + ylab("-log 10 p-value")

p.df.anc <- p.df.anc + 
  theme(axis.text.x = element_text(color = "black", size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title  = element_text(size = 20, face = "bold"),
        legend.position = "none") +
  ylim(0, 8.5)

p.df.anc









# ----------------------- medium-aged gene set 
# settings: 28 genes, exclude MHC, use FDR correction, date: 22/05/18

sum.medage     <- fread("data/raw/FUMA/medage_FUMA_gene2func82287/summary.txt", h = T)

geneTab.medage <- read.table("data/raw/FUMA/medage_FUMA_gene2func82287/geneTable.txt", h = T)
geneIDs.medage <- read.table("data/raw/FUMA/medage_FUMA_gene2func82287/geneIDs.txt", h = T)

GS.medage      <- fread("data/raw/FUMA/medage_FUMA_gene2func82287/GS.txt", h = T)

# for supplementary materials
# write_xlsx(GS.medage, "output/supp_mat_new/sup_mat_06b.xlsx")

GS.medage.GWAS <- GS.medage[GS.medage$Category == "GWAScatalog",]

# none              


gtex.30.log.medage  <- read.table("data/raw/FUMA/medage_FUMA_gene2func82287/gtex_v8_ts_general_avg_log2TPM_exp.txt", h = T)
gtex.30.norm.medage <- read.table("data/raw/FUMA/medage_FUMA_gene2func82287/gtex_v8_ts_general_avg_normTPM_exp.txt", h = T)
gtex.30.DEG.medage  <- fread("data/raw/FUMA/medage_FUMA_gene2func82287/gtex_v8_ts_general_DEG.txt", h = T)

tissueEnrich30medage <- gtex.30.DEG.medage[gtex.30.DEG.medage$Category == "DEG.up"  & gtex.30.DEG.medage$adjP < 0.05, ]

# bladder & blood vessel


# ..... continue working with "gtex.30.DEG.medage"


gtex30.medage.up <- subset(gtex.30.DEG.medage, grepl("DEG.up", gtex.30.DEG.medage$Category)) 
gtex30.medage.up$significant <- "n"
gtex30.medage.up$significant[gtex30.medage.up$adjP <= 0.05]  <- "y"

gtex30.medage.up$log10pval <- -log10(gtex30.medage.up$p)


df.medage <- data.frame(c(gtex30.medage.up[, c(2, 8, 9)]))    # select relevant columns

df.medage$GeneSet <- str_replace_all(df.medage$GeneSet, "_T", " t")
df.medage$GeneSet <- str_replace_all(df.medage$GeneSet, "_V", " v")
df.medage$GeneSet <- str_replace_all(df.medage$GeneSet, "_G", " g")
df.medage$GeneSet <- str_replace_all(df.medage$GeneSet, "_I", " i")
df.medage$GeneSet <- str_replace_all(df.medage$GeneSet, "_Uteri", "")

df.medage <- df.medage %>%                                    # order descending by mean
  arrange(desc(log10pval))

p.df.medage <- ggplot(df.medage, aes(x = GeneSet, y = log10pval, fill = significant)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#999999", "#8D6EFF"),   #E242AC #06945A
                    breaks = c("n", "y")) +
  theme_classic() 

p.df.medage <- p.df.medage + 
  aes(x = fct_inorder(GeneSet)) + xlab("Tissue") + ylab("-log 10 p-value")

p.df.medage <- p.df.medage + 
  theme(axis.text.x = element_text(color = "black", size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title  = element_text(size = 20, face = "bold"),
        legend.position = "none") +
  ylim(0, 8.5)

p.df.medage





# ----------------------- modern gene set 
# settings: 98 genes, exclude MHC, use FDR correction, date: 21/05/26

sum.mod <- fread("data/raw/FUMA/mod_FUMA_gene2func82298/summary.txt", h = T)

geneTab.mod <- read.table("data/raw/FUMA/mod_FUMA_gene2func82298/geneTable.txt", h = T)
geneIDs.mod <- read.table("data/raw/FUMA/mod_FUMA_gene2func82298/geneIDs.txt", h = T)

GS.mod <- fread("data/raw/FUMA/mod_FUMA_gene2func82298/GS.txt", h = T)

# for supplementary materials
# write_xlsx(GS.mod, "output/supp_mat_new/sup_mat_06c.xlsx")

GS.mod.GWAS <- GS.mod[GS.mod$Category == "GWAScatalog",]

# 3 with teeth things              


gtex.30.log.mod <- read.table("data/raw/FUMA/mod_FUMA_gene2func82298/gtex_v8_ts_general_avg_log2TPM_exp.txt", h = T)
gtex.30.norm.mod <- read.table("data/raw/FUMA/mod_FUMA_gene2func82298/gtex_v8_ts_general_avg_normTPM_exp.txt", h = T)
gtex.30.DEG.mod <- fread("data/raw/FUMA/mod_FUMA_gene2func82298/gtex_v8_ts_general_DEG.txt", h = T)


tissueEnrich30mod <- gtex.30.DEG.mod[gtex.30.DEG.mod$Category == "DEG.up"  & gtex.30.DEG.mod$adjP < 0.05, ]

# brain and muscle


# ..... continue working with "gtex.30.DEG.mod"



gtex30.mod.up <- subset(gtex.30.DEG.mod, grepl("DEG.up", gtex.30.DEG.mod$Category)) 
gtex30.mod.up$significant <- "n"
gtex30.mod.up$significant[gtex30.mod.up$adjP <= 0.05]  <- "y"

gtex30.mod.up$log10pval <- -log10(gtex30.mod.up$p)


df.mod <- data.frame(c(gtex30.mod.up[, c(2, 8, 9)]))    # select relevant columns

df.mod$GeneSet <- str_replace_all(df.mod$GeneSet, "_T", " t")
df.mod$GeneSet <- str_replace_all(df.mod$GeneSet, "_V", " v")
df.mod$GeneSet <- str_replace_all(df.mod$GeneSet, "_G", " g")
df.mod$GeneSet <- str_replace_all(df.mod$GeneSet, "_I", " i")
df.mod$GeneSet <- str_replace_all(df.mod$GeneSet, "_Uteri", "")

df.mod <- df.mod %>%                                    # order descending by mean
  arrange(desc(log10pval))

p.df.mod <- ggplot(df.mod, aes(x = GeneSet, y = log10pval, fill = significant)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#999999", "#49C9B7"),   #E242AC #06945A
                    breaks = c("n", "y")) +
  theme_classic() 

p.df.mod <- p.df.mod + 
  aes(x = fct_inorder(GeneSet)) + xlab("Tissue") + ylab("-log 10 p-value")

p.df.mod <- p.df.mod + 
  theme(axis.text.x = element_text(color = "black", size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title  = element_text(size = 20, face = "bold"),
        legend.position = "none") +
  ylim(0, 8.5)

p.df.mod







# ..... combine p.df.nv and p.df.v 
p.df.medage.mod <- plot_grid(p.df.medage, p.df.mod,
                      nrow = 2,
                      labels = c('B', 'C'),
                      label_x = .94,
                      label_size = 28) # 12 * 20

p.df.medage.mod


# ggsave(filename = "output/FUMA_tissueEnrich01_b.pdf", p.df.medage.mod,
#        width = 10, height = 5, units = "in", device= 'pdf')





# .....combine p.df.anc and p.df.medage and p.df.mod 

p.df.medagec <- p.df.medage + 
  theme(axis.title.x = element_blank())

p.df.ancc <- p.df.anc + 
  theme(axis.title.x = element_blank())


p.df.ancmedagemod <- plot_grid(p.df.ancc, p.df.medagec, p.df.mod,
                      nrow = 3,
                      labels = c('B', 'C', 'D'),
                      label_x = .94,
                      label_size = 28) # 12 * 20

p.df.ancmedagemod

# ggsave(filename = "output/FUMA_tissueEnrich02.pdf", p.df.ancmedagemod,
#       width = 10, height = 10, units = "in", device= 'pdf')


#### END OF SCRIPT ####



