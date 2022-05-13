############################
# FUMA analysis prep
############################

# prepare old and new gene list for submission (this mm and the mm the phylostratigraphy is based on have the same dimensions and genes, just saying :))
M1C_mm <- M1C_tai$mm
OTgenesPS <-M1C_mm[,1:2]

# new gene set
yPS <- c(4, 5, 6, 7, 12)
youngOT <- OTgenesPS[OTgenesPS$Phylostratum %in% yPS, ]  
youngOT$index <- c(1:21)
youngOT$GeneID <- toupper(youngOT$GeneID)

# old gene set
oPS <- c(1, 2, 3)
oldOT <- OTgenesPS[OTgenesPS$Phylostratum %in% oPS, ]  
oldOT$index <- c(1:121)
oldOT$GeneID <- toupper(oldOT$GeneID)




###############################
# FUMA Main OT analysis results 
###############################

# ----------------------- OLD gene set 
# settings: 121 genes, exclude MHC, use FDR correction, date: 21/07/01

summary_o <- fread(paste0(BASE, "analyses/FUMA/FUMA_oldOT61601/summary.txt"), h=T)

geneTab_o <- read.table(paste0(BASE, "analyses/FUMA/FUMA_oldOT61601/geneTable.txt"), h=T)
geneIDs_o <- read.table(paste0(BASE, "analyses/FUMA/FUMA_oldOT61601/geneIDs.txt"), h=T)

GS_o <- fread(paste0(BASE, "analyses/FUMA/FUMA_oldOT61601/GS.txt"), h=T)
GS_oGWAS <- GS_o[GS_o$Category=="GWAScatalog",]

# Systolic blood preasure: 1/2        
# Extreme obesity: 2/2              


gtex_54tissues_log_o <- read.table(paste0(BASE, "analyses/FUMA/FUMA_oldOT61601/gtex_v8_ts_avg_log2TPM_exp.txt"), h=T)
gtex_54tissues_norm_o <- read.table(paste0(BASE, "analyses/FUMA/FUMA_oldOT61601/gtex_v8_ts_avg_normTPM_exp.txt"), h=T)
gtex_54tissues_DEG_o <- fread(paste0(BASE, "analyses/FUMA/FUMA_oldOT61601/gtex_v8_ts_DEG.txt"), h=T)

gtex_30tissues_log_o <- read.table(paste0(BASE, "analyses/FUMA/FUMA_oldOT61601/gtex_v8_ts_general_avg_log2TPM_exp.txt"), h=T)
gtex_30tissues_norm_o <- read.table(paste0(BASE, "analyses/FUMA/FUMA_oldOT61601/gtex_v8_ts_general_avg_normTPM_exp.txt"), h=T)
gtex_30tissues_DEG_o <- fread(paste0(BASE, "analyses/FUMA/FUMA_oldOT61601/gtex_v8_ts_general_DEG.txt"), h=T)


tissueEnrich30old <- gtex_30tissues_DEG_o[gtex_30tissues_DEG_o$Category == "DEG.up"  & gtex_30tissues_DEG_o$adjP < 0.05, ]

# ..... continue working with "gtex_30tissues_DEG_o"

gtex30_Oup <- subset(gtex_30tissues_DEG_o, grepl("DEG.up", gtex_30tissues_DEG_o$Category)) 
gtex30_Oup$significant <- "n"
gtex30_Oup$significant[gtex30_Oup$adjP <= 0.05]  <- "y"

gtex30_Oup$log10pval <- -log10(gtex30_Oup$p)


df_o <- data.frame(c(gtex30_Oup[, c(2, 8, 9)]))

df_o$GeneSet <- str_replace_all(df_o$GeneSet, "_T", " t")
df_o$GeneSet <- str_replace_all(df_o$GeneSet, "_V", " v")
df_o$GeneSet <- str_replace_all(df_o$GeneSet, "_G", " g")
df_o$GeneSet <- str_replace_all(df_o$GeneSet, "_I", " i")
df_o$GeneSet <- str_replace_all(df_o$GeneSet, "_Uteri", "")

df_o <- df_o %>%           # order descending by mean
  arrange(desc(log10pval))

tplot_o <- ggplot(df_o, aes(x = GeneSet, y = log10pval, fill = significant)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values=c("#999999", "#49C9B7"),    #56B4E9 #F2B342
                    breaks=c("n", "y")) +
  theme_classic() 
tplot_o <- tplot_o + 
  aes(x = fct_inorder(GeneSet)) + xlab("Tissue") + ylab("-log 10 p-value")
tplot_o <- tplot_o + 
  theme(axis.text.x = element_text(color="black", size=12, angle=45, hjust = 1),
        axis.text.y = element_text(size=15),
        axis.title = element_text(size=20, face="bold"),
        legend.position = "none") 

tplot_o




# ----------------------- YOUNG gene set 
# settings: 20 genes, exclude MHC, use FDR correction, date: 21/05/26

summary_y <- fread(paste0(BASE, "analyses/FUMA/FUMA_youngOT61602/summary.txt"), h=T)

geneTab_y <- read.table(paste0(BASE, "analyses/FUMA/FUMA_youngOT61602/geneTable.txt"), h=T)
geneIDs_y <- read.table(paste0(BASE, "analyses/FUMA/FUMA_youngOT61602/geneIDs.txt"), h=T)

GS_y <- fread(paste0(BASE, "analyses/FUMA/FUMA_youngOT61602/GS.txt"), h=T)
GS_yGWAS <- GS_y[GS_y$Category=="GWAScatalog",]

# not enriched in any GWAS



gtex_54tissues_log_y <- read.table(paste0(BASE, "analyses/FUMA/FUMA_youngOT61602/gtex_v8_ts_avg_log2TPM_exp.txt"), h=T)
gtex_54tissues_norm_y <- read.table(paste0(BASE, "analyses/FUMA/FUMA_youngOT61602/gtex_v8_ts_avg_normTPM_exp.txt"), h=T)
gtex_54tissues_DEG_y <- fread(paste0(BASE, "analyses/FUMA/FUMA_youngOT61602/gtex_v8_ts_DEG.txt"), h=T)


gtex_30tissues_log_y <- read.table(paste0(BASE, "analyses/FUMA/FUMA_youngOT61602/gtex_v8_ts_general_avg_log2TPM_exp.txt"), h=T)
gtex_30tissues_norm_y <- read.table(paste0(BASE, "analyses/FUMA/FUMA_youngOT61602/gtex_v8_ts_general_avg_normTPM_exp.txt"), h=T)
gtex_30tissues_DEG_y <- fread(paste0(BASE, "analyses/FUMA/FUMA_youngOT61602/gtex_v8_ts_general_DEG.txt"), h=T)


tissueEnrich30young <- gtex_30tissues_DEG_y[gtex_30tissues_DEG_y$Category == "DEG.up" & gtex_30tissues_DEG_y$adjP < 0.05, ]


# ..... continue working with "gtex_30tissues_DEG_y"

gtex30_Yup <- subset(gtex_30tissues_DEG_y, grepl("DEG.up", gtex_30tissues_DEG_y$Category)) 
gtex30_Yup$significant <- "n"
gtex30_Yup$significant[gtex30_Yup$adjP <= 0.05]  <- "y"

gtex30_Yup$log10pval <- -log10(gtex30_Yup$p)


df_y <- data.frame(c(gtex30_Yup[, c(2, 8, 9)]))

df_y$GeneSet <- str_replace_all(df_y$GeneSet, "_T", " t")
df_y$GeneSet <- str_replace_all(df_y$GeneSet, "_V", " v")
df_y$GeneSet <- str_replace_all(df_y$GeneSet, "_G", " g")
df_y$GeneSet <- str_replace_all(df_y$GeneSet, "_I", " i")
df_y$GeneSet <- str_replace_all(df_y$GeneSet, "_Uteri", "")

df_y <- df_y %>%           # order descending by mean
  arrange(desc(log10pval))

tplot_y <- ggplot(df_y, aes(x = GeneSet, y = log10pval, fill = significant)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values=c("#999999", "#8D6EFF"),   #E242AC #06945A
                    breaks=c("n", "y")) +
  theme_classic() 
tplot_y <- tplot_y + 
  aes(x = fct_inorder(GeneSet)) + xlab("Tissue") + ylab("-log 10 p-value")
tplot_y <- tplot_y + 
  theme(axis.text.x = element_text(color="black", size=12, angle=45, hjust = 1),
        axis.text.y = element_text(size=15),
        axis.title = element_text(size=20, face="bold"),
        legend.position = "none") +
  ylim(0, 7.08)

tplot_y




# ..... combine tplot_o with tplot_y from above
tplot_oy <- plot_grid(tplot_o, tplot_y,
                      nrow = 2,
                      labels = c('B', 'C'),
                      label_x = .94,
                      label_size = 28) # 12 * 20

tplot_oy


ggsave(filename = paste0(BASE, "RProject/OT-bd-analyses_loc/output/tissueEnrich_d.pdf"), tplot_oy,
       width = 10, height = 10, units = "in", device='pdf')




###################################
# FUMA Specificity analysis results 
###################################

## ................ Insulin

# ----------------------- OLD gene set 
# settings: 111 genes, exclude MHC, use FDR correction, date: 21/08/24

sum_oINS <- fread(paste0(BASE, "analyses/FUMA/FUMA_oldINS64577/summary.txt"), h=T)

geneTab_oINS <- read.table(paste0(BASE, "analyses/FUMA/FUMA_oldINS64577/geneTable.txt"), h=T)
geneIDs_oINS <- read.table(paste0(BASE, "analyses/FUMA/FUMA_oldINS64577/geneIDs.txt"), h=T)

GS_oINS <- fread(paste0(BASE, "analyses/FUMA/FUMA_oldINS64577/GS.txt"), h=T)
GS_oGWASins <- GS_oINS[GS_oINS$Category=="GWAScatalog",]

# no GWAS enrichment             


gtex_30tissues_log_oINS <- read.table(paste0(BASE, "analyses/FUMA/FUMA_oldINS64577/gtex_v8_ts_general_avg_log2TPM_exp.txt"), h=T)
gtex_30tissues_norm_oINS <- read.table(paste0(BASE, "analyses/FUMA/FUMA_oldINS64577/gtex_v8_ts_general_avg_normTPM_exp.txt"), h=T)
gtex_30tissues_DEG_oINS <- fread(paste0(BASE, "analyses/FUMA/FUMA_oldINS64577/gtex_v8_ts_general_DEG.txt"), h=T)


tissueEnrich30oldINS <- gtex_30tissues_DEG_oINS[gtex_30tissues_DEG_oINS$Category == "DEG.up"  & gtex_30tissues_DEG_oINS$adjP < 0.05, ]

# ..... continue working with "gtex_30tissues_DEG_oINS"

gtex30_OupINS <- subset(gtex_30tissues_DEG_oINS, grepl("DEG.up", gtex_30tissues_DEG_oINS$Category)) 
gtex30_OupINS$significant <- "n"
gtex30_OupINS$significant[gtex30_OupINS$adjP <= 0.05]  <- "y"

gtex30_OupINS$log10pval <- -log10(gtex30_OupINS$p)


df_oINS <- data.frame(c(gtex30_OupINS[, c(2, 8, 9)]))  # fetch columns of interest

df_oINS$GeneSet <- str_replace_all(df_oINS$GeneSet, "_T", " t")
df_oINS$GeneSet <- str_replace_all(df_oINS$GeneSet, "_V", " v")
df_oINS$GeneSet <- str_replace_all(df_oINS$GeneSet, "_G", " g")
df_oINS$GeneSet <- str_replace_all(df_oINS$GeneSet, "_I", " i")
df_oINS$GeneSet <- str_replace_all(df_oINS$GeneSet, "_Uteri", "")

df_oINS <- df_oINS %>%           # order descending by mean
  arrange(desc(log10pval))

tplot_oINS <- ggplot(df_oINS, aes(x = GeneSet, y = log10pval, fill = significant)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values=c("#999999", "#56B4E9"),
                    breaks=c("n", "y")) +
  theme_classic() 
tplot_oINS <- tplot_oINS + 
  aes(x = fct_inorder(GeneSet)) + xlab("Tissue") + ylab("-log 10 p-value")
tplot_oINS <- tplot_oINS + 
  theme(axis.text.x = element_text(color="black", size=12, angle=45, hjust = 1),
        axis.text.y = element_text(size=15),
        axis.title = element_text(size=20, face="bold"),
        legend.position = "none") 

tplot_oINS



# ----------------------- YOUNG gene set 
# settings: 15 genes, exclude MHC, use FDR correction, date: 21/08/24

sum_yINS <- fread(paste0(BASE, "analyses/FUMA/FUMA_youngINS64578/summary.txt"), h=T)

geneTab_yINS <- read.table(paste0(BASE, "analyses/FUMA/FUMA_youngINS64578/geneTable.txt"), h=T)
geneIDs_yINS <- read.table(paste0(BASE, "analyses/FUMA/FUMA_youngINS64578/geneIDs.txt"), h=T)

GS_yINS <- fread(paste0(BASE, "analyses/FUMA/FUMA_youngINS64578/GS.txt"), h=T)
GS_yGWASins <- GS_yINS[GS_yINS$Category=="GWAScatalog",]

# not enriched in any GWAS


gtex_30tissues_log_yINS <- read.table(paste0(BASE, "analyses/FUMA/FUMA_youngINS64578/gtex_v8_ts_general_avg_log2TPM_exp.txt"), h=T)
gtex_30tissues_norm_yINS <- read.table(paste0(BASE, "analyses/FUMA/FUMA_youngINS64578/gtex_v8_ts_general_avg_normTPM_exp.txt"), h=T)
gtex_30tissues_DEG_yINS <- fread(paste0(BASE, "analyses/FUMA/FUMA_youngINS64578/gtex_v8_ts_general_DEG.txt"), h=T)


tissueEnrich30youngINS <- gtex_30tissues_DEG_yINS[gtex_30tissues_DEG_yINS$Category == "DEG.up" & gtex_30tissues_DEG_yINS$adjP < 0.05, ]


# ..... continue working with "gtex_30tissues_DEG_yINS"

gtex30_YupINS <- subset(gtex_30tissues_DEG_yINS, grepl("DEG.up", gtex_30tissues_DEG_yINS$Category)) 
gtex30_YupINS$significant <- "n"
gtex30_YupINS$significant[gtex30_YupINS$adjP <= 0.05]  <- "y"

gtex30_YupINS$log10pval <- -log10(gtex30_YupINS$p)


df_yINS <- data.frame(c(gtex30_YupINS[, c(2, 8, 9)]))

df_yINS$GeneSet <- str_replace_all(df_yINS$GeneSet, "_T", " t")
df_yINS$GeneSet <- str_replace_all(df_yINS$GeneSet, "_V", " v")
df_yINS$GeneSet <- str_replace_all(df_yINS$GeneSet, "_G", " g")
df_yINS$GeneSet <- str_replace_all(df_yINS$GeneSet, "_I", " i")
df_yINS$GeneSet <- str_replace_all(df_yINS$GeneSet, "_Uteri", "")

df_yINS <- df_yINS %>%           # order descending by mean
  arrange(desc(log10pval))

tplot_yINS <- ggplot(df_yINS, aes(x = GeneSet, y = log10pval, fill = significant)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values=c("#999999", "#E242AC"),
                    breaks=c("n", "y")) +
  theme_classic() 
tplot_yINS <- tplot_yINS + 
  aes(x = fct_inorder(GeneSet)) + xlab("Tissue") + ylab("-log 10 p-value")
tplot_yINS <- tplot_yINS + 
  theme(axis.text.x = element_text(color="black", size=12, angle=45, hjust = 1),
        axis.text.y = element_text(size=15),
        axis.title = element_text(size=20, face="bold"),
        legend.position = "none") +
  ylim(0, 5.23)

tplot_yINS




# ..... combine tplot_oINS with tplot_yINS from above
tplot_oyINS <- plot_grid(tplot_oINS, tplot_yINS,
                         nrow = 2,
                         labels = c('B', 'C'),
                         label_x = .94,
                         label_size = 28) # 12 * 20

tplot_oyINS


ggsave(filename = paste0(BASE, "RProject/OT-bd-analyses_loc/output/tissueEnrich_INS.pdf"), tplot_oyINS,
       width = 10, height = 10, units = "in", device='pdf')






## ................ Estrogen

# ----------------------- OLD gene set 
# settings: 81 genes (one gene NA), exclude MHC, use FDR correction, date: 21/08/24

sum_oES <- fread(paste0(BASE, "analyses/FUMA/FUMA_oldES64579/summary.txt"), h=T)

geneTab_oES <- read.table(paste0(BASE, "analyses/FUMA/FUMA_oldES64579/geneTable.txt"), h=T)
geneIDs_oES <- read.table(paste0(BASE, "analyses/FUMA/FUMA_oldES64579/geneIDs.txt"), h=T)

GS_oES <- fread(paste0(BASE, "analyses/FUMA/FUMA_oldES64579/GS.txt"), h=T)
GS_oGWASes <- GS_oES[GS_oES$Category=="GWAScatalog",]

# no GWAS enrichment             


gtex_30tissues_log_oES <- read.table(paste0(BASE, "analyses/FUMA/FUMA_oldES64579/gtex_v8_ts_general_avg_log2TPM_exp.txt"), h=T)
gtex_30tissues_norm_oES <- read.table(paste0(BASE, "analyses/FUMA/FUMA_oldES64579/gtex_v8_ts_general_avg_normTPM_exp.txt"), h=T)
gtex_30tissues_DEG_oES <- fread(paste0(BASE, "analyses/FUMA/FUMA_oldES64579/gtex_v8_ts_general_DEG.txt"), h=T)


tissueEnrich30oldES <- gtex_30tissues_DEG_oES[gtex_30tissues_DEG_oES$Category == "DEG.up"  & gtex_30tissues_DEG_oES$adjP < 0.05, ]

# ..... continue working with "gtex_30tissues_DEG_oES"

gtex30_OupES <- subset(gtex_30tissues_DEG_oES, grepl("DEG.up", gtex_30tissues_DEG_oES$Category)) 
gtex30_OupES$significant <- "n"
gtex30_OupES$significant[gtex30_OupES$adjP <= 0.05]  <- "y"

gtex30_OupES$log10pval <- -log10(gtex30_OupES$p)


df_oES <- data.frame(c(gtex30_OupES[, c(2, 8, 9)]))  # fetch columns of interest

df_oES$GeneSet <- str_replace_all(df_oES$GeneSet, "_T", " t")
df_oES$GeneSet <- str_replace_all(df_oES$GeneSet, "_V", " v")
df_oES$GeneSet <- str_replace_all(df_oES$GeneSet, "_G", " g")
df_oES$GeneSet <- str_replace_all(df_oES$GeneSet, "_I", " i")
df_oES$GeneSet <- str_replace_all(df_oES$GeneSet, "_Uteri", "")

df_oES <- df_oES %>%           # order descending by mean
  arrange(desc(log10pval))

tplot_oES <- ggplot(df_oES, aes(x = GeneSet, y = log10pval, fill = significant)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values=c("#999999", "#56B4E9"),
                    breaks=c("n", "y")) +
  theme_classic() 
tplot_oES <- tplot_oES + 
  aes(x = fct_inorder(GeneSet)) + xlab("Tissue") + ylab("-log 10 p-value")
tplot_oES <- tplot_oES + 
  theme(axis.text.x = element_text(color="black", size=12, angle=45, hjust = 1),
        axis.text.y = element_text(size=15),
        axis.title = element_text(size=20, face="bold"),
        legend.position = "none") +
  ylim(0, 6.52)

tplot_oES



# ----------------------- YOUNG gene set 
# settings: 15 genes, exclude MHC, use FDR correction, date: 21/08/24

sum_yES <- fread(paste0(BASE, "analyses/FUMA/FUMA_youngES64580/summary.txt"), h=T)

geneTab_yES <- read.table(paste0(BASE, "analyses/FUMA/FUMA_youngES64580/geneTable.txt"), h=T)
geneIDs_yES <- read.table(paste0(BASE, "analyses/FUMA/FUMA_youngES64580/geneIDs.txt"), h=T)

GS_yES <- fread(paste0(BASE, "analyses/FUMA/FUMA_youngES64580/GS.txt"), h=T)
GS_yGWASes <- GS_yES[GS_yES$Category=="GWAScatalog",]

# not enriched in any GWAS


gtex_30tissues_log_yES <- read.table(paste0(BASE, "analyses/FUMA/FUMA_youngES64580/gtex_v8_ts_general_avg_log2TPM_exp.txt"), h=T)
gtex_30tissues_norm_yES <- read.table(paste0(BASE, "analyses/FUMA/FUMA_youngES64580/gtex_v8_ts_general_avg_normTPM_exp.txt"), h=T)
gtex_30tissues_DEG_yES <- fread(paste0(BASE, "analyses/FUMA/FUMA_youngES64580/gtex_v8_ts_general_DEG.txt"), h=T)


tissueEnrich30youngES <- gtex_30tissues_DEG_yES[gtex_30tissues_DEG_yES$Category == "DEG.up" & gtex_30tissues_DEG_yES$adjP < 0.05, ]


# ..... continue working with "gtex_30tissues_DEG_yES"

gtex30_YupES <- subset(gtex_30tissues_DEG_yES, grepl("DEG.up", gtex_30tissues_DEG_yES$Category)) 
gtex30_YupES$significant <- "n"
gtex30_YupES$significant[gtex30_YupES$adjP <= 0.05]  <- "y"

gtex30_YupES$log10pval <- -log10(gtex30_YupES$p)


df_yES <- data.frame(c(gtex30_YupES[, c(2, 8, 9)]))

df_yES$GeneSet <- str_replace_all(df_yES$GeneSet, "_T", " t")
df_yES$GeneSet <- str_replace_all(df_yES$GeneSet, "_V", " v")
df_yES$GeneSet <- str_replace_all(df_yES$GeneSet, "_G", " g")
df_yES$GeneSet <- str_replace_all(df_yES$GeneSet, "_I", " i")
df_yES$GeneSet <- str_replace_all(df_yES$GeneSet, "_Uteri", "")

df_yES <- df_yES %>%           # order descending by mean
  arrange(desc(log10pval))

tplot_yES <- ggplot(df_yES, aes(x = GeneSet, y = log10pval, fill = significant)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values=c("#999999", "#E242AC"),
                    breaks=c("n", "y")) +
  theme_classic() 

tplot_yES <- tplot_yES + 
  aes(x = fct_inorder(GeneSet)) + xlab("Tissue") + ylab("-log 10 p-value")

tplot_yES <- tplot_yES + 
  theme(axis.text.x = element_text(color="black", size=12, angle=45, hjust = 1),
        axis.text.y = element_text(size=15),
        axis.title = element_text(size=20, face="bold"),
        legend.position = "none") 

tplot_yES




# ..... combine tplot_oES with tplot_yES from above
tplot_oyES <- plot_grid(tplot_oES, tplot_yES,
                         nrow = 2,
                         labels = c('B', 'C'),
                         label_x = .94,
                         label_size = 28) # 12 * 20

tplot_oyES


ggsave(filename = paste0(BASE, "RProject/OT-bd-analyses_loc/output/tissueEnrich_ES.pdf"), tplot_oyES,
       width = 10, height = 10, units = "in", device='pdf')






## ................ Relaxin

# ----------------------- OLD gene set 
# settings: 80 genes, exclude MHC, use FDR correction, date: 21/08/24

sum_oRX <- fread(paste0(BASE, "analyses/FUMA/FUMA_oldRX64587/summary.txt"), h=T)

geneTab_oRX <- read.table(paste0(BASE, "analyses/FUMA/FUMA_oldRX64587/geneTable.txt"), h=T)
geneIDs_oRX <- read.table(paste0(BASE, "analyses/FUMA/FUMA_oldRX64587/geneIDs.txt"), h=T)

GS_oRX <- fread(paste0(BASE, "analyses/FUMA/FUMA_oldRX64587/GS.txt"), h=T)
GS_oGWASrx <- GS_oRX[GS_oRX$Category=="GWAScatalog",]

# quite a few enrichments            


gtex_30tissues_log_oRX <- read.table(paste0(BASE, "analyses/FUMA/FUMA_oldRX64587/gtex_v8_ts_general_avg_log2TPM_exp.txt"), h=T)
gtex_30tissues_norm_oRX <- read.table(paste0(BASE, "analyses/FUMA/FUMA_oldRX64587/gtex_v8_ts_general_avg_normTPM_exp.txt"), h=T)
gtex_30tissues_DEG_oRX <- fread(paste0(BASE, "analyses/FUMA/FUMA_oldRX64587/gtex_v8_ts_general_DEG.txt"), h=T)


tissueEnrich30oldRX <- gtex_30tissues_DEG_oRX[gtex_30tissues_DEG_oRX$Category == "DEG.up" & gtex_30tissues_DEG_oRX$adjP < 0.05, ]

# ..... continue working with "gtex_30tissues_DEG_oRX"

gtex30_OupRX <- subset(gtex_30tissues_DEG_oRX, grepl("DEG.up", gtex_30tissues_DEG_oRX$Category)) 
gtex30_OupRX$significant <- "n"
gtex30_OupRX$significant[gtex30_OupRX$adjP <= 0.05]  <- "y"

gtex30_OupRX$log10pval <- -log10(gtex30_OupRX$p)


df_oRX <- data.frame(c(gtex30_OupRX[, c(2, 8, 9)]))  # fetch columns of interest

df_oRX$GeneSet <- str_replace_all(df_oRX$GeneSet, "_T", " t")
df_oRX$GeneSet <- str_replace_all(df_oRX$GeneSet, "_V", " v")
df_oRX$GeneSet <- str_replace_all(df_oRX$GeneSet, "_G", " g")
df_oRX$GeneSet <- str_replace_all(df_oRX$GeneSet, "_I", " i")
df_oRX$GeneSet <- str_replace_all(df_oRX$GeneSet, "_Uteri", "")

df_oRX <- df_oRX %>%           # order descending by mean
  arrange(desc(log10pval))

tplot_oRX <- ggplot(df_oRX, aes(x = GeneSet, y = log10pval, fill = significant)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values=c("#999999", "#56B4E9"),
                    breaks=c("n", "y")) +
  theme_classic() 

tplot_oRX <- tplot_oRX + 
  aes(x = fct_inorder(GeneSet)) + xlab("Tissue") + ylab("-log 10 p-value")

tplot_oRX <- tplot_oRX + 
  theme(axis.text.x = element_text(color="black", size=12, angle=45, hjust = 1),
        axis.text.y = element_text(size=15),
        axis.title = element_text(size=20, face="bold"),
        legend.position = "none") #+
 # ylim(0, 6.52)

tplot_oRX



# ----------------------- YOUNG gene set 
# settings: 37 genes, exclude MHC, use FDR correction, date: 21/08/24

sum_yRX <- fread(paste0(BASE, "analyses/FUMA/FUMA_youngRX64588/summary.txt"), h=T)

geneTab_yRX <- read.table(paste0(BASE, "analyses/FUMA/FUMA_youngRX64588/geneTable.txt"), h=T)
geneIDs_yRX <- read.table(paste0(BASE, "analyses/FUMA/FUMA_youngRX64588/geneIDs.txt"), h=T)

GS_yRX <- fread(paste0(BASE, "analyses/FUMA/FUMA_youngRX64588/GS.txt"), h=T)
GS_yGWASrx <- GS_yRX[GS_yRX$Category=="GWAScatalog",]

# not enriched in any GWAS


gtex_30tissues_log_yRX <- read.table(paste0(BASE, "analyses/FUMA/FUMA_youngRX64588/gtex_v8_ts_general_avg_log2TPM_exp.txt"), h=T)
gtex_30tissues_norm_yRX <- read.table(paste0(BASE, "analyses/FUMA/FUMA_youngRX64588/gtex_v8_ts_general_avg_normTPM_exp.txt"), h=T)
gtex_30tissues_DEG_yRX <- fread(paste0(BASE, "analyses/FUMA/FUMA_youngRX64588/gtex_v8_ts_general_DEG.txt"), h=T)


tissueEnrich30youngRX <- gtex_30tissues_DEG_yRX[gtex_30tissues_DEG_yRX$Category == "DEG.up" & gtex_30tissues_DEG_yRX$adjP < 0.05, ]


# ..... continue working with "gtex_30tissues_DEG_yRX"

gtex30_YupRX <- subset(gtex_30tissues_DEG_yRX, grepl("DEG.up", gtex_30tissues_DEG_yRX$Category)) 
gtex30_YupRX$significant <- "n"
gtex30_YupRX$significant[gtex30_YupRX$adjP <= 0.05]  <- "y"

gtex30_YupRX$log10pval <- -log10(gtex30_YupRX$p)


df_yRX <- data.frame(c(gtex30_YupRX[, c(2, 8, 9)]))

df_yRX$GeneSet <- str_replace_all(df_yRX$GeneSet, "_T", " t")
df_yRX$GeneSet <- str_replace_all(df_yRX$GeneSet, "_V", " v")
df_yRX$GeneSet <- str_replace_all(df_yRX$GeneSet, "_G", " g")
df_yRX$GeneSet <- str_replace_all(df_yRX$GeneSet, "_I", " i")
df_yRX$GeneSet <- str_replace_all(df_yRX$GeneSet, "_Uteri", "")

df_yRX <- df_yRX %>%           # order descending by mean
  arrange(desc(log10pval))

tplot_yRX <- ggplot(df_yRX, aes(x = GeneSet, y = log10pval, fill = significant)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values=c("#999999", "#E242AC"),
                    breaks=c("n", "y")) +
  theme_classic() 

tplot_yRX <- tplot_yRX + 
  aes(x = fct_inorder(GeneSet)) + xlab("Tissue") + ylab("-log 10 p-value")

tplot_yRX <- tplot_yRX + 
  theme(axis.text.x = element_text(color="black", size=12, angle=45, hjust = 1),
        axis.text.y = element_text(size=15),
        axis.title = element_text(size=20, face="bold"),
        legend.position = "none") +
  ylim(0, 1.93)

tplot_yRX




# ..... combine tplot_oRX with tplot_yRX from above
tplot_oyRX <- plot_grid(tplot_oRX, tplot_yRX,
                        nrow = 2,
                        labels = c('B', 'C'),
                        label_x = .94,
                        label_size = 28) # 12 * 20

tplot_oyRX


ggsave(filename = paste0(BASE, "RProject/OT-bd-analyses_loc/output/tissueEnrich_RX.pdf"), tplot_oyRX,
       width = 10, height = 10, units = "in", device='pdf')




## grid plot with the three specificity pathways

# alternative ylim: 6.52

tplot_oINS.2 <- tplot_oINS + ylim(0, 7.08) + theme(axis.text.x = element_text(size=8), 
                                                   axis.title = element_text(size=15, face="bold"))

tplot_yINS.2 <- tplot_yINS + ylim(0, 7.08) + theme(axis.text.x = element_text(size=8),
                                                   axis.title = element_text(size=15, face="bold"))

tplot_oES.2 <- tplot_oES + ylim(0, 7.08) + theme(axis.text.x = element_text(size=8),
                                                 axis.title = element_text(size=15, face="bold"))

tplot_yES.2 <- tplot_yES + ylim(0, 7.08) + theme(axis.text.x = element_text(size=8),
                                                 axis.title = element_text(size=15, face="bold"))

tplot_oRX.2 <- tplot_oRX + ylim(0, 7.08) + theme(axis.text.x = element_text(size=8),
                                                 axis.title = element_text(size=15, face="bold"))

tplot_yRX.2 <- tplot_yRX + ylim(0, 7.08) + theme(axis.text.x = element_text(size=8),
                                                 axis.title = element_text(size=15, face="bold"))


tplot_ALLspec <- plot_grid(tplot_oINS.2, tplot_oES.2, tplot_oRX.2, tplot_yINS.2, tplot_yES.2, tplot_yRX.2,
                           ncol = 3,
                           nrow = 2,
                           labels = c('Insulin signaling pathway', 'Estrogen signaling pathway', 'Relaxin signaling pathway', '', '', ''),
                           label_x = 0.1,
                           label_size = 12) 

tplot_ALLspec



## grid plot with all four pathways

tplot_o.2 <- tplot_o + ylim(0, 7.08) + theme(axis.text.x = element_text(size=8),
                                                axis.title = element_text(size=15, face="bold"))

tplot_y.2 <- tplot_y + ylim(0, 7.08) + theme(axis.text.x = element_text(size=8),
                                                axis.title = element_text(size=15, face="bold"))

tplot_ALLspecOT <- plot_grid(tplot_o.2, tplot_oINS.2, tplot_oES.2, tplot_oRX.2, tplot_y.2, tplot_yINS.2, tplot_yES.2, tplot_yRX.2,
                             ncol = 4,
                             nrow = 2,
                             labels = c('Oxytocin signaling pathway', 'Insulin signaling pathway', 'Estrogen signaling pathway', 'Relaxin signaling pathway', '', '', '', ''),
                             label_x = 0.1,
                             label_size = 12) 
tplot_ALLspecOT


ggsave(filename = paste0(BASE, "RProject/OT-bd-analyses_loc/output/tissueEnrich_OTINSESRX.pdf"), tplot_ALLspecOT,
       width = 18, height = 10, units = "in", device='pdf')














