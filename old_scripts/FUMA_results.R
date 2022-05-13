############################
# FUMA analysis results 
############################

# ----------------------- OLD gene set 
# settings: 118 genes, exclude MHC, use FDR correction, date: 21/05/26

summary_o <- fread(paste0(BASE, "analyses/FUMA/FUMA_59444_oldOT_exMHC_FDR/summary.txt"), h=T)

geneTab_o <- read.table(paste0(BASE, "analyses/FUMA/FUMA_59444_oldOT_exMHC_FDR/geneTable.txt"), h=T)
geneIDs_o <- read.table(paste0(BASE, "analyses/FUMA/FUMA_59444_oldOT_exMHC_FDR/geneIDs.txt"), h=T)

GS_o <- fread(paste0(BASE, "analyses/FUMA/FUMA_59444_oldOT_exMHC_FDR/GS.txt"), h=T)
GS_oGWAS <- GS_o[GS_o$Category=="GWAScatalog",]

# Bodymass and weight: 4/9         
# Teeth/dentures: 2/9              
# Blood pressure: 1/9               
# Pubertal anthropometrics: 1/9
# Hand grip strength: 1/9               


gtex_54tissues_log_o <- read.table(paste0(BASE, "analyses/FUMA/FUMA_59444_oldOT_exMHC_FDR/gtex_v8_ts_avg_log2TPM_exp.txt"), h=T)
gtex_54tissues_norm_o <- read.table(paste0(BASE, "analyses/FUMA/FUMA_59444_oldOT_exMHC_FDR/gtex_v8_ts_avg_normTPM_exp.txt"), h=T)
gtex_54tissues_DEG_o <- fread(paste0(BASE, "analyses/FUMA/FUMA_59444_oldOT_exMHC_FDR/gtex_v8_ts_DEG.txt"), h=T)

gtex_30tissues_log_o <- read.table(paste0(BASE, "analyses/FUMA/FUMA_59444_oldOT_exMHC_FDR/gtex_v8_ts_general_avg_log2TPM_exp.txt"), h=T)
gtex_30tissues_norm_o <- read.table(paste0(BASE, "analyses/FUMA/FUMA_59444_oldOT_exMHC_FDR/gtex_v8_ts_general_avg_normTPM_exp.txt"), h=T)
gtex_30tissues_DEG_o <- fread(paste0(BASE, "analyses/FUMA/FUMA_59444_oldOT_exMHC_FDR/gtex_v8_ts_general_DEG.txt"), h=T)



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
  scale_fill_manual(values=c("#999999", "#56B4E9"),
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

summary_y <- fread(paste0(BASE, "analyses/FUMA/FUMA_59445_youngOT_exMHC_FDR/summary.txt"), h=T)

geneTab_y <- read.table(paste0(BASE, "analyses/FUMA/FUMA_59445_youngOT_exMHC_FDR/geneTable.txt"), h=T)
geneIDs_y <- read.table(paste0(BASE, "analyses/FUMA/FUMA_59445_youngOT_exMHC_FDR/geneIDs.txt"), h=T)

GS_y <- fread(paste0(BASE, "analyses/FUMA/FUMA_59445_youngOT_exMHC_FDR/GS.txt"), h=T)

gtex_54tissues_log_y <- read.table(paste0(BASE, "analyses/FUMA/FUMA_59445_youngOT_exMHC_FDR/gtex_v8_ts_avg_log2TPM_exp.txt"), h=T)
gtex_54tissues_norm_y <- read.table(paste0(BASE, "analyses/FUMA/FUMA_59445_youngOT_exMHC_FDR/gtex_v8_ts_avg_normTPM_exp.txt"), h=T)
gtex_54tissues_DEG_y <- fread(paste0(BASE, "analyses/FUMA/FUMA_59445_youngOT_exMHC_FDR/gtex_v8_ts_DEG.txt"), h=T)


gtex_30tissues_log_y <- read.table(paste0(BASE, "analyses/FUMA/FUMA_59445_youngOT_exMHC_FDR/gtex_v8_ts_general_avg_log2TPM_exp.txt"), h=T)
gtex_30tissues_norm_y <- read.table(paste0(BASE, "analyses/FUMA/FUMA_59445_youngOT_exMHC_FDR/gtex_v8_ts_general_avg_normTPM_exp.txt"), h=T)
gtex_30tissues_DEG_y <- fread(paste0(BASE, "analyses/FUMA/FUMA_59445_youngOT_exMHC_FDR/gtex_v8_ts_general_DEG.txt"), h=T)



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
  scale_fill_manual(values=c("#999999", "#E242AC"),
                    breaks=c("n", "y")) +
  theme_classic() 
tplot_y <- tplot_y + 
  aes(x = fct_inorder(GeneSet)) + xlab("Tissue") + ylab("-log 10 p-value")
tplot_y <- tplot_y + 
  theme(axis.text.x = element_text(color="black", size=12, angle=45, hjust = 1),
        axis.text.y = element_text(size=15),
        axis.title = element_text(size=20, face="bold"),
        legend.position = "none") 

tplot_y




# ..... combine tplot_o with tplot_y from above
tplot_oy <- plot_grid(tplot_o, tplot_y,
                      nrow = 2,
                      labels = c('B', 'C'),
                      label_x = .94,
                      label_size = 28) # 12 * 20

tplot_oy


ggsave(filename = paste0(BASE, "RProject/OT-bd-analyses_loc/output/tissueEnrich.pdf"), tplot_oy,
       width = 10, height = 10, units = "in", device='pdf')








                       