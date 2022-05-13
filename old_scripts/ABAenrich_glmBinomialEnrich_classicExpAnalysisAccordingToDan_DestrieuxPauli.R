

########################################################
########################################################
########################################################
# load atlases
# INFO: Parcellations on the AHBA data were performed with the abagen toolbox (https://abagen.readthedocs.io/en/stable/index.html).
#       Atlases were fetched via nilearn (http://nilearn.github.io/modules/generated/nilearn.datasets.fetch_atlas_destrieux_2009.html#nilearn.datasets.fetch_atlas_destrieux_2009, 
#                                         http://nilearn.github.io/modules/generated/nilearn.datasets.fetch_atlas_pauli_2017.html#nilearn.datasets.fetch_atlas_pauli_2017)
#       View python/jupyter lab script "abagen_destrieux_pauli.ipynb" for script and detailed information. 

# ABAEnrichment: https://bioconductor.org/packages/release/bioc/vignettes/ABAEnrichment/inst/doc/ABAEnrichment.html#schematic-1-hypergeometric-test-and-fwer-calculation

# Papers: Pauli et al., 2018, https://www.nature.com/articles/sdata201863
#         Destrieux et al., 2010, https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2937159/
#         Abraham et al., 2014, https://www.frontiersin.org/articles/10.3389/fninf.2014.00014/full
#         Grote et al., 2016, https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5048072/ 



AHBAdestrieux <- read_csv(paste0(BASE, "/data/processed/exp_AHBAdestrieux.csv"))

# annotate with region labels (file is downloaded together with the atlas map/image in Jupyter Lab)
destrieux_labels <- read_csv("C:/Users/alina/nilearn_data/destrieux_2009/destrieux2009_rois_labels_lateralized.csv")
destrieux_labels <- destrieux_labels[ ! destrieux_labels$index %in% c(0, 42, 117), ] # remove background, and medial left and right medial wall because these are missing in the expression data (rows 42 and 117)

AHBAdestrieux_a <- dplyr::full_join(destrieux_labels, AHBAdestrieux, by=c("index"="label"))

AHBAdestrieux_L <- AHBAdestrieux_a[grepl("^L", AHBAdestrieux_a$name), ]   # keep only left hemisphere


## Pauli subcortical

AHBApauli <- read_csv(paste0(BASE, "/data/processed/exp_AHBA_sc_pauli.csv"))
dim(AHBApauli)
AHBApauli$label <- c(0:15)

# annotate with region labels (file is downloaded together with the atlas map/image in Jupyter Lab)
pauli_labels <- read.table("C:/Users/alina/nilearn_data/pauli_2017/labels.txt")
pauli_labels$node_info <-  c("Putamen", "Caudate nucleus", "Nucleus accumbens", "Extended amygdala", "Globus pallidus, external", "Globus pallidus, internal", "Substantia nigra, pars compacta", 
                             "Red nucleus", "Substantia nigra, pars reticulata", "Parabrachial pigmented nucleus", "Ventral tegmental area", "Ventral pallidum", "Habenular nuclei", "Hypothalamus", 
                             "Mammillary nucleus", "Subthalamic nucleus")

AHBApauli_a <- dplyr::full_join(pauli_labels, AHBApauli, by=c("V1"="label"))
names(AHBApauli_a)[names(AHBApauli_a) == "V1"] <- "ids"
names(AHBApauli_a)[names(AHBApauli_a) == "V2"] <- "regions"




# --------------------------------------------------- #


# vector with younger phylostrata
yPS <- c("4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19")

# load phylo map
HomoSapiens.PhyloMap <- read_excel(paste0(BASE, "RProject/OT-bd-analyses_loc/MBE_2008_Homo_Sapiens_PhyloMap.xls"), 
                                   sheet = 1, skip = 1)
HomoSapiens.PhyloMap <- HomoSapiens.PhyloMap[ , c("Phylostratum","Gene_ID")]          # 22845 genes


# get gene names for all genes in the genome but only those that are available in the phylo map (NOT limited to OT pathway genes)
genensid <- ensembldb::select(EnsDb.Hsapiens.v79, 
                              keys= HomoSapiens.PhyloMap$Gene_ID, 
                              keytype = "GENEID", 
                              columns = c("SYMBOL","GENEID"))
genensid <- genensid %>% distinct(SYMBOL, .keep_all = TRUE)                          

# annotate list of gene names with phylostrata
genensid_ps <- dplyr::inner_join(HomoSapiens.PhyloMap, genensid, by = c("Gene_ID" = "GENEID")) # 18844 genes



# subset into more modern genes list
yGenes <- genensid_ps[genensid_ps$Phylostratum %in% yPS, ] # keep only more modern phylostrata, 5771 genes

# subset into older genes list
oGenes <- genensid_ps[!genensid_ps$Phylostratum %in% yPS, ] # keep only older phylostrata, 13073 genes



# ...... subset modern OT gene names ref list

M1C_mm <- M1C_tai$mm                       # I'm using this df as reference for OT pathway genes and not "kop" because it is the set of genes that we use in the main 
M1C_mm$GeneID <- toupper(M1C_mm$GeneID)    # myTAI analyses. For consistency it might be best to stick with the same list of genes in as many instances/analyses as possible.

yOT <- yGenes[yGenes$Gene_ID %in% M1C_mm$GeneID, ] # keep only more modern OT genes, 20 genes

# ...... subset older OT gene names ref list
oOT <- oGenes[oGenes$Gene_ID %in% M1C_mm$GeneID, ] # keep only older OT genes, 118 genes






# subset abagen Destrieux CORTICAL expression matrix into young genes
# ....... only left hemisphere
ahbaDes_yc <- AHBAdestrieux_L[, names(AHBAdestrieux_L) %in% c("index", "name", yGenes$SYMBOL)] # extract subset of more modern genes with cortical expression
dim(ahbaDes_yc)                                     # 3227 genes, 74 regions

# extract young OT subset
OTdes_yc2 <- ahbaDes_yc[, names(ahbaDes_yc) %in% c("index", "name", yOT$SYMBOL)] # extract subset of modern OT genes with cortical expression
dim(OTdes_yc2) # 19 modern genes, 74 cortical regions


# ....... both hemispheres 
ahbaDes_yc_lr <- AHBAdestrieux_a[, names(AHBAdestrieux_a) %in% c("index", "name", yGenes$SYMBOL)] # extract subset of more modern genes with cortical expression
dim(ahbaDes_yc_lr)                                     # 3227 genes, 148 regions

# extract young OT subset
OTdes_yc_lr <- ahbaDes_yc_lr[, names(ahbaDes_yc_lr) %in% c("index", "name", yOT$SYMBOL)] # extract subset of modern OT genes with cortical expression
dim(OTdes_yc_lr) # 19 modern genes, 148 cortical regions





# subset abagen SUBCORTICAL expression matrix into young genes
ahbaDes_ysc <- AHBApauli_a[, names(AHBApauli_a) %in% c("ids", "region", "node_info", yGenes$SYMBOL)] # extract subset of more modern genes with cortical expression
dim(ahbaDes_ysc)                                     # 3227 genes, 16 regions

# extract young OT subset
OTdes_ysc2 <- ahbaDes_ysc[, names(ahbaDes_ysc) %in% c("ids", "region", "node_info", yOT$SYMBOL)] # extract subset of modern OT genes with cortical expression
dim(OTdes_ysc2) # 19 modern genes, 16 subcortical regions





# CORTICAL overexpression analysis (Destrieux)

# create extra column for lateralization
OTdes_yc_lr <- data.frame(OTdes_yc_lr)
OTdes_yc_lr$hemisphere <- "L"
OTdes_yc_lr$hemisphere[75:148] <- "R"

OTdes_yc_lr$name  <- gsub('L ', '', OTdes_yc_lr$name)
OTdes_yc_lr$name  <- gsub('R ', '', OTdes_yc_lr$name)


# average expression across all genes in c OT df
OTdes_yc_lr$mean_OTy <- rowSums(OTdes_yc_lr[,3:21])/19 
OT_yc_m <- OTdes_yc_lr[,c(1:2,22,23)]

library(powerAnalysis)

mean <- mean(OT_yc_m$mean_OTy)             
regions <- unique(OT_yc_m$name)
pvals <- numeric()
tStatistic <- numeric()
cohensD <- numeric()

for (i in 1:length(regions)){
  Temp <- t.test(OT_yc_m[OT_yc_m$name==regions[i], "mean_OTy"], 
                 mu = mean, alternative = "two.sided")
  CohD <- ES.t.one(m=Temp$estimate,
                   sd=sd(OT_yc_m[OT_yc_m$name==regions[i],"mean_OTy"]),
                   mu=mean,
                   alternative = "two.sided")
  pvals <- c(pvals,Temp$p.value)
  tStatistic <- c(tStatistic,Temp$statistic)
  cohensD <- c(cohensD,CohD$d)
}

results <- data.frame(Region=regions,P_unadjusted=pvals,
                      t_stat=tStatistic,CohD=cohensD)
results$P_fdr <- p.adjust(results$P_unadjusted, method="fdr")

results$stat_sign <- ""
results$stat_sign[results$P_fdr <= 0.05]  <- "*"
results



# SUBCORTICAL overexpression analysis (Pauli)

# not possible because the data is not lateralized, t.test needs minimum N = 2, non-lateralization leads to N = 1







#### ENRICHMENT ANALYSES

# ---------------- manual with binomial generalized linear model
### subcortical

AHBApauli_x <- AHBApauli_a[, -c(1:2)]                                                      
AHBApauli_x <- AHBApauli_x %>% remove_rownames %>% column_to_rownames(var="node_info")
AHBApauli_x_t <- data.frame(t(AHBApauli_x))
AHBApauli_x_t <- rownames_to_column(AHBApauli_x_t)         # transpose and create a new column with gene symbols


cnames_p <- data.frame(AHBApauli_x_t$rowname)                # prepare matrix with gene symbols and new column to categorize genes as "OT young" = 1 or "not" = 0
cnames_p$OTgs <- 0
names(cnames_p)[names(cnames_p) == "AHBApauli_x_t.rowname"] <- "GeneSymbols"

# mark young OT genes with "1"
OTy_genesym <- yOT$SYMBOL

for (i in OTy_genesym) {
  cnames_p[cnames_p$GeneSymbols==i,"OTgs"] <- 1  
}

# combine
AHBApauli_x_ot <- dplyr::inner_join(cnames_p, AHBApauli_x_t, by=c("GeneSymbols"="rowname"))


glm_list_p <- list()
glm_coefs_pauli <- list()
for (i in 3:18) {
  tmp=as.formula(paste0("OTgs~",colnames(AHBApauli_x_ot)[i]))            # define formula in an object
  glm_test_p <- glm(tmp, family = binomial, data = AHBApauli_x_ot)
  sum_glm_p <- summary(glm_test_p)
  glm_list_p[[i]] <- sum_glm_p                                           # store sum stats in a list
  glm_coefs_pauli[[i]] <- glm_list_p[[i]]$coefficients                    # store coefs of sum stats in a list
  glm_coefs_p <- data.frame(do.call(rbind, glm_coefs_pauli))   
  glm_coefs_pp <- rownames_to_column(glm_coefs_p)                                                             # combine into one df
}

glmCoefs_p <- glm_coefs_pp[!grepl("^X.Intercept.", glm_coefs_pp$rowname), ]
View(glmCoefs_p)


pvals_p <- glmCoefs_p$Pr...z..



glmCoefs_p$pFDR <- p.adjust(pvals_p, method = "fdr", n = length(pvals_p))

glmCoefs_p$stat_sign = ""
glmCoefs_p[glmCoefs_p$pFDR <= .05,"stat_sign"] <- "*"
glmCoefs_p[glmCoefs_p$pFDR <= .01,"stat_sign"] <- "**"
glmCoefs_p[glmCoefs_p$pFDR <= .001,"stat_sign"] <- "***"


# so not significant anywhere?
# but FUMA said its enriched in the brain?





### cortical


## Destrieux cortical


AHBAdes_L <- AHBAdestrieux_L[, -1]                                                      
AHBAdes_L <- AHBAdes_L %>% remove_rownames %>% column_to_rownames(var="name")
AHBAdes_L_t <- data.frame(t(AHBAdes_L))
AHBAdes_L_t <- rownames_to_column(AHBAdes_L_t)         # transpose and create a new column with gene symbols


cnames <- data.frame(AHBAdes_L_t$rowname)              # prepare matrix with gene symbols and new column to categorize genes as "OT young" = 1 or "not" = 0
cnames$OTgs <- 0
names(cnames)[names(cnames) == "AHBAdes_L_t.rowname"] <- "GeneSymbols"

# mark young OT genes with "1"
OTy_genesym <- yOT$SYMBOL

for (i in OTy_genesym) {
  cnames[cnames$GeneSymbols==i,"OTgs"] <- 1  
}

# combine
AHBAdes_L_ot <- dplyr::inner_join(cnames, AHBAdes_L_t, by=c("GeneSymbols"="rowname"))



glm_list7 <- list()
glm_coefs_l_3 <- list()
for (i in 3:76) {
  tmp=as.formula(paste0("OTgs~",colnames(AHBAdes_L_ot)[i]))            # define formula in an object
  glm_test1 <- glm(tmp, family = binomial, data = AHBAdes_L_ot)
  sum_glm1 <- summary(glm_test1)
  glm_list7[[i]] <- sum_glm1                                           # store sum stats in a list
  glm_coefs_l_3[[i]] <- glm_list7[[i]]$coefficients                    # store coefs of sum stats in a list
  glm_coefs4 <- data.frame(do.call(rbind, glm_coefs_l_3))   
  glm_coefs5 <- rownames_to_column(glm_coefs4)                                                             # combine into one df
}

View(glmCoefs2)

glmCoefs2 <- glm_coefs5[!grepl("^X.Intercept.", glm_coefs5$rowname), ]

pvals2 <- glmCoefs2$Pr...z..



glmCoefs2$pFDR <- p.adjust(pvals2, method = "fdr", n = length(pvals2))

glmCoefs2$stat_sign = ""
glmCoefs2[glmCoefs2$pFDR <= .05,"stat_sign"] <- "*"
glmCoefs2[glmCoefs2$pFDR <= .01,"stat_sign"] <- "**"
glmCoefs2[glmCoefs2$pFDR <= .001,"stat_sign"] <- "***"



# generate 100 random "young" gene sets with 19 random genes out of transposed "abagen subcortical expression for young genes" - subset
fun_ysc1 <- function(i) {                                # function to generate random gene set of Ngene = 19
  sample_n(ahba_ysc_t_s, 19) 
}

erg_ysc <- lapply(1:100,fun_ysc1) 





# so not significant anywhere?
# but FUMA said its enriched in the brain?




# # ---------------- enrichment of young OT pathway gene set in the brain with ABAEnrichment

library(ABAEnrichment)
gene_ids = OTy_genesym
input_hyper = data.frame(gene_ids, is_candidate=1)
head(input_hyper)

res_adult = aba_enrich(input_hyper, dataset='adult', gene_len = T, n_randsets=10000) # default for random permutations: 1000
View(res_adult$results)

res <- res_adult$results

pattern <- c("Left|left")

res_lhemi <- res[grep(pattern, res$structure),]

# no enrichment either...?

##########################################