#########################
# setting up working environment 
#########################

rm(list=ls()) # delete all objects in the workspace
gc(reset=T) # resest memory (especially useful when working with large data sets)

options(stringsAsFactors=F) # disables automatic conversion of char. strings into factors

Sys.setenv(LANG = "en")

#install.packages("installr")

#library(installr)
# updateR()

BASE="C:/01Dokumente/07Erwerbstaetigkeit/01PhD_UiO/projects/first_author/oxytocinHealth/RQ4/"

#########################
# setting up working environment 
#########################


######################################################################################################
######################################################################################################
######################################################################################################


cols_metat <- t(cols_meta)
new_hd_cols <- c("row_num", cols_meta$structure_acronym)
colnames(ex_m) <- new_hd_cols 

ex_m_t2 <- t(ex_m)

ex_m_meta <- merge(rows_meta, ex_m, by.x = "row_num", by.y = "row_num") # merge expression matrix with gene meta info. there is a problem here because some of the columns will have the same name... how to fix this?
ex_m_meta <- ex_m_meta[, -c(1:2, 5)]

#Gender <- c(c("", "", ""), cols_meta$gender)

test <- merge(phylo, ex_m_meta, by.x="Gene_ID", by.y="ensembl_gene_id") # Ngenes = 19256
test2 <- merge(OT_pathway_genes, test, by.x="OT_pathway_genes", by.y="gene_symbol") # Ngenes = 146

Age <- c(c("", "", ""), grepl("8 pcw", cols_meta$age))
donor_id <- c(c("", "", ""), cols_meta$donor_id) # continue with rbind()
test3 <- rbind(Age, test2)
test4 <- rbind(donor_id, test3)

test4t <- t(test4)

pren = c("8 pcw|9 pcw|12 pcw|13 pcw|16 pcw|17 pcw|19 pcw|21 pcw|24 pcw|25 pcw|26 pcw|35 pcw|37 pcw")
prental <- subset(test4t, grepl(pren, test4t$V2))




######################################################################################################
######################################################################################################
######################################################################################################


pat = c("^M1C$", "^DFC$", "^VFC$", "^OFC$", "^S1C$", "^IPC$", "^A1C$", "^STC$", "^ITC$", "^V1C$", "^MFC$", "^HIP$", 
            "^STR$","^AMY$", "^MD$", "^CBC$")

testx <- data.frame(row_num = rep(NA, 52376))

for (i in pat) {
  temp = rowMeans(prenatal_df[, grep(pattern=i, colnames(prenatal_df))]) # 16 new columns with 52376 rows
  testx[ , i] <- temp
}
testx <- testx[, -c(1)]

View(testx)

######################################################################################################
######################################################################################################
######################################################################################################


pat = c("^M1C$", "^DFC$", "^VFC$", "^OFC$", "^S1C$", "^IPC$", "^A1C$", "^STC$", "^ITC$", "^V1C$", "^MFC$", "^HIP$", 
        "^STR$","^AMY$", "^MD$", "^CBC$")

testx <- vector("list", 16)

for (i in pat) {
  testx[[i]] = rowMeans(prenatal_df[, grep(pattern=i, colnames(prenatal_df))]) # 16 new columns with 52376 rows
}
testx <- do.call(cbind, testx)

View(testx)



######################################################################################################
######################################################################################################
######################################################################################################


# ----------------------------------------------- data structuring for prenatal ... loop this somehow or put it in a function to apply to the other ontogenetic stages....

prenatal = NULL
pre = c("8pcw", "9pcw", "12pcw", "13pcw", "16pcw", "17pcw", "19pcw", "21pcw", "24pcw", "25pcw", "26pcw", "35pcw", "37pcw")
for (i in pre) {
  pre_temp <- ex_m_col[which(ex_m_col$age == i), ]
  prenatal = rbind(prenatal,pre_temp)
}

hd_pre <- prenatal$structure_acronym
prenatal <- prenatal[, -c(1:8)]
prenatal_t <- t(prenatal)
prenatal_df <- as.data.frame(prenatal_t)
colnames(prenatal_df) <- hd_pre 

# calculate row means for every row of a specific column, i.e., all cols named "A1C"
pat = c("M1C", "DFC", "VFC", "OFC", "S1C", "IPC", "A1C", "STC", "ITC", "V1C", "MFC", "HIP", 
        "STR","AMY", "MD", "CBC")
prenat_str_m <- data.frame(row_num = 1:52376)
for (i in pat) {
  temp = rowMeans(prenatal_df[, grep(pattern=paste0("^",i,"$"), colnames(prenatal_df))]) # 16 new columns with 52376 rows
  prenat_str_m[ , i] <- temp
}

prenat_meta <- merge(rows_meta, prenat_str_m, by.x="row_num", by.y="row_num", all=T)
prenat_meta <- prenat_meta[,-c(1:2, 5)]
prenat_m_phy <- merge(phylo, prenat_meta, by.x="Gene_ID", by.y="ensembl_gene_id", all=T)
prenat_mpp <- merge(OT_pathway_genes, prenat_m_phy, by.x="OT_pathway_genes", by.y="gene_symbol", all=T)


# ----------------------------------------------- data structuring for infancy

infancy = NULL
inf = c("4mos", "10mos", "1yrs")
for (i in inf) {
  inf_temp <- ex_m_col[which(ex_m_col$age == i), ]
  infancy = rbind(infancy,inf_temp)
}

hd_inf <- infancy$structure_acronym
infancy <- infancy[, -c(1:8)]
infancy_t <- t(infancy)
infancy_df <- as.data.frame(infancy_t)
colnames(infancy_df) <- hd_inf 

pat = c("M1C", "DFC", "VFC", "OFC", "S1C", "IPC", "A1C", "STC", "ITC", "V1C", "MFC", "HIP", 
        "STR","AMY", "MD", "CBC")
infan_str_m <- data.frame(row_num = 1:52376)
for (i in pat) {
  temp = rowMeans(infancy_df[, grep(pattern=paste0("^",i,"$"), colnames(infancy_df))]) # 16 new columns with 52376 rows
  infan_str_m[ , i] <- temp
}

infan_meta <- merge(rows_meta, infan_str_m, by.x="row_num", by.y="row_num", all=T)
infan_meta <- infan_meta[,-c(1:2, 5)]
infan_m_phy <- merge(phylo, infan_meta, by.x="Gene_ID", by.y="ensembl_gene_id", all=T)
infan_mpp <- merge(OT_pathway_genes, infan_m_phy, by.x="OT_pathway_genes", by.y="gene_symbol", all=T)



# ----------------------------------------------- data structuring for childhood

childhood = NULL
chi = c("2yrs", "3yrs", "4yrs", "8yrs", "11yrs")
for (i in chi) {
  chi_temp <- ex_m_col[which(ex_m_col$age == i), ]
  childhood = rbind(childhood,chi_temp)
}

hd_chi <- childhood$structure_acronym
childhood <- childhood[, -c(1:8)]
childhood_t <- t(childhood)
childhood_df <- as.data.frame(childhood_t)
colnames(childhood_df) <- hd_chi 

pat = c("M1C", "DFC", "VFC", "OFC", "S1C", "IPC", "A1C", "STC", "ITC", "V1C", "MFC", "HIP", 
        "STR","AMY", "MD", "CBC")
child_str_m <- data.frame(row_num = 1:52376)
for (i in pat) {
  temp = rowMeans(childhood_df[, grep(pattern=paste0("^",i,"$"), colnames(childhood_df))]) # 16 new columns with 52376 rows
  child_str_m[ , i] <- temp
}

child_meta <- merge(rows_meta, child_str_m, by.x="row_num", by.y="row_num", all=T)
child_meta <- child_meta[,-c(1:2, 5)]
child_m_phy <- merge(phylo, child_meta, by.x="Gene_ID", by.y="ensembl_gene_id", all=T)
child_mpp <- merge(OT_pathway_genes, child_m_phy, by.x="OT_pathway_genes", by.y="gene_symbol", all=T)


# ----------------------------------------------- data structuring for adolescence

adolescence = NULL
ado = c("13yrs", "15yrs", "18yrs", "19yrs")
for (i in ado) {
  ado_temp <- ex_m_col[which(ex_m_col$age == i), ]
  adolescence = rbind(adolescence,ado_temp)
}

hd_ado <- adolescence$structure_acronym
adolescence <- adolescence[, -c(1:8)]
adolescence_t <- t(adolescence)
adolescence_df <- as.data.frame(adolescence_t)
colnames(adolescence_df) <- hd_ado 

pat = c("M1C", "DFC", "VFC", "OFC", "S1C", "IPC", "A1C", "STC", "ITC", "V1C", "MFC", "HIP", 
        "STR","AMY", "MD", "CBC")
adols_str_m <- data.frame(row_num = 1:52376)
for (i in pat) {
  temp = rowMeans(adolescence_df[, grep(pattern=paste0("^",i,"$"), colnames(adolescence_df))]) # 16 new columns with 52376 rows
  adols_str_m[ , i] <- temp
}

adols_meta <- merge(rows_meta, adols_str_m, by.x="row_num", by.y="row_num", all=T)
adols_meta <- adols_meta[,-c(1:2, 5)]
adols_m_phy <- merge(phylo, adols_meta, by.x="Gene_ID", by.y="ensembl_gene_id", all=T)
adols_mpp <- merge(OT_pathway_genes, adols_m_phy, by.x="OT_pathway_genes", by.y="gene_symbol", all=T)


# ----------------------------------------------- data structuring for adulthood

adulthood = NULL
adu = c("21yrs", "23yrs", "30yrs", "36yrs", "37yrs", "40yrs")
for (i in adu) {
  adu_temp <- ex_m_col[which(ex_m_col$age == i), ]
  adulthood = rbind(adulthood,adu_temp)
}

hd_adu <- adulthood$structure_acronym
adulthood <- adulthood[, -c(1:8)]
adulthood_t <- t(adulthood)
adulthood_df <- as.data.frame(adulthood_t)
colnames(adulthood_df) <- hd_adu 

pat = c("M1C", "DFC", "VFC", "OFC", "S1C", "IPC", "A1C", "STC", "ITC", "V1C", "MFC", "HIP", 
        "STR","AMY", "MD", "CBC")
adult_str_m <- data.frame(row_num = 1:52376)
for (i in pat) {
  temp = rowMeans(adulthood_df[, grep(pattern=paste0("^",i,"$"), colnames(adulthood_df))]) # 16 new columns with 52376 rows
  adult_str_m[ , i] <- temp
}

adult_meta <- merge(rows_meta, adult_str_m, by.x="row_num", by.y="row_num", all=T)
adult_meta <- adult_meta[,-c(1:2, 5)]
adult_m_phy <- merge(phylo, adult_meta, by.x="Gene_ID", by.y="ensembl_gene_id", all=T)
adult_mpp <- merge(OT_pathway_genes, adult_m_phy, by.x="OT_pathway_genes", by.y="gene_symbol", all=T)





######################################################################################################
######################################################################################################
######################################################################################################





# retrieve expression data
exp <- get_expression(structure_ids=c("Allen:10163"),
                      gene_ids=kop, 
                      dataset='5_stages') # Get expression values

# make it look "palatable"
exp_l <- data.frame(matrix(unlist(exp), 
                           nrow=length(exp), byrow=T),
                    use.names = T)  # Unlists without row and column names

# rename rows
rownames(exp_l) <- c("prenatal", 
                     "infant", 
                     "child", 
                     "adolescent", 
                     "adult")

# transpose and rename cols
exp_s1 <- exp$age_category_1 # or: exp_s1 <- exp[["age_category_1"]]
cn <- colnames(exp_s1)
colnames(exp_l) <- cn
exp_l <- as.data.frame(exp_l)
exp_l_t <- as.data.frame(t(exp_l))
exp_l_t <- rownames_to_column(exp_l_t, var = "gene")
exp_l_t <- head(exp_l_t, -1)

############################################
# Alternative, slightly shorter code:

exp <- get_expression(structure_ids=c("Allen:10163"),
                      gene_ids=kop, 
                      dataset='5_stages') 

new <- as.data.frame(do.call(rbind, exp)) 
rownames(new) <- NULL
rownames(new) <- c("Prenatal", "Infant", "Child", "Adolescent", "Adult")
new_t <- as.data.frame(t(new))
new_t <- rownames_to_column(new_t, var = "gene")
############################################






# get gene ensemble ID w/ matchting gene name/symbol
ens_gene_dev <- ensembldb::select(EnsDb.Hsapiens.v79, keys= exp_l_t$gene, keytype = "SYMBOL", columns = c("SYMBOL","GENEID")) 

# removes duplicates?
ens_gene_dev <- ens_gene_dev %>% distinct(SYMBOL, .keep_all = TRUE)

# rename SYMBOL column into "gene"
names(ens_gene_dev)[names(ens_gene_dev) == "SYMBOL"] <- "gene"   # alternative: df %>% rename(Gene_ID = GENEID) !! requires package tidyverse !!

ex_t_f <- dplyr::inner_join(ens_gene_dev, exp_l_t, by = "gene")   # alternative: merge, all=T

ex_t_f <- dplyr::select(ex_t_f, -gene)                           # alternative: df <- df[, -c(col_num)]

names(ex_t_f)[names(ex_t_f) == "GENEID"] <- "Gene_ID"            # alternative: df %>% rename(Gene_ID = GENEID) !! requires package tidyverse !!

############################################
# Alternative, slightly shorter code:
ens_gene_dev_al <- as.data.frame(ensembldb::select(EnsDb.Hsapiens.v79, 
                                  keys= new_t$gene, 
                                  keytype = "SYMBOL", 
                                  columns = c("SYMBOL","GENEID"))) 
ens_gene_dev_al <- ens_gene_dev_al %>% distinct(SYMBOL, .keep_all = TRUE)
# this should work, but it doesn't.. why?                rename(ens_gene_dev_al,  Gene_ID = SYMBOL)
names(ens_gene_dev_al)[names(ens_gene_dev_al) == "SYMBOL"] <- "gene"
new_t_f <- merge(ens_gene_dev_al, new_t, by="gene", all=T)
new_t_f <- new_t_f[, -c(1)]
names(new_t_f)[names(new_t_f) == "GENEID"] <- "Gene_ID" 
############################################







HomoSapiens.PhyloMap <- read_excel("D:/eigene_dateien/02Dokumente/07Erwerbstaetigkeit/01PhD_UiO/projects/first_author/oxytocinHealth/RQ4/github_dan/MBE_2008_Homo_Sapiens_PhyloMap.xls", 
                              sheet = 1, skip = 1)

HomoSapiens.PhyloMap <- HomoSapiens.PhyloMap[ , c("Phylostratum","Gene_ID")]

mm <- MatchMap(HomoSapiens.PhyloMap, ex_t_f)

############################################
hs_phylomap <- as.data.frame(read_excel("D:/eigene_dateien/02Dokumente/07Erwerbstaetigkeit/01PhD_UiO/projects/first_author/oxytocinHealth/RQ4/github_dan/analysis/MBE_2008_Homo_Sapiens_PhyloMap_allgenes.xls", skip=1))
ex_phylo <- merge(hs_phylomap, new_t_f, by="Gene_ID", all=F)
ex_phylo <- ex_phylo[, -c(2, 4)]
ex_phylo <- ex_phylo %>% distinct(Gene_ID, .keep_all = TRUE) # if i do it this way, i have one more row?
ex_phylo = ex_phylo %>% relocate(Phylostratum)
ex_phylo <- ex_phylo[order(ex_phylo$Phylostratum),]
names(ex_phylo)[names(ex_phylo) == "Gene_ID"] <- "GeneID" 
############################################

# my method ends up with one additional row.. why?

# optional: Expressed() function? What would be the cut-offs?



# DONE! Now the dataset is prepared for further analysis.

# ==================================== QC and pre-analysis tests
# flatline test to test for deviation from phylotrancriptomic pattern 
# myTAI documentation p22
# (great, what does this tell me after all? What does it mean or indicate if
# this dataset deviates from the pattern? low quality, exciting event, etc.?

flatline <- FlatLineTest(ExpressionSet = mm,
                         permutations = 1000,
                         plotHistogram = F,
                         runs = 100,
                         parallel = F
                         )
summary(flatline)


# function performing quality checks (still gotta understand this)
# I don't know if this makes sense
# myTAI documentation p61

Rep_quality <- PlotReplicateQuality(ExpressionSet = mm,
                                    nrep=5, # this must equal the number of stages?
                                    legend.pos = "topright")
Rep_quality






## ==================================== PLOT frequency distribution of genes per PS ##
# myTAI documentation p 42
# this does not plot frequencies for PS 10 and 12 for some reason?
# this is somewhat the same as Dan did with the barplot in the beginning, but withoout the annotation of the PS and PS 10 and 12 are missing as per usual

rel_freq_dis <- PlotDistribution(PhyloExpressionSet = mm,
                                 legendName="PS",
                                 as.ratio=T,
                                 use.only.map = F)
rel_freq_dis

abs_fre_dis <- PlotDistribution(PhyloExpressionSet = mm,
                                legendName="PS",
                                as.ratio=F,
                                use.only.map = F)
abs_fre_dis


# QUESTION
# myTAI documentation p 62
# difference to PlotDistribution function?
# --> this allows to specify for a subset of genes while the other function uses the whole exp-phylo map. Since we used the subset of 
# interest as our base exp-phylo map, the two functions yield the same output if we do not specify further within the OT pathway gene set.

# PROBLEM: doesnt plot PS 10 and 12 again...., had to add dummy rows to make it plot until 10 and 12. other solutions?
ps8 = as.numeric(c(8, 0, 0, 0, 0, 0, 0))
mm_tem <- rbind(mm, ps8)
ps9 = as.numeric(c(9, 0, 0, 0, 0, 0, 0))
mm_tem <- rbind(mm_tem, ps9)
ps11 = as.numeric(c(11, 0, 0, 0, 0, 0, 0))
mm_tem <- rbind(mm_tem, ps11)
PS_dis <- PlotSelectedAgeDistr(ExpressionSet = mm_tem,
                               gene.set = c("ENSG00000180914", "ENSG00000101405", "ENSG00000004468"),
                               legendName = "PS",
                               as.ratio = F,
                               use.only.map = F)
# ENSG00000180914 = OXTR, ENSG00000101405 = OXT, ENSG00000004468 = CD38
PS_dis








# ==================================== MAIN TAI FUNCTIONS 
## compute TAI w/o plot
# myTAI documentation p 86
TAI <- TAI(PhyloExpressionSet = mm)
TAI

# Plot TAI signature, different phylostrata merged in one graph line
# myTAI documentation p63
## Function for the TAI without plot: TAI()

# QUESTION: why does it generate NaNs?
# for the FlatLineTest(), there is a similar issue, but it's been explained: In case there are extreme outlier 
# expression values stored in the dataset (PhyloExpressionSet or DivergenceExpressionSet), the internal fitdist 
# function that is based on the bootMatrix output might return a warning: "In densfun(x, parm[1], parm[2], ...) : NaNs were produced" 
# which indicates that permutation results caused by extreme outlier expression values that could not be fitted accordingly. This 
# warning will not be printed out when the corresponding outlier values are extracted from the dataset.

TAI_plot <- PlotSignature(ExpressionSet= mm, 
                          measure="TAI", 
                          TestStatistic = "FlatLineTest", 
                          p.value=F)
TAI_plot

# the higher the TAI, the younger the transcriptome
# transcriptome is very young during prental, meaning more younger genes contribute to the transcriptome
# while more ancient gene contribute to the mRNA profile in adults?
# isn't that a bit weird since some genes in the 2nd phylostratum are extremely highly expressed during the prenatal stage and the
# genes from the 10th and 12th PS are almost not at all expressed during prenatal but have higher expression values in later stages?
# maybe, we need to exclude some expression outliers?


## Plot pattern - TAI of a PhyloExpSet ##
# myTAI documenation p55
# QUESTION
# what is the difference to the PlotSignature function? it looks less pretty, but apart from that...?
# with the PlotSignature(), one can choose between TAI, TDI and TPI, but apart from that the input parameters are almos identical

#pattern <- PlotPattern(ExpressionSet = mm,
                       #TestStatistic = "FlatLineTest")

#pattern

# Plot TAI signature, phylostrata displayed separately
# myTAI documentation p39

TAIci_ps <- PlotContribution(ExpressionSet = ex_phylo, legendName = "PS")
TAIci_ps

## data only for PlotContribution function
# myTAI documentation p70

PS_TAI_contrib <- pTAI(ExpressionSet = ex_phylo)
PS_TAI_contrib

## partial TAI for each PS, i.e., the overall TAI of a brain region is distributed onto the 9 PS. I guess, the PS that had most effect 
# is assigned the highest partial value?
# myTAI documentation p69

TAI_to_ps <- pStrata(ExpressionSet = ex_phylo)
TAI_to_ps

## partial TAI per gene, NO PLOTTING
# myTAI documenation p68

TAIpar <- pMatrix(ExpressionSet = ex_phylo)
TAIpar <- data.frame(TAIpar)
TAIpar$GeneID <- ex_phylo$GeneID
View(TAIpar)

## ... but if you want to plot:
# again, very weird plot, maybe don't use it or prettify it?
TAIpar_bp <- boxplot(pMatrix(ExpressionSet = ex_phylo))
TAIpar_bp




## ==================================== Early Conservation Test ##
# The early conservation model predicts that the highest developmental constraints occur at the beginning of embryogenesis.
# This could indicate that embryos of different species progressively diverge from one another during ontogeny.
# Alternative hypothesis is the hourglass model, highest contraint happens during mid-embryogensis
# Biological constraints are factors which make populations resistant to evolutionary change. One proposed definition of 
# constraint is "A property of a trait that, although possibly adaptive in the environment in which it originally evolved, 
# acts to place limits on the production of new phenotypic variants." P-values < .05 = TAI follows early conservation like pattern. 

## there is also a test for the reductive hourglass model! ReductiveHourglassTest()
## there is also a test for the reverse hourglass score! reversehourglassScore()
## there is also a test for the reverse hourglass model! ReverseHourglassTest()

#RECT <- EarlyConservationTest(ExpressionSet = mm, modules= list(early= 1:2, mid=3:4, late=5))
#RECT # p = 1, TAI does not follow an early conservation like pattern. What does this mean?






# ==================================== Plot Phylostratum Enrichment ##
# myTAI documentation p44 and p 18
# This Phylostratum or Divergence Stratum enrichment analysis is motivated by Sestak and Domazet-
# Loso (2015) who perform Phylostratum or Divergence Stratum enrichment analyses to correlate
# organ evolution with the origin of organ specific genes.

# --> this is just an addition to the bar plots Dan did in the beginning, proving that OT genes are indeed enriched in the first PS

# I think for this analysis we need a Expression Set with *all* genes which is compared to our OT pathway testset

# get full expression set with all genes available**
  rows_meta <- read.csv("D:/eigene_dateien/02Dokumente/07Erwerbstaetigkeit/01PhD_UiO/projects/first_author/oxytocinHealth/RQ4/analyses/data/raw/brainspanatlas/rows_metadata.csv")
  all_geneids <- rows_meta$ensembl_gene_id
  # available entries: 52376 (might contain duplicates) 

  # get expression data from all ontogenetic stages for each brain region for as many genes as possible
  exp_all <- get_expression(structure_ids=c("Allen:10163"),
                         gene_ids=all_geneids, 
                         dataset='5_stages') 
  ex_all_df <- as.data.frame(do.call(rbind, exp_all)) 
  rownames(ex_all_df) <- NULL
  rownames(ex_all_df) <- c("Prenatal", "Infant", "Child", "Adolescent", "Adult")
  ex_all_t <- as.data.frame(t(ex_all_df))
  ex_all_t <- rownames_to_column(ex_all_t, var = "Gene_ID")
  ex_all_t <- ex_all_t %>% distinct(Gene_ID, .keep_all = TRUE)
  # available genes: 17259
  
  # get PS data
  HomoSapiens.PhyloMap <- read_excel("D:/eigene_dateien/02Dokumente/07Erwerbstaetigkeit/01PhD_UiO/projects/first_author/oxytocinHealth/RQ4/github_dan/MBE_2008_Homo_Sapiens_PhyloMap.xls", 
                                     sheet = 1, skip = 1)
  HomoSapiens.PhyloMap <- HomoSapiens.PhyloMap[ , c("Phylostratum","Gene_ID")]
  # available genes: 22845

  # merge expression and phylogenetic data
  mm_all <- MatchMap(HomoSapiens.PhyloMap, ex_all_t)

# ** ... how many genes are available in that data set originally? in the downloadable version, there are over 52.000, here, I get about 17.000 (after matching
# and all the other steps: Ngenes = 15858)



enrich <- EnrichmentTest(ExpressionSet = mm_all,
                         test.set = OT_subset,
                         use.only.map = F,
                         measure = "log-foldchange",
                         complete.bg = F) # not sure about this
enrich$p.values

enrich_p <- PlotEnrichment(ExpressionSet = mm_all, 
                         test.set = OT_subset, 
                         legendName = "PS",
                         measure = "log-foldchange",
                         p.adjust.method = "fdr",
                         use.only.map = F, 
                         complete.bg = F,
                         plot.bars=T)
                           
enrich_p$p.values
enrich_p$p_fdr <-  p.adjust(enrich_p$p.values, method = "fdr", n = 19)
enrich_p$p_fdr




res_m_test <- data.frame(enrich_p$enrichment.matrix)

res_m_test$pvals_fdr <- enrich_p$p.values
res_m_test <- data.frame(res_m_test)
# the plot looks interesting generally, but a bit weird... find out why and how to change



## Plot Gene Set expression profiles ##
# myTAI documentation p47
# (very messy and weird, don't use this or make it prettier) 

#ex_prof <- PlotGeneSet(ExpressionSet = ex_phylo,
                      # gene.set = genes,
                      # get.subset = F,
                      # use.only.map = F,
                      # plot.legend = F)



## compare and PLOT PS group differences ##
# plot looks weird, I don't know about this...

#Group_diffs_ex <- PlotBarRE(ExpressionSet = mm, 
                           # Groups = list(c(1:3), c(4:6), c(7,12)),
                           # ratio = T,
                           # p.adjust.method = "fdr")






## ==================================== mean expression values of gene set of each PS per ontogenetic stage ##
# myTAI documentation p 51
# apparently it works with one group only?
#mm1to7 <- mm[!grepl("^10$", mm$Phylostratum), ]
#mm1to7 <- mm1to7[!grepl("^12$", mm1to7$Phylostratum), ]
            
means <- PlotMeans(ExpressionSet = mm,
                   Groups = list(c(1:7, 10, 12)),
                   
                   legendName="PS")
means


# 1/2 generally, genes belonging to the 2nd PS have the highest mean expression values across all ontogenetic stages, with mRNA during
# the prenatal stage containing the highest values .....

# an issue that comes to my mind here is that the sample sizes of the different PS are not the same, in one PS, there are
# 20 genes while in another, there is 1 or even zero - is this relevant?



## plots exp. levels of each PS per ontogenetic stage (if you will) ##
# myTAI documentation p 35
# why do these functions always turn the 10th and 12th PS into NAs?

# 2/2 ... as becomes furhter apparent from these violin plots. We further see that genes from the 2nd and 3rd PS have expression values 
# that are characterized by great variability (long stretched violins). 

# Side note: when the values are not log-transformed, the mRNA from the 2nd PS during
# prenatal literally shoots through the roof. Is this maybe an outlier?

violinPlot <- PlotCategoryExpr(ExpressionSet = mm,
                             legendName = "PS",
                             test.stat = T,
                             distr.type = "violin",
                             log.expr = T,
                             gene.set = mm$GeneID,
                             type= "category-centered")
violinPlot







## ==================================== shenanigans


# for the TPI, we would need a polymorphism expression set....

test <- bootMatrix(ExpressionSet = mm, permutations=10)
test
# shouldn't this be random and changing?




if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("biomaRt")

library(biomaRt)
listMarts()
ensembl <- useMart("ensembl")
listAttributes(ensembl)
ds_ensembl <- listDatasets(ensembl)
# biomaRt ensembl has only sapiens and macaca, not neanderthalensis ir denisova



# relative expression levels of each PS per ontogenetic stage averaged across all brain regions

M1C_re <- M1C$re
M1C_re_dat <- M1C_re$data
colnames(M1C_re_dat)[3] <- "M1C"

DFC_re <- DFC$re
DFC_re_dat <- DFC_re$data
colnames(DFC_re_dat)[3] <- "DFC"

VFC_re <- VFC$re
VFC_re_dat <- VFC_re$data
colnames(VFC_re_dat)[3] <- "VFC"

OFC_re <- OFC$re
OFC_re_dat <- OFC_re$data
colnames(OFC_re_dat)[3] <- "OFC"

S1C_re <- S1C$re
S1C_re_dat <- S1C_re$data
colnames(S1C_re_dat)[3] <- "S1C"

IPC_re <- IPC$re
IPC_re_dat <- IPC_re$data
colnames(IPC_re_dat)[3] <- "IPC"

A1C_re <- A1C$re
A1C_re_dat <- A1C_re$data
colnames(A1C_re_dat)[3] <- "A1C"

STC_re <- STC$re
STC_re_dat <- STC_re$data
colnames(STC_re_dat)[3] <- "STC"

ITC_re <- ITC$re
ITC_re_dat <- ITC_re$data
colnames(ITC_re_dat)[3] <- "ITC"

V1C_re <- V1C$re
V1C_re_dat <- V1C_re$data
colnames(V1C_re_dat)[3] <- "V1C"

MFC_re <- MFC$re
MFC_re_dat <- MFC_re$data
colnames(MFC_re_dat)[3] <- "MFC"

HIP_re <- HIP$re
HIP_re_dat <- HIP_re$data
colnames(HIP_re_dat)[3] <- "HIP"

STR_re <- STR$re
STR_re_dat <- STR_re$data
colnames(STR_re_dat)[3] <- "STR"

AMY_re <- AMY$re
AMY_re_dat <- AMY_re$data
colnames(AMY_re_dat)[3] <- "AMY"

MD_re <- MD$re
MD_re_dat <- MD_re$data
colnames(MD_re_dat)[3] <- "MD"

CBC_re <- CBC$re
CBC_re_dat <- CBC_re$data
colnames(CBC_re_dat)[3] <- "CBC"


re_all_struct <- M1C_re_dat %>%
  left_join(DFC_re_dat, by=c('age', 'stage')) %>%
  left_join(VFC_re_dat, by=c('age', 'stage')) %>%
  left_join(OFC_re_dat, by=c('age', 'stage')) %>%
  left_join(S1C_re_dat, by=c('age', 'stage')) %>%
  left_join(IPC_re_dat, by=c('age', 'stage')) %>%
  left_join(A1C_re_dat, by=c('age', 'stage')) %>%
  left_join(STC_re_dat, by=c('age', 'stage')) %>%
  left_join(ITC_re_dat, by=c('age', 'stage')) %>%
  left_join(V1C_re_dat, by=c('age', 'stage')) %>%
  left_join(MFC_re_dat, by=c('age', 'stage')) %>%
  left_join(HIP_re_dat, by=c('age', 'stage')) %>%
  left_join(STR_re_dat, by=c('age', 'stage')) %>%
  left_join(AMY_re_dat, by=c('age', 'stage')) %>%
  left_join(MD_re_dat, by=c('age', 'stage')) %>%
  left_join(CBC_re_dat, by=c('age', 'stage')) 


re_all_struct$sum = rowSums(re_all_struct[,c("M1C", "DFC", "VFC", "OFC", "S1C", "IPC", "A1C", "STC", "ITC", "V1C", "MFC", "HIP", "STR", "AMY", "MD", "CBC")])
re_all_struct$all_struct_mean <- re_all_struct$sum / 16
re_all_struct <- re_all_struct[, -c(3:19)]
names(re_all_struct)[names(re_all_struct) == "age"] <- "Phylostratum"
re_all_struct$Phylostratum <- paste("PS", re_all_struct$Phylostratum, sep="")

re_all_struct$stage <- as.factor(re_all_struct$stage)
re_all_struct$Phylostratum <- factor(re_all_struct$Phylostratum, levels=c("PS1", "PS2", "PS3", "PS4", "PS5", "PS6","PS7", "PS10", "PS12"))


p_re <- ggplot(re_all_struct, aes(x= stage, y=all_struct_mean, group=Phylostratum, colour=Phylostratum)) +
  geom_line(size = 2) +
  theme_classic() +
  ylab("Relative brain-wide mRNA intensities") + 
  xlab("\nOntogenetic stage") +
  theme(axis.title=element_text(size=20,face="bold"),
        axis.text.y = element_text(size=15),
        axis.text.x = element_text(size = 15, angle=40, hjust=1, vjust = 1),
        legend.title = element_text(size=16),
        legend.text = element_text(size=13)) +
  viridis::scale_color_viridis(option = 'mako', 
                               discrete = T)

p_re







exp <- get_expression(structure_ids=c("Allen:10163"),
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

HomoSapiens.PhyloMap <- read_excel("D:/eigene_dateien/02Dokumente/07Erwerbstaetigkeit/01PhD_UiO/projects/first_author/oxytocinHealth/RQ4/github_dan/MBE_2008_Homo_Sapiens_PhyloMap.xls", 
                                   sheet = 1, skip = 1)

HomoSapiens.PhyloMap <- HomoSapiens.PhyloMap[ , c("Phylostratum","Gene_ID")]

mm <- MatchMap(HomoSapiens.PhyloMap, ex_t_f)

RECT <- EarlyConservationTest(ExpressionSet = PhyloExpressionSetExample, 
                              modules= list(early= 1:2, mid=3:5, late=6:7))

RECT # p = 1, TAI does not follow an early conservation like pattern. What does this mean?
data(PhyloExpressionSetExample)
View(PhyloExpressionSetExample)
# This dataset covers 7 developmental stages of Arabidopsis thaliana EMBRYO development, so logically, the conservation
# test and hourglass test can be computed (p. 31)






# first way to complicated version copied from source code

#age_table <- table(M1C_enrich$mm_all[, 1])
#nPS <- length(age_table)
#UpRegulated <- which(M1C_enrich$res_m[, 2] >= 0)
#DownRegulated <- which(M1C_enrich$res_m[, 2] < 0)


#bar.colors <- vector("character", nPS)
#bar.colors[UpRegulated] <- "#FF0564"
#bar.colors[DownRegulated] <- "#3D477B"

#ylim.range_alt <-  range(floor(min(M1C_enrich$res_m[, 2])), ceiling(max(M1C_enrich$res_m[, 2])))

#barPlotFoldChanges_M1C <- graphics::barplot(M1C_enrich$res_m[, 2], beside = FALSE, col = bar.colors, lwd = 2, 
#                                          space = 1, names.arg = paste0("PS", names(age_table)), 
#                                          ylim = ylim.range_alt, border = NA, xlab = "Phylostrata", cex.lab = 1.5, 
#                                          las=1)

# how can I limit the length of the abline??
#graphics::abline(h = 0, col = "black", lty = 1, lwd = 2)
#graphics::legend("topright", 
#                 legend = c("Over-represented", "Under-represented"), 
#                 fill = c("#FF0564", "#3D477B"), 
#                 bty = "n", cex = 1,
#                 border = NA)
#pValNames <- vector("character", nPS)
#pValNames <- rep("", nPS)

#enrichment.p_valsM1C <- stats::p.adjust(M1C_enrich$res_pvals, 
#                                       method = "fdr", n = nPS)

#pValNames[which(enrichment.p_valsM1C <= 0.05)] <- "*"
#pValNames[which(enrichment.p_valsM1C <= 0.005)] <- "**"
#pValNames[which(enrichment.p_valsM1C <= 5e-04)] <- "***"
#pValNames[which(is.na(pValNames))] <- ""
#graphics::text(x = barPlotFoldChanges_M1C, 
#              y = ylim.range_alt[2] - .1, labels = pValNames, cex = 1)



writeLines('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', con = "~/.Renviron")
Sys.which("make")
install.packages("jsonlite", type = "source")
# worked!



## ==================================== NSS score

library(data.table)
hs_vs_ns_align <- fread("D:/eigene_dateien/02Dokumente/07Erwerbstaetigkeit/01PhD_UiO/projects/first_author/oxytocinHealth/RQ4/selective_sweep/data/ntSssZScorePMVar.bed")
hs_vs_ns_align <- data.frame(hs_vs_ns_align)
names(hs_vs_ns_align)[names(hs_vs_ns_align) == "V1"] <- "chr"
names(hs_vs_ns_align)[names(hs_vs_ns_align) == "V2"] <- "bp_start"
names(hs_vs_ns_align)[names(hs_vs_ns_align) == "V3"] <- "bp_end"
names(hs_vs_ns_align)[names(hs_vs_ns_align) == "V4"] <- "nss_score"

# create working copy
hsns_align <- hs_vs_ns_align

# check chromosome availability
unique(hsns_align$chr) # 1-22, and X



## ================================= dN/dS values from Dumas et al., 2021 (working draft)

library(data.table)
library(readr)
library(ggplot2)
library(tidyverse)


# load file generated by https://genevo.pasteur.fr/, entered the OT pathway geneset as gene list

dnds_primates_OT <- fread(paste0(BASE, "dNdS-ratio_TDI/genevo_2021-04-28.tsv"))
View(dnds_primates_OT)

ids_OT <- data.frame(dnds_primates_OT$`Entrez ID`)

# generally, all the genes from the OT pathway geneset in HUMANS compared to a common primate ancestor have a dN/dS ratio < 1, indicating evolutionary constraint and high conservation.
# For direct comparison against the neanderthals, we can use the NSS score. 




# Basic box plot OLD VERSION w/o whisker length adjustment

subset_hq <- dnds_primates_OT[, c(2, 6, 9, 12, 15, 18, 21, 24, 27)]

colnames(subset_hq) <- c("Gene",
                         "H. sapiens", 
                         "H. neanderthalensis", 
                         "H. denisova", 
                         "P. troglodytes", 
                         "G. gorilla", 
                         "P. abelii", 
                         "M. mulatta", 
                         "C. jacchus")

subset_hq <- subset_hq %>% pivot_longer(c("H. sapiens", 
                                          "H. neanderthalensis", 
                                          "H. denisova", 
                                          "P. troglodytes", 
                                          "G. gorilla",
                                          "P. abelii",
                                          "M. mulatta",
                                          "C. jacchus"),
                                        names_to = "Species", values_to = "value")


subset_hq$Species <- factor(subset_hq$Species , levels = c("H. sapiens", 
                                                           "H. neanderthalensis", 
                                                           "H. denisova", 
                                                           "P. troglodytes", 
                                                           "G. gorilla",
                                                           "P. abelii",
                                                           "M. mulatta",
                                                           "C. jacchus"))


subset_hq$Gene <- as.factor(subset_hq$Gene)

p_hq <- ggplot(subset_hq, aes(x = value, y = Species, fill = Species)) +
  geom_boxplot( show.legend = FALSE, alpha=0.5) 
p_hq <- p_hq + geom_point(aes(fill=Species), size = 5, colour="black", shape = 21, position = position_jitterdodge()) 
p_hq <- p_hq +  ylab("") +  xlab("dN/dS ratio") 
p_hq <- p_hq + theme(axis.title=element_text(size=20,face="bold"),
                     axis.text.y = element_text(size=15),
                     axis.text.x = element_text(size = 15),
                     legend.title = element_text(size=16)) 
p_hq <- p_hq + viridis::scale_fill_viridis(option = 'viridis', discrete = T)

p_hq <- p_hq + theme(legend.position = "none")

p_hq

ggsave(filename = paste0(BASE, "RProject/OT-bd-analyses/output/dnds_OTgeneset_hq.pdf"), p_hq,
       width = 15, height = 8, units = "in", device='pdf')



##########################

p_lq_2 <- ggplot(subset_lq, aes(x = value, y = Species, fill = Species)) 
p_lq_2 <- p_lq_2 +  
  stat_summary(fun.data=f, geom="boxplot", alpha=0.5, show.legend = FALSE, width=.75) +
  stat_summary(fun = o, geom="point") 
p_lq_2 <- p_lq_2 + scale_y_discrete(limits = rev(levels(subset_lq$Species)))
p_lq_2 <- p_lq_2 + geom_point(aes(fill=Species), size = 4, colour="Black", shape = 21, position = position_jitterdodge()) 
p_lq_2 <- p_lq_2 +  ylab("") +  xlab("dN/dS ratio") 
p_lq_2 <- p_lq_2 + theme(axis.text.y= element_blank(),         # removes y axis text
                     axis.text.x = element_blank(),
                     axis.ticks.y = element_blank(),
                     axis.title.x = element_blank(),
                     axis.title.y = element_blank(),
                     axis.ticks.x = element_blank()) 


p_lq_2 <- p_lq_2 + viridis::scale_fill_viridis(option = 'rocket', discrete = T)

p_lq_2 <- p_lq_2 + theme(legend.position = "none")

p_lq_2


#install.packages("phytools")
library(phytools)
## Loading required package: maps
#BiocManager::install("ggtree")
library(ggtree)
library(gridExtra)
library(grid)
library(lattice)




# w/o x axis
#gg_tr_2 <- ggtree(treee) + geom_tiplab(align=TRUE) +
#  scale_x_continuous(expand=expansion(1.5)) 
#gg_tr_2



lq_tree <- gg_tr + p_lq
lq_tree
ggsave(filename = paste0(BASE, "RProject/OT-bd-analyses/output/lq_tree.pdf"), lq_tree,
       width = 15, height = 8, units = "in", device='pdf')


## OR

tree6 <- grid.arrange(gg_tr, p_lq_2, ncol=2, widths = 1:2) #, 
             #heights=unit(c(), c("in", "mm")))
ggsave(filename = paste0(BASE, "RProject/OT-bd-analyses/output/tree6.pdf"), tree6,
       width = 15, height = 8, units = "in", device='pdf')


## OR

pushViewport(viewport(layout = grid.layout(1, 2)))
print(gg_tr, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p_lq, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
gg_tr
p_lq


## OR

ggarrange(gg_tr, p_lq_2, ncol = 2)
gg_tr_2 <- ggtree(treee) + geom_tiplab(align=TRUE) +
  scale_x_continuous(expand=expansion(1.5)) +
  scale_y_continuous(expand=expansion(.2))
gg_tr_2

# ways to generate a phylo tree
rtree(treee)
ggtree(treee) + geom_tiplab(align=T) +
scale_x_continuous(expand=expansion(1.5))
plot(treee,no.margin=T,edge.width=2)
roundPhylogram(treee)
plotTree(treee)



butt <- p_hq 
butt

bum <- gg_tr 
bum

peach <- grid.arrange(bum, butt,
                      layout_matrix = rbind(c(1, 2, 2, 2, 2),
                                            c(1, 2, 2, 2, 2)))

apple <- grid.arrange(bum, butt, layout_matrix = rbind(c(1, 1, 2, 2, 2, 2, 2,2,2),
                                                       c(1, 1, 2, 2, 2, 2, 2,2,2)))


ggsave(filename = paste0(BASE, "RProject/OT-bd-analyses/output/apple3.pdf"), apple,
       width = 20, height = 8, units = "in", device='pdf')



#tree_2 <- ape::read.tree(text='(C._jacchus,(M._mulatta,(P._abelii,(G._gorilla,(P._troglodytes,(H._denisova,(H._neanderthalensis,H._sapiens)))))));') 

#ape::write.nexus(tree_2, file='E:/Other_data/test.nex')
#tree_21 <- ape::read.tree("E:/Other_data/test.nex")


#browseVignettes("ggtree")
#library(ggstance)
#pt_test <- facet_plot(p_tree, panel = 'Barplot', data = subset_lq, 
#                 geom = geom_boxploth, 
#                 mapping = aes(x = value, y = Species, fill = Species))





# %%%%%%%%%%%% old tree %%%%%%%%%%%%%%% #

# create phylogenetic template tree
# template/orientation, generated with https://www.megasoftware.net/: 1. download software, 2. create custom tree with option "user tree"
tree_t <- fread(paste0(BASE, "RProject/OT-bd-analyses_loc/Data/phylo_tree_primates"), header = F)
View(tree_t)
# generate tree accordingly
text.string<-
  "(C. jacchus,(M. mulatta,(P. abelii,(G. gorilla,(P. troglodytes,(H. denisova,(H. neanderthalensis,H. sapiens)))))));"
treee<-read.tree(text=text.string)
ape::write.nexus(treee, file=paste0(BASE, "RProject/OT-bd-analyses_loc/output/primate_pt.nex"))




############################### 
###############################


supp_tab1 <- read_csv(paste0(BASE, "dNdS-ratio_TDI/Dumas_2020/supplementary_material_publishedPaper/Supplemental_Tables/Supplemental_Table_S1.csv"), skip=1, col_names = T)
View(supp_tab1)

#
supp_tab3 <- read_csv(paste0(BASE, "dNdS-ratio_TDI/Dumas_2020/supplementary_material_publishedPaper/Supplemental_Tables/Supplemental_Table_S3.csv"), skip=1, col_names = T)
supp_tab3$asterisk <- ""
supp_tab3$asterisk[supp_tab3$bonferroni == "TRUE"]  <- "*"
View(supp_tab3)

# supplementary table for gene set with top SCG. These genes were more associated with brain disease, particularly intellectual disability and autism. 
# ... ok, cool, but which genes with which disease?
supp_tab4 <- read_csv(paste0(BASE, "dNdS-ratio_TDI/Dumas_2020/supplementary_material_publishedPaper/Supplemental_Tables/Supplemental_Table_S4.csv"), skip=1, col_names = T)
View(supp_tab4)
ids_supp4 <- data.frame(supp_tab4$EntrezId)
ids_supp4$diseases <- supp_tab4$Diseases


# compare OT gene set with supp4 gene set and check disease relevance
ids_supp4OT <- merge(ids_supp4, ids_OT, by.x="supp_tab4.EntrezId", by.y="dnds_primates_OT..Entrez.ID.")
View(ids_supp4OT)

# --> 10 matches



names <- data.frame(colnames(temp2))








temp2$CD38 <- as.numeric(temp2$CD38)
CD38_mean <- mean(temp2$CD38)
CD38_regions <- unique(temp2$region)
CD38_pvals <- numeric()
CD38_t <- numeric()
CD38_cohd <- numeric()

for (i in 1:length(CD38_regions)){
  Temp <- t.test(temp2[temp2$region==CD38_regions[i],"CD38"], 
                 mu = CD38_mean, alternative = "two.sided")
  CohD <- ES.t.one(m=Temp$estimate,
                   sd=sd(temp2[temp2$region==CD38_regions[i],"CD38"]),
                   mu=CD38_mean,
                   alternative = "two.sided")
  CD38_pvals <- c(CD38_pvals,Temp$p.value)
  CD38_t <- c(CD38_t,Temp$statistic)
  CD38_cohd <- c(CD38_cohd,CohD$d)
}

CD38_results <- data.frame(Region=CD38_regions,P_unadjusted=CD38_pvals,
                           t_stat=CD38_t,CohD=CD38_cohd)
CD38_results$P_fdr <- p.adjust(CD38_results$P_unadjusted, method="fdr")

CD38_results 

#ASB3_sum <- summarySE(all_genes_l, measurevar="ASB3_mRNA", 
#groupvars=c("area","region"), na.rm = T)

#ASB3_sum <- merge(ASB3_sum, ASB3_results[,c("Area", "t_stat", "CohD", "P_unadjusted", "P_fdr")], by.x = "area", by.y = "Area")
CD38_results$stat_sign <- ""
CD38_results$stat_sign[CD38_results$P_fdr <= 0.05]  <- "*"



#################################################
#################################################
#################################################



# Prepare variables

M1Csgs <- M1C$sgs
M1Csgs$Region <- "M1C"

DFCsgs <- DFC$sgs
DFCsgs$Region <- "DFC"

VFCsgs <- VFC$sgs
VFCsgs$Region <- "VFC"

OFCsgs <- OFC$sgs
OFCsgs$Region <- "OFC"

S1Csgs <- S1C$sgs
S1Csgs$Region <- "S1C"

IPCsgs <- IPC$sgs
IPCsgs$Region <- "IPC"

A1Csgs <- A1C$sgs
A1Csgs$Region <- "A1C"

STCsgs <- STC$sgs
STCsgs$Region <- "STC"

ITCsgs <- ITC$sgs
ITCsgs$Region <- "ITC"

V1Csgs <- V1C$sgs
V1Csgs$Region <- "V1C"

MFCsgs <- MFC$sgs
MFCsgs$Region <- "MFC"

HIPsgs <- HIP$sgs
HIPsgs$Region <- "HIP"

STRsgs <- STR$sgs
STRsgs$Region <- "STR"

AMYsgs <- AMY$sgs
AMYsgs$Region <- "AMY"

MDsgs <- MD$sgs
MDsgs$Region <- "MD"

CBCsgs <- CBC$sgs
CBCsgs$Region <- "CBC"


bind <- M1Csgs %>%
  rbind(DFCsgs) %>%
  rbind(VFCsgs) %>%
  rbind(OFCsgs) %>%
  rbind(S1Csgs) %>%
  rbind(IPCsgs) %>%
  rbind(A1Csgs) %>%
  rbind(STCsgs) %>%
  rbind(ITCsgs) %>%
  rbind(V1Csgs) %>%
  rbind(MFCsgs) %>%
  rbind(HIPsgs) %>%
  rbind(STRsgs) %>%
  rbind(AMYsgs) %>%
  rbind(MDsgs) %>%
  rbind(CBCsgs)

bind_subset <- bind %>% pivot_longer(c("Prenatal", 
                                     "Infant", 
                                     "Child", 
                                     "Adolescent", 
                                     "Adult"),
                                   names_to = "Stages", values_to = "mRNA")
bind_subset$Stages <- factor(bind_subset$Stages , levels = c("Prenatal", 
                                                             "Infant", 
                                                             "Child", 
                                                             "Adolescent", 
                                                             "Adult"))
bind_subset$Region <- as.factor(bind_subset$Region)

View(sgs_subset)

px <- ggplot(bind_subset, aes(x = Stages, y = mRNA, group=Region, colour=Region)) +
  geom_line(size=1.5) + 
  theme_classic() 
px <- px + xlab("\nOntogenetic stages") + ylab("OXTR mRNA intensity\n")
px <- px + theme(axis.title.x = element_text(size=20, face="bold"),
                 axis.text.x = element_text(size=15, angle=40, hjust=1, vjust = 1),
                 axis.title.y = element_text(size=20, face="bold"),
                 axis.text.y = element_text(size=15),
                 legend.position = c(.94, .5))
px <- px + viridis::scale_colour_viridis(option = 'mako', discrete = T)
#px <- px + geom_smooth(method="loess")

px <- px + stat_summary(aes(y = mRNA,group=1), fun=mean, colour="#F6044F", geom="line",group=1, size=3)
px






oxtr_rib <- ggplot(bind_oxtr, aes(x=as.numeric(factor(Stages)), y = mRNA, fill=Region)) +
  geom_area(lwd=1, 
            linetype=1, 
            position='stack', 
            colour='gray90', 
            alpha=.8,  
            na.rm=TRUE) +
  scale_x_discrete(labels = levels(bind_oxtr$Stages)) #+ 
  #geom_dl(aes(label = Region), list("top.points", cex = .6)) + 
  #guides(fill = FALSE)
oxtr_rib <- oxtr_rib + viridis::scale_fill_viridis(option = 'rocket', discrete = T) 
oxtr_rib <- oxtr_rib + xlab("Stages") + ylab("OXTR mRNA intensities") 
oxtr_rib <- oxtr_rib + theme(legend.position = "right",
                             axis.title.y = element_text(size=20, face="bold"))
oxtr_rib <- oxtr_rib + theme_minimal()
oxtr_rib




####################################
####################################
####################################

# >>>>>>>>>>>>>>> abagen <<<<<<<<<<<<<<<<<<

## IMPORTANT INFO ##
# This part of the script is using the toolbox 'abagen' (https://github.com/rmarkello/abagen), which is written in PYTHON.
# Do not attempt to run in R, it will not work. Run in python shell in the terminal or in Jupyter Lab/Notebook.
# Use 'abagen_analysis.ipynb' for detailed information, pre-requisites, dependencies, etc.
# Most importantly, read https://abagen.readthedocs.io/en/stable/ carefully for optimal use. 


# NOT RUN {
# initiate abagen 
import abagen

# fetch atlas/parcellation
atlas = abagen.fetch_desikan_killiany()

# get expression data and do the parcellation 
# --> IF data is already downloaded to local machine, specify the directory it is stored in with " data_dir='your/local/path/microarray/' " parameter
# --> otherwise, the function will download the data per default
expression_main_mis, report_main_mis = abagen.get_expression_data(atlas['image'], atlas['info'],
                                                                  missing = 'interpolate',       # fill missing values
                                                                  norm_matched=False,            # necessary for the missing parameter to work
                                                                  return_report=True)            # returns report about what was done

# check if everything worked
print(expression_main_mis)

# export to csv
expression_main.to_csv('C:/01Dokumente/07Erwerbstaetigkeit/01PhD_UiO/projects/first_author/oxytocinHealth/RQ4/data/processed/expression_main_mis.csv', 
                        index = False, header=True)

print(report_main_mis)
# }

#############################
##### end python script #####
#############################

exp_main_mis <- read_csv("C:/01Dokumente/07Erwerbstaetigkeit/01PhD_UiO/projects/first_author/oxytocinHealth/RQ4/data/processed/expression_main_mis.csv")
dim(exp_main_mis)



## M1C

re1<-M1C$re
re1dat <- re1$data
 
re1dat$group <- re1dat$age
re1dat$group <- as.factor(re1dat$group)

pattern_old = c("1", "2", "3")
re1dat$group <- str_replace_all(re1dat$group, pattern_old, "old")

pattern_young = c("4|5|6|7|10|12")
re1dat$group <- str_replace_all(re1dat$group, pattern_young, "young")



pattern_old = c("1", "2", "3")
pattern_you = c("4", "5", "6", "7", "10", "12")



# create new df

## ----------------------------- STR

m_pre_o <- with(STRre_dat, mean(expr[age %in% pattern_old & stage == "Prenatal"])) 
m_pre_y <- with(STRre_dat, mean(expr[age %in% pattern_you & stage == "Prenatal"]))
m_inf_o <- with(STRre_dat, mean(expr[age %in% pattern_old & stage == "Infant"])) 
m_inf_y <- with(STRre_dat, mean(expr[age %in% pattern_you & stage == "Infant"]))
m_chi_o <- with(STRre_dat, mean(expr[age %in% pattern_old & stage == "Child"])) 
m_chi_y <- with(STRre_dat, mean(expr[age %in% pattern_you & stage == "Child"]))
m_ado_o <- with(STRre_dat, mean(expr[age %in% pattern_old & stage == "Adolescent"])) 
m_ado_y <- with(STRre_dat, mean(expr[age %in% pattern_you & stage == "Adolescent"]))
m_adu_o <- with(STRre_dat, mean(expr[age %in% pattern_old & stage == "Adult"])) 
m_adu_y <- with(STRre_dat, mean(expr[age %in% pattern_you & stage == "Adult"]))

exp_str <- m_pre_o %>%
  rbind(m_pre_y) %>%
  rbind(m_inf_o) %>%
  rbind(m_inf_y) %>%
  rbind(m_chi_o) %>%
  rbind(m_chi_y) %>%
  rbind(m_ado_o) %>%
  rbind(m_ado_y) %>%
  rbind(m_adu_o) %>%
  rbind(m_adu_y)

exp_str <- data.frame(exp_str)
exp_str$PSclass <- rep(c("Old", "Young"), 5)
exp_str$stage <- c(rep("1. Prenatal", 2), rep("2. Infant", 2), rep("3. Child",2), rep("4. Adolescent",2), rep("5. Adult",2))

STRre_p <- 
  exp_str %>% 
  group_by(stage) %>% 
  mutate(sd = sd(exp_str)) %>% 
  ggplot(aes(x = stage, y = exp_str, fill=PSclass)) + 
  geom_bar(stat="identity", alpha=0.75, position=position_dodge()) 
STRre_p <- STRre_p + ylim(0, 1.25) + ggtitle("S1C") + theme_minimal()
STRre_p <- STRre_p + ylab("Relative expression intensity\n") + xlab("Ontogenetic stage") +
  labs(fill = "Phylostratum \nage category") 
STRre_p <- STRre_p + theme(axis.title.x = element_text(size=17, face="bold"),
                           axis.text.x = element_text(size=11, angle=40), 
                           axis.title.y = element_text(size=17, face="bold"))
STRre_p <- STRre_p + scale_fill_manual(values=c("#20244C", "#A7DBFF"))
STRre_p <- STRre_p + geom_errorbar(aes(ymin=exp_str, ymax=exp_str+sd), width=.2, colour="orange", 
                                   position=position_dodge(.9))
STRre_p






# dk cortical atlas: https://www.sciencedirect.com/science/article/abs/pii/S1053811906000437 

# aseg subcortical atlas: https://pubmed.ncbi.nlm.nih.gov/11832223/ 



# include/add spaces in the raw data manually.. is there a better way of doing this? removing the spaces in the data included in ggseg?
dk_info2 <- read_csv(paste0(BASE, "RProject/OT-bd-analyses_loc/Data/atlas-desikankilliany_spaces.csv"))
names(dk_info2)[names(dk_info2) == "label"] <- "region"



### create cortical and subcortical dfs

rs <- c("pallidum", "thalamus proper", "caudate", "putamen", "accumbensarea", "hippocampus", "amygdala", "brain stem") # subcortical regions

# -------------- cortical
dk_info2_c <- dk_info2[!dk_info2$region %in% rs,] # remove subcortical regions

exp_m_anno_c <- dplyr::inner_join(dk_info2_c, exp_main_mis, by = c("id" = "label")) # annotate expression matrix
dim(exp_m_anno_c) # 15633 genes, 68 regions (multilateral)

# -------------- subcortical
dk_info2_sc <- dk_info2[dk_info2$region %in% rs,] # keep subcortical regions
dk_info2_sc <- dk_info2_sc[!dk_info2_sc$region=="accumbensarea",]

exp_m_anno_sc <- dplyr::inner_join(dk_info2_sc, exp_main_mis, by = c("id" = "label")) # annotate expression matrix
dim(exp_m_anno_sc) # 15633 genes, 13 regions (multilateral)




### create subsets with all OT pathway genes, modern genes and ancient genes for cortical and subcortical

modPS <- c("4", "5", "6", "7", "10", "12")

## oT pathway gene subset *ALL* PS
#exp_ma_OTall <- exp_m_anno2[, names(exp_m_anno2) %in% c("id", "region", "hemisphere", "structure", kop)] # select OT subset
#dim(exp_ma_OTall) # 136 genes, 17 skipped


## *MODERN* OT Pathway gene subset (PS4-12)
modern_OT <- M1C_mm[M1C_mm$Phylostratum %in% modPS, ] # keep only modern phylostrata
modern_OT$GeneID <- toupper(modern_OT$GeneID)


# annotate with gene name
ens_gene_mod <- ensembldb::select(EnsDb.Hsapiens.v79, 
                                  keys= modern_OT$GeneID, 
                                  keytype = "GENEID", 
                                  columns = c("SYMBOL","GENEID"))
ens_gene_mod <- ens_gene_mod %>% distinct(SYMBOL, .keep_all = TRUE) 
geneids_mod <- ens_gene_mod$SYMBOL  # reference list with gene names

# cortical
OTmod_c <- exp_m_anno_c[, names(exp_m_anno_c) %in% c("id", "region", "hemisphere", "structure", geneids_mod)] 
dim(OTmod_c) # 19 modern genes, 68 cortical regions

# subcortical
OTmod_sc <- exp_m_anno_sc[, names(exp_m_anno_sc) %in% c("id", "region", "hemisphere", "structure", geneids_mod)] 
dim(OTmod_sc) # 19 modern genes, 13 subcortical regions



## *ANCIENT* OT Pathway gene subset (PS1-3)
ancient_OT <- M1C_mm[!M1C_mm$Phylostratum %in% modPS, ] # remove modern phylostrata
ancient_OT$GeneID <- toupper(ancient_OT$GeneID)

# annotate with gene name
ens_gene_old <- ensembldb::select(EnsDb.Hsapiens.v79, 
                                  keys= ancient_OT$GeneID, 
                                  keytype = "GENEID", 
                                  columns = c("SYMBOL","GENEID"))
ens_gene_old <- ens_gene_old %>% distinct(SYMBOL, .keep_all = TRUE) 
geneids_old <- ens_gene_old$SYMBOL  # reference list with gene names

# cortical
OTold_c <- exp_m_anno_c[, names(exp_m_anno_c) %in% c("id", "region", "hemisphere", "structure", geneids_old)] 
dim(OTold_c) # 107 ancient genes, 68 cortical regions

# subcortical
OTold_sc <-    exp_m_anno_sc[, names(exp_m_anno_sc) %in% c("id", "region", "hemisphere", "structure", geneids_old)] 
dim(OTold_sc) # 107 ancient genes, 13 subcortical regions





# ------------ test for OXTR
# cortical modern
oxtr_oldc <- ggseg(.data=OTmod_c, mapping=aes(fill=OXTR), colour="black", size=.7, position="stacked", hemisphere="left" ) +
  theme_void()
oxtr_oldc <- oxtr_oldc + viridis::scale_fill_viridis(option = 'turbo', discrete = F)
oxtr_oldc

# subcortical modern
oxtr_modsc <- ggseg(.data=OTmod_sc, atlas="aseg", mapping=aes(fill=OXTR), colour="black", size=.7, position="stacked") +
  theme_void()
oxtr_modsc <- oxtr_modsc + viridis::scale_fill_viridis(option = 'turbo', discrete = F)
oxtr_modsc




# test for OXT
oxt <- ggseg(.data=exp_ma_OT2, mapping=aes(fill=OXT), colour="gray50", size=1, position="stacked",hemisphere="left" ) +
  theme_void()
oxt <- oxt + viridis::scale_fill_viridis(option = 'turbo', discrete = F)
oxt

# test for CD38
cd38 <- ggseg(.data=exp_ma_OT2, mapping=aes(fill=CD38), colour="gray50", size=1, position="stacked",hemisphere="left" ) +
  theme_void()
cd38 <- cd38 + viridis::scale_fill_viridis(option = 'turbo', discrete = F)
cd38





# default

plot(aseg)




# additional information on atlases

remotes::install_github("LCBC-UiO/ggsegExtra", build_vignettes = F)
library(ggsegExtra)
ggseg_atlas_repos()

ggseg$labels
plot(aseg)
help(package = ggseg)



ggseg
c("3rd ventricle", "4th ventricle", "CC anterior", "CC central", "cc mid anterior", 
  "CC mid posterior", "CC posterior", "cerebellum cortex", "cerebellum white matter", 
  "lateral ventricle", "ventral DC")




  
  
c_exp <- function(gene) {
  gene <- enquo(gene)
  cortical <- ggseg(.data=OTmod_c, mapping=aes(fill=!!gene), colour="black", size=.7, position="stacked", hemisphere="left" ) +
    theme_void()
  cortical <- cortical + viridis::scale_fill_viridis(option = 'turbo') + 
    labs(title="OXTR cortical expression distribution", subtitle="lateral | medial", fill="mRNA intensity") +
    theme(plot.title = element_text(hjust = .5),
          plot.subtitle = element_text(hjust = .5))
  value <- list(                                    
    cortical = cortical
  ) # Create a list of output objects
  attr(value, "class") <- "c_exp"
  value
}

cacnb1_c <- c_exp(CACNB1)
cacnb2_c <- c_exp(CACNB2)
cacnb3_c <- c_exp(CACNB3)
cacnb4_c <- c_exp(CACNB4)
oxtr_c <- c_exp(OXTR)
elk1_c <- c_exp(ELK1)
nppa_c <- c_exp(NPPA)
nfatc2_c <- c_exp(NFATC2)
nfatc3_c <- c_exp(NFATC3)
nfatc4_c <- c_exp(NFATC4)
pik3r5_c <- c_exp(PIK3R5)
oxt_c <- c_exp(OXT)
cd38_c <- c_exp(CD38)
fos_c <- c_exp(FOS)
cdkn1a_c <- c_exp(CDKN1A)
cacng1_c <- c_exp(CACNG1)
cacng2_c <- c_exp(CACNG2)
cacng3_c <- c_exp(CACNG3)
cacng6_c <- c_exp(CACNG6)
cacng7_c <- c_exp(CACNG7)

sc_exp <- function(gene) {
  gene <- enquo(gene)
  
  subcortical <- ggseg(.data=OTmod_sc, atlas="aseg", mapping=aes(fill=!!gene), colour="black", size=.7) +
    theme_void()
  
  subcortical <- subcortical + viridis::scale_fill_viridis(option = 'turbo', discrete = F) + 
    labs(title="OXTR subcortical expression", subtitle = "axial | sagittal", fill="mRNA intensity") 
  
  subcortical <- subcortical + theme(plot.title = element_text(hjust=0.5),
                                     plot.subtitle = element_text(hjust=0.515),
                                     legend.text=element_text(size=10)) 
  
  value <- list(                                    
    subcortical = subcortical
  ) # Create a list of output objects
  attr(value, "class") <- "sc_exp"
  value
}

cacnb1_sc <- sc_exp(CACNB1)
cacnb2_sc <- sc_exp(CACNB2)
cacnb3_sc <- sc_exp(CACNB3)
cacnb4_sc <- sc_exp(CACNB4)
oxtr_sc <- sc_exp(OXTR)
elk1_sc <- sc_exp(ELK1)
nppa_sc <- sc_exp(NPPA)
nfatc2_sc <- sc_exp(NFATC2)
nfatc3_sc <- sc_exp(NFATC3)
nfatc4_sc <- sc_exp(NFATC4)
pik3r5_sc <- sc_exp(PIK3R5)
oxt_sc <- sc_exp(OXT)
cd38_sc <- sc_exp(CD38)
fos_sc <- sc_exp(FOS)
cdkn1a_sc <- sc_exp(CDKN1A)
cacng1_sc <- sc_exp(CACNG1)
cacng2_sc <- sc_exp(CACNG2)
cacng3_sc <- sc_exp(CACNG3)
cacng6_sc <- sc_exp(CACNG6)
cacng7_sc <- sc_exp(CACNG7)


  
  

  

test <- res_modc
test[[20]]

names(test) <- c("CACNB1", "CACNB2", "CACNB3", "CACNB4", "CACNG1", "CACNG2", "CACNG3", "CACNG6", "CACNG7", "CD38", "CDKN1A", 
                 "ELK1", "FOS", "NFATC3", "NFATC4", "NPPA", "OXT", "OXTR", "PIK3R5")






fig <- OT_yc[, -c(12, 15, 26)]
fig$mean <- rowSums(fig[,5:23])/19

figp <- ggseg(.data=fig, mapping=aes_string(fill="mean"), colour="white", size=1, position="stacked", hemisphere="left" ) +
  theme_void()
figp <- figp + viridis::scale_fill_viridis(option = 'turbo') + 
  labs(title="Average cortical expression", subtitle="lateral | medial", fill="mRNA intensity") +
  theme(plot.title = element_text(hjust = .5, size=16),
        plot.subtitle = element_text(hjust = .5, size=14),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12))
figp





# get expression data on all structures for all genes in the Phylo Map
# get structure ids for the adult brain
ontology_AHBA <- read_csv(paste0(BASE, "data/allenbrainatlas/allenbrainatlas_structureontology.csv")) # for structure ontology reference
allen_ids <- data.frame(ontology_AHBA$id)
allen_ids$ontology_AHBA.id <- paste0("Allen:", allen_ids$ontology_AHBA.id)
feed <- allen_ids
feed <- feed[-25,]
feed <- data.frame(feed)
feed$feed <- as.character(feed$feed)

# get data
full_ahba <- get_expression(structure_ids=feed$feed,
                              gene_ids=genensid$SYMBOL, 
                              dataset='adult')

full_ahba <- as.data.frame(full_ahba)

full_ahba_t <- as.data.frame(t(full_ahba))

full_ahba_t <- rownames_to_column(full_ahba_t, var = "gene")



##################################################

# young OT pathway genes cortical expression averaged

# ...... all genes
yc_means <- OTmod_c
yc_means$sum = rowSums(yc_means[,list_g_yc])
yc_means$mean <- yc_means$sum / 19

genset_yc <- ggseg(.data=yc_means, mapping=aes_string(fill="mean"), colour="white", size=1, position="stacked", hemisphere="left" ) +
  theme_void()
genset_yc <- genset_yc + viridis::scale_fill_viridis(option = 'turbo') + 
  labs(title="Average cortical expression", subtitle="lateral | medial", fill="mRNA intensity") +
  theme(plot.title = element_text(hjust = .5, size=16),
        plot.subtitle = element_text(hjust = .5, size=14),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12))
genset_yc




# young OT pathway genes subcortical expression averaged

# ...... all genes
ysc_means <- OTmod_sc
ysc_means$sum = rowSums(ysc_means[,list_g_ysc])
ysc_means$mean <- ysc_means$sum / 19

genset_ysc <- ggseg(.data=ysc_means, atlas="aseg", mapping=aes_string(fill="mean"), colour="white", size=1) +
  theme_void()

genset_ysc <- genset_ysc + viridis::scale_fill_viridis(option = 'turbo', discrete = F) + 
  labs(title="Average subcortical expression", subtitle = "axial | sagittal", fill="mRNA intensity") 

genset_ysc <- genset_ysc + theme(plot.title = element_text(hjust=0.5, size=16),
                                 plot.subtitle = element_text(hjust=0.515, size=14),
                                 legend.title = element_text(size=14),
                                 legend.text = element_text(size=12)) 

genset_ysc



##################################################



# https://github.com/LCBC-UiO/ggseg/
# https://lcbc-uio.github.io/ggsegExtra/
# https://lcbc-uio.github.io/ggseg/articles/ggseg.html 


# load parcellated AHBA data with DK parcellation (this data should be available on your machine per default when you used the
# DK atlas for parcellation in the abagen toolbox from before)
exp_main_mis <- read_csv("C:/01Dokumente/07Erwerbstaetigkeit/01PhD_UiO/projects/first_author/oxytocinHealth/RQ4/data/processed/expression_main_mis.csv")
dim(exp_main_mis) #83 15634


# dk CORTICAL atlas: https://www.sciencedirect.com/science/article/abs/pii/S1053811906000437 

# aseg SUBCORTICAL atlas: https://pubmed.ncbi.nlm.nih.gov/11832223/ , https://freesurfer.net/fswiki/SubcorticalSegmentation 


# load additional DK labelling information
# include/add spaces in the raw data manually.. is there a better way of doing this? removing the spaces in the data included in ggseg?
dk_info <- read_csv(paste0(BASE, "RProject/OT-bd-analyses_loc/Data/atlas-desikankilliany_spaces.csv"))
names(dk_info)[names(dk_info) == "label"] <- "region"



### create cortical and subcortical dfs

rs <- c("pallidum", "thalamus proper", "caudate", "putamen", "accumbens area", "hippocampus", "amygdala", "brain stem") # subcortical regions

# -------------- cortical
dk_info_c <- dk_info[!dk_info$region %in% rs,]                                       # remove subcortical regions

exp_m_anno_c <- dplyr::inner_join(dk_info_c, exp_main_mis, by = c("id" = "label"))   # annotate expression matrix
dim(exp_m_anno_c) # 15633 genes, 68 regions (multilateral)

# -------------- subcortical
dk_info_sc <- dk_info[dk_info$region %in% rs,]                                       # keep subcortical regions
#dk_info_sc <- dk_info_sc[!dk_info_sc$region=="accumbensarea",]

exp_m_anno_sc <- dplyr::inner_join(dk_info_sc, exp_main_mis, by = c("id" = "label")) # annotate expression matrix
dim(exp_m_anno_sc) # 15633 genes, 15 regions (multilateral)





### create subsets with (all OT pathway genes,) modern genes and ancient genes for cortical and subcortical regions

M1C_mm <- M1C_tai$mm

modPS <- c("4", "5", "6", "7", "10", "12")

## oT pathway gene subset *ALL* PS
#exp_ma_OTall <- exp_m_anno2[, names(exp_m_anno2) %in% c("id", "region", "hemisphere", "structure", kop)] # select OT subset
#dim(exp_ma_OTall) # 136 genes, 17 skipped


## *MODERN* OT Pathway gene subset (PS4-12)
# subset OT pathway genes into modern genes only
modern_OT <- M1C_mm[M1C_mm$Phylostratum %in% modPS, ] # keep only modern phylostrata
modern_OT$GeneID <- toupper(modern_OT$GeneID)


# annotate ensembl genes IDs with gene name/symbol
ens_gene_mod <- ensembldb::select(EnsDb.Hsapiens.v79, 
                                  keys= modern_OT$GeneID, 
                                  keytype = "GENEID", 
                                  columns = c("SYMBOL","GENEID"))
ens_gene_mod <- ens_gene_mod %>% distinct(SYMBOL, .keep_all = TRUE) 
geneids_mod <- ens_gene_mod$SYMBOL  # reference list with gene names

# cortical
OTmod_c <- exp_m_anno_c[, names(exp_m_anno_c) %in% c("id", "region", "hemisphere", "structure", geneids_mod)] # extract subset of modern genes with cortical expression
dim(OTmod_c) # 19 modern genes, 68 cortical regions

# subcortical
OTmod_sc <- exp_m_anno_sc[, names(exp_m_anno_sc) %in% c("id", "region", "hemisphere", "structure", geneids_mod)] 
dim(OTmod_sc) # 19 modern genes, 15 subcortical regions





## *ANCIENT* OT Pathway gene subset (PS1-3)
ancient_OT <- M1C_mm[!M1C_mm$Phylostratum %in% modPS, ] # remove modern phylostrata
ancient_OT$GeneID <- toupper(ancient_OT$GeneID)

# annotate with gene name
ens_gene_old <- ensembldb::select(EnsDb.Hsapiens.v79, 
                                  keys= ancient_OT$GeneID, 
                                  keytype = "GENEID", 
                                  columns = c("SYMBOL","GENEID"))
ens_gene_old <- ens_gene_old %>% distinct(SYMBOL, .keep_all = TRUE) 
geneids_old <- ens_gene_old$SYMBOL  # reference list with gene names

# cortical
OTold_c <- exp_m_anno_c[, names(exp_m_anno_c) %in% c("id", "region", "hemisphere", "structure", geneids_old)] 
dim(OTold_c) # 107 ancient genes, 68 cortical regions

# subcortical
OTold_sc <-    exp_m_anno_sc[, names(exp_m_anno_sc) %in% c("id", "region", "hemisphere", "structure", geneids_old)] 
dim(OTold_sc) # 107 ancient genes, 15 subcortical regions











## ------------ ancient OT pathway genes


## cortical expression

list_g_oc <- colnames(OTold_c[5:111])
list_g_oc  

exp_oc <- function(gene) {
  oc <- ggseg(.data=OTold_c, mapping=aes_string(fill=gene), colour="white", size=1, position="stacked", hemisphere="left" ) +
    theme_void()
  
  oc <- oc + viridis::scale_fill_viridis(option = 'turbo') + 
    labs(title=paste(gene, "cortical expression"), subtitle="lateral | medial", fill="mRNA intensity") +
    theme(plot.title = element_text(hjust = .5, size=16),
          plot.subtitle = element_text(hjust = .5, size=14),
          legend.title = element_text(size=14),
          legend.text = element_text(size=12))
  
  value <- list(                                    
    oc = oc
  ) # Create a list of output objects
  attr(value, "class") <- "exp_oc"
  value
}


res_oc <- sapply(list_g_oc, exp_oc, USE.NAMES = T)
res_oc[["MYLK.oc"]]

## plotting 
#ogene_ccom <- plot_grid(plotlist = res_oc, ncol=4, scale = 0.9)
#ogene_ccom

# save plot
#ggsave(filename = paste0(BASE, "RProject/OT-bd-analyses_loc/output/ogene_ccom.pdf"), ogene_ccom,
#       width = 35, height = 165, units = "in", device='pdf')



## subcortical expression

list_g_osc <- colnames(OTold_sc[5:111])
list_g_osc  

exp_osc <- function(gene) {
  osc <- ggseg(.data=OTold_sc, atlas="aseg", mapping=aes_string(fill=gene), colour="white", size=1) +
    theme_void()
  
  osc <- osc + viridis::scale_fill_viridis(option = 'turbo', discrete = F) + 
    labs(title=paste(gene, "subcortical expression"), subtitle = "axial | sagittal", fill="mRNA intensity") 
  
  osc <- osc + theme(plot.title = element_text(hjust=0.5, size=16),
                     plot.subtitle = element_text(hjust=0.515, size=14),
                     legend.title = element_text(size=14),
                     legend.text = element_text(size=12)) 
  
  value <- list(                                    
    osc = osc
  ) # Create a list of output objects
  attr(value, "class") <- "exp_osc"
  value
}

res_osc <- sapply(list_g_osc, exp_osc, USE.NAMES = T)
res_osc[["MYLK.osc"]]



## plotting 
#ogene_sccom <- plot_grid(plotlist = res_osc, ncol=4, scale = 0.9)
#ogene_sccom


# save plot
#ggsave(filename = paste0(BASE, "RProject/OT-bd-analyses_loc/output/ogene_sccom.pdf"), ogene_sccom,
#       width = 35, height = 165, units = "in", device='pdf')






# function to concatenate left and right hemi region vals into one joint region value for each geneset average, so each region is unique and there is no L and R

test_fun <- function(colname) {

e_list <- list()
regions <- unique(temp81$region)

  for (i in 1:length(regions)) {
  mean <- mean(temp81[temp81$region==regions[i], colname])
  e_list[[i]] <- mean 
  }

vals_stack = do.call(rbind, e_list)
colnames(vals_stack) <- colname

vals_stack
#value <- list(
#  vals_stack = vals_stack) # Create a list of output objects
#attr(value, "class") <- "test_fun"
#value

}

names_col <- colnames(temp81[,5:105])

res <- lapply(names, test_fun)
exp_uniqueR <- do.call(cbind, res)
exp_uniqueR <- cbind(regions, exp_uniqueR)
  


# NEXT UP: compare the columns to make statements like: 5% of the random genesets had a higher exression in region x than the OT geneset, 80% had a higher
# expression in region y than OT geneset, and so on and so forth.


col_nam <- colnames(temp81_L[5:105])
  
temp81_L <- temp81[!temp81$hemisphere=="R",]

temp81_LL <- temp81_L %>% pivot_longer(all_of(col_nam),
                                        names_to = "gensets", values_to = "ex_value")


p81 <- ggplot(temp81_LL, aes(x=region, y=ex_value)) + 
  geom_point(shape=21, position=position_dodge2(width = .5)) +
  geom_point(data=temp81_LL[temp81_LL$gensets=="mean_OTy",], aes(x=region, y=ex_value), colour="red", size=4)
p81 <- p81 + theme(axis.text.x = element_text(angle=35, hjust=1))
p81 <- p81 + 
  geom_hline(yintercept = median(temp81_LL$ex_value), color="blue") +
  geom_hline(yintercept = mean(temp81_LL$ex_value), color="darkgreen") +
  geom_ribbon(aes(y = mean(ex_value), 
                  ymin = mean(ex_value) - sd(ex_value), 
                  ymax = mean(ex_value) + sd(ex_value)), 
              alpha = .2) 
p81


























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







regions_a2 <- unique(des_ad_info_l$name)

lista22 <- list()
for (i in 1:length(regions)) {
  
  mean_test <- desALL_infoL[desALL_infoL$name==regions_a2[i],"meanALL"]
  lista22[[i]] <- mean_test
}



desALL_infoL[desALL_infoL$name==,]


View(lista22)
unique(meanS_dp$mean_allstruc - dp_sd_mean)















###### install atlas

ggseg_atlas_repos("Desterieux")
des_repo <- ggseg_atlas_repos("desterieux", ignore.case = TRUE)
if (!requireNamespace("Desterieux", quietly = TRUE)) {
  install_ggseg_atlas(repo = des_repo$repo, source = des_repo$source)
}


library(ggsegDesterieux)

# inspect atlas
ggseg(atlas = ggsegDesterieux::desterieux, 
      mapping=aes(fill = region)) +
  scale_fill_brain("desterieux", package="ggsegDesterieux")

#####

# both hemispheres

reg_labels <- data.frame(brain_labels(desterieux))
reg_labels

dat <- desterieux$data
reg_regions <- unique(dat$region)
reg_regions


################

des_ad_reordered <- data.frame(des_ad_info_regnavn[match(reg_labels, des_ad_info_regnavn$region),])

des_ad_reordered$mean <- rowSums(des_ad_reordered[,3:8])/6

des_ad_reor_n <- des_ad_reordered[, c(2, 9)]





order <- unique(dat$label)


df3 <- data.frame(df2[match(order, df2$region),])


df3$region <- gsub('lh_', '', df3$region) 
df3$region <- gsub('rh_', '', df3$region)

df3$region <- gsub('_', ' ', df3$region) 
df3$mean <- as.numeric(df3$mean)

new.row <- data.frame(region="lh_Medial_wall", mean=NA, stringsAsFactors=F)
df2 <- des_ad_reor_n
df2[42,] <- c("lh_Medial_wall", "")
df2[117,] <- c("rh_Medial_wall", "")


# plotting


ggseg(.data = df2,
      atlas = desterieux)


ggplot(df3) +
  geom_brain(atlas = desterieux,
             fill = mean) #+
 # scale_fill_viridis_c(option = "cividis", direction = -1) +
 # theme_void() 


# adjust labeling of input data to order of template atlas labels




# left hemisphere

# data preparation
resdat_left <- results_dpL
resdat_left$Region <- gsub('_', ' ', resdat_left$Region)
resdat_left_c <- resdat_left[-c(75:90),]

resdat_left_ro <- data.frame(resdat_left_c[match(reg_regions, resdat_left_c$Region),])

names(resdat_left_ro)[names(resdat_left_ro) == "Region"] <- "region"

# plotting

p987 <- ggseg(.data = resdat_left_ro,
               atlas = desterieux,
               hemisphere = "left", 
               mapping = aes(fill = t_stat),
               colour = "black",
               size = .1, )
p987 <- p987 + scale_fill_viridis_c(option = "mako", direction = -1) +
  theme_void() +
  labs(fill='t statistic') 

p987






#des_all_donors_ALL$meanALL <- ml_dada           
#desALL <- des_all_donors_ALL[, c(1, 8)]


#desALL_info <- dplyr::right_join(des_info_s, desALL, by=c("index"="label"))                # annotate with brain region labels
#desALL_infoL <- subset(desALL_info, grepl("^L ", desALL_info$name))                        # keep only left hemisphere and remove right
#desALL_infoL$name <- gsub('L ', '', desALL_infoL$name)                                     # remove "L " in front of each region
#desALL_infoL <- data.frame(desALL_infoL)



#pauALL <- data.frame(ml_pada)     
#pauALL <- cbind(pau_all_donors_ALL[, 1], pauALL)
#pauALL_info <- dplyr::right_join(pau_info, pauALL, by=c("V1"="label"))                     # annotate with brain region labels
#pauALL_info <- data.frame(pauALL_info)


sd(desALL_infoL[desALL_infoL$name=="S_precentral-inf-part","meanALL"])


STR_pce <- STR_tai$pCE
STR_pce_dat <- STR_pce$data
STR_pce_dat


ppx <- ggplot(data=STR_pce_dat, aes(x=Stage, y=value)) +
  geom_point(aes(colour=PS), size=5, alpha=.6, position = position_jitterdodge(dodge.width = 1, jitter.width = 1, jitter.height = 50))
ppx

ppy <- ggplot(data=CBCdat, aes(y=CBC, x=Stage)) +
  geom_line()
ppy
CBC_tai









exp <- get_expression(structure_ids=c("Allen:10163"),
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
                     p.value = FALSE,
                     xlab          = "Developmental stage", 
                     ylab          = "TAI" )
ps

ps_dat <- ps$data                                               # Plotted

PlotContribution(mm,                                      # Plotted
                       legendName = "PS")

PlotRE(mm,                                                 # Weirdly and funnily, if the "Groups" parameter is NOT defined correctly, as is the case here, one
            Groups = list(c(1:7)),                              # gets access to the raw data. If, however, two groups are defined as its supposed to be, the
            legendName = "PS")                                  # raw data is not provided in the object "re".

 PlotMeans(mm,                                             # Plotted as average. Weirdly, even though the group has been limited to PS 1-7, it plots all PS?
                Groups = list(c(1:7)),          
                legendName = "PS",
                adjust.range = TRUE)

#pb <- PlotBarRE(ExpressionSet = mm,                             
#                Groups        = list(c(1:3), c(4:7, 10, 12)),           # PS 1-3 = "old", PS 4-12 = "young" (myTAI documentation p. 32/33)?
#                p.adjust.method = "fdr")

PlotCategoryExpr(mm,
                        legendName = "PS",
                        test.stat = TRUE,
                        type = "stage-centered",
                        distr.type = "violin",
                        log.expr = T)
View(mm)




M1C_pmm <- M1C_tai$pmm_l
DFC_pmm <- DFC_tai$pmm_l
VFC_pmm <- VFC_tai$pmm_l
OFC_pmm <- OFC_tai$pmm_l
S1C_pmm <- S1C_tai$pmm_l
IPC_pmm <- IPC_tai$pmm_l
A1C_pmm <- A1C_tai$pmm_l
STC_pmm <- STC_tai$pmm_l
ITC_pmm <- ITC_tai$pmm_l
V1C_pmm <- V1C_tai$pmm_l
MFC_pmm <- MFC_tai$pmm_l
HIP_pmm <- HIP_tai$pmm_l
STR_pmm <- STR_tai$pmm_l
AMY_pmm <- AMY_tai$pmm_l
MD_pmm <- MD_tai$pmm_l
CBC_pmm <- CBC_tai$pmm_l

colnames(M1C_pmm)[6] <- "M1C"
colnames(DFC_pmm)[6] <- "DFC"
colnames(VFC_pmm)[6] <- "VFC"
colnames(OFC_pmm)[6] <- "OFC"
colnames(S1C_pmm)[6] <- "S1C"
colnames(IPC_pmm)[6] <- "IPC"
colnames(A1C_pmm)[6] <- "A1C"
colnames(STC_pmm)[6] <- "STC"
colnames(ITC_pmm)[6] <- "ITC"
colnames(V1C_pmm)[6] <- "V1C"
colnames(MFC_pmm)[6] <- "MFC"
colnames(HIP_pmm)[6] <- "HIP"
colnames(STR_pmm)[6] <- "STR"
colnames(AMY_pmm)[6] <- "AMY"
colnames(MD_pmm)[6] <- "MD"
colnames(CBC_pmm)[6] <- "CBC"


pmm_all_struct <- M1C_pmm %>%
  left_join(DFC_pmm, by=c('Phylostratum', 'GeneID', 'index', 'Phylostratum class', 'Stages')) %>%
  left_join(VFC_pmm, by=c('Phylostratum', 'GeneID', 'index', 'Phylostratum class', 'Stages')) %>%
  left_join(OFC_pmm, by=c('Phylostratum', 'GeneID', 'index', 'Phylostratum class', 'Stages')) %>%
  left_join(S1C_pmm, by=c('Phylostratum', 'GeneID', 'index', 'Phylostratum class', 'Stages')) %>%
  left_join(IPC_pmm, by=c('Phylostratum', 'GeneID', 'index', 'Phylostratum class', 'Stages')) %>%
  left_join(A1C_pmm, by=c('Phylostratum', 'GeneID', 'index', 'Phylostratum class', 'Stages')) %>%
  left_join(STC_pmm, by=c('Phylostratum', 'GeneID', 'index', 'Phylostratum class', 'Stages')) %>%
  left_join(ITC_pmm, by=c('Phylostratum', 'GeneID', 'index', 'Phylostratum class', 'Stages')) %>%
  left_join(V1C_pmm, by=c('Phylostratum', 'GeneID', 'index', 'Phylostratum class', 'Stages')) %>%
  left_join(MFC_pmm, by=c('Phylostratum', 'GeneID', 'index', 'Phylostratum class', 'Stages')) %>%
  left_join(HIP_pmm, by=c('Phylostratum', 'GeneID', 'index', 'Phylostratum class', 'Stages')) %>%
  left_join(STR_pmm, by=c('Phylostratum', 'GeneID', 'index', 'Phylostratum class', 'Stages')) %>%
  left_join(AMY_pmm, by=c('Phylostratum', 'GeneID', 'index', 'Phylostratum class', 'Stages')) %>%
  left_join(MD_pmm, by=c('Phylostratum', 'GeneID', 'index', 'Phylostratum class', 'Stages')) %>%
  left_join(CBC_pmm, by=c('Phylostratum', 'GeneID', 'index', 'Phylostratum class', 'Stages')) 


pmm_mean_as = data.frame(rowMeans(pmm_all_struct[,6:21]))
pmm_means_as = cbind(pmm_all_struct[,1:5], pmm_mean_as)
colnames(pmm_means_as)[6] <- "pmm_means_as"

pmm_means_as$Phylostratum <- paste("PS", pmm_means_as$Phylostratum, sep="")   # turn column entries from "1" to "PS1", "5" to "PS5", etc.

pmm_means_as$Stages <- as.factor(pmm_means_as$Stages)
pmm_means_as$Phylostratum <- factor(pmm_means_as$Phylostratum, levels=c("PS1", "PS2", "PS3", "PS4", "PS5", "PS6","PS7", "PS10", "PS12"))
pmm_means_as$Stages <- factor(pmm_means_as$Stages , levels = c("Prenatal", 
                                                               "Infant",
                                                               "Child",
                                                               "Adolescent",
                                                               "Adult")) 


# for log transformed values
pmm_llog <- pmm_means_as
pmm_llog$pmm_means_asLOG <- log(pmm_llog$pmm_means_as)
pmm_llog$pmm_means_asLOG2 <- log2(pmm_llog$pmm_means_as)
pmm_llog$pmm_means_asLOG10 <- log10(pmm_llog$pmm_means_as)


# plotting no log
pmm_p <- ggplot(pmm_means_as, aes(x=Stages, y=pmm_means_as, fill=Phylostratum.class, colour=Phylostratum.class)) + theme_classic() +
  geom_jitter(mapping = aes(colour=Phylostratum.class), size=5, alpha=.5) +
  viridis::scale_colour_viridis(option = "mako", discrete=T, begin=.2, end=.8) #+ ylim(-15, 1) 
pmm_p <- pmm_p + xlab("\nOntogenetic stages") + ylab("Partial TAIs (log)") +
  theme(axis.title = element_text(size=25, face="bold"),
        axis.text = element_text(size=15),
        axis.text.x = element_text(angle=30, vjust=.75, hjust=.8))
pmm_p


# plotting LOG
pmm_pLOG <- ggplot(pmm_llog, aes(x=Stages, y=pmm_means_asLOG, fill=Phylostratum.class, colour=Phylostratum.class)) + theme_classic() +
  geom_jitter(mapping = aes(colour=Phylostratum.class), size=5, alpha=.5) +
  viridis::scale_colour_viridis(option = "mako", discrete=T, begin=.2, end=.8) #+ ylim(-15, 1) 
pmm_pLOG <- pmm_pLOG + xlab("\nOntogenetic stages") + ylab("Partial TAIs (log)") +
  theme(axis.title = element_text(size=25, face="bold"),
        axis.text = element_text(size=15),
        axis.text.x = element_text(angle=30, vjust=.75, hjust=.8))
pmm_pLOG

# plotting LOG 2
pmm_pLOG2 <- ggplot(pmm_llog, aes(x=Stages, y=pmm_means_asLOG2, fill=Phylostratum.class, colour=Phylostratum.class)) + theme_classic() +
  geom_jitter(mapping = aes(colour=Phylostratum.class), size=5, alpha=.5) +
  viridis::scale_colour_viridis(option = "mako", discrete=T, begin=.2, end=.8) #+ ylim(-15, 1) 
pmm_pLOG2 <- pmm_pLOG2 + xlab("\nOntogenetic stages") + ylab("Partial TAIs (log)") +
  theme(axis.title = element_text(size=25, face="bold"),
        axis.text = element_text(size=15),
        axis.text.x = element_text(angle=30, vjust=.75, hjust=.8))
pmm_pLOG2


# plotting LOG 10
pmm_pLOG10 <- ggplot(pmm_llog, aes(x=Stages, y=pmm_means_asLOG10, fill=Phylostratum.class, colour=Phylostratum.class)) + theme_classic() +
  geom_jitter(mapping = aes(colour=Phylostratum.class), size=5, alpha=.5) +
  viridis::scale_colour_viridis(option = "mako", discrete=T, begin=.2, end=.8) #+ ylim(-15, 1) 
pmm_pLOG10 <- pmm_pLOG10 + xlab("\nOntogenetic stages") + ylab("Partial TAIs (log)") +
  theme(axis.title = element_text(size=25, face="bold"),
        axis.text = element_text(size=15),
        axis.text.x = element_text(angle=30, vjust=.75, hjust=.8))
pmm_pLOG10






desD0<- read_csv(paste0(BASE, "/data/processed/abagen_DesPauli_rm/AHBAdesD0.csv"))

desD1 <- read_csv(paste0(BASE, "/data/processed/abagen_DesPauli_rm/AHBAdesD1.csv"))







value1 <- abs(rnorm(26))*2
data <- data.frame(
  x=LETTERS[1:26], 
  value1=value1, 
  value2=value1+1+rnorm(26, sd=1) 
)

View(data)





ggseg_atlas_repos()


###### fetch DK atlas map

ggseg_atlas_repos("DKT")
ggsegDKT <- ggseg_atlas_repos("DKT", ignore.case = TRUE)
if (!requireNamespace("DKT", quietly = TRUE)) {
  install_ggseg_atlas(repo = ggsegDKT$repo, source = ggsegDKT$source, force=T)
}


library(ggsegDKT)

# inspect atlas
ggseg(atlas = ggsegDKT::dkt, 
      mapping=aes(fill = region)) +
  scale_fill_brain("dkt", package="ggsegDKT")

#####



lpo <- ggsegDefaultExtra::hcpa_3d
dat <- lpo$ggseg_3d
View(dat[[1]])


someData %>%
  group_by(groups) %>%
  ggplot() +
  geom_brain(atlas = dk, 
             position = position_brain(hemi ~ side),
             aes(fill = p)) +
  facet_wrap(~groups)


dat <- dk$data





#### ALL genes
ens_gene_ALL <- ensembldb::select(EnsDb.Hsapiens.v79, 
                                  keys    = HomoSapiens.PhyloMap$Gene_ID, 
                                  keytype = "GENEID", 
                                  columns = c("SYMBOL","GENEID"))

ens_gene_ALL <- ens_gene_ALL %>% distinct(SYMBOL, .keep_all = TRUE)
mm_ALL       <- dplyr::inner_join(HomoSapiens.PhyloMap, ens_gene_ALL, by=c("Gene_ID"="GENEID"))


# calculate relative frequencies for ALL genes
tabALL <- table(mm_ALL$Phylostratum)/length(mm_ALL$Phylostratum)
o<- data.frame(table(mm_ALL$Phylostratum))
o <- o[-c(13:19),]



q<-table(mm_hyp$Phylostratum)

c1 <- c(94,19,5,4,1,7,4,0,0,1,0,3)    # absolute freqs
oqc <- data.frame(labs=labs, counts_abs_ALL=d1, counts_abs_OT=c1)
oqc$rel <- oqc$counts_abs_OT/oqc$counts_abs_ALL




NSS <- fread("C:/01Dokumente/07Erwerbstaetigkeit/01PhD_UiO/projects/first_author/oxytocinHealth/RQ4/RProject/OT-bd-analyses_loc/Data/ntSssZScorePMVar.bed")

colnames(NSS) <- c("Chr", "BPstart", "BPend", "NSSscore")








# get PS data
HomoSapiens.PhyloMap <- read_excel(paste0(BASE, "RProject/OT-bd-analyses_loc/MBE_2008_Homo_Sapiens_PhyloMap.xls"), 
                                   sheet = 1, skip = 1)
HomoSapiens.PhyloMap <- HomoSapiens.PhyloMap[ , c("Phylostratum","Gene_ID")]
# available genes: 22845

# get expression data from all ontogenetic stages for each brain region for as many genes as possible (ens gene ids in the phylo map as reference)
exp_all <- get_expression(structure_ids=c("Allen:10163"),
                          gene_ids=HomoSapiens.PhyloMap$Gene_ID, 
                          dataset='5_stages') 

ex_all_df <- as.data.frame(do.call(rbind, exp_all)) 

rownames(ex_all_df) <- c("Prenatal", "Infant", "Child", "Adolescent", "Adult")

ex_all_t <- as.data.frame(t(ex_all_df))
ex_all_t <- rownames_to_column(ex_all_t, var = "Gene_ID")
ex_all_t <- ex_all_t %>% distinct(Gene_ID, .keep_all = TRUE)                    # remove duplicates?
# available genes: 15858


# merge expression and phylogenetic data
mm_all <- MatchMap(HomoSapiens.PhyloMap, ex_all_t)
# available genes: 15858

# I end up with about 15858 after all matching and distincting steps. 

enrich_p_dat <- PlotEnrichment(ExpressionSet = mm_all, 
                               test.set = OT_subset, 
                               legendName = "PS",
                               measure = "log-foldchange",
                               p.adjust.method = "fdr",
                               use.only.map = F, 
                               complete.bg = F,
                               plot.bars = F)   

res_m <- enrich_p_dat$enrichment.matrix

res_pvals <- enrich_p_dat$p.values

res_pFdr <-  p.adjust(enrich_p_dat$p.values, method = "fdr", n = 19) 

# thresholds copied from PlotEnrichment() source code
M1C_res_m$star <- ""
M1C_res_m$star[M1C_res_m$pFdr <= .05]  <- "*"
M1C_res_m$star[M1C_res_m$pFdr <= .005]  <- "**"
M1C_res_m$star[M1C_res_m$pFdr <= 5e-04] <- "***"










PlotEnrichment()

age.table <- table(mm_all[, 1])
nPS <- length(age.table)
FactorBackgroundSet <- factor(mm_all[, 1], levels = names(age.table))


MatchedGeneIDs <- stats::na.omit(match(tolower(unlist(OT_subset)), 
                                       tolower(unlist(mm_all[, 2]))))
age.distr.test.set <- mm_all[MatchedGeneIDs, 1:2]


FactorTestSet <- factor(age.distr.test.set[, 1], levels = names(age.table))

N_ij <- cbind(table(FactorBackgroundSet), table(FactorTestSet))
get.contingency.tbl(N_ij, 
                    index)


enrichment.p_vals <- sapply(1:3, function(index) stats::fisher.test(cv_sim)$p.value, simplify = "array")


cv<-data.frame(V1=mm_absDF$counts_rel_ALL, V2=mm_absDF$counts_rel_OT)
cv_sim <- cv[1:3,]
stats::fisher.test(cv, simulate.p.value = T)









# calculate relative frequencies for ALL genes
rel_fALL <- data.frame(table(mm_all$Phylostratum)/length(mm_all$Phylostratum))
abs_fALL <- data.frame(table(mm_all$Phylostratum))
rel_fALL_r <- data.frame(round((rel_fALL$Freq*100), digits = 1))

d1 <- abs_fALL$Freq[1:12]               # absolute freqs
d2 <- rel_fALL$Freq[1:12]               # relative freqs







struc_ont <- read_csv(paste0(BASE, "/data/allenbrainatlas_dat/allenbrainatlas_structureontology_adultbrain.csv"))






########################################
########################################

tai_devALL <-
  function(region) 
  {
    exp <- get_expression(structure_ids=c(region),
                          gene_ids=exp_all_t_dIDs, 
                          dataset='5_stages') # Get expression values
    
    
    exp_df <- as.data.frame(do.call(rbind, exp)) 
    rownames(exp_df) <- NULL
    rownames(exp_df) <- c("Prenatal", "Infant", "Child", "Adolescent", "Adult")
    exp_t <- as.data.frame(t(exp_df))
    exp_t <- rownames_to_column(exp_t, var = "gene")                         
    
    ens_gene_dev <- ensembldb::select(EnsDb.Hsapiens.v79, 
                                      keys= exp_t$gene, 
                                      keytype = "GENEID", 
                                      columns = c("SYMBOL","GENEID"))
    
    ens_gene_dev <- ens_gene_dev %>% distinct(SYMBOL, .keep_all = TRUE)
    
    names(ens_gene_dev)[names(ens_gene_dev) == "SYMBOL"] <- "gene"
    
    ex_t_f <- dplyr::full_join(ens_gene_dev, exp_t, by = c("GENEID"="gene"))
    
    ex_t_f <- dplyr::select(ex_t_f, -gene)
    
    names(ex_t_f)[names(ex_t_f) == "GENEID"] <- "Gene_ID"
    
    HomoSapiens.PhyloMap <- read_excel(paste0(BASE, "RProject/OT-bd-analyses_loc/MBE_2008_Homo_Sapiens_PhyloMap.xls"), 
                                       sheet = 1, skip = 1)
    
    HomoSapiens.PhyloMap <- HomoSapiens.PhyloMap[ , c("Phylostratum","Gene_ID")]
    
    mm <- MatchMap(HomoSapiens.PhyloMap, ex_t_f)
    
    
    ps <- PlotSignature( ExpressionSet = mm,
                         measure       = "TAI", 
                         TestStatistic = "FlatLineTest",
                         p.value = FALSE,
                         xlab          = "Developmental stage", 
                         ylab          = "TAI" )
    
    ps_dat <- ps$data                                               # Plotted
    
    
    
    value <- list(
      mm = mm,
      ps = ps,
      ps_dat = ps_dat
    ) # Create a list of output objects
    attr(value, "class") <- "tai_devALL"
    value
  }



## get main data

get_name(10163)
M1C_taiALL <- tai_devALL("Allen:10163") 

get_name(10173)
DFC_taiALL <- tai_devALL("Allen:10173")

get_name(10185)
VFC_taiALL <- tai_devALL("Allen:10185")

get_name(10194)
OFC_taiALL <- tai_devALL("Allen:10194")

get_name(10209)
S1C_taiALL <- tai_devALL("Allen:10209")

get_name(10225)
IPC_taiALL <- tai_devALL("Allen:10225")

get_name(10236)
A1C_taiALL <- tai_devALL("Allen:10236")

get_name(10243)
STC_taiALL <- tai_devALL("Allen:10243")

get_name(10252)
ITC_taiALL <- tai_devALL("Allen:10252")

get_name(10269)
V1C_taiALL <- tai_devALL("Allen:10269")

get_name(10278)
MFC_taiALL <- tai_devALL("Allen:10278")

get_name(10294)
HIP_taiALL <- tai_devALL("Allen:10294")

get_name(10333)
STR_taiALL <- tai_devALL("Allen:10333")

get_name(10361)
AMY_taiALL <- tai_devALL("Allen:10361")

get_name(10398)
MD_taiALL <- tai_devALL("Allen:10398")

get_name(10657)
CBC_taiALL <- tai_devALL("Allen:10657")




## -------------- PLOTTING

## ================================= Combine individual TAI plots 

# Prepare variables

M1CdatALL <- M1C_taiALL$ps_dat
colnames(M1CdatALL)[2] <- "M1C"

DFCdatALL <- DFC_taiALL$ps_dat
colnames(DFCdatALL)[2] <- "DFC"

VFCdatALL <- VFC_taiALL$ps_dat
colnames(VFCdatALL)[2] <- "VFC"

OFCdatALL <- OFC_taiALL$ps_dat
colnames(OFCdatALL)[2] <- "OFC"

S1CdatALL <- S1C_taiALL$ps_dat
colnames(S1CdatALL)[2] <- "S1C"

IPCdatALL <- IPC_taiALL$ps_dat
colnames(IPCdatALL)[2] <- "IPC"

A1CdatALL <- A1C_taiALL$ps_dat
colnames(A1CdatALL)[2] <- "A1C"

STCdatALL <- STC_taiALL$ps_dat
colnames(STCdatALL)[2] <- "STC"

ITCdatALL <- ITC_taiALL$ps_dat
colnames(ITCdatALL)[2] <- "ITC"

V1CdatALL <- V1C_taiALL$ps_dat
colnames(V1CdatALL)[2] <- "V1C"

MFCdatALL <- MFC_taiALL$ps_dat
colnames(MFCdatALL)[2] <- "MFC"

HIPdatALL <- HIP_taiALL$ps_dat
colnames(HIPdatALL)[2] <- "HIP"

STRdatALL <- STR_taiALL$ps_dat
colnames(STRdatALL)[2] <- "STR"

AMYdatALL <- AMY_taiALL$ps_dat
colnames(AMYdatALL)[2] <- "AMY"

MDdatALL <- MD_taiALL$ps_dat
colnames(MDdatALL)[2] <- "MD"

CBCdatALL <- CBC_taiALL$ps_dat
colnames(CBCdatALL)[2] <- "CBC"

jtALL <- M1CdatALL %>%
  left_join(DFCdatALL, by='Stage') %>%
  left_join(VFCdatALL, by='Stage') %>%
  left_join(OFCdatALL, by='Stage') %>%
  left_join(S1CdatALL, by='Stage') %>%
  left_join(IPCdatALL, by='Stage') %>%
  left_join(A1CdatALL, by='Stage') %>%
  left_join(STCdatALL, by='Stage') %>%
  left_join(ITCdatALL, by='Stage') %>%
  left_join(V1CdatALL, by='Stage') %>%
  left_join(MFCdatALL, by='Stage') %>%
  left_join(HIPdatALL, by='Stage') %>%
  left_join(STRdatALL, by='Stage') %>%
  left_join(AMYdatALL, by='Stage') %>%
  left_join(MDdatALL, by='Stage') %>%
  left_join(CBCdatALL, by='Stage')

jtALL$Stage <- as.factor(jtALL$Stage)

level_order <- c("Prenatal", "Infant",
                 "Child", "Adolescent", 
                 "Adult")
jtlALL <- gather(jtALL, 
                 brain_region, 
                 TAI, 
                 M1C:CBC, 
                 factor_key=TRUE)

jtlALL$cat <- "All protein-coding genes"
jtlALL$Stage <- factor(jtlALL$Stage , levels = c("Prenatal", "Infant",
                                           "Child", "Adolescent", 
                                           "Adult"))

# single plot
p <- ggplot() +
  geom_line(jtlALL, mapping=aes(x= Stage, y=TAI, group=brain_region, color=cat), size = 1.5) +
  theme_bw()

p <- p + geom_line(jtl, mapping= aes(x= Stage, y=TAI, group=brain_region, color=cat), size= 1.5) +
  theme_bw() +
  scale_color_manual(values = c("grey50", "#FAA100", "#81C800", "#FF0564", "#104C7E"),
                     labels = c("All protein-coding genes", "Adolescence drop", "Childhood peak", "Combined", "Majority pattern"))  

p <- p + 
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(face="bold", size=25),
        axis.text.x  = element_text(angle=20, vjust=.5, hjust=0.5, size=21),
        axis.title.y = element_text(face="bold", size=24),
        axis.text.y = element_text(size=21),
        axis.line = element_line(color = 'black'))

p <- p + labs(color = "TAI trajectories") + xlab("\nOntogenetic stage") + ylab("Transcriptome Age Index\n")

p <- p + facet_zoom(ylim=c(1.4, 1.75)) 

p <- p +
  annotate('text', x = 0.75, y = 1.692, label = 'STR', colour= "#FF0564", size=7) +
  annotate('text', x = 0.75, y = 1.627, label = 'V1C', colour= "#81C800", size=7) +
  annotate('text', x = 0.75, y = 1.645, label = 'MD', colour= "#81C800", size=7) +
  annotate('text', x = 0.75, y = 1.72, label = 'AMY', colour= "#FAA100", size=7) +
  annotate('text', x = 0.75, y = 1.67, label = 'M1C', colour= "#FAA100", size=7) +
  annotate('text', x = 0.75, y = 1.636, label = 'IPC', colour= "#FAA100", size=7) 

  
p <- p + theme(legend.position = c(.45, .8),
               legend.title = element_text(size=20, face="bold"),
               legend.text = element_text(size=20))

pb <- ggplot_build(p)
pb$data[[3]][1, 'alpha'] <- 0
pb$data[[4]][1, 'alpha'] <- 0
pb$data[[5]][1, 'alpha'] <- 0
pb$data[[6]][1, 'alpha'] <- 0
pb$data[[7]][1, 'alpha'] <- 0
pb$data[[8]][1, 'alpha'] <- 0
pg <- ggplot_gtable(pb)
plot(pg)

ggsave(filename = paste0(BASE, "RProject/OT-bd-analyses_loc/output/TAIcomOT_ALL4.pdf"), pg,
       width = 15, height = 10, units = "in", device='pdf')


View(pb$data)



jtl$cat <- "Majority pattern"

jtl[jtl$brain_region %in% c('V1C', 'MD'),]$cat <- "Childhood peak"

jtl[jtl$brain_region %in% c('M1C', 'IPC', 'AMY'),]$cat <- "BAdolescence drop"

jtl[jtl$brain_region %in% c('STR'),]$cat <- "Combined"

jtl$Stage <- factor(jtl$Stage , levels = c("Prenatal", "Infant",
                                           "Child", "Adolescent", 
                                           "Adult"))








library(ggforce)






M1C_taiALL$ps
DFC_taiALL$ps
VFC_taiALL$ps
OFC_taiALL$ps

S1C_taiALL$ps
IPC_taiALL$ps
A1C_taiALL$ps
STC_taiALL$ps

ITC_taiALL$ps
V1C_taiALL$ps
MFC_taiALL$ps
HIP_taiALL$ps   ### weird TAI

STR_taiALL$ps   ### peak at adulthood
AMY_taiALL$ps   ### peak at adulthood
MD_taiALL$ps
CBC_taiALL$ps   ###  low TAIs






M1C_mm <- M1C_tai$mm
DFC_mm <- DFC_tai$mm
VFC_mm <- VFC_tai$mm
V1C_mm <- V1C_tai$mm
MFC_mm <- MFC_tai$mm


# ------------------------------------------------------ 210624



3217.093*1 + 1825.311*2 + 45.1848*3 + 60.78174*4 + 9.95168000*5 + 24.77607*6 + 84.18666*7 + 0.34655267*10+ 13.37*12


test <- "(t5:0.89,((t4:0.59,t1:0.37):0.34,(t2:0.03,t3:0.67):0.9):0.04);" 
test2 <- read.tree(text=test)




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

ggsave(filename = paste0(BASE, "RProject/OT-bd-analyses_loc/output/testb5.jpg"), testxb,
       width = 46, height = 20, units = "in", device='jpg')







treetext = "(((ADH2:0.11[&&NHX:S=], ADH1:0.11[&&NHX:S=Euteleostomi]):

0.05 [&&NHX:S=primates:D=Y:B=100],ADHY:

0.1[&&NHX:S=nematode],ADHX:0.12 [&&NHX:S=insect]):

0.1[&&NHX:S=metazoa:D=N],


(ADH4:0.09[&&NHX:S=yeast],

ADH3:0.13[&&NHX:S=yeast], ADH2:0.12[&&NHX:S=yeast],

ADH1:0.11[&&NHX:S=yeast]):0.1[&&NHX:S=Fungi])[&&NHX:D=CellularOrganisms];"

tree <- read.nhx(textConnection(treetext))
ggtree(tree) + geom_tiplab() + 
  geom_label(aes(x=branch, label=S), fill='lightgreen') + 
  geom_label(aes(label=D), fill='steelblue') + 
  geom_text(aes(label=B), hjust=-.5)





install.packages("treeio")
library(treeio)



# ------------------------------------------------------ 210628

otgene <- fread(paste0(BASE, "analyses/NSS_selective_sweep/oxytocin.bed"), h=T)

otreg <- fread(paste0(BASE, "analyses/NSS_selective_sweep/oxytocinreg.bed"), h=T)



# ------------------------------------------------------ 210630


# NSS score file can be downloaded at : https://genome.ucsc.edu/Neandertal/, "Selective Sweep Scan (S) on Neandertal vs. Human Polymorphisms", hg19

NTSS <- fread(paste0(BASE, "analyses/NSS_selective_sweep/data/ntSssZScorePMVar.bed"), h=T)
# 10.973.733 rows, 4 columns (rows are SNPs)


# retrieve gene list from KEGG
# NOTE: the genes belonging to a pathway are updated constantly, new genes are added on a regular basis, thus the
# list you retrieve might vary slightly from the list we initially retrieved

# first fetch information from KEGG to fetch label of the pathway we're looking for
kg.hsa=kegg.gsets("hsa") # retrieve gene set data from KEGG for Homo SApiens
kg.hsa.sets <- kg.hsa$kg.sets # fetch gene sets and search for "oxytocin signaling pathway"
# ----> pathway label: hsa04921

# get gene names
OTpthwy <- keggGet("hsa04921")[[1]]$GENE    
OTpthwy_nms_base <- test[seq(2,length(OTpthwy),2)]
OTpthwy_nms <- gsub("\\;.*","",namesodd)
OTpthwy_df <- data.frame(OTpthwy_nms)

OTensIDs <- ensembldb::select(EnsDb.Hsapiens.v79,                   # fetch ensemble gene IDs
                                  keys= OTpthwy_df$OTpthwy_nms, 
                                  keytype = "SYMBOL", 
                                  columns = c("SYMBOL","GENEID"))

OTensIDs <- OTensIDs %>% distinct(SYMBOL, .keep_all = TRUE)  




hg19 <- fread(paste0(BASE, "analyses/NSS_selective_sweep/hg19.fa/hg19.fa"), h=T)

DS <- read_csv(paste0(BASE, "/analyses/differential_stability/hawrylycz_2015/NIHMS731485-supplement-4.csv"))

DStop50 <- DS[1:8674,]
DSlower50 <- DS[8675:17348,]


DSOTtop50 <- dplyr::inner_join(OTpthwy_df, DStop50, by=c("OTpthwy_nms"="Gene"))
DSOTlower50 <- dplyr::inner_join(OTpthwy_df, DSlower50, by=c("OTpthwy_nms"="Gene"))

OTgenes.out <- dplyr::anti_join(OTpthwy_df, DS, by=c("OTpthwy_nms"="Gene"))
# no data available for 14 genes
# of 140 genes, 97 were among the top 50% genes with highest DS values, 43 were among the genes with lower DS values
# of the genes for which data was available, 69.29% have high DS values and are thus consistently and reliabley expressed across brain regions independet of donor characteristics. 



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

hyp_ex2 <- get_expression(structure_ids=left_hyp,      # fetch expression data
                         gene_ids=OTpthwy_df$OTpthwy_nms, 
                         dataset='adult')

hyp_ex2 <- as.data.frame(hyp_ex2)                       # prepare data: transpose, create column with gene labels

hyp_ex_t2 <- as.data.frame(t(hyp_ex2))

hyp_ex_t2 <- rownames_to_column(hyp_ex_t2, var = "gene")

hyp_ex_t2$gene

ens_gene_hyp2 <- ensembldb::select(EnsDb.Hsapiens.v79,                   # fetch ensemble gene names
                                  keys= hyp_ex_t2$gene, 
                                  keytype = "SYMBOL", 
                                  columns = c("SYMBOL","GENEID"))

ens_gene_hyp2 <- ens_gene_hyp2 %>% distinct(SYMBOL, .keep_all = TRUE)     # remove duplicates (?)
names(ens_gene_hyp2)[names(ens_gene_hyp2) == "SYMBOL"] <- "gene"
hyp_ex_t_f2 <- dplyr::full_join(ens_gene_hyp2, hyp_ex_t2, by = "gene")     # join with expression data
hyp_ex_t_f2 <- dplyr::select(hyp_ex_t_f2, -gene)

names(hyp_ex_t_f2)[names(hyp_ex_t_f2) == "GENEID"] <- "Gene_ID"

HomoSapiens.PhyloMap <- read_excel(paste0(BASE, "RProject/OT-bd-analyses_loc/MBE_2008_Homo_Sapiens_PhyloMap.xls"),  # load phylo map
                                   sheet = 1, skip = 1)

HomoSapiens.PhyloMap <- HomoSapiens.PhyloMap[ , c("Phylostratum","Gene_ID")]

counts <- table(HomoSapiens.PhyloMap$Phylostratum)

barplot(counts, main="Phylostratum levels",
        xlab="Phylostratum")

mm_hyp2 <- MatchMap(HomoSapiens.PhyloMap, hyp_ex_t_f2)



##################

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





download.file( url      = "http://www.biomedcentral.com/content/supplementary/1741-7007-8-66-s1.xls", 
               destfile = "BMCBiology_2010_Homo_Sapiens_PhyloMap.xls" )

HomoSapiensPhyloMap.BMCBiology <- read_excel(paste0(BASE, "RProject/OT-bd-analyses_loc/12915_2009_362_MOESM1_ESM.xls"), 
                                                    sheet = 1, skip = 3, col_names = FALSE)

colnames(HomoSapiensPhyloMap.BMCBiology)[1:2] <- c("Gene_ID","Phylostratum")

# format Phylostratigraphic Map for use with myTAI
HomoSapiens.PhyloMap2010 <- HomoSapiensPhyloMap.BMCBiology[ , c("Phylostratum","Gene_ID")]

# have a look at the final format
head(HomoSapiens.PhyloMap2010)
View(HomoSapiens.PhyloMap2010)


HomoSapiensPhyloMap.BMCGenomics <- read_excel(paste0(BASE, "RProject/OT-bd-analyses_loc/12864_2012_4867_MOESM1_ESM.xlsx"), sheet = 2)

HomoSapiens.PhyloMap2013 <- HomoSapiensPhyloMap.BMCGenomics[ , c(2,1)]

# have a look at the final format
View(HomoSapiens.PhyloMap2013)

mm_hyp2013 <- MatchMap(HomoSapiens.PhyloMap2013, hyp_ex_t_f)


# calculate relative frequencies for OT genes
rel_fOT <- table(mm_hyp2013$`Oldest Phylostratum`)/length(mm_hyp2013$`Oldest Phylostratum`)
rel_fOT_r <- round((rel_fOT*100), digits = 1)
abs_fOT <- table(mm_hyp2013$`Oldest Phylostratum`)











#### ALL genes
# get expression data from all ontogenetic stages for each brain region for as many genes as possible (ens gene ids in the phylo map as reference)
exp_allx <- get_expression(structure_ids=left_hyp,      # fetch expression data
                          gene_ids=HomoSapiens.PhyloMap$Gene_ID, # I'm using the gene ids for genes that are avail in the phylo map because these are the ones we need and will use in the end. No need to call genes for which there is no phylogenetic data avail anyways.
                          dataset='adult')

exp_allx <- as.data.frame(exp_allx)                       # prepare data: transpose, create column with gene labels

exp_all_tx <- as.data.frame(t(exp_allx))

exp_all_tx <- rownames_to_column(exp_all_tx, var = "gene")

exp_all_tx <- exp_all_tx %>% distinct(gene, .keep_all = T)

# merge expression and phylogenetic data
mmALLEx <- MatchMap(HomoSapiens.PhyloMap, exp_all_tx)
# final available genes: 14957




# calculate relative frequencies for ALL genes
rel_fALLx <- data.frame(table(mmALLEx$Phylostratum)/length(mmALLEx$Phylostratum))
rel_fALL_rx <- data.frame(round((rel_fALL$Freq*100), digits = 1))
abs_fALLx <- data.frame(table(mmALLEx$Phylostratum))


d1 <- abs_fALLx$Freq[1:12]               # absolute freqs
d2 <- rel_fALLx$Freq[1:12]               # relative freqs
d3 <- rel_fALL_rx$round..rel_fALL.Freq...100...digits...1.[1:12]  # relative freqs, rounded

dim(exp_allx)









exp <- get_expression(structure_ids=c("Allen:10163"),
                      gene_ids=OTpw.genes, 
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


HomoSapiensPhyloMap.BMCGenomics <- read_excel(paste0(BASE, "RProject/OT-bd-analyses_loc/12864_2012_4867_MOESM1_ESM.xlsx"), sheet = 2)

HomoSapiens.PhyloMap2013 <- HomoSapiensPhyloMap.BMCGenomics[ , c(2,1)]

names(HomoSapiens.PhyloMap2013)[names(HomoSapiens.PhyloMap2013) == "Oldest Phylostratum"] <- "Phylostratum"
names(HomoSapiens.PhyloMap2013)[names(HomoSapiens.PhyloMap2013) == "Ensembl Gene ID"] <- "EnsemblGeneID"


mm <- MatchMap(HomoSapiens.PhyloMap2013, ex_t_f)

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



pmm <- data.frame(pMatrix(mm))
pmm <- rownames_to_column(pmm)
pmm <- pmm %>% 
  mutate(index = 1:142)
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


expMD <- get_expression(structure_ids=c("Allen:10398"),
                      gene_ids=OTpw.genes, 
                      dataset='5_stages')          # Get expression values


expMD_df <- as.data.frame(do.call(rbind, expMD)) 
rownames(expMD_df) <- NULL
rownames(expMD_df) <- c("Prenatal", "Infant", "Child", "Adolescent", "Adult")
expMD_t <- as.data.frame(t(expMD_df))
expMD_t <- rownames_to_column(expMD_t, var = "gene")                         

ens_gene_devMD <- ensembldb::select(EnsDb.Hsapiens.v79, 
                                  keys= expMD_t$gene, 
                                  keytype = "SYMBOL", 
                                  columns = c("SYMBOL","GENEID"))

ens_gene_devMD <- ens_gene_devMD %>% distinct(SYMBOL, .keep_all = TRUE)

names(ens_gene_devMD)[names(ens_gene_devMD) == "SYMBOL"] <- "gene"

exMD_t_f <- dplyr::full_join(ens_gene_devMD, expMD_t, by = "gene")

exMD_t_f <- dplyr::select(exMD_t_f, -gene)

names(exMD_t_f)[names(exMD_t_f) == "GENEID"] <- "Gene_ID"     

mmMD5stages <- MatchMap(HomoSapiens.PhyloMap2013, exMD_t_f)







expefwee <- get_expression(structure_ids = c("Allen:10163"),
                           gene_ids      = HomoSapiens.PhyloMap2013$EnsemblGeneID, 
                           dataset       = '5_stages') # Get expression values

exp_allv <- get_expression(structure_ids = c("Allen:10389"),      # AHBA structure ontology ID for the diencephalon. The other IDs for the hypothalamus did not work. 
                           gene_ids      = HomoSapiens.PhyloMap2013$EnsemblGeneID, 
                           dataset       = '5_stages')


# ---------------------------------------------------- 210714
library(data.table)
library(readr)

blastoxt02 <- read_csv(paste0(BASE, "analyses/phylogeny_synteny/OXT_EXSHX8Y6013-Alignment-Descriptions.csv"), col_names = T)


OXTfasta <- fread(paste0(BASE, "analyses/phylogeny_synteny/OXT_GenSeq/ncbi_dataset/data/gene.fna"))


# ---------------------------------------------------- 210824

### .......... insulin signaling pathway: 137 genes, endocrine system

# get gene names
INSpthwy <- keggGet("hsa04910")[[1]]$GENE    
INSpthwy_nms_base <- INSpthwy[seq(2,length(INSpthwy),2)]
INSpthwy_nms <- gsub("\\;.*","",INSpthwy_nms_base)


INSexpDie <- get_expression(structure_ids=c("Allen:10389"),      # AHBA structure ontology ID for the diencephalon. The other IDs for the hypothalamus did not work. 
                            gene_ids=INSpthwy_nms, 
                            dataset='5_stages')   
# --> no exp. data for 10 genes

INSexpDie_df <- as.data.frame(do.call(rbind, INSexpDie)) 
rownames(INSexpDie_df) <- NULL
rownames(INSexpDie_df) <- c("Prenatal", "Infant", "Child", "Adolescent", "Adult")
INSexpDie_t <- as.data.frame(t(INSexpDie_df))
INSexpDie_t <- rownames_to_column(INSexpDie_t, var = "gene")                         

ens_gene_devInsDie <- ensembldb::select(EnsDb.Hsapiens.v79, 
                                        keys= INSexpDie_t$gene, 
                                        keytype = "SYMBOL", 
                                        columns = c("SYMBOL","GENEID"))

ens_gene_devInsDie <- ens_gene_devInsDie %>% distinct(SYMBOL, .keep_all = TRUE)

names(ens_gene_devInsDie)[names(ens_gene_devInsDie) == "SYMBOL"] <- "gene"

INSexpDie_t_f <- dplyr::full_join(ens_gene_devInsDie, INSexpDie_t, by = "gene")

INSexpDie_t_f <- dplyr::select(INSexpDie_t_f, -gene)

names(INSexpDie_t_f)[names(INSexpDie_t_f) == "GENEID"] <- "Gene_ID"     

mmInsDie5stages <- MatchMap(HomoSapiens.PhyloMap2013, INSexpDie_t_f) # 126 genes available, i.e., one gene (ENSG00000275176) discarded
 

# ---- get supp. table with genes and phylostrata
mmInsDie5stages_reduced <- mmInsDie5stages[,1:2] 
mmInsDie5stages_reduced$GeneID <- toupper(mmInsDie5stages_reduced$GeneID)

supp_matINS <- dplyr::left_join(mmInsDie5stages_reduced, ens_gene_devInsDie, by=c("GeneID"="GENEID"))
names(supp_matINS)[names(supp_matINS) == "gene"] <- "GeneLabel"
supp_matINS <- supp_matINS %>%       # order ascending by PS
  arrange((Phylostratum))


col1_INS <- NA
col2_INS <- NA
supp_matINS.2 <- dplyr::anti_join(data.frame(INSpthwy_nms), supp_matINS, by= c("INSpthwy_nms"="GeneLabel"))
supp_matINS.22 <- cbind(col1_INS, col2_INS, supp_matINS.2)
names(supp_matINS.22)[names(supp_matINS.22) == "col1_INS"] <- "Phylostratum"
names(supp_matINS.22)[names(supp_matINS.22) == "col2_INS"] <- "GeneID"
names(supp_matINS.22)[names(supp_matINS.22) == "INSpthwy_nms"] <- "GeneLabel"

supp_matINS_fin <- rbind(supp_matINS, supp_matINS.22)

library("writexl")
write_xlsx(supp_matINS_fin, paste0(BASE, "supp_matINSpthwy.xlsx"))






### .......... estrogen signaling pathway: 138 genes, endocrine system

# get gene names
ESpthwy <- keggGet("hsa04915")[[1]]$GENE    
ESpthwy_nms_base <- ESpthwy[seq(2,length(ESpthwy),2)]
ESpthwy_nms <- gsub("\\;.*","",ESpthwy_nms_base)


ESexpDie <- get_expression(structure_ids = c("Allen:10389"),      # AHBA structure ontology ID for the diencephalon. The other IDs for the hypothalamus did not work. 
                           gene_id       = ESpthwy_nms, 
                           dataset       = '5_stages')   
# --> no exp. data for 9 genes

ESexpDie_df <- as.data.frame(do.call(rbind, ESexpDie)) 
rownames(ESexpDie_df) <- NULL
rownames(ESexpDie_df) <- c("Prenatal", "Infant", "Child", "Adolescent", "Adult")
ESexpDie_t <- as.data.frame(t(ESexpDie_df))
ESexpDie_t <- rownames_to_column(ESexpDie_t, var = "gene")                         

ens_gene_devESDie <- ensembldb::select(EnsDb.Hsapiens.v79, 
                                       keys= ESexpDie_t$gene, 
                                       keytype = "SYMBOL", 
                                       columns = c("SYMBOL","GENEID"))

ens_gene_devESDie <- ens_gene_devESDie %>% distinct(SYMBOL, .keep_all = TRUE)

names(ens_gene_devESDie)[names(ens_gene_devESDie) == "SYMBOL"] <- "gene"

ESexpDie_t_f <- dplyr::full_join(ens_gene_devESDie, ESexpDie_t, by = "gene")

ESexpDie_t_f <- dplyr::select(ESexpDie_t_f, -gene)

names(ESexpDie_t_f)[names(ESexpDie_t_f) == "GENEID"] <- "Gene_ID"     

mmESDie5stages <- MatchMap(HomoSapiens.PhyloMap2013, ESexpDie_t_f) # 129 genes available, i.e., no genes discarded


# ---- get supp. table with genes and phylostrata
mmESDie5stages_reduced <- mmESDie5stages[,1:2] 
mmESDie5stages_reduced$GeneID <- toupper(mmESDie5stages_reduced$GeneID)

supp_matES <- dplyr::left_join(mmESDie5stages_reduced, ens_gene_devESDie, by=c("GeneID"="GENEID"))
names(supp_matES)[names(supp_matES) == "gene"] <- "GeneLabel"
supp_matES <- supp_matES %>%       # order ascending by PS
  arrange((Phylostratum))


col1_ES <- NA
col2_ES <- NA
supp_matES.2 <- dplyr::anti_join(data.frame(ESpthwy_nms), supp_matES, by= c("ESpthwy_nms"="GeneLabel"))
supp_matES.22 <- cbind(col1_ES, col2_ES, supp_matES.2)
names(supp_matES.22)[names(supp_matES.22) == "col1_ES"] <- "Phylostratum"
names(supp_matES.22)[names(supp_matES.22) == "col2_ES"] <- "GeneID"
names(supp_matES.22)[names(supp_matES.22) == "ESpthwy_nms"] <- "GeneLabel"

supp_matES_fin <- rbind(supp_matES, supp_matES.22)

library("writexl")
write_xlsx(supp_matES_fin, paste0(BASE, "supp_matESpthwy.xlsx"))





### .......... relaxin signaling pathway: 129 genes, endocrine system

# get gene names
RXpthwy <- keggGet("hsa04926")[[1]]$GENE    
RXpthwy_nms_base <- RXpthwy[seq(2,length(RXpthwy),2)]
RXpthwy_nms <- gsub("\\;.*","",RXpthwy_nms_base)


RXexpDie <- get_expression(structure_ids=c("Allen:10389"),      # AHBA structure ontology ID for the diencephalon. The other IDs for the hypothalamus did not work. 
                           gene_ids=RXpthwy_nms, 
                           dataset='5_stages')   
# --> no exp. data for 12 genes

RXexpDie_df <- as.data.frame(do.call(rbind, RXexpDie)) 
rownames(RXexpDie_df) <- NULL
rownames(RXexpDie_df) <- c("Prenatal", "Infant", "Child", "Adolescent", "Adult")
RXexpDie_df_t <- as.data.frame(t(RXexpDie_df))
RXexpDie_df_t <- rownames_to_column(RXexpDie_df_t, var = "gene")                         

ens_gene_devRXDie <- ensembldb::select(EnsDb.Hsapiens.v79, 
                                       keys= RXexpDie_df_t$gene, 
                                       keytype = "SYMBOL", 
                                       columns = c("SYMBOL","GENEID"))

ens_gene_devRXDie <- ens_gene_devRXDie %>% distinct(SYMBOL, .keep_all = TRUE)

names(ens_gene_devRXDie)[names(ens_gene_devRXDie) == "SYMBOL"] <- "gene"

RXexpDie_t_f <- dplyr::full_join(ens_gene_devRXDie, RXexpDie_df_t, by = "gene")

RXexpDie_t_f <- dplyr::select(RXexpDie_t_f, -gene)

names(RXexpDie_t_f)[names(RXexpDie_t_f) == "GENEID"] <- "Gene_ID"     

mmRXDie5stages <- MatchMap(HomoSapiens.PhyloMap2013, RXexpDie_t_f) # 117 genes available, no genes discarded


# ---- get supp. table with genes and phylostrata
mmRXDie5stages_reduced <- mmRXDie5stages[,1:2] 
mmRXDie5stages_reduced$GeneID <- toupper(mmRXDie5stages_reduced$GeneID)

supp_matRX <- dplyr::left_join(mmRXDie5stages_reduced, ens_gene_devRXDie, by=c("GeneID"="GENEID"))
names(supp_matRX)[names(supp_matRX) == "gene"] <- "GeneLabel"
supp_matRX <- supp_matRX %>%       # order ascending by PS
  arrange((Phylostratum))


col1_RX <- NA
col2_RX <- NA
supp_matRX.2 <- dplyr::anti_join(data.frame(RXpthwy_nms), supp_matRX, by= c("RXpthwy_nms"="GeneLabel"))
supp_matRX.22 <- cbind(col1_RX, col2_RX, supp_matRX.2)
names(supp_matRX.22)[names(supp_matRX.22) == "col1_RX"] <- "Phylostratum"
names(supp_matRX.22)[names(supp_matRX.22) == "col2_RX"] <- "GeneID"
names(supp_matRX.22)[names(supp_matRX.22) == "RXpthwy_nms"] <- "GeneLabel"

supp_matRX_fin <- rbind(supp_matRX, supp_matRX.22)

#library("writexl")
write_xlsx(supp_matRX_fin, paste0(BASE, "supp_matRXpthwy.xlsx"))







# ------------------------------------------------------ 210901


# get list of all protein-coding genes

library(biomaRt)  
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = 'www.ensembl.org')
genes <- biomaRt::getBM(attributes = c("external_gene_name", "chromosome_name", "transcript_biotype"), 
                        filters = c("transcript_biotype", "chromosome_name"), values = list("protein_coding", c(1:22)), mart = mart)
genenames <- data.frame(genes$external_gene_name)

genes1to1000 <- data.frame(genenames[1:1000,])
genes1001to2000 <- data.frame(genenames[1001:2000,])
genes2001to3000 <- data.frame(genenames[2001:3000,])
genes3001to4000 <- data.frame(genenames[3001:4000,])
genes4001to5000 <- data.frame(genenames[4001:5000,])

genes5001to6000 <- data.frame(genenames[5001:6000,])
genes6001to7000 <- data.frame(genenames[6001:7000,])
genes7001to8000 <- data.frame(genenames[7001:8000,])
genes8001to9000 <- data.frame(genenames[8001:9000,])
genes9001to10000 <- data.frame(genenames[9001:10000,])

genes10001to11000 <- data.frame(genenames[10001:11000,])
genes11001to12000 <- data.frame(genenames[11001:12000,])
genes12001to13000 <- data.frame(genenames[12001:13000,])
genes13001to14000 <- data.frame(genenames[13001:14000,])
genes14001to15000 <- data.frame(genenames[14001:15000,])

genes15001to16000 <- data.frame(genenames[15001:16000,])
genes16001to17000 <- data.frame(genenames[16001:17000,])
genes17001to18000 <- data.frame(genenames[17001:18000,])
genes18001to18445 <- data.frame(genenames[18001:18445,])




genes1to1000_t <- t(genes1to1000)
genes1001to2000_t <- t(genes1001to2000)
genes2001to3000_t <- t(genes2001to3000)
genes3001to4000_t <- t(genes3001to4000)
genes4001to5000_t <- t(genes4001to5000)

genes5001to6000_t <- t(genes5001to6000)
genes6001to7000_t <- t(genes6001to7000)
genes7001to8000_t <- t(genes7001to8000)
genes8001to9000_t <- t(genes8001to9000)
genes9001to10000_t <- t(genes9001to10000)

genes10001to11000_t <- t(genes10001to11000)
genes11001to12000_t <- t(genes11001to12000)
genes12001to13000_t <- t(genes12001to13000)
genes13001to14000_t <- t(genes13001to14000)
genes14001to15000_t <- t(genes14001to15000)

genes15001to16000_t <- t(genes15001to16000)
genes16001to17000_t <- t(genes16001to17000)
genes17001to18000_t <- t(genes17001to18000)
genes18001to18445_t <- t(genes18001to18445)




write.csv(genes1to1000_t, paste0(BASE, "genes1to1000.csv"), row.names = F, quote=F)
write.csv(genes1001to2000_t, paste0(BASE, "genes1001to2000.csv"), row.names = F, quote=F)
write.csv(genes2001to3000_t, paste0(BASE, "genes2001to3000.csv"), row.names = F, quote=F)
write.csv(genes3001to4000_t, paste0(BASE, "genes3001to4000.csv"), row.names = F, quote=F)
write.csv(genes4001to5000_t, paste0(BASE, "genes4001to5000.csv"), row.names = F, quote=F)

write.csv(genes5001to6000_t, paste0(BASE, "genes5001to6000.csv"), row.names = F, quote=F)
write.csv(genes6001to7000_t, paste0(BASE, "genes6001to7000.csv"), row.names = F, quote=F)
write.csv(genes7001to8000_t, paste0(BASE, "genes7001to8000.csv"), row.names = F, quote=F)
write.csv(genes8001to9000_t, paste0(BASE, "genes8001to9000.csv"), row.names = F, quote=F)
write.csv(genes9001to10000_t, paste0(BASE, "genes9001to10000.csv"), row.names = F, quote=F)

write.csv(genes10001to11000_t, paste0(BASE, "genes10001to11000.csv"), row.names = F, quote=F)
write.csv(genes11001to12000_t, paste0(BASE, "genes11001to12000.csv"), row.names = F, quote=F)
write.csv(genes12001to13000_t, paste0(BASE, "genes12001to13000.csv"), row.names = F, quote=F)
write.csv(genes13001to14000_t, paste0(BASE, "genes13001to14000.csv"), row.names = F, quote=F)
write.csv(genes14001to15000_t, paste0(BASE, "genes14001to15000.csv"), row.names = F, quote=F)

write.csv(genes15001to16000_t, paste0(BASE, "genes15001to16000.csv"), row.names = F, quote=F)
write.csv(genes16001to17000_t, paste0(BASE, "genes16001to17000.csv"), row.names = F, quote=F)
write.csv(genes17001to18000_t, paste0(BASE, "genes17001to18000.csv"), row.names = F, quote=F)
write.csv(genes18001to18445_t, paste0(BASE, "genes18001to18445.csv"), row.names = F, quote=F)



# ------------------------------------------------------ 210902

test <- fread(paste0(BASE, "RProject/OT-bd-analyses_loc/Data/dNdS/genevo_2021-09-01_1to1000.tsv"))

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
empty_l5 <- list()
newdf <- NULL

for (i in df_list1) {
  df2 <- data.frame(i[, c(2, 7, 10, 13, 16, 19, 22, 25, 28)])
  #empty_l5[[length(empty_l5)+1]] <- df2                                  # if you want to create a list and from that list
  #assign(paste0("set", length(empty_l5)), df2) # saves to environment    # 'export' data frames into the environment
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

## functions for whisker length adjustment

o <- function(x) {
  subset(x, x == max(x) | x == min(x))
}

f <- function(x) {
  r <- quantile(x, probs = c(0.00, 0.25, 0.5, 0.75, 1))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

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
  viridis::scale_fill_viridis(option = 'magma', discrete = T, begin = .2, end = .85) +
  viridis::scale_colour_viridis(option = 'magma', discrete = T, begin = .2, end = .85)
pALL_mq <- pALL_mq + 
  geom_vline(xintercept = 1, linetype = "dashed", size = .75, color="gray50")

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
ggsave(filename = paste0(BASE, "RProject/OT-bd-analyses_loc/output/dndsALL_tree_mq2.pdf"), dnds_ALL_tree_mq,
       width = 22, height = 8, units = "in", device='pdf')




# ---------------------------------------- 211006

test <- g_means_lt_sc_L_asc[g_means_lt_sc_L_asc$mean1 < 0.4515017, g_means_lt_sc_L_asc$expression]

g_means_lt_sc_L_asc[2, 3]

low_exp <- g_means_lt_sc_L_asc[g_means_lt_sc_L_asc$expression < 0.4515017, c(8:9)]

hi_exp <- g_means_lt_sc_L_asc[g_means_lt_sc_L_asc$expression > 0.4515017, c(8:9)]





(1*-1) + (1*3) - (2*1) 
21

# ---------------------------------------- 211020

# Enable this universe
options(repos = c(
  ggseg = 'https://ggseg.r-universe.dev',
  CRAN = 'https://cloud.r-project.org'))

# Install some packages
install.packages('ggsegDefaultExtra')

# install.packages("remotes")
remotes::install_github("LCBC-UiO/ggsegDefaultExtra")





p127 <-  ggplot() + theme_void() +
  geom_brain(data = DKresdat_l, 
             atlas = dk, 
             colour = "white",
             position = position_brain(side + hemi ~ .),
             aes(fill = t_stat)) +
  viridis::scale_fill_viridis(option="mako", discrete=F) +
  labs(fill='t statistic') 

p127

p127 <- p127 + ylab("\nlateral   |   medial\n") + labs(title="Cortical expression") +
  theme(plot.title = element_text(hjust=.5, size=25, face="bold"),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_text(size=23, angle=90, face="bold"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.title = element_text(size=23),
        legend.text = element_text(size=16),
        legend.key.size = unit(1, 'cm'))

p127


plot(dkextra) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 7)) +
  guides(fill = guide_legend(ncol = 3))


View(dkextra$data)

ggplot() +
  geom_brain(atlas = dk, hemi = "left", side = "medial")

install.packages("ggseg")


plot(dkextra) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 7)) +
  guides(fill = guide_legend(ncol = 3))




# ---------------------------------------- 211021

# Enable this universe
options(repos = c(
  ggseg = 'https://ggseg.r-universe.dev',
  CRAN = 'https://cloud.r-project.org'))

# Install some packages
install.packages('ggsegDefaultExtra')

library(ggsegDefaultExtra)
library(ggseg)
library(ggplot2)

plot(dkextra) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 7)) +
  guides(fill = guide_legend(ncol = 3))


# ---------------------------------------- 211026

jcr19 <- read.csv("D:/eigene_dateien/02Dokumente/07Erwerbstaetigkeit/01PhD_UiO/courses/01UiO_phd_programme/HS2021-22/SV9107_international-publishing/exam/JCR2019.csv",
                  header = F)

names(jcr19) <- jcr19[2,]

jcr19 <-jcr19[-c(1:2),]

jcr19_2 = jcr19

jcr19$`Total Cites` <- str_replace_all(jcr19$`Total Cites`, ",", "")
jcr19_2$`Total Cites` <- str_replace_all(jcr19_2$`Total Cites`, ",", ".")

jcr19$`Total Cites` <- as.numeric(jcr19$`Total Cites`)
jcr19_2$`Total Cites` <- as.numeric(jcr19_2$`Total Cites`)


sum(jcr19$`Total Cites`, na.rm=T)
sum(jcr19_2$`Total Cites`, na.rm=T)








pee <- ggseg(DKresdat_l, 
        atlas = dkextra, 
        colour = "white",
        position = "dispersed",
        hemisphere = "left",
        mapping = aes(fill = t_stat)) +
  viridis::scale_fill_viridis(option = "mako", discrete = F) +
  labs(fill='t statistic') 

pee


pee2 <- pee + ylab("\ninferior | superior                             \n") + 
  labs(title="Cortical expression", subtitle = "\nlateral   |   medial\n") +
  theme(plot.title = element_text(hjust=.5, size=25, face="bold"),
        plot.subtitle = element_text(hjust=.5, size=18, face="bold"),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title.y = element_text(size=18, angle=90, face="bold"),
        axis.title.x = element_blank(),
        legend.title = element_text(size=23),
        legend.text = element_text(size=16),
        legend.key.size = unit(1, 'cm'),
        text = element_text(family="Arial"))

pee2





library(extrafont)
loadfonts(device = "win")


p127 <-  ggplot() + theme_void() +
  geom_brain(data = DKresdat_l, 
             atlas = dk, 
             colour = "white",
             position = position_brain(side + hemi ~ .),
             aes(fill = t_stat)) +
  viridis::scale_fill_viridis(option="mako", discrete=F) +
  labs(fill='t statistic') 




someData <- data.frame(
  region = c("transverse temporal", "insula",
             "precentral","superior parietal"),
  p = sample(seq(0,.5,.001), 4),
  stringsAsFactors = FALSE
)

ggseg(someData, atlas = dk, mapping = aes(fill = p))

DKresdat_l_new <- DKresdat_l[-1,]
DKresdat_l_new2 <- DKresdat_l_new[, c(1,4)]

ggseg(someData, atlas = dk, mapping = aes(fill = p))



# ---------------------------------------- 220506

temp2 = c("001OXT-c-carcharias","002OXTR-p-marinus", "003GNAQ-c-carcharias", "004HRAS-p-marinus", "005KRAS-c-carcharias", 
          "006NRAS-c-carcharias", "007RAF1-c-carcharias", "008MAP2K1-c-carcharias", "009MAP2K2-c-carcharias",
          "010MAPK1-c-carcharias", "013PLA2G4A-p-marinus", "018PLA2G4F-c-carcharias", "019PTGS2-p-marinus", "020MAP2K5-p-marinus",
          "022JUN-c-carcharias", "023FOS-c-carcharias", "024MEF2C-c-carcharias", "025CCND1-c-carcharias", "026ELK1-c-carcharias", 
          "028RYR2-c-carcharias", "029RYR3-c-carcharias", "030CD38-c-carcharias", "032KCNJ2-c-carcharias", "033KCNJ12-c-carcharias",
          "035KCNJ4-c-carcharias", "037PLCB1-c-carcharias", "040PLCB4-c-carcharias", "041PRKCA-c-carcharias", "042PRKCB-c-carcharias", 
          "044EEF2K-c-carcharias", "045EEF2-c-carcharias", "046CACNA1C-c-carcharias", "047CACNA1D-c-carcharias", "050CACNB1-p-marinus",
          "051CACNB2-c-carcharias", "053CACNB4-c-carcharias", "054CACNA2D1-c-carcharias", "056CACNA2D3-p-marinus", "057CACNA2D4-p-marinus",
          "058CACNG1-c-carcharias", "059CACNG2-c-carcharias", "060CACNG3-c-carcharias", "061CACNG4-c-carcharias", "062CACNG5-c-carcharias", 
          "066ITPR1-c-carcharias", "067ITPR2-c-carcharias", "070CALM2-c-carcharias", "073CALML6-c-carcharias", "075CALML4-c-carcharias",
          "076PPP3CA-c-carcharias", "077PPP3CB-p-marinus", "078PPP3CC-c-carcharias", "079PPP3R1-c-carcharias", "081NFATC1-c-carcharias",
          "082NFATC2-c-carcharias", "083NFATC3-c-carcharias", "087CAMK2-c-carcharias", "088PRKAA1-c-carcharias", "089PRKAA2-c-carcharias", 
          "090PRKAB1-c-carcharias", "092PRKAG1-c-carcharias", "094PRKAG2-c-carcharias", "095CAMK1D-c-carcharias", "096CAMK1G-c-carcharias", 
          "097CAMK1-c-carcharias", "098CAMK2A-c-carcharias", "099CAMK2D-c-carcharias", "101CAMK2G-c-carcharias",
          "102CAMK4-c-carcharias", "104GUCY1A2-p-marinus", "106GUCY1B1-c-carcharias", "108NPR1-c-carcharias", "112MYLK3-c-carcharias", 
          "116MYL9-c-carcharias", "119RHOA-c-carcharias", "120ROCK1-c-carcharias", "121ROCK2-c-carcharias", "123PPP1CB-c-carcharias", 
          "124PPP1CC-p-marinus", "125PPP1R12A-p-marinus", "129ADCY1-c-carcharias", "130ADCY2-c-carcharias", "133ADCY5-c-carcharias",
          "135ADCY7-c-carcharias", "136ADCY8-c-carcharias", "137ADCY9-p-marinus", "138PRKACA-c-carcharias", "139PRKACB-c-carcharias", 
          "141GNAI1-c-carcharias", "146PIK3R5-p-marinus", "148SRC-c-carcharias", "149KCNJ3-c-carcharias", "150KCNJ6-c-carcharias", 
          "152KCNJ5-c-carcharias", "153EGFR-c-carcharias")

list.1 <- list()
for (i in temp2) {
  list.1[[i]] = read_csv(paste0(path, "data/raw/BLASTp/all/", i, ".csv")) # read in dfs
  colnames(list.1) <- temp3
  names(list.1[[i]]) <- make.names(names(list.1[[i]]), unique = TRUE)
  
  # df = as.data.frame(list.1)
  #assign(i, df) #make new df
}




# ---------------------------------------- 220509

path = "D:/eigene_dateien/02Dokumente/07Erwerbstaetigkeit/01PhD_UiO/projects/first_author/OTpathway_evolution/RProject/OT-bd-analyses_loc/"
temp <- list.files(path = paste0(path, "data/raw/BLASTp/all/"), full.names = T)
temp3 <- list.files(path = paste0(path, "data/raw/BLASTp/all/"), full.names = F)


L1 <- list()
for (i in 1:length(temp)) {
 L1[[i]] = data.frame(read_csv(temp[i])) # read in csvs into list
 names(L1)[i] <- paste(temp3[i])
 names(L1[[i]]) <- make.names(names(L1[[i]]), unique = TRUE)
 #L1 <- tools::file_path_sans_ext(L1) 
}

# list2env(L1 ,.GlobalEnv)   # unlists L1 list with dataframes and 

pat = data.frame(c("Ciona savignyi", "Ciona intestinalis", "Branchiostoma lanceolatum", "Branchiostoma floridae", "Strongylocentrotus purpuratus", 
        "Asterias rubens", "Saccoglossus kowalevskii", "Xenoturbella bocki", "Symsagittifera roscoffensis", "Caenorhabditis elegans", 
        "Priapulus caudatus", "Octopus vulgaris", "Hydra vulgaris", "Acropora millepora", "Clytia hemisphaerica", "Ephydatia muelleri", 
        "Amphimedon queenslandica", "Salpingoeca rosetta", "Monosiga brevicollis", "Fusarium culmorum", "Aspergillus nidulans", 
        "Saccharomyces cerevisiae", "Dictyostelium discoideum", "Spironucleus salmonicida", "Methanolobus zinderi", "Escherichia coli",
        "Caulobacter crescentus", "Caulobacter vibrioides"))
pat$c..Ciona.savignyi....Ciona.intestinalis....Branchiostoma.lanceolatum... <- as.character(pat$c..Ciona.savignyi....Ciona.intestinalis....Branchiostoma.lanceolatum...)


# loops through ONE ref gene csv file and returns list with max scores, accession lengths, 
# protein length of ref gene for each species, and results


### ----------------------------------

# -------------------------------------------- template
df.temp <- data.frame(1:nrow(pat))
L.temp <- list()

for (k in 1:nrow(pat)) {
  L.temp[[k]] <- data.frame(L1[[0000]][grep(pat[k,], L1[[0000]][["Scientific.Name"]]), ])
  df.temp[k,] <- max((data.frame(L1[[0000]][grep(pat[k,], L1[[0000]][["Scientific.Name"]]), ]))$Max.Score)
  df.temp[k, 2] <- first(L.temp[[k]][["Acc..Len"]][L.temp[[k]][["Max.Score"]] == max(L.temp[[k]][["Max.Score"]])])
  df.temp[, 3] <- rep(0000, 28)
  df.temp[, 4] <- df.temp[, 1] / ((df.temp[, 3]) + df.temp[, 2])
}

### ----------------------------------



# OXT
df.001oxt <- data.frame(1:nrow(pat))
L.001oxt  <- list()

for (k in 1:nrow(pat)) {
  L.001oxt[[k]] <- data.frame(L1[[1]][grep(pat[k,], L1[[1]][["Scientific.Name"]]), ])
  df.001oxt[k,] <- max((data.frame(L1[[1]][grep(pat[k,], L1[[1]][["Scientific.Name"]]), ]))$Max.Score)
  df.001oxt[k, 2] <- first(L.001oxt[[k]][["Acc..Len"]][L.001oxt[[k]][["Max.Score"]] == max(L.001oxt[[k]][["Max.Score"]])])
  df.001oxt[, 3] <- rep(126, 28)
  df.001oxt[, 4] <- df.001oxt[, 1] / ((df.001oxt[, 3]) + df.001oxt[, 2])
}


# OXTR
df.002oxtr <- data.frame(1:nrow(pat))
L.002oxtr  <- list()

for (k in 1:nrow(pat)) {
  L.002oxtr[[k]] <- data.frame(L1[[2]][grep(pat[k,], L1[[2]][["Scientific.Name"]]), ])
  df.002oxtr[k,] <- max((data.frame(L1[[2]][grep(pat[k,], L1[[2]][["Scientific.Name"]]), ]))$Max.Score)
  df.002oxtr[k, 2] <- first(L.002oxtr[[k]][["Acc..Len"]][L.002oxtr[[k]][["Max.Score"]] == max(L.002oxtr[[k]][["Max.Score"]])])
  df.002oxtr[, 3] <- rep(431, 28)
  df.002oxtr[, 4] <- df.002oxtr[, 1] / ((df.002oxtr[, 3]) + df.002oxtr[, 2])
}


# GNAQ
df.003gnaq <- data.frame(1:nrow(pat))
L.003gnaq  <- list()

for (k in 1:nrow(pat)) {
  L.003gnaq[[k]] <- data.frame(L1[[3]][grep(pat[k,], L1[[3]][["Scientific.Name"]]), ])
  df.003gnaq[k,] <- max((data.frame(L1[[3]][grep(pat[k,], L1[[3]][["Scientific.Name"]]), ]))$Max.Score)
  df.003gnaq[k, 2] <- first(L.003gnaq[[k]][["Acc..Len"]][L.003gnaq[[k]][["Max.Score"]] == max(L.003gnaq[[k]][["Max.Score"]])])
  df.003gnaq[, 3] <- rep(359, 28)
  df.003gnaq[, 4] <- df.003gnaq[, 1] / ((df.003gnaq[, 3]) + df.003gnaq[, 2])
}


# HRAS
df.004hras <- data.frame(1:nrow(pat))
L.004hras <- list()

for (k in 1:nrow(pat)) {
  L.004hras[[k]] <- data.frame(L1[[4]][grep(pat[k,], L1[[4]][["Scientific.Name"]]), ])
  df.004hras[k,] <- max((data.frame(L1[[4]][grep(pat[k,], L1[[4]][["Scientific.Name"]]), ]))$Max.Score)
  df.004hras[k, 2] <- first(L.004hras[[k]][["Acc..Len"]][L.004hras[[k]][["Max.Score"]] == max(L.004hras[[k]][["Max.Score"]])])
  df.004hras[, 3] <- rep(188, 28)
  df.004hras[, 4] <- df.004hras[, 1] / ((df.004hras[, 3]) + df.004hras[, 2])
}

# KRAS
df.005kras <- data.frame(1:nrow(pat))
L.005kras <- list()

for (k in 1:nrow(pat)) {
  L.005kras[[k]] <- data.frame(L1[[5]][grep(pat[k,], L1[[5]][["Scientific.Name"]]), ])
  df.005kras[k,] <- max((data.frame(L1[[5]][grep(pat[k,], L1[[5]][["Scientific.Name"]]), ]))$Max.Score)
  df.005kras[k, 2] <- first(L.005kras[[k]][["Acc..Len"]][L.005kras[[k]][["Max.Score"]] == max(L.005kras[[k]][["Max.Score"]])])
  df.005kras[, 3] <- rep(207, 28)
  df.005kras[, 4] <- df.005kras[, 1] / ((df.005kras[, 3]) + df.005kras[, 2])
}


# NRAS
df.006nras <- data.frame(1:nrow(pat))
L.006nras <- list()

for (k in 1:nrow(pat)) {
  L.006nras[[k]] <- data.frame(L1[[6]][grep(pat[k,], L1[[6]][["Scientific.Name"]]), ])
  df.006nras[k,] <- max((data.frame(L1[[6]][grep(pat[k,], L1[[6]][["Scientific.Name"]]), ]))$Max.Score)
  df.006nras[k, 2] <- first(L.006nras[[k]][["Acc..Len"]][L.006nras[[k]][["Max.Score"]] == max(L.006nras[[k]][["Max.Score"]])])
  df.006nras[, 3] <- rep(189, 28)
  df.006nras[, 4] <- df.006nras[, 1] / ((df.006nras[, 3]) + df.006nras[, 2])
}


# RAF1
df.raf1 <- data.frame(1:nrow(pat))
L.raf1 <- list()

for (k in 1:nrow(pat)) {
  L.raf1[[k]] <- data.frame(L1[[7]][grep(pat[k,], L1[[7]][["Scientific.Name"]]), ])
  df.raf1[k,] <- max((data.frame(L1[[7]][grep(pat[k,], L1[[7]][["Scientific.Name"]]), ]))$Max.Score)
  df.raf1[k, 2] <- first(L.raf1[[k]][["Acc..Len"]][L.raf1[[k]][["Max.Score"]] == max(L.raf1[[k]][["Max.Score"]])])
  df.raf1[, 3] <- rep(650, 28)
  df.raf1[, 4] <- df.raf1[, 1] / ((df.raf1[, 3]) + df.raf1[, 2])
}


# MAP2K1
df.008map2k1 <- data.frame(1:nrow(pat))
L.008map2k1 <- list()

for (k in 1:nrow(pat)) {
  L.008map2k1[[k]] <- data.frame(L1[[8]][grep(pat[k,], L1[[8]][["Scientific.Name"]]), ])
  df.008map2k1[k,] <- max((data.frame(L1[[8]][grep(pat[k,], L1[[8]][["Scientific.Name"]]), ]))$Max.Score)
  df.008map2k1[k, 2] <- first(L.008map2k1[[k]][["Acc..Len"]][L.008map2k1[[k]][["Max.Score"]] == max(L.008map2k1[[k]][["Max.Score"]])])
  df.008map2k1[, 3] <- rep(393, 28)
  df.008map2k1[, 4] <- df.008map2k1[, 1] / ((df.008map2k1[, 3]) + df.008map2k1[, 2])
}


# MAP2K2
df.009map2k2 <- data.frame(1:nrow(pat))
L.009map2k2 <- list()

for (k in 1:nrow(pat)) {
  L.009map2k2[[k]] <- data.frame(L1[[9]][grep(pat[k,], L1[[9]][["Scientific.Name"]]), ])
  df.009map2k2[k,] <- max((data.frame(L1[[9]][grep(pat[k,], L1[[9]][["Scientific.Name"]]), ]))$Max.Score)
  df.009map2k2[k, 2] <- first(L.009map2k2[[k]][["Acc..Len"]][L.009map2k2[[k]][["Max.Score"]] == max(L.009map2k2[[k]][["Max.Score"]])])
  df.009map2k2[, 3] <- rep(400, 28)
  df.009map2k2[, 4] <- df.009map2k2[, 1] / ((df.009map2k2[, 3]) + df.009map2k2[, 2])
}


# MAPK1
df.010mapk1 <- data.frame(1:nrow(pat))
L.010mapk1 <- list()

for (k in 1:nrow(pat)) {
  L.010mapk1[[k]] <- data.frame(L1[[10]][grep(pat[k,], L1[[10]][["Scientific.Name"]]), ])
  df.010mapk1[k,] <- max((data.frame(L1[[10]][grep(pat[k,], L1[[10]][["Scientific.Name"]]), ]))$Max.Score)
  df.010mapk1[k, 2] <- first(L.010mapk1[[k]][["Acc..Len"]][L.010mapk1[[k]][["Max.Score"]] == max(L.010mapk1[[k]][["Max.Score"]])])
  df.010mapk1[, 3] <- rep(364, 28)
  df.010mapk1[, 4] <- df.010mapk1[, 1] / ((df.010mapk1[, 3]) + df.010mapk1[, 2])
}


# PLA2G4A
df.013pla2g4a <- data.frame(1:nrow(pat))
L.013pla2g4a <- list()

for (k in 1:nrow(pat)) {
  L.013pla2g4a[[k]] <- data.frame(L1[[11]][grep(pat[k,], L1[[11]][["Scientific.Name"]]), ])
  df.013pla2g4a[k,] <- max((data.frame(L1[[11]][grep(pat[k,], L1[[11]][["Scientific.Name"]]), ]))$Max.Score)
  df.013pla2g4a[k, 2] <- first(L.013pla2g4a[[k]][["Acc..Len"]][L.013pla2g4a[[k]][["Max.Score"]] == max(L.013pla2g4a[[k]][["Max.Score"]])])
  df.013pla2g4a[, 3] <- rep(780, 28)
  df.013pla2g4a[, 4] <- df.013pla2g4a[, 1] / ((df.013pla2g4a[, 3]) + df.013pla2g4a[, 2])
}


# PLA2G4F
df.018pla2g4f <- data.frame(1:nrow(pat))
L.018pla2g4f <- list()

for (k in 1:nrow(pat)) {
  L.018pla2g4f[[k]] <- data.frame(L1[[12]][grep(pat[k,], L1[[12]][["Scientific.Name"]]), ])
  df.018pla2g4f[k,] <- max((data.frame(L1[[12]][grep(pat[k,], L1[[12]][["Scientific.Name"]]), ]))$Max.Score)
  df.018pla2g4f[k, 2] <- first(L.018pla2g4f[[k]][["Acc..Len"]][L.018pla2g4f[[k]][["Max.Score"]] == max(L.018pla2g4f[[k]][["Max.Score"]])])
  df.018pla2g4f[, 3] <- rep(856, 28)
  df.018pla2g4f[, 4] <- df.018pla2g4f[, 1] / ((df.018pla2g4f[, 3]) + df.018pla2g4f[, 2])
}


# PTGS2
df.019ptgs2 <- data.frame(1:nrow(pat))
L.019ptgs2 <- list()

for (k in 1:nrow(pat)) {
  L.019ptgs2[[k]] <- data.frame(L1[[13]][grep(pat[k,], L1[[13]][["Scientific.Name"]]), ])
  df.019ptgs2[k,] <- max((data.frame(L1[[13]][grep(pat[k,], L1[[13]][["Scientific.Name"]]), ]))$Max.Score)
  df.019ptgs2[k, 2] <- first(L.019ptgs2[[k]][["Acc..Len"]][L.019ptgs2[[k]][["Max.Score"]] == max(L.019ptgs2[[k]][["Max.Score"]])])
  df.019ptgs2[, 3] <- rep(731, 28)
  df.019ptgs2[, 4] <- df.019ptgs2[, 1] / ((df.019ptgs2[, 3]) + df.019ptgs2[, 2])
}


# MAP2K5
df.020map2k5 <- data.frame(1:nrow(pat))
L.020map2k5 <- list()

for (k in 1:nrow(pat)) {
  L.020map2k5[[k]] <- data.frame(L1[[14]][grep(pat[k,], L1[[14]][["Scientific.Name"]]), ])
  df.020map2k5[k,] <- max((data.frame(L1[[14]][grep(pat[k,], L1[[14]][["Scientific.Name"]]), ]))$Max.Score)
  df.020map2k5[k, 2] <- first(L.020map2k5[[k]][["Acc..Len"]][L.020map2k5[[k]][["Max.Score"]] == max(L.020map2k5[[k]][["Max.Score"]])])
  df.020map2k5[, 3] <- rep(444, 28)
  df.020map2k5[, 4] <- df.020map2k5[, 1] / ((df.020map2k5[, 3]) + df.020map2k5[, 2])
}


# JUN
df.022jun <- data.frame(1:nrow(pat))
L.022jun <- list()

for (k in 1:nrow(pat)) {
  L.022jun[[k]] <- data.frame(L1[[15]][grep(pat[k,], L1[[15]][["Scientific.Name"]]), ])
  df.022jun[k,] <- max((data.frame(L1[[15]][grep(pat[k,], L1[[15]][["Scientific.Name"]]), ]))$Max.Score)
  df.022jun[k, 2] <- first(L.022jun[[k]][["Acc..Len"]][L.022jun[[k]][["Max.Score"]] == max(L.022jun[[k]][["Max.Score"]])])
  df.022jun[, 3] <- rep(318, 28)
  df.022jun[, 4] <- df.022jun[, 1] / ((df.022jun[, 3]) + df.022jun[, 2])
}


# FOS
df.023fos <- data.frame(1:nrow(pat))
L.023fos <- list()

for (k in 1:nrow(pat)) {
  L.023fos[[k]] <- data.frame(L1[[16]][grep(pat[k,], L1[[16]][["Scientific.Name"]]), ])
  df.023fos[k,] <- max((data.frame(L1[[16]][grep(pat[k,], L1[[16]][["Scientific.Name"]]), ]))$Max.Score)
  df.023fos[k, 2] <- first(L.023fos[[k]][["Acc..Len"]][L.023fos[[k]][["Max.Score"]] == max(L.023fos[[k]][["Max.Score"]])])
  df.023fos[, 3] <- rep(301, 28)
  df.023fos[, 4] <- df.023fos[, 1] / ((df.023fos[, 3]) + df.023fos[, 2])
}


# 024MEF2C
df.024mef2c <- data.frame(1:nrow(pat))
L.024mef2c <- list()

for (k in 1:nrow(pat)) {
  L.024mef2c[[k]] <- data.frame(L1[[17]][grep(pat[k,], L1[[17]][["Scientific.Name"]]), ])
  df.024mef2c[k,] <- max((data.frame(L1[[17]][grep(pat[k,], L1[[17]][["Scientific.Name"]]), ]))$Max.Score)
  df.024mef2c[k, 2] <- first(L.024mef2c[[k]][["Acc..Len"]][L.024mef2c[[k]][["Max.Score"]] == max(L.024mef2c[[k]][["Max.Score"]])])
  df.024mef2c[, 3] <- rep(465, 28)
  df.024mef2c[, 4] <- df.024mef2c[, 1] / ((df.024mef2c[, 3]) + df.024mef2c[, 2])
}


# 025CCND1
df.025ccnd1 <- data.frame(1:nrow(pat))
L.025ccnd1 <- list()

for (k in 1:nrow(pat)) {
  L.025ccnd1[[k]] <- data.frame(L1[[18]][grep(pat[k,], L1[[18]][["Scientific.Name"]]), ])
  df.025ccnd1[k,] <- max((data.frame(L1[[18]][grep(pat[k,], L1[[18]][["Scientific.Name"]]), ]))$Max.Score)
  df.025ccnd1[k, 2] <- first(L.025ccnd1[[k]][["Acc..Len"]][L.025ccnd1[[k]][["Max.Score"]] == max(L.025ccnd1[[k]][["Max.Score"]])])
  df.025ccnd1[, 3] <- rep(292, 28)
  df.025ccnd1[, 4] <- df.025ccnd1[, 1] / ((df.025ccnd1[, 3]) + df.025ccnd1[, 2])
}


# 026ELK1
df.026elk1 <- data.frame(1:nrow(pat))
L.026elk1 <- list()

for (k in 1:nrow(pat)) {
  L.026elk1[[k]] <- data.frame(L1[[19]][grep(pat[k,], L1[[19]][["Scientific.Name"]]), ])
  df.026elk1[k,] <- max((data.frame(L1[[19]][grep(pat[k,], L1[[19]][["Scientific.Name"]]), ]))$Max.Score)
  df.026elk1[k, 2] <- first(L.026elk1[[k]][["Acc..Len"]][L.026elk1[[k]][["Max.Score"]] == max(L.026elk1[[k]][["Max.Score"]])])
  df.026elk1[, 3] <- rep(372, 28)
  df.026elk1[, 4] <- df.026elk1[, 1] / ((df.026elk1[, 3]) + df.026elk1[, 2])
}


# 028RYR2
df.028ryr2 <- data.frame(1:nrow(pat))
L.028ryr2 <- list()

for (k in 1:nrow(pat)) {
  L.028ryr2[[k]] <- data.frame(L1[[20]][grep(pat[k,], L1[[20]][["Scientific.Name"]]), ])
  df.028ryr2[k,] <- max((data.frame(L1[[20]][grep(pat[k,], L1[[20]][["Scientific.Name"]]), ]))$Max.Score)
  df.028ryr2[k, 2] <- first(L.028ryr2[[k]][["Acc..Len"]][L.028ryr2[[k]][["Max.Score"]] == max(L.028ryr2[[k]][["Max.Score"]])])
  df.028ryr2[, 3] <- rep(4953, 28)
  df.028ryr2[, 4] <- df.028ryr2[, 1] / ((df.028ryr2[, 3]) + df.028ryr2[, 2])
}


# 029RYR3
df.029RYR3 <- data.frame(1:nrow(pat))
L.029RYR3 <- list()

for (k in 1:nrow(pat)) {
  L.029RYR3[[k]] <- data.frame(L1[[21]][grep(pat[k,], L1[[21]][["Scientific.Name"]]), ])
  df.029RYR3[k,] <- max((data.frame(L1[[21]][grep(pat[k,], L1[[21]][["Scientific.Name"]]), ]))$Max.Score)
  df.029RYR3[k, 2] <- first(L.029RYR3[[k]][["Acc..Len"]][L.029RYR3[[k]][["Max.Score"]] == max(L.029RYR3[[k]][["Max.Score"]])])
  df.029RYR3[, 3] <- rep(4881, 28)
  df.029RYR3[, 4] <- df.029RYR3[, 1] / ((df.029RYR3[, 3]) + df.029RYR3[, 2])
}


# 030CD38
df.030CD38 <- data.frame(1:nrow(pat))
L.030CD38 <- list()

for (k in 1:nrow(pat)) {
  L.030CD38[[k]] <- data.frame(L1[[22]][grep(pat[k,], L1[[22]][["Scientific.Name"]]), ])
  df.030CD38[k,] <- max((data.frame(L1[[22]][grep(pat[k,], L1[[22]][["Scientific.Name"]]), ]))$Max.Score)
  df.030CD38[k, 2] <- first(L.030CD38[[k]][["Acc..Len"]][L.030CD38[[k]][["Max.Score"]] == max(L.030CD38[[k]][["Max.Score"]])])
  df.030CD38[, 3] <- rep(291, 28)
  df.030CD38[, 4] <- df.030CD38[, 1] / ((df.030CD38[, 3]) + df.030CD38[, 2])
}


# -------------------------------------------- 032KCNJ2
df.032KCNJ2 <- data.frame(1:nrow(pat))
L.032KCNJ2 <- list()

for (k in 1:nrow(pat)) {
  L.032KCNJ2[[k]] <- data.frame(L1[[23]][grep(pat[k,], L1[[23]][["Scientific.Name"]]), ])
  df.032KCNJ2[k,] <- max((data.frame(L1[[23]][grep(pat[k,], L1[[23]][["Scientific.Name"]]), ]))$Max.Score)
  df.032KCNJ2[k, 2] <- first(L.032KCNJ2[[k]][["Acc..Len"]][L.032KCNJ2[[k]][["Max.Score"]] == max(L.032KCNJ2[[k]][["Max.Score"]])])
  df.032KCNJ2[, 3] <- rep(426, 28)
  df.032KCNJ2[, 4] <- df.032KCNJ2[, 1] / ((df.032KCNJ2[, 3]) + df.032KCNJ2[, 2])
}


# -------------------------------------------- 035KCNJ4
df.035KCNJ4 <- data.frame(1:nrow(pat))
L.035KCNJ4 <- list()

for (k in 1:nrow(pat)) {
  L.035KCNJ4[[k]] <- data.frame(L1[[24]][grep(pat[k,], L1[[24]][["Scientific.Name"]]), ])
  df.035KCNJ4[k,] <- max((data.frame(L1[[24]][grep(pat[k,], L1[[24]][["Scientific.Name"]]), ]))$Max.Score)
  df.035KCNJ4[k, 2] <- first(L.035KCNJ4[[k]][["Acc..Len"]][L.035KCNJ4[[k]][["Max.Score"]] == max(L.035KCNJ4[[k]][["Max.Score"]])])
  df.035KCNJ4[, 3] <- rep(426, 28)
  df.035KCNJ4[, 4] <- df.035KCNJ4[, 1] / ((df.035KCNJ4[, 3]) + df.035KCNJ4[, 2])
}


# -------------------------------------------- 037PLCB1
df.037PLCB1 <- data.frame(1:nrow(pat))
L.037PLCB1 <- list()

for (k in 1:nrow(pat)) {
  L.037PLCB1[[k]] <- data.frame(L1[[25]][grep(pat[k,], L1[[25]][["Scientific.Name"]]), ])
  df.037PLCB1[k,] <- max((data.frame(L1[[25]][grep(pat[k,], L1[[25]][["Scientific.Name"]]), ]))$Max.Score)
  df.037PLCB1[k, 2] <- first(L.037PLCB1[[k]][["Acc..Len"]][L.037PLCB1[[k]][["Max.Score"]] == max(L.037PLCB1[[k]][["Max.Score"]])])
  df.037PLCB1[, 3] <- rep(1260, 28)
  df.037PLCB1[, 4] <- df.037PLCB1[, 1] / ((df.037PLCB1[, 3]) + df.037PLCB1[, 2])
}


# -------------------------------------------- 040PLCB4
df.040PLCB4 <- data.frame(1:nrow(pat))
L.040PLCB4 <- list()

for (k in 1:nrow(pat)) {
  L.040PLCB4[[k]] <- data.frame(L1[[26]][grep(pat[k,], L1[[26]][["Scientific.Name"]]), ])
  df.040PLCB4[k,] <- max((data.frame(L1[[26]][grep(pat[k,], L1[[26]][["Scientific.Name"]]), ]))$Max.Score)
  df.040PLCB4[k, 2] <- first(L.040PLCB4[[k]][["Acc..Len"]][L.040PLCB4[[k]][["Max.Score"]] == max(L.040PLCB4[[k]][["Max.Score"]])])
  df.040PLCB4[, 3] <- rep(1207, 28)
  df.040PLCB4[, 4] <- df.040PLCB4[, 1] / ((df.040PLCB4[, 3]) + df.040PLCB4[, 2])
}


# -------------------------------------------- 041PRKCA
df.041PRKCA <- data.frame(1:nrow(pat))
L.041PRKCA <- list()

for (k in 1:nrow(pat)) {
  L.041PRKCA[[k]] <- data.frame(L1[[27]][grep(pat[k,], L1[[27]][["Scientific.Name"]]), ])
  df.041PRKCA[k,] <- max((data.frame(L1[[27]][grep(pat[k,], L1[[27]][["Scientific.Name"]]), ]))$Max.Score)
  df.041PRKCA[k, 2] <- first(L.041PRKCA[[k]][["Acc..Len"]][L.041PRKCA[[k]][["Max.Score"]] == max(L.041PRKCA[[k]][["Max.Score"]])])
  df.041PRKCA[, 3] <- rep(664, 28)
  df.041PRKCA[, 4] <- df.041PRKCA[, 1] / ((df.041PRKCA[, 3]) + df.041PRKCA[, 2])
}


# -------------------------------------------- 042PRKCB
df.042PRKCB <- data.frame(1:nrow(pat))
L.042PRKCB <- list()

for (k in 1:nrow(pat)) {
  L.042PRKCB[[k]] <- data.frame(L1[[28]][grep(pat[k,], L1[[28]][["Scientific.Name"]]), ])
  df.042PRKCB[k,] <- max((data.frame(L1[[28]][grep(pat[k,], L1[[28]][["Scientific.Name"]]), ]))$Max.Score)
  df.042PRKCB[k, 2] <- first(L.042PRKCB[[k]][["Acc..Len"]][L.042PRKCB[[k]][["Max.Score"]] == max(L.042PRKCB[[k]][["Max.Score"]])])
  df.042PRKCB[, 3] <- rep(671, 28)
  df.042PRKCB[, 4] <- df.042PRKCB[, 1] / ((df.042PRKCB[, 3]) + df.042PRKCB[, 2])
}


# -------------------------------------------- 044EEF2K
df.044EEF2K <- data.frame(1:nrow(pat))
L.044EEF2K <- list()

for (k in 1:nrow(pat)) {
  L.044EEF2K[[k]] <- data.frame(L1[[29]][grep(pat[k,], L1[[29]][["Scientific.Name"]]), ])
  df.044EEF2K[k,] <- max((data.frame(L1[[29]][grep(pat[k,], L1[[29]][["Scientific.Name"]]), ]))$Max.Score)
  df.044EEF2K[k, 2] <- first(L.044EEF2K[[k]][["Acc..Len"]][L.044EEF2K[[k]][["Max.Score"]] == max(L.044EEF2K[[k]][["Max.Score"]])])
  df.044EEF2K[, 3] <- rep(715, 28)
  df.044EEF2K[, 4] <- df.044EEF2K[, 1] / ((df.044EEF2K[, 3]) + df.044EEF2K[, 2])
}


# -------------------------------------------- 045EEF2
df.045EEF2 <- data.frame(1:nrow(pat))
L.045EEF2 <- list()

for (k in 1:nrow(pat)) {
  L.045EEF2[[k]] <- data.frame(L1[[30]][grep(pat[k,], L1[[30]][["Scientific.Name"]]), ])
  df.045EEF2[k,] <- max((data.frame(L1[[30]][grep(pat[k,], L1[[30]][["Scientific.Name"]]), ]))$Max.Score)
  df.045EEF2[k, 2] <- first(L.045EEF2[[k]][["Acc..Len"]][L.045EEF2[[k]][["Max.Score"]] == max(L.045EEF2[[k]][["Max.Score"]])])
  df.045EEF2[, 3] <- rep(858, 28)
  df.045EEF2[, 4] <- df.045EEF2[, 1] / ((df.045EEF2[, 3]) + df.045EEF2[, 2])
}


# -------------------------------------------- 046CACNA1C
df.046CACNA1C <- data.frame(1:nrow(pat))
L.046CACNA1C <- list()

for (k in 1:nrow(pat)) {
  L.046CACNA1C[[k]] <- data.frame(L1[[31]][grep(pat[k,], L1[[31]][["Scientific.Name"]]), ])
  df.046CACNA1C[k,] <- max((data.frame(L1[[31]][grep(pat[k,], L1[[31]][["Scientific.Name"]]), ]))$Max.Score)
  df.046CACNA1C[k, 2] <- first(L.046CACNA1C[[k]][["Acc..Len"]][L.046CACNA1C[[k]][["Max.Score"]] == max(L.046CACNA1C[[k]][["Max.Score"]])])
  df.046CACNA1C[, 3] <- rep(2200, 28)
  df.046CACNA1C[, 4] <- df.046CACNA1C[, 1] / ((df.046CACNA1C[, 3]) + df.046CACNA1C[, 2])
}


# -------------------------------------------- 047CACNA1D
df.047CACNA1D <- data.frame(1:nrow(pat))
L.047CACNA1D <- list()

for (k in 1:nrow(pat)) {
  L.047CACNA1D[[k]] <- data.frame(L1[[32]][grep(pat[k,], L1[[32]][["Scientific.Name"]]), ])
  df.047CACNA1D[k,] <- max((data.frame(L1[[32]][grep(pat[k,], L1[[32]][["Scientific.Name"]]), ]))$Max.Score)
  df.047CACNA1D[k, 2] <- first(L.047CACNA1D[[k]][["Acc..Len"]][L.047CACNA1D[[k]][["Max.Score"]] == max(L.047CACNA1D[[k]][["Max.Score"]])])
  df.047CACNA1D[, 3] <- rep(2130, 28)
  df.047CACNA1D[, 4] <- df.047CACNA1D[, 1] / ((df.047CACNA1D[, 3]) + df.047CACNA1D[, 2])
}


# -------------------------------------------- 050CACNB1
df.050CACNB1 <- data.frame(1:nrow(pat))
L.050CACNB1 <- list()

for (k in 1:nrow(pat)) {
  L.050CACNB1[[k]] <- data.frame(L1[[33]][grep(pat[k,], L1[[33]][["Scientific.Name"]]), ])
  df.050CACNB1[k,] <- max((data.frame(L1[[33]][grep(pat[k,], L1[[33]][["Scientific.Name"]]), ]))$Max.Score)
  df.050CACNB1[k, 2] <- first(L.050CACNB1[[k]][["Acc..Len"]][L.050CACNB1[[k]][["Max.Score"]] == max(L.050CACNB1[[k]][["Max.Score"]])])
  df.050CACNB1[, 3] <- rep(800, 28)
  df.050CACNB1[, 4] <- df.050CACNB1[, 1] / ((df.050CACNB1[, 3]) + df.050CACNB1[, 2])
}


# -------------------------------------------- 051CACNB2
df.051CACNB2 <- data.frame(1:nrow(pat))
L.051CACNB2 <- list()

for (k in 1:nrow(pat)) {
  L.051CACNB2[[k]] <- data.frame(L1[[34]][grep(pat[k,], L1[[34]][["Scientific.Name"]]), ])
  df.051CACNB2[k,] <- max((data.frame(L1[[34]][grep(pat[k,], L1[[34]][["Scientific.Name"]]), ]))$Max.Score)
  df.051CACNB2[k, 2] <- first(L.051CACNB2[[k]][["Acc..Len"]][L.051CACNB2[[k]][["Max.Score"]] == max(L.051CACNB2[[k]][["Max.Score"]])])
  df.051CACNB2[, 3] <- rep(613, 28)
  df.051CACNB2[, 4] <- df.051CACNB2[, 1] / ((df.051CACNB2[, 3]) + df.051CACNB2[, 2])
}


# -------------------------------------------- 053CACNB4
df.053CACNB4 <- data.frame(1:nrow(pat))
L.053CACNB4 <- list()

for (k in 1:nrow(pat)) {
  L.053CACNB4[[k]] <- data.frame(L1[[35]][grep(pat[k,], L1[[35]][["Scientific.Name"]]), ])
  df.053CACNB4[k,] <- max((data.frame(L1[[35]][grep(pat[k,], L1[[35]][["Scientific.Name"]]), ]))$Max.Score)
  df.053CACNB4[k, 2] <- first(L.053CACNB4[[k]][["Acc..Len"]][L.053CACNB4[[k]][["Max.Score"]] == max(L.053CACNB4[[k]][["Max.Score"]])])
  df.053CACNB4[, 3] <- rep(476, 28)
  df.053CACNB4[, 4] <- df.053CACNB4[, 1] / ((df.053CACNB4[, 3]) + df.053CACNB4[, 2])
}


# -------------------------------------------- 054CACNA2D1
df.054CACNA2D1 <- data.frame(1:nrow(pat))
L.054CACNA2D1 <- list()

for (k in 1:nrow(pat)) {
  L.054CACNA2D1[[k]] <- data.frame(L1[[36]][grep(pat[k,], L1[[36]][["Scientific.Name"]]), ])
  df.054CACNA2D1[k,] <- max((data.frame(L1[[36]][grep(pat[k,], L1[[36]][["Scientific.Name"]]), ]))$Max.Score)
  df.054CACNA2D1[k, 2] <- first(L.054CACNA2D1[[k]][["Acc..Len"]][L.054CACNA2D1[[k]][["Max.Score"]] == max(L.054CACNA2D1[[k]][["Max.Score"]])])
  df.054CACNA2D1[, 3] <- rep(1086, 28)
  df.054CACNA2D1[, 4] <- df.054CACNA2D1[, 1] / ((df.054CACNA2D1[, 3]) + df.054CACNA2D1[, 2])
}


# -------------------------------------------- 056CACNA2D3
df.056CACNA2D3 <- data.frame(1:nrow(pat))
L.056CACNA2D3 <- list()

for (k in 1:nrow(pat)) {
  L.056CACNA2D3[[k]] <- data.frame(L1[[37]][grep(pat[k,], L1[[37]][["Scientific.Name"]]), ])
  df.056CACNA2D3[k,] <- max((data.frame(L1[[37]][grep(pat[k,], L1[[37]][["Scientific.Name"]]), ]))$Max.Score)
  df.056CACNA2D3[k, 2] <- first(L.056CACNA2D3[[k]][["Acc..Len"]][L.056CACNA2D3[[k]][["Max.Score"]] == max(L.056CACNA2D3[[k]][["Max.Score"]])])
  df.056CACNA2D3[, 3] <- rep(1085, 28)
  df.056CACNA2D3[, 4] <- df.056CACNA2D3[, 1] / ((df.056CACNA2D3[, 3]) + df.056CACNA2D3[, 2])
}


# -------------------------------------------- 057CACNA2D4
df.057CACNA2D4 <- data.frame(1:nrow(pat))
L.057CACNA2D4 <- list()

for (k in 1:nrow(pat)) {
  L.057CACNA2D4[[k]] <- data.frame(L1[[38]][grep(pat[k,], L1[[38]][["Scientific.Name"]]), ])
  df.057CACNA2D4[k,] <- max((data.frame(L1[[38]][grep(pat[k,], L1[[38]][["Scientific.Name"]]), ]))$Max.Score)
  df.057CACNA2D4[k, 2] <- first(L.057CACNA2D4[[k]][["Acc..Len"]][L.057CACNA2D4[[k]][["Max.Score"]] == max(L.057CACNA2D4[[k]][["Max.Score"]])])
  df.057CACNA2D4[, 3] <- rep(1195, 28)
  df.057CACNA2D4[, 4] <- df.057CACNA2D4[, 1] / ((df.057CACNA2D4[, 3]) + df.057CACNA2D4[, 2])
}


# -------------------------------------------- 058CACNG1
df.058CACNG1 <- data.frame(1:nrow(pat))
L.058CACNG1 <- list()

for (k in 1:nrow(pat)) {
  L.058CACNG1[[k]] <- data.frame(L1[[39]][grep(pat[k,], L1[[39]][["Scientific.Name"]]), ])
  df.058CACNG1[k,] <- max((data.frame(L1[[39]][grep(pat[k,], L1[[39]][["Scientific.Name"]]), ]))$Max.Score)
  df.058CACNG1[k, 2] <- first(L.058CACNG1[[k]][["Acc..Len"]][L.058CACNG1[[k]][["Max.Score"]] == max(L.058CACNG1[[k]][["Max.Score"]])])
  df.058CACNG1[, 3] <- rep(227, 28)
  df.058CACNG1[, 4] <- df.058CACNG1[, 1] / ((df.058CACNG1[, 3]) + df.058CACNG1[, 2])
}


# -------------------------------------------- 059CACNG2
df.059CACNG2 <- data.frame(1:nrow(pat))
L.059CACNG2 <- list()

for (k in 1:nrow(pat)) {
  L.059CACNG2[[k]] <- data.frame(L1[[40]][grep(pat[k,], L1[[40]][["Scientific.Name"]]), ])
  df.059CACNG2[k,] <- max((data.frame(L1[[40]][grep(pat[k,], L1[[40]][["Scientific.Name"]]), ]))$Max.Score)
  df.059CACNG2[k, 2] <- first(L.059CACNG2[[k]][["Acc..Len"]][L.059CACNG2[[k]][["Max.Score"]] == max(L.059CACNG2[[k]][["Max.Score"]])])
  df.059CACNG2[, 3] <- rep(322, 28)
  df.059CACNG2[, 4] <- df.059CACNG2[, 1] / ((df.059CACNG2[, 3]) + df.059CACNG2[, 2])
}





# ...... now for all CSV files in one loop?
















