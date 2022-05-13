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

# annotate with phylostrata
genensid_ps <- dplyr::inner_join(HomoSapiens.PhyloMap, genensid, by = c("Gene_ID" = "GENEID")) # 18844 genes


# subset more modern genes
yGenes <- genensid_ps[genensid_ps$Phylostratum %in% yPS, ] # keep only more modern phylostrata, 5771 genes

# subset older genes
oGenes <- genensid_ps[!genensid_ps$Phylostratum %in% yPS, ] # keep only older phylostrata, 13073 genes


# ...... subset modern OT genes
yOT <- yGenes[yGenes$SYMBOL %in% kop, ] # keep only more modern OT genes, 23 genes

# ...... subset older OT genes
oOT <- oGenes[oGenes$SYMBOL %in% kop, ] # keep only older OT genes, 123 genes






##############################################################
##############################################################
################### cortical expression ######################
##############################################################
##############################################################



# subset abagen cortical expression matrix into young genes
ahba_yc <- exp_m_anno_c[, names(exp_m_anno_c) %in% c("id", "region", "hemisphere", "structure", yGenes$SYMBOL)] # extract subset of more modern genes with cortical expression
dim(ahba_yc)                                     # 3227 genes
ahba_yc_t <- data.frame(t(ahba_yc))
ahba_yc_t_s <- ahba_yc_t[-c(1:4),]


# extract young OT subset
OT_yc <- ahba_yc[, names(ahba_yc) %in% c("id", "region", "hemisphere", "structure", yOT$SYMBOL)] # extract subset of modern OT genes with cortical expression
dim(OT_yc) # 21 modern genes, 68 cortical regions

# average expression across all genes in OT df
OT_yc[, 5:25] <- OT_yc[,5:25] %>% mutate_if(is.character,as.numeric)
OT_yc$mean_OTy <- rowSums(OT_yc[,5:25])/21 
OT_yc_m <- OT_yc[,c(1:4,26)]


# generate 100 random "young" gene sets with 21 random genes out of transposed "abagen cortical expression for young genes" - subset
fun_yc1 <- function(i) {                                # function to generate random gene set of Ngene = 19
  sample_n(ahba_yc_t_s, 21) 
}

erg <- lapply(1:100,fun_yc1)                            # execute function 100 times


l_yc1 <- list() 
l_yc2 <- list()
l_yc3 <- list()
for (i in 1:100) {                                      
  temp1 <- rbind(ahba_yc_t[(1:4),], erg[[i]])           # sub-function to add the 3 descriptive rows back at the beginning of each new gene set (df)
  l_yc1[[length(l_yc1)+1]] = temp1                      # store results in a list
  temp2 <- t(l_yc1[[i]])                                # sub-function to transpose each gene set (df) back to its original format
  temp3 <- assign(paste0("rgsc_",i), temp2)              # store results matrices in the **environment**
  l_yc2[[length(l_yc2)+1]] = temp2                      # store results as vectors (?) in a list
  l_yc3[[i]] <- data.frame(temp3)                       # store results as dfs in a list
  }

# check random gene set and compare matrix structure with OT_yc
View(rgs_100)


# average expression across all genes in each random df
# for some reason the columns are not numeric, so first convert them to numeric
list1 <- list()
for (i in l_yc3) {
temp4 <- i[,5:25] %>% mutate_if(is.character,as.numeric)
list1[[length(list1)+1]] <- temp4
}

# calculate average expression across genes for each random gene set
list2 <- list() 
for (i in list1) {
  temp5 <- rowSums(i)/21
  list2[[length(list2)+1]] <- temp5
}

# add 4 descriptive cols, rename col storing average expression vals, and generate final **environment** dataframes
list3 <- list()
for (i in 1:100){
 # temp6 <- cbind(ahba_yc[,1:4], list2[[i]])   # append 4 info cols
 # colnames(temp6)[5] = paste0("rgs_", i)      # change colname
 # assign(paste0("rgs_m_",i), temp6)           # make new dfs
 # list3[[length(list3)+1]] = temp6            # store everything in a list just in case
 
   temp7  <- do.call(cbind, list2)       
   temp8 <- cbind(OT_yc_m, temp7)
   colnames(temp8)[6:105] <- paste0("rgs_", colnames(temp8[, 6:105]))
   }



# ----------------------------- OT young cortical

yc_ex_fun <- function(colname) {
mean <- mean(temp8[[colname]])             
regions <- unique(temp8$region)
pvals <- numeric()
tStatistic <- numeric()
cohensD <- numeric()

for (i in 1:length(regions)){
  Temp <- t.test(temp8[temp8$region==regions[i],colname], 
                 mu = mean, alternative = "two.sided")
  CohD <- ES.t.one(m=Temp$estimate,
                   sd=sd(temp8[temp8$region==regions[i],colname]),
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

value <- list(
  results = results
) # Create a list of output objects
attr(value, "class") <- "yc_ex_fun"
value

}


names <- colnames(temp8[5:105])

seika_yrgs <- sapply(names, yc_ex_fun, USE.NAMES = T)

# plot OTyc only
p_fun <- function(colname) {
pOTyc <- ggseg(.data=temp8, mapping=aes_string(fill=colname), colour="white", size=1, position="stacked", hemisphere="left" ) +
  theme_void()
pOTyc <- pOTyc + viridis::scale_fill_viridis(option = 'turbo') + 
  labs(title=paste0(colname, "average cortical expression"), subtitle="lateral | medial", fill="mRNA intensity") +
  theme(plot.title = element_text(hjust = .5, size=16),
        plot.subtitle = element_text(hjust = .5, size=14),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12))
value <- list(                                    
  pOTyc = pOTyc
) # Create a list of output objects
attr(value, "class") <- "p_fun"
value
}

p_res <- sapply(names, p_fun, USE.NAMES = T)


## NEXT UP: compare results



## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> OLDER PHYLOSTRATA



# subset abagen cortical expression matrix into older genes
ahba_oc <- exp_m_anno_c[, names(exp_m_anno_c) %in% c("id", "region", "hemisphere", "structure", oGenes$SYMBOL)] # extract subset of older genes with cortical expression
dim(ahba_oc)
ahba_oc_t <- data.frame(t(ahba_oc))
ahba_oc_t_s <- ahba_oc_t[-c(1:4),]



##############################################################
##############################################################
################# subcortical expression #####################
##############################################################
##############################################################




# subset abagen subcortical expression matrix into young genes
ahba_ysc <- exp_m_anno_sc[, names(exp_m_anno_sc) %in% c("id", "region", "hemisphere", "structure", yGenes$SYMBOL)] # extract subset of more modern genes with subcortical expression
dim(ahba_ysc)                                     # 3223 genes
ahba_ysc_t <- data.frame(t(ahba_ysc))
ahba_ysc_t_s <- ahba_ysc_t[-c(1:4),]


# extract young OT subset
OT_ysc <- ahba_ysc[, names(ahba_ysc) %in% c("id", "region", "hemisphere", "structure", yOT$SYMBOL)] # extract subset of modern OT genes with cortical expression
dim(OT_ysc) # 21 modern genes, 15 cortical regions

# average expression across all genes in OT df
OT_ysc[, 5:25] <- OT_ysc[,5:25] %>% mutate_if(is.character,as.numeric)
OT_ysc$mean_OTy <- rowSums(OT_ysc[,5:25])/21 
OT_ysc_m <- OT_ysc[,c(1:4,26)]


# generate 100 random "young" gene sets with 21 random genes out of transposed "abagen subcortical expression for young genes" - subset
fun_ysc1 <- function(i) {                                # function to generate random gene set of Ngene = 21
  sample_n(ahba_ysc_t_s, 21) 
}

erg_ysc <- lapply(1:100,fun_ysc1)                            # execute function 100 times


l_ysc1 <- list() 
l_ysc2 <- list()
l_ysc3 <- list()
for (i in 1:100) {                                      
  temp12 <- rbind(ahba_ysc_t[(1:4),], erg_ysc[[i]])           # sub-function to add the 3 descriptive rows back at the beginning of each new gene set (df)
  l_ysc1[[length(l_ysc1)+1]] = temp12                        # store results in a list
  temp22 <- t(l_ysc1[[i]])                                   # sub-function to transpose each gene set (df) back to its original format
  temp32 <- assign(paste0("rgssc_",i), temp22)                 # store results matrices in the **environment**
  l_ysc2[[length(l_ysc2)+1]] = temp22                        # store results as vectors (?) in a list
  l_ysc3[[i]] <- data.frame(temp32)                          # store results as dfs in a list
}

# check random gene set and compare matrix structure with OT_yc
View(rgssc_100)


# average expression across all genes in each random df
# for some reason the columns are not numeric, so first convert them to numeric
list12 <- list()
for (i in l_ysc3) {
  temp42 <- i[,5:25] %>% mutate_if(is.character,as.numeric)
  list12[[length(list12)+1]] <- temp42
}

# calculate average expression across genes for each random gene set
list22 <- list() 
for (i in list12) {
  temp52 <- rowSums(i)/21
  list22[[length(list22)+1]] <- temp52
}

# add 4 descriptive cols, rename col storing average expression vals, and generate final **environment** dataframes
list32 <- list()
for (i in 1:100){
  # temp6 <- cbind(ahba_yc[,1:4], list2[[i]])   # append 4 info cols
  # colnames(temp6)[5] = paste0("rgs_", i)      # change colname
  # assign(paste0("rgs_m_",i), temp6)           # make new dfs
  # list3[[length(list3)+1]] = temp6            # store everything in a list just in case
  
  temp72  <- do.call(cbind, list22)       
  temp82 <- cbind(OT_ysc_m, temp72)
  colnames(temp82)[6:105] <- paste0("rgssc_", colnames(temp82[, 6:105]))
}



# ----------------------------- OT young subcortical

temp82_nobs <- temp82[-15,]

ysc_ex_fun <- function(colname) {
  mean <- mean(temp82_nobs[[colname]])             
  regions <- unique(temp82_nobs$region)
  pvals <- numeric()
  tStatistic <- numeric()
  cohensD <- numeric()
  
  for (i in 1:length(regions)){
    Temp <- t.test(temp82_nobs[temp82_nobs$region==regions[i],colname], 
                   mu = mean, alternative = "two.sided")
    CohD <- ES.t.one(m=Temp$estimate,
                     sd=sd(temp82_nobs[temp82_nobs$region==regions[i],colname]),
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
  
  value <- list(
    results = results
  ) # Create a list of output objects
  attr(value, "class") <- "ysc_ex_fun"
  value
  
}


names_ysc <- colnames(temp82[5:105])

seika_yscrgs <- sapply(names_ysc, ysc_ex_fun, USE.NAMES = T)



















# plot OTyc only
#p_fun <- function(colname) {
#  pOTyc <- ggseg(.data=temp8, mapping=aes_string(fill=colname), colour="white", size=1, position="stacked", hemisphere="left" ) +
#    theme_void()
#  pOTyc <- pOTyc + viridis::scale_fill_viridis(option = 'turbo') + 
#    labs(title=paste0(colname, "average cortical expression"), subtitle="lateral | medial", fill="mRNA intensity") +
#    theme(plot.title = element_text(hjust = .5, size=16),
#          plot.subtitle = element_text(hjust = .5, size=14),
#          legend.title = element_text(size=14),
#          legend.text = element_text(size=12))
#  value <- list(                                    
#    pOTyc = pOTyc
#  ) # Create a list of output objects
#  attr(value, "class") <- "p_fun"
#  value
#}

#p_res <- sapply(names, p_fun, USE.NAMES = T)


## NEXT UP: compare results





























# >>>>>>>>>>>>>>>>>>>>>< junk and trial and error <<<<<<<<<<<<<<<<<<<<<<<<<<< #

#toffee3 <- list()
#for (i in 1:100){
#  nougat <- cbind(ahba_yc[,1:4], tatin5[[i]]) 
#  toffee3[[length(toffee3)+1]] = nougat
#}



list1 <- list()
for (j in 1:100) {
  cherry <- assign(paste0("che_",j), test) 
  list1[[i]] <- data.frame(cherry)
}




testfun <- function(colname) {
  mean(temp8[[colname]])     
}

testfun("mean_OTy")



mean(temp8[["mean_OTy"]])


mean(temp8$mean_OTy)







rgs_33 <- data.frame(rgs_33)
rgs_33[,5:25] <- sapply(rgs_33[,5:25], as.numeric)
rgs_33$mean <- rowSums(rgs_33[,5:25])/21



rgs_100 <- data.frame(rgs_100)
rgs_100[,5:25] <- sapply(rgs_100[,5:25], as.numeric)
rgs_100$mean <- rowSums(rgs_100[,5:25])/21





