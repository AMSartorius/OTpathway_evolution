


# annotate original data frame

jtl$cat <- "Majority pattern"

jtl[jtl$brain_region %in% c('V1C', 'MD'),]$cat <- "Childhood peak"

jtl[jtl$brain_region %in% c('M1C', 'IPC', 'AMY'),]$cat <- "BAdolescence drop"

jtl[jtl$brain_region %in% c('STR'),]$cat <- "Combined"

jtl$Stage <- factor(jtl$Stage , levels = c("Prenatal", "Infant",
                                           "Child", "Adolescent", 
                                           "Adult"))




# ------------- install package(s)

library(ggforce)



# ------------- analyses

tai_devALL <-
  function(region) 
  {
    exp <- get_expression(structure_ids=c(region),
                          gene_ids=HomoSapiens.PhyloMap2013$EnsemblGeneID, 
                          dataset='5_stages') # Get expression values
    
    
    exp_df <- as.data.frame(do.call(rbind, exp)) 
    rownames(exp_df) <- NULL
    rownames(exp_df) <- c("Prenatal", "Infant", "Child", "Adolescent", "Adult")
    exp_t <- as.data.frame(t(exp_df))
    exp_t <- rownames_to_column(exp_t, var = "gene")                         
    
    names(exp_t)[names(exp_t) == "gene"] <- "Gene_ID"
    
    mm <- MatchMap(HomoSapiens.PhyloMap2013, exp_t)
    
    
    ps <- PlotSignature( ExpressionSet = mm,
                         measure       = "TAI", 
                         TestStatistic = "FlatLineTest",
                         p.value       = FALSE,
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
  geom_line(jtlALL, mapping=aes(x= Stage, y=TAI, group=brain_region, color=cat), size = 1) +
  theme_bw()

p <- p + geom_line(jtl, mapping= aes(x= Stage, y=TAI, group=brain_region, color=cat), size= 1) +
  theme_bw() +
  scale_color_manual(values = c("grey50", "#FAA100", "#81C800", "#FF0564", "#104C7E"),
                     labels = c("All protein-coding genes", "Adolescence drop", "Childhood peak", "Combined", "Majority pattern"))  

p <- p + 
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(face="bold", size=23),
        axis.text.x  = element_text(angle=18, vjust=.5, hjust=0.5, size=21),
        axis.title.y = element_text(face="bold", size=21),
        axis.text.y = element_text(size=21),
        axis.line = element_line(color = 'black'))

p <- p + labs(color = "TAI trajectories") + xlab("\nOntogenetic stage") + ylab("Transcriptome Age Index\n")

p <- p + facet_zoom(ylim=c(1.4, 1.75)) 

p <- p +
  annotate('text', x = 0.75, y = 1.683, label = 'STR', colour= "#FF0564", size=7) +
  annotate('text', x = 0.75, y = 1.618, label = 'V1C', colour= "#81C800", size=7) +
  annotate('text', x = 0.75, y = 1.636, label = 'MD', colour= "#81C800", size=7) +
  annotate('text', x = 0.75, y = 1.711, label = 'AMY', colour= "#FAA100", size=7) +
  annotate('text', x = 0.75, y = 1.661, label = 'M1C', colour= "#FAA100", size=7) +
  annotate('text', x = 0.75, y = 1.627, label = 'IPC', colour= "#FAA100", size=7) 


p <- p + theme(legend.position = c(.45, .8),
               legend.title = element_text(size=20, face="bold"),
               legend.text = element_text(size=20)) +
  guides(color = guide_legend(override.aes = list(size=5)))

pb <- ggplot_build(p)
pb$data[[3]][1, 'alpha'] <- 0
pb$data[[4]][1, 'alpha'] <- 0
pb$data[[5]][1, 'alpha'] <- 0
pb$data[[6]][1, 'alpha'] <- 0
pb$data[[7]][1, 'alpha'] <- 0
pb$data[[8]][1, 'alpha'] <- 0
pb$data[[2]][81:160,'size'] <- 1.5
pg <- ggplot_gtable(pb)
plot(pg)

ggsave(filename = paste0(BASE, "RProject/OT-bd-analyses_loc/output/TAIcomOT_ALL2013b.pdf"), pg,
       width = 15, height = 10, units = "in", device='pdf')



