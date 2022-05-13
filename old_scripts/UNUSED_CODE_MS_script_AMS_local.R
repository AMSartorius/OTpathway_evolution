## ================================= Plot OXTR/OXT/CD38 expression across ontogenetic stages per brain region

# ISSUES: I'd like the lines to be smooth, but I can't seem to get geom_smooth to work...

# ----------------------------- OXTR

M1C_oxtr <- M1C_tai$gs_oxtr
M1C_oxtr$Region <- "M1C"

DFC_oxtr <- DFC_tai$gs_oxtr
DFC_oxtr$Region <- "DFC"

VFC_oxtr <- VFC_tai$gs_oxtr
VFC_oxtr$Region <- "VFC"

OFC_oxtr <- OFC_tai$gs_oxtr
OFC_oxtr$Region <- "OFC"

S1C_oxtr <- S1C_tai$gs_oxtr
S1C_oxtr$Region <- "S1C"

IPC_oxtr <- IPC_tai$gs_oxtr
IPC_oxtr$Region <- "IPC"

A1C_oxtr <- A1C_tai$gs_oxtr
A1C_oxtr$Region <- "A1C"

STC_oxtr <- STC_tai$gs_oxtr
STC_oxtr$Region <- "STC"

ITC_oxtr <- ITC_tai$gs_oxtr
ITC_oxtr$Region <- "ITC"

V1C_oxtr <- V1C_tai$gs_oxtr
V1C_oxtr$Region <- "V1C"

MFC_oxtr <- MFC_tai$gs_oxtr
MFC_oxtr$Region <- "MFC"

HIP_oxtr <- HIP_tai$gs_oxtr
HIP_oxtr$Region <- "HIP"

STR_oxtr <- STR_tai$gs_oxtr
STR_oxtr$Region <- "STR"

AMY_oxtr <- AMY_tai$gs_oxtr
AMY_oxtr$Region <- "AMY"

MD_oxtr <- MD_tai$gs_oxtr
MD_oxtr$Region <- "MD"

CBC_oxtr <- CBC_tai$gs_oxtr
CBC_oxtr$Region <- "CBC"

bind <- M1C_oxtr %>%
  rbind(DFC_oxtr) %>%
  rbind(VFC_oxtr) %>%
  rbind(OFC_oxtr) %>%
  rbind(S1C_oxtr) %>%
  rbind(IPC_oxtr) %>%
  rbind(A1C_oxtr) %>%
  rbind(STC_oxtr) %>%
  rbind(ITC_oxtr) %>%
  rbind(V1C_oxtr) %>%
  rbind(MFC_oxtr) %>%
  rbind(HIP_oxtr) %>%
  rbind(STR_oxtr) %>%
  rbind(AMY_oxtr) %>%
  rbind(MD_oxtr) %>%
  rbind(CBC_oxtr)

bind_oxtr <- bind %>% pivot_longer(c("Prenatal", 
                                     "Infant", 
                                     "Child", 
                                     "Adolescent", 
                                     "Adult"),
                                   names_to = "Stages", values_to = "mRNA")

bind_oxtr$Stages <- factor(bind_oxtr$Stages , levels = c("Prenatal", 
                                                         "Infant", 
                                                         "Child", 
                                                         "Adolescent", 
                                                         "Adult"))
bind_oxtr$Region <- as.factor(bind_oxtr$Region)

psgsOXTR <- ggplot(bind_oxtr, aes(x = Stages, y = mRNA, group=Region, colour=Region)) +
  geom_line(size=1.5) + 
  theme_minimal() 
psgsOXTR <- psgsOXTR + xlab("\n") + ylab("mRNA intensity\n") + ggtitle("OXTR")
psgsOXTR <- psgsOXTR + theme(axis.title.x = element_text(size=20, face="bold"),
                             axis.text.x = element_text(size=15, angle=40, hjust=1, vjust = 1),
                             axis.title.y = element_text(size=20, face="bold"),
                             axis.text.y = element_text(size=15),
                             plot.title = element_text(size=20))
psgsOXTR <- psgsOXTR + viridis::scale_colour_viridis(option = 'mako', discrete = T)
psgsOXTR <- psgsOXTR + stat_summary(aes(y = mRNA,group=1), fun=mean, colour="#F6044F", geom="line",group=1, size=3)
psgsOXTR <- psgsOXTR + theme(legend.position = "none")
psgsOXTR <- psgsOXTR + ylim(0,8.2)
psgsOXTR

#ggsave(filename = paste0(BASE, "RProject/OT-bd-analyses_loc/output/OXTR_expression.pdf"), psgsOXTR,
#       width = 10, height = 8, units = "in", device='pdf')



# as ribbon plots because why not

## looks cool but probably doesn't work in this context


#oxtr_rib <- ggplot(bind_oxtr, aes(x=as.numeric(factor(Stages)), y = mRNA, fill=Region)) +
#  geom_area(lwd=1, 
#            linetype=1, 
#            position='stack', 
#            colour='gray90', 
#            alpha=.8,  
#            na.rm=TRUE) +
#  scale_x_discrete(labels = levels(bind_oxtr$Stages)) 

#oxtr_rib <- oxtr_rib + viridis::scale_fill_viridis(option = 'rocket', discrete = T) 
#oxtr_rib <- oxtr_rib + xlab("Stages") + ylab("OXTR mRNA intensities") 
#oxtr_rib <- oxtr_rib + theme(legend.position = "right",
#                             axis.title.y = element_text(size=20, face="bold"))
#oxtr_rib <- oxtr_rib + theme_minimal()
#oxtr_rib


# ----------------------------- OXT

M1C_oxt <- M1C_tai$gs_oxt
M1C_oxt$Region <- "M1C"

DFC_oxt <- DFC_tai$gs_oxt
DFC_oxt$Region <- "DFC"

VFC_oxt <- VFC_tai$gs_oxt
VFC_oxt$Region <- "VFC"

OFC_oxt <- OFC_tai$gs_oxt
OFC_oxt$Region <- "OFC"

S1C_oxt <- S1C_tai$gs_oxt
S1C_oxt$Region <- "S1C"

IPC_oxt <- IPC_tai$gs_oxt
IPC_oxt$Region <- "IPC"

A1C_oxt <- A1C_tai$gs_oxt
A1C_oxt$Region <- "A1C"

STC_oxt <- STC_tai$gs_oxt
STC_oxt$Region <- "STC"

ITC_oxt <- ITC_tai$gs_oxt
ITC_oxt$Region <- "ITC"

V1C_oxt <- V1C_tai$gs_oxt
V1C_oxt$Region <- "V1C"

MFC_oxt <- MFC_tai$gs_oxt
MFC_oxt$Region <- "MFC"

HIP_oxt <- HIP_tai$gs_oxt
HIP_oxt$Region <- "HIP"

STR_oxt <- STR_tai$gs_oxt
STR_oxt$Region <- "STR"

AMY_oxt <- AMY_tai$gs_oxt
AMY_oxt$Region <- "AMY"

MD_oxt <- MD_tai$gs_oxt
MD_oxt$Region <- "MD"

CBC_oxt <- CBC_tai$gs_oxt
CBC_oxt$Region <- "CBC"

bind_oxt <- M1C_oxt %>%
  rbind(DFC_oxt) %>%
  rbind(VFC_oxt) %>%
  rbind(OFC_oxt) %>%
  rbind(S1C_oxt) %>%
  rbind(IPC_oxt) %>%
  rbind(A1C_oxt) %>%
  rbind(STC_oxt) %>%
  rbind(ITC_oxt) %>%
  rbind(V1C_oxt) %>%
  rbind(MFC_oxt) %>%
  rbind(HIP_oxt) %>%
  rbind(STR_oxt) %>%
  rbind(AMY_oxt) %>%
  rbind(MD_oxt) %>%
  rbind(CBC_oxt)

bind_oxt <- bind_oxt %>% pivot_longer(c("Prenatal", 
                                        "Infant", 
                                        "Child", 
                                        "Adolescent", 
                                        "Adult"),
                                      names_to = "Stages", values_to = "mRNA")

bind_oxt$Stages <- factor(bind_oxt$Stages , levels = c("Prenatal", 
                                                       "Infant", 
                                                       "Child", 
                                                       "Adolescent", 
                                                       "Adult"))
bind_oxt$Region <- as.factor(bind_oxt$Region)

psgsOXT <- ggplot(bind_oxt, aes(x = Stages, y = mRNA, group=Region, colour=Region)) +
  geom_line(size=1.5) + 
  theme_minimal()
psgsOXT <- psgsOXT + xlab("\nOntogenetic stages") + ylab("") + ggtitle("OXT")
psgsOXT <- psgsOXT + theme(axis.title.x = element_text(size=20, face="bold"),
                           axis.text.x = element_text(size=15, angle=40, hjust=1, vjust = 1),
                           axis.title.y = element_text(size=20, face="bold"),
                           axis.text.y = element_blank(),
                           plot.title = element_text(size=20))
psgsOXT <- psgsOXT + viridis::scale_colour_viridis(option = 'mako', discrete = T)
psgsOXT <- psgsOXT + stat_summary(aes(y = mRNA,group=1), fun=mean, colour="#F6044F", geom="line",group=1, size=3)
psgsOXT <- psgsOXT + theme(legend.position = "none") + ylim(0,8.2)
psgsOXT

#ggsave(filename = paste0(BASE, "RProject/OT-bd-analyses_loc/output/OXT_expression.pdf"), psgsOXT,
#       width = 10, height = 8, units = "in", device='pdf')


# ----------------------------- CD38

M1C_cd38 <- M1C_tai$gs_cd38
M1C_cd38$Region <- "M1C"

DFC_cd38 <- DFC_tai$gs_cd38
DFC_cd38$Region <- "DFC"

VFC_cd38 <- VFC_tai$gs_cd38
VFC_cd38$Region <- "VFC"

OFC_cd38 <- OFC_tai$gs_cd38
OFC_cd38$Region <- "OFC"

S1C_cd38 <- S1C_tai$gs_cd38
S1C_cd38$Region <- "S1C"

IPC_cd38 <- IPC_tai$gs_cd38
IPC_cd38$Region <- "IPC"

A1C_cd38 <- A1C_tai$gs_cd38
A1C_cd38$Region <- "A1C"

STC_cd38 <- STC_tai$gs_cd38
STC_cd38$Region <- "STC"

ITC_cd38 <- ITC_tai$gs_cd38
ITC_cd38$Region <- "ITC"

V1C_cd38 <- V1C_tai$gs_cd38
V1C_cd38$Region <- "V1C"

MFC_cd38 <- MFC_tai$gs_cd38
MFC_cd38$Region <- "MFC"

HIP_cd38 <- HIP_tai$gs_cd38
HIP_cd38$Region <- "HIP"

STR_cd38 <- STR_tai$gs_cd38
STR_cd38$Region <- "STR"

AMY_cd38 <- AMY_tai$gs_cd38
AMY_cd38$Region <- "AMY"

MD_cd38 <- MD_tai$gs_cd38
MD_cd38$Region <- "MD"

CBC_cd38 <- CBC_tai$gs_cd38
CBC_cd38$Region <- "CBC"

bind_cd38 <- M1C_cd38 %>%
  rbind(DFC_cd38) %>%
  rbind(VFC_cd38) %>%
  rbind(OFC_cd38) %>%
  rbind(S1C_cd38) %>%
  rbind(IPC_cd38) %>%
  rbind(A1C_cd38) %>%
  rbind(STC_cd38) %>%
  rbind(ITC_cd38) %>%
  rbind(V1C_cd38) %>%
  rbind(MFC_cd38) %>%
  rbind(HIP_cd38) %>%
  rbind(STR_cd38) %>%
  rbind(AMY_cd38) %>%
  rbind(MD_cd38) %>%
  rbind(CBC_cd38)

bind_cd38 <- bind_cd38 %>% pivot_longer(c("Prenatal", 
                                          "Infant", 
                                          "Child", 
                                          "Adolescent", 
                                          "Adult"),
                                        names_to = "Stages", values_to = "mRNA")

bind_cd38$Stages <- factor(bind_cd38$Stages , levels = c("Prenatal", 
                                                         "Infant", 
                                                         "Child", 
                                                         "Adolescent", 
                                                         "Adult"))
bind_cd38$Region <- as.factor(bind_cd38$Region)

psgscd38 <- ggplot(bind_cd38, aes(x = Stages, y = mRNA, group=Region, colour=Region)) +
  geom_line(size=1.5) + 
  theme_minimal() 
psgscd38 <- psgscd38 + xlab("\n") + ylab("") +  ggtitle("CD38")
psgscd38 <- psgscd38 + theme(axis.title.x = element_text(size=20, face="bold"),
                             axis.text.x = element_text(size=15, angle=40, hjust=1, vjust = 1),
                             axis.title.y = element_text(size=20, face="bold"),
                             axis.text.y = element_blank(),
                             legend.position = "right", 
                             plot.title = element_text(size=20))
psgscd38 <- psgscd38 + viridis::scale_colour_viridis(option = 'mako', discrete = T)

psgscd38 <- psgscd38 + stat_summary(aes(y = mRNA,group=1), fun=mean, colour="#F6044F", geom="line",group=1, size=3)
psgscd38 <- psgscd38 + ylim(0,8.2)
psgscd38

#ggsave(filename = paste0(BASE, "RProject/OT-bd-analyses_loc/output/CD38_expression.pdf"), psgscd38,
#  width = 10, height = 8, units = "in", device='pdf')




# ------------------------------------------- combine OXTR, OXT, CD38 plots

grid_psgs <- plot_grid(psgsOXTR, psgsOXT, psgscd38,
                       nrow = 1)
grid_psgs


ggsave(filename = paste0(BASE, "RProject/OT-bd-analyses_loc/output/OT_leadGenes_expression_AHBA.pdf"), grid_psgs,
       width = 18, height = 8, units = "in", device='pdf')




###################################################
###################################################
###################################################



## ---------- COMBINED *relative* expression levels of each PS per ontogenetic stage averaged across all brain regions

colnames(M1Cre_dat)[3] <- "M1C"
colnames(DFCre_dat)[3] <- "DFC"
colnames(VFCre_dat)[3] <- "VFC"
colnames(OFCre_dat)[3] <- "OFC"
colnames(S1Cre_dat)[3] <- "S1C"
colnames(IPCre_dat)[3] <- "IPC"
colnames(A1Cre_dat)[3] <- "A1C"
colnames(STCre_dat)[3] <- "STC"
colnames(ITCre_dat)[3] <- "ITC"
colnames(V1Cre_dat)[3] <- "V1C"
colnames(MFCre_dat)[3] <- "MFC"
colnames(HIPre_dat)[3] <- "HIP"
colnames(STRre_dat)[3] <- "STR"
colnames(AMYre_dat)[3] <- "AMY"
colnames(MDre_dat)[3] <- "MD"
colnames(CBCre_dat)[3] <- "CBC"


re_all_struct <- M1Cre_dat %>%
  left_join(DFCre_dat, by=c('age', 'stage')) %>%
  left_join(VFCre_dat, by=c('age', 'stage')) %>%
  left_join(OFCre_dat, by=c('age', 'stage')) %>%
  left_join(S1Cre_dat, by=c('age', 'stage')) %>%
  left_join(IPCre_dat, by=c('age', 'stage')) %>%
  left_join(A1Cre_dat, by=c('age', 'stage')) %>%
  left_join(STCre_dat, by=c('age', 'stage')) %>%
  left_join(ITCre_dat, by=c('age', 'stage')) %>%
  left_join(V1Cre_dat, by=c('age', 'stage')) %>%
  left_join(MFCre_dat, by=c('age', 'stage')) %>%
  left_join(HIPre_dat, by=c('age', 'stage')) %>%
  left_join(STRre_dat, by=c('age', 'stage')) %>%
  left_join(AMYre_dat, by=c('age', 'stage')) %>%
  left_join(MDre_dat, by=c('age', 'stage')) %>%
  left_join(CBCre_dat, by=c('age', 'stage')) 


re_all_struct$sum = rowSums(re_all_struct[,c("M1C", "DFC", "VFC", "OFC", "S1C", "IPC", "A1C", "STC", "ITC", "V1C", "MFC", "HIP", "STR", "AMY", "MD", "CBC")])
re_all_struct$all_struct_mean <- re_all_struct$sum / 16
re_all_struct <- re_all_struct[, -c(3:19)]
names(re_all_struct)[names(re_all_struct) == "age"] <- "Phylostratum"
re_all_struct$Phylostratum <- paste("PS", re_all_struct$Phylostratum, sep="")

re_all_struct$stage <- as.factor(re_all_struct$stage)
re_all_struct$Phylostratum <- factor(re_all_struct$Phylostratum, levels=c("PS1", "PS2", "PS3", "PS4", "PS5", "PS6","PS7", "PS10", "PS12"))

p_re <- ggplot(re_all_struct, aes(x= stage, y=all_struct_mean, group=Phylostratum, colour=Phylostratum)) +
  geom_line(size = 2) +
  theme_classic() 
p_re <- p_re + ylab("Relative brain-wide mRNA intensities") + 
  xlab("\nOntogenetic stage") 
p_re <- p_re + theme(axis.title=element_text(size=20,face="bold"),
                     axis.text.y = element_text(size=15),
                     axis.text.x = element_text(size = 15, angle=40, hjust=1, vjust = 1),
                     legend.title = element_text(size=16),
                     legend.text = element_text(size=13)) +
  viridis::scale_color_viridis(option = 'mako', 
                               discrete = T)
p_re

ggsave(filename = paste0(BASE, "output/PlotRE_all_struct.pdf"), p_re,
       width = 10, height = 8, units = "in", device='pdf')





###################################################
###################################################
###################################################




## ================================= MEAN expression levels of each PS per ontogenetic stage averaged across all brain regions

M1C_pm<- M1C_tai$pm
M1C_pm_dat <- M1C_pm$data
colnames(M1C_pm_dat)[3] <- "M1C"

DFC_pm <- DFC_tai$pm
DFC_pm_dat <- DFC_pm$data
colnames(DFC_pm_dat)[3] <- "DFC"

VFC_pm <- VFC_tai$pm
VFC_pm_dat <- VFC_pm$data
colnames(VFC_pm_dat)[3] <- "VFC"

OFC_pm <- OFC_tai$pm
OFC_pm_dat <- OFC_pm$data
colnames(OFC_pm_dat)[3] <- "OFC"

S1C_pm <- S1C_tai$pm
S1C_pm_dat <- S1C_pm$data
colnames(S1C_pm_dat)[3] <- "S1C"

IPC_pm <- IPC_tai$pm
IPC_pm_dat <- IPC_pm$data
colnames(IPC_pm_dat)[3] <- "IPC"

A1C_pm <- A1C_tai$pm
A1C_pm_dat <- A1C_pm$data
colnames(A1C_pm_dat)[3] <- "A1C"

STC_pm <- STC_tai$pm
STC_pm_dat <- STC_pm$data
colnames(STC_pm_dat)[3] <- "STC"

ITC_pm <- ITC_tai$pm
ITC_pm_dat <- ITC_pm$data
colnames(ITC_pm_dat)[3] <- "ITC"

V1C_pm <- V1C_tai$pm
V1C_pm_dat <- V1C_pm$data
colnames(V1C_pm_dat)[3] <- "V1C"

MFC_pm <- MFC_tai$pm
MFC_pm_dat <- MFC_pm$data
colnames(MFC_pm_dat)[3] <- "MFC"

HIP_pm <- HIP_tai$pm
HIP_pm_dat <- HIP_pm$data
colnames(HIP_pm_dat)[3] <- "HIP"

STR_pm <- STR_tai$pm
STR_pm_dat <- STR_pm$data
colnames(STR_pm_dat)[3] <- "STR"

AMY_pm <- AMY_tai$pm
AMY_pm_dat <- AMY_pm$data
colnames(AMY_pm_dat)[3] <- "AMY"

MD_pm <- MD_tai$pm
MD_pm_dat <- MD_pm$data
colnames(MD_pm_dat)[3] <- "MD"

CBC_pm <- CBC_tai$pm
CBC_pm_dat <- CBC_pm$data
colnames(CBC_pm_dat)[3] <- "CBC"


pm_all_struct <- M1C_pm_dat %>%
  left_join(DFC_pm_dat, by=c('age', 'stage')) %>%
  left_join(VFC_pm_dat, by=c('age', 'stage')) %>%
  left_join(OFC_pm_dat, by=c('age', 'stage')) %>%
  left_join(S1C_pm_dat, by=c('age', 'stage')) %>%
  left_join(IPC_pm_dat, by=c('age', 'stage')) %>%
  left_join(A1C_pm_dat, by=c('age', 'stage')) %>%
  left_join(STC_pm_dat, by=c('age', 'stage')) %>%
  left_join(ITC_pm_dat, by=c('age', 'stage')) %>%
  left_join(V1C_pm_dat, by=c('age', 'stage')) %>%
  left_join(MFC_pm_dat, by=c('age', 'stage')) %>%
  left_join(HIP_pm_dat, by=c('age', 'stage')) %>%
  left_join(STR_pm_dat, by=c('age', 'stage')) %>%
  left_join(AMY_pm_dat, by=c('age', 'stage')) %>%
  left_join(MD_pm_dat, by=c('age', 'stage')) %>%
  left_join(CBC_pm_dat, by=c('age', 'stage')) 


pm_all_struct$sum = rowSums(pm_all_struct[,c("M1C", "DFC", "VFC", "OFC", "S1C", "IPC", "A1C", "STC", "ITC", "V1C", "MFC", "HIP", "STR", "AMY", "MD", "CBC")])
pm_all_struct$all_struct_mean <- pm_all_struct$sum / 16
pm_all_struct <- pm_all_struct[, -c(3:19)]
names(pm_all_struct)[names(pm_all_struct) == "age"] <- "Phylostratum"
pm_all_struct$Phylostratum <- paste("PS", pm_all_struct$Phylostratum, sep="")

pm_all_struct$stage <- as.factor(pm_all_struct$stage)
pm_all_struct$Phylostratum <- factor(pm_all_struct$Phylostratum, levels=c("PS1", "PS2", "PS3", "PS4", "PS5", "PS6","PS7", "PS10", "PS12"))

p_pm <- ggplot(pm_all_struct, aes(x= stage, y=all_struct_mean, group=Phylostratum, colour=Phylostratum)) +
  geom_line(size = 2) +
  theme_classic() 
p_pm <- p_pm + ylab("Brain-wide mean mRNA intensities") + 
  xlab("\nOntogenetic stage") 
p_pm <- p_pm +
  theme(axis.title=element_text(size=20,face="bold"),
        axis.text.y = element_text(size=15),
        axis.text.x = element_text(size = 15, angle=40, hjust=1, vjust = 1),
        legend.title = element_text(size=16),
        legend.text = element_text(size=13)) +
  viridis::scale_color_viridis(option = 'cividis', 
                               discrete = T)
p_pm

ggsave(filename = paste0(BASE, "output/PlotMean_all_struct.pdf"), p_pm,
       width = 10, height = 8, units = "in", device='pdf')









## ================================= STR in depth visualization
# ... because it is the only region that behaves against the majority pattern at several stages?

pco <- STR_tai$pc
pco_dat <- pco$data

pd <- ggplot(pco_dat, aes(x= stage, y=par_value, group=PS, colour=PS)) +
  geom_line(size = 2) +
  theme_classic() 
pd <- pd + 
  annotate('text', x = 0.75, y = 0.3980758, label = 'PS1') +
  annotate('text', x = 5.25, y = 1.3031204, label = 'PS2') +
  annotate('text', x = 0.75, y = 1.5669891, label = "PS3") +
  annotate('text', x = 5.25, y = 1.3730324, label = "PS4") +
  annotate('text', x = 0.75, y = 1.5958877, label = "PS5") +
  annotate('text', x = 0.75, y = 1.6328868, label = "PS6") +
  annotate('text', x = 0.75, y = 1.6822660, label = "PS7") +
  annotate('text', x = 5.3, y = 1.4837371, label = "PS10") +
  annotate('text', x = 0.75, y = 1.72, label = "PS12") 
pd

pd <- pd  + theme(axis.title.x = element_text(face="bold", size=16),
                  axis.text.x  = element_text(angle=40, hjust=1, 
                                              vjust = 1, size=12),
                  axis.title.y = element_text(face="bold", size=16),
                  axis.text.y = element_text(size=12)) 
pd <- pd + ylab("Transcriptome age index (cumulative)") +
  xlab("Ontogenetic stage")+ 
  labs(color = "Phylostratum") +
  # plots the PS from 1 - 9, skips 10 and 12 even though OT genes are not present in PS 8 + 9???
  viridis::scale_color_viridis(option = 'plasma', 
                               discrete = T, 
                               labels = c("PS1", "PS2", "PS3", "PS4", "PS5", "PS6", "PS7", "PS10", "PS12")) +
  theme(legend.position = c(.12, .44))

pd

ggsave(filename = paste0(BASE, "RProject/OT-bd-analyses_loc/output/TAI_cumulative_STR.pdf"), pd,
       width = 8, height = 8, units = "in", device='pdf')

# as ribbon plot --- this is not working? 
#pco_dat$par_value <- as.numeric(pco_dat$par_value)
#str_rib <- ggplot(pco_dat, aes(x=stage, y = par_value, fill=PS, group=PS)) +
#  geom_area(lwd=1, 
#            linetype=1, 
#            position='stack', 
#            colour='gray90', 
#            alpha=.8,  
#            na.rm=TRUE) +
#  scale_x_discrete(labels = levels(pco_dat$stage)) #+ 
#geom_dl(aes(label = Region), list("top.points", cex = .6)) + 
#guides(fill = FALSE)
#str_rib <- str_rib + viridis::scale_fill_viridis(option = 'rocket', discrete = T) 
#str_rib <- str_rib + xlab("Ontogenetic stage") + ylab("Transcriptome Age Index (cumulative)") 
#str_rib <- str_rib + theme(legend.position = "right",
#                             axis.title.y = element_text(size=20, face="bold"))
#str_rib <- str_rib + theme_minimal()
#str_rib





