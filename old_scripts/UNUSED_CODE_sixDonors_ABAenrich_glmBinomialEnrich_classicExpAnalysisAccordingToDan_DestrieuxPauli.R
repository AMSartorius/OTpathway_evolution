
### >>>>>>>>>>>>>>>>>>> 1 DESTRIEUX cortical

### DO NOT RUN ### 

{
  
  # convert to long format
  #des_ad_info_l <- des_ad_info %>% pivot_longer(c("mean9861", "mean10021", "mean12876", "mean14380", "mean15496", "mean15697"),        # transform into long format, the way t.test likes it :)
  #                                              names_to = "donors", values_to = "expression")
  #des_ad_info_l <- data.frame(des_ad_info_l)
  
  
  # analysis for both hemispheres
  
  #mean_des <- mean(des_ad_info_l$expression)             
  #regions_des <- unique(des_ad_info_l$name)
  #pvals_des <- numeric()
  #tStatistic_des <- numeric()
  #cohensD_des <- numeric()
  
  #for (i in 1:length(regions_des)){
  #    Temp <- t.test(des_ad_info_l[des_ad_info_l$name==regions_des[i],"expression"], 
  #                   mu = mean_des, alternative = "two.sided")
  #    CohD <- ES.t.one(m=Temp$estimate,
  #                     sd=sd(des_ad_info_l[des_ad_info_l$name==regions_des[i],"expression"]),
  #                     mu=mean_des,
  #                     alternative = "two.sided")
  #    pvals_des <- c(pvals_des,Temp$p.value)
  #    tStatistic_des <- c(tStatistic_des,Temp$statistic)
  #    cohensD_des <- c(cohensD_des,CohD$d)
  #  }
  
  #results_des <- data.frame(Region=regions_des,P_unadjusted=pvals_des,
  #                        t_stat=tStatistic_des,CohD=cohensD_des)
  #results_des$P_fdr <- p.adjust(results_des$P_unadjusted, method="fdr", n = length(results_des$P_unadjusted))
  
  #results_des$stat_sign <- ""
  #results_des$stat_sign[results_des$P_fdr <= 0.05]  <- "*"
  
  #results_des  
  
  
}

#########




### >>>>>>>>>>>>>>>>>>> 2 PAULI subcortical

### DO NOT RUN ### 

{
  
  # convert to long format
  #pau_ad_info_l <- pau_ad_info %>% pivot_longer(c("mean9861", "mean10021", "mean12876", "mean14380", "mean15496", "mean15697"),        # transform into long format, the way t.test likes it :)
  #                                               names_to = "donors", values_to = "expression")
  #pau_ad_info_l <- data.frame(pau_ad_info_l)
  
  
  # analysis
  
  #mean_pau <- mean(pau_ad_info_l$expression)             
  #regions_pau <- unique(pau_ad_info_l$node_info)
  #pvals_pau <- numeric()
  #tStatistic_pau <- numeric()
  #cohensD_pau <- numeric()
  
  #for (i in 1:length(regions_pau)){
  #  Temp <- t.test(pau_ad_info_l[pau_ad_info_l$node_info==regions_pau[i],"expression"], 
  #                 mu = mean_pau, alternative = "two.sided")
  #  CohD <- ES.t.one(m=Temp$estimate,
  #                   sd=sd(pau_ad_info_l[pau_ad_info_l$node_info==regions_pau[i],"expression"]),
  #                   mu=mean_pau,
  #                   alternative = "two.sided")
  #  pvals_pau <- c(pvals_pau,Temp$p.value)
  #  tStatistic_pau <- c(tStatistic_pau,Temp$statistic)
  #  cohensD_pau <- c(cohensD_pau,CohD$d)
  #}
  
  #results_pau <- data.frame(Region=regions_pau,P_unadjusted=pvals_pau,
  #                      t_stat=tStatistic_pau,CohD=cohensD_pau)
  #results_pau$P_fdr <- p.adjust(results_pau$P_unadjusted, method="fdr", n = length(results_pau$P_unadjusted))
  
  #results_pau$stat_sign <- ""
  #results_pau$stat_sign[results_pau$P_fdr <= 0.05]  <- "*"
  
  #results_pau  
  
}

######





### >>>>>>>>>>>>>>>>>>> 2 cortical and subcortical combined, both hemispheres (in the original script, only left hemi is used)


### DO NOT RUN ###

{
  
  
  ## both hemispheres
  
  #pau_ad_info_noV2 <- pau_ad_info[,-2]
  #names(pau_ad_info_noV2)[names(pau_ad_info_noV2) == "node_info"] <- "name"
  #names(pau_ad_info_noV2)[names(pau_ad_info_noV2) == "V1"] <- "index"
  
  #des_pau <- rbind(des_ad_info, pau_ad_info_noV2)
  
  #des_pauL <- des_pau %>% pivot_longer(c("mean9861", "mean10021", "mean12876", "mean14380", "mean15496", "mean15697"),        # transform into long format, the way t.test likes it :)
  #                                               names_to = "donors", values_to = "expression")
  
  #des_pauL <- data.frame(des_pauL)
  #names(des_pauL)[names(des_pauL) == "name"] <- "region"
  
  
  # analysis
  
  #mean_dp <- mean(des_pauL$expression)             
  #regions_dp <- unique(des_pauL$region)
  #pvals_dp <- numeric()
  #tStatistic_dp <- numeric()
  #cohensD_dp <- numeric()
  #estimate_dp <- numeric()
  
  #for (i in 1:length(regions_dp)){
  #  Temp_dp <- t.test(des_pauL[des_pauL$region==regions_dp[i],"expression"], 
  #                 mu = mean_dp, alternative = "two.sided")
  #  CohD_dp <- ES.t.one(m=Temp_dp$estimate,
  #                   sd=sd(des_pauL[des_pauL$region==regions_dp[i],"expression"]),
  #                   mu=mean_dp,
  #                   alternative = "two.sided")
  #  pvals_dp <- c(pvals_dp,Temp_dp$p.value)
  #  tStatistic_dp <- c(tStatistic_dp,Temp_dp$statistic)
  #  cohensD_dp <- c(cohensD_dp,CohD_dp$d)
  #  estimate_dp <- c(estimate_dp, Temp_dp$estimate)
  #}
  
  #results_dp <- data.frame(Region=regions_dp,P_unadjusted=pvals_dp,
  #                          t_stat=tStatistic_dp,CohD=cohensD_dp)
  #results_dp$P_fdr <- p.adjust(results_dp$P_unadjusted, method="fdr")
  
  #results_dp$stat_sign <- ""
  #results_dp$stat_sign[results_dp$P_fdr <= 0.05]  <- "*"
  
  #results_dp  
  
}

##########



# probably not mention this one for now #

################################################################
## OPTION 3: Brain region enrichment analysis with ABAEnrichment
################################################################

# simple enrichment analysis: test user-defined genes (here: young OT genset) for expression enrichment in different human brain regions.


# # ---------------- enrichment of young OT pathway gene set in the brain with ABAEnrichment

#library(ABAEnrichment)
#gene_ids = yOT$SYMBOL
#input_hyper = data.frame(gene_ids, is_candidate=1)
#head(input_hyper)

#res_adult = aba_enrich(input_hyper, dataset='adult', gene_len = T, n_randsets=10000) # default for random permutations: 1000

#res <- res_adult$results # assign and view

#pattern <- c("Left|left")  # define pattern to fetch left hemisphere

#res_lhemi <- res[grep(pattern, res$structure),]  # keep only left hemisphere

# no enrichment either...?








### >>>>>>>>>>>>>>>>>>> Analysis 1: Visualization of subcortical expression with barplot (not implemented, violin rain plots used instead)


# --------------------- barplot with jitter

#region_bpj <- unique(pauli_allR_L$region)
#sd_bpj <- numeric()

# get SD for each region 
#for (i in 1:length(region_bpj)){
#  sd_temp =sd(pauli_allR_L[pauli_allR_L$region==region_bpj[i],"expression"])
#  sd_bpj = c(sd_bpj, sd_temp)
#}

#df_sd <- data.frame(sd=sd_bpj)

#df_sum <- data.frame(cbind(pauli_allR$region, pauli_allR$mean2, df_sd$sd))
#names(df_sum)[names(df_sum) == "X1"] <- "region"
#names(df_sum)[names(df_sum) == "X2"] <- "expression"
#names(df_sum)[names(df_sum) == "X3"] <- "sd"
#df_sum$expression <- as.numeric(df_sum$expression)
#df_sum$sd <- as.numeric(df_sum$sd)



#bp_j <-ggplot(pauli_allR_L, aes(x=factor(region), y=expression, fill=region, color=region)) + theme_classic() +
#  viridis::scale_fill_viridis(option="plasma", discrete=T, end=.8) +
#  viridis::scale_color_viridis(option="plasma", discrete=T, end=.8) +

#  geom_bar(stat = "identity", data = df_sum, alpha=.3) +

#  geom_jitter( position = position_jitter(0.2), size=4, shape=21, color="black") + 

#  geom_errorbar(data=df_sum, aes(ymin = expression-sd, ymax = expression+sd), width = 0.2) +

#  theme(legend.position="none",
#        axis.text.x = element_text(angle=50, vjust=1, hjust=1, size=13),
#        axis.title = element_text(size=20, face="bold")) +
#  ylab("mRNA intensity") + xlab("\nSubcortical brain regions")

#bp_j <- bp_j + geom_label_repel(data=(pauli_allR_L %>% dplyr::filter(genes %in% c("OXTR"))), 
#                                aes(label     = genes), 
#                                nudge_x       = 0, 
#                                nudge_y       = 0,
#                                label.padding =.25,
#                                colour        = "black",
#                                fill          = alpha(c("black"),0.25),
#                                size          = 3) 
#bp_j <- bp_j + geom_label_repel(data=(pauli_allR_L %>% dplyr::filter(genes %in% c("CD38"))), 
#                                aes(label     = genes), 
#                                nudge_x       = 0, 
#                                nudge_y       = 0,
#                                label.padding = .25,
#                                colour        = "gray45",
#                                fill          = alpha(c("gray33"),0.25),
#                                size          = 3)
#bp_j <- bp_j + geom_label_repel(data=(pauli_allR_L %>% dplyr::filter(genes %in% c("OXT"))), 
#                                aes(label     = genes), 
#                                nudge_x       = 0, 
#                                nudge_y       = 0,
#                                label.padding = .25,
#                                colour        = "gray66",
#                                fill          = alpha(c("gray66"),0.25),
#                                size          = 3)

#bp_j <- bp_j + geom_hline(yintercept=mean(pauli_allR_L$mean2), linetype="dashed", size=.75)

#bp_j

#ggsave(filename = paste0(BASE, "RProject/OT-bd-analyses_loc/output/subcortical_expression_analysis1_barplot2.pdf"), bp_j,
#       width = 20, height = 10, units = "in", device='pdf')








#######################

# ........... analysis 1 figure like figure 1 in Dans nat com paper

#hlinesL <- c(unique(meanS_dpL$mean_allstruc - SD_mean_dpL), unique(meanS_dpL$mean_allstruc), unique(meanS_dpL$mean_allstruc + SD_mean_dpL))

#p2 <- ggplot(meanS_dpL, aes(x=regions_dpL, y=estimate_dpL)) + theme_classic() +
#  geom_hline(yintercept=hlinesL,  linetype=c("dashed", "solid", "dashed"), color=c("gray50", "black", "gray50"), size=c(1, 2, 1))  
#p2 <- p2 +
#  geom_point(size=6, aes(colour = regions_dpL)) +
#  geom_errorbar(aes(ymin = estimate_dpL - sd, ymax = estimate_dpL + sd, color=factor(regions_dpL), width = .5)) +
#  viridis::scale_colour_viridis(option="viridis", discrete = T)
#p2 <- p2 + ylab("mRNA intensity\n") + xlab("\nBrain regions")
#p2 <- p2 + theme(axis.text.x = element_text(angle=50, vjust=1, hjust=1, size=12),
#                axis.title.x = element_text(size=25, face="bold"),
#                axis.title.y = element_text(size=22, face="bold"),
#                axis.text.y = element_text(size=12),
#                legend.position = "none")  
#p2 

# important: this way of annotation only work if you keep the width and height of the pdf file EXACTLY the same

#pa <- grid.arrange(p2, right = textGrob("+1SD\n\n\n\nMean\n\n\n\n-1SD", vjust =-.23, gp=gpar(fontsize=15,font=3)))

#ggsave(filename = paste0(BASE, "RProject/OT-bd-analyses_loc/output/test1004.pdf"), pa,
#       width = 24, height = 9.5, units = "in", device='pdf')




##########################

# ........... horizontal cortical plotting of t-values for analysis 1

# horizontal

#p987 <- ggseg(.data = resdat_left_ro,
#              atlas = desterieux,
#              hemisphere = "left", 
#              mapping = aes(fill = t_stat),
#              colour = "black",
#              size = .1, 
#              position = "dispersed")
#p987 <- p987 + scale_fill_viridis_c(option = "mako", direction = -1) +
#  theme_void() +
#  labs(fill='t statistic') 
#p987 <- p987 + labs(title= "Left hemisphere", subtitle="lateral | medial\n") +
#  theme(plot.title = element_text(hjust=.5, size=25, face="bold"),
#        plot.subtitle = element_text(hjust=.5, size=18))

#p987

#ggsave(filename = paste0(BASE, "RProject/OT-bd-analyses_loc/output/t_statAna1_cortical.pdf"), p987,
#       width = 10, height = 4, units = "in", device='pdf')






#########################

# ..... horizontal version of the half violin plots for visualization of the subcortical expression of analysis 1

# plot specific data prep
# order by mean expression, DESCENDING
pauli_aR_L_ord <- pauli_aR_L %>% 
  arrange(desc(mean2))

# manually set levels to make sure they are plotted in descending order on the x-axis
order <- unique(pauli_aR_L_ord$region)
pauli_aR_L_ord$region <- factor(pauli_aR_L_ord$region, levels = order)


# separate values for OXT, OXTR, and CD38
pauli_aR_L_cd38 <- subset(pauli_aR_L_ord, grepl("^CD38$", pauli_aR_L_ord$genes)) 
pauli_aR_L_oxt <- subset(pauli_aR_L_ord, grepl("^OXT$", pauli_aR_L_ord$genes)) 
pauli_aR_L_oxtr <- subset(pauli_aR_L_ord, grepl("^OXTR$", pauli_aR_L_ord$genes)) 
pauli_aR_L_ot <- rbind(pauli_aR_L_cd38, pauli_aR_L_oxt, pauli_aR_L_oxtr)     # combine
pauli_aR_L_ot <- pauli_aR_L_ot %>%                                               # it should still be in the right order, but just to be on the very safe side
  arrange(desc(mean2))


# ... and remove them from base jitter (they are kept for the violin plot to ensure appropriate curving of the plots)
pat_ot <- c("OXT|OXTR|CD38")
pauli_aR_L_short <- subset(pauli_aR_L_ord, !grepl(pat_ot, pauli_aR_L_ord$genes))


# plotting
#g1 <- 
#  ggplot() +
#  geom_flat_violin(pauli_all_R_L_ord, mapping = aes(x = factor(region), y = expression,                         # use full dataset with all genes for violin plots
#                                                    fill = factor(region), color = factor(region)), 
#                   position = position_nudge(x = .2, y = 0), trim = F, alpha = .1, scale = "width") +
#  theme_classic() +
#  theme(axis.title.x = element_text(size = 22, face="bold"),
#        axis.title.y = element_text(size = 22, face="bold"),
#        axis.text.x = element_text(size=15, angle=50, vjust=1, hjust=1),
#        axis.text.y = element_text(size=15),
#        legend.position="none") +
#  ylab("mRNA intensity") + xlab("\nSubcortical brain regions")

#g1 <- g1 + geom_point(pauli_allR_L_short, mapping = aes(x= factor(region), y = expression,                     # use dataset w/o OT genes for base jitter
#                                                        fill = factor(region), color = factor(region)), 
#                      position = position_jitter(width = .15), size = 4, alpha = .6)

#g1 <- g1 + geom_point(pauli_allR_L_ot,    mapping = aes(x = factor(region), y = expression),                   # use dataset with only OT genes for highlight
#                      position = position_jitter(width = .15), size = 4, alpha = .5) 

#g1 <- g1 + geom_point(pauli_all_R_L_ord, mapping = aes(x = factor(region), y = mean2),                         # add mean expression points
#                    position = position_nudge(x = 0.3), size = 7) 

#g1 <- g1 + geom_hline(yintercept=mean(pauli_all_R_L_ord$mean2), linetype="dashed", size=.75)                   # add mean line of mean expression points

#g1

#ggsave(filename = paste0(BASE, "RProject/OT-bd-analyses_loc/output/subcortical_expression_analysis1_7.pdf"), g1,
#       width = 30, height = 10, units = "in", device='pdf')





