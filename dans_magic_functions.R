

# plot absolute OT genes per PS in horizontal bar plot

mm_hyp %>% 
  group_by(Phylostratum) %>%
  summarise(no_rows = length(Phylostratum))

plyr::count(mm_hyp$Phylostratum) # Counts the numbers, add zeros

# Manually extract counts and add zeros for phylostrata 9, 10, 11

c1 <- c(94,19,5,4,1,7,4,0,0,1,0,3) 
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

mmhypdf <- data.frame(x=c1, counts=labs)

mmhypdf$counts <- factor(mmhypdf$counts, levels = mmhypdf$counts) 


# TO DO: Plot each bar in a different colour...
# Update: I think its fine this way. Or else it gets too colourful.And the bars are easy to discriminate as is.

pplot <- ggplot(mmhypdf, 
                # keep all aesthetics in one place
                aes(x = x, y = counts, label = x)) +
  # replacement of geom_bar(stat = "identity")
  geom_col(fill = "#56B4E9") +
  # avoid overlap of text and bar to make text visible as bar and text have the same colour 
  geom_text(nudge_x = 5, size = 8) +
  theme_classic() +
  xlab("Number of genes \n per phylostratum ") + ylab("\nPhylostratum")

pplot <- pplot + theme(axis.text.x = element_text(size = 20),
                       axis.text.y = element_text(size = 20),
                       axis.title.x = element_text(size=22, face="bold"),
                       axis.title.y = element_text(size=22, face="bold"))

pplot       





# A function to modify plots

vert_axis_plot <-
  function(object) 
  {
    plot <- object + theme_classic()
    plot <- plot + theme(axis.text.x = 
                           element_text(color="white", 
                                        size=12, angle=45, hjust = 1),
                         axis.text.y = element_text(size = 12 )) +
      xlab("") + ylab("")
    plot
    
    plot <- plot + scale_y_continuous(breaks=seq(0.5,3,0.5))
    value <- list(
      plot = plot
    ) # Create a list of output objects
    attr(value, "class") <- "vert_axis_plot"
    value
  }

# Another function to modify plots

dis_plot <-
  function(object) 
  {
    q <- ggplot_build(object)
    q$data[[2]]$size <- 2
    q$data[[2]]$colour <- "#56B4E9"
    q$data[[1]]$fill <- "#56B4E9"
    q$data[[1]]$alpha <- 0.12
    q <- ggplot_gtable(q)
    plot(q)
    object <- as.ggplot(function() plot(q))
    plot
    value <- list(
      plot = plot
    ) # Create a list of output objects
    attr(value, "class") <- "dis_plot"
    value
  }



# ----------------------------------------------- M1C

#M1C_pm <- M1C$pm + ggtitle("M1C")
#M1C_ps <- M1C$ps + ggtitle("M1C")


#M1C_ps <- M1C_ps + theme_classic()
#M1C_ps <- M1C_ps + theme(axis.text.x = 
#                           element_text(color="white", 
#                                        size=12, angle=45, hjust = 1),
#                         axis.text.y = element_text(size = 12 )) +
#                        xlab("") + ylab("")
#M1C_ps

#M1C_ps <- M1C_ps + scale_y_continuous(breaks=seq(0.5,3,0.5))
#M1C_ps <- M1C_ps + theme(plot.title = element_text(size=22, face="bold.italic"),
#                        axis.text.y = element_text(size=22))
#M1C_ps

#q <- ggplot_build(M1C_ps)
#q$data[[2]]$size <- 2
#q$data[[2]]$colour <- "#56B4E9"
#q$data[[1]]$fill <- "#56B4E9"
#q$data[[1]]$alpha <- 0.12
#q <- ggplot_gtable(q)
#plot(q)
#M1C_ps <- as.ggplot(function() plot(q))
#M1C_ps

# ----------------------------------------------- DFC

#DFC_pm <- DFC$pm + ggtitle("DFC") 

#DFC_ps <- DFC$ps + ggtitle("DFC") 
#DFC_ps <- DFC_ps + theme_classic()
#DFC_ps <- DFC_ps + theme(axis.text.x = 
#                           element_text(color="white", 
#                                        size=12, angle=45, hjust = 1),
#                         axis.text.y = element_text(size = 12, color = "white")) +
#                        xlab("") + ylab("")

#DFC_ps <- DFC_ps + scale_y_continuous(breaks=seq(0.5,3,0.5))
#DFC_ps <- DFC_ps + theme(plot.title = element_text(size=22, face="bold.italic"))
#DFC_ps

#q <- ggplot_build(DFC_ps)
#q$data[[2]]$size <- 2
#q$data[[2]]$colour <- "#56B4E9"
#q$data[[1]]$fill <- "#56B4E9"
#q$data[[1]]$alpha <- 0.12
#q <- ggplot_gtable(q)
#plot(q)
#DFC_ps <- as.ggplot(function() plot(q))
#DFC_ps


# ----------------------------------------------- VFC

#VFC_ps <- VFC$ps + ggtitle("VFC")
#VFC_ps <- VFC_ps + theme_classic()
#VFC_ps <- VFC_ps + theme(axis.text.x = 
#                           element_text(color="white", 
#                                        size=12, angle=45, hjust = 1),
#                         axis.text.y = element_text(size = 12, color = "white")) +
#                        xlab("") + ylab("")

#VFC_ps <- VFC_ps + scale_y_continuous(breaks=seq(0.5,3,0.5))
#VFC_ps <- VFC_ps + theme(plot.title = element_text(size=22, face="bold.italic"))
#VFC_ps

#q <- ggplot_build(VFC_ps)
#q$data[[2]]$size <- 2
#q$data[[2]]$colour <- "#56B4E9"
#q$data[[1]]$fill <- "#56B4E9"
#q$data[[1]]$alpha <- 0.12
#q <- ggplot_gtable(q)
#plot(q)
#VFC_ps <- as.ggplot(function() plot(q))
#VFC_ps

#VFC_pm <- VFC$pm + ggtitle("VFC") 


# ----------------------------------------------- OFC

#OFC_pm <- OFC$pm + ggtitle("OFC") 

#OFC_ps <- OFC$ps + ggtitle("OFC") 
#OFC_ps <- OFC_ps + theme_classic()
#OFC_ps <- OFC_ps + theme(axis.text.x = 
#                           element_text(color="white", 
#                                        size=12, angle=45, hjust = 1),
#                         axis.text.y = element_text(size = 12, color = "white")) +
#                        xlab("") + ylab("")

#OFC_ps <- OFC_ps + scale_y_continuous(breaks=seq(0.5,3,0.5))
#OFC_ps <- OFC_ps + theme(plot.title = element_text(size=22, face="bold.italic"))
#OFC_ps

#q <- ggplot_build(OFC_ps)
#q$data[[2]]$size <- 2
#q$data[[2]]$colour <- "#56B4E9"
#q$data[[1]]$fill <- "#56B4E9"
#q$data[[1]]$alpha <- 0.12
#q <- ggplot_gtable(q)
#plot(q)
#OFC_ps <- as.ggplot(function() plot(q))
#OFC_ps


# ----------------------------------------------- S1C

#S1C_pm <- S1C$pm + ggtitle("S1C")

#S1C_ps <- S1C$ps + ggtitle("S1C")
#S1C_ps <- S1C_ps + theme_classic()
#S1C_ps <- S1C_ps + theme(axis.text.x = 
#                           element_text(color="white", 
#                                        size=12, angle=45, hjust = 1),
#                         axis.text.y = element_text(size = 12 )) +
#                       xlab("") + ylab("")
#S1C_ps

#S1C_ps <- S1C_ps + scale_y_continuous(breaks=seq(0.5,3,0.5))
#S1C_ps <- S1C_ps + theme(plot.title = element_text(size=22, face="bold.italic"),
#                         axis.text.y = element_text(size=22))
#S1C_ps

#q <- ggplot_build(S1C_ps)
#q$data[[2]]$size <- 2
#q$data[[2]]$colour <- "#56B4E9"
#q$data[[1]]$fill <- "#56B4E9"
#q$data[[1]]$alpha <- 0.12
#q <- ggplot_gtable(q)
#plot(q)
#S1C_ps <- as.ggplot(function() plot(q))
#S1C_ps


# ----------------------------------------------- IPC

#IPC_pm <- IPC$pm + ggtitle("IPC")

#IPC_ps <- IPC$ps + ggtitle("IPC")
#IPC_ps <- IPC_ps + theme_classic()
#IPC_ps <- IPC_ps + theme(axis.text.x = 
#                           element_text(color="white", 
#                                        size=12, angle=45, hjust = 1),
#                         axis.text.y = element_text(size = 12, color = "white")) +
#                        xlab("") + ylab("")

#IPC_ps <- IPC_ps + scale_y_continuous(breaks=seq(0.5,3,0.5))
#IPC_ps <- IPC_ps + theme(plot.title = element_text(size=22, face="bold.italic"))
#IPC_ps

#q <- ggplot_build(IPC_ps)
#q$data[[2]]$size <- 2
#q$data[[2]]$colour <- "#56B4E9"
#q$data[[1]]$fill <- "#56B4E9"
#q$data[[1]]$alpha <- 0.12
#q <- ggplot_gtable(q)
#plot(q)
#IPC_ps <- as.ggplot(function() plot(q))
#IPC_ps


# ----------------------------------------------- A1C

#A1C_pm <- A1C$pm + ggtitle("A1C")

#A1C_ps <- A1C$ps + ggtitle("A1C")
#A1C_ps <- A1C_ps + theme_classic()
#A1C_ps <- A1C_ps + theme(axis.text.x = 
#                           element_text(color="white", 
#                                        size=12, angle=45, hjust = 1),
#                         axis.text.y = element_text(size = 12, color = "white")) +
#                        xlab("") + ylab("")

#A1C_ps <- A1C_ps + scale_y_continuous(breaks=seq(0.5,3,0.5))
#A1C_ps <- A1C_ps + theme(plot.title = element_text(size=22, face="bold.italic"))
#A1C_ps

#q <- ggplot_build(A1C_ps)
#q$data[[2]]$size <- 2
#q$data[[2]]$colour <- "#56B4E9"
#q$data[[1]]$fill <- "#56B4E9"
#q$data[[1]]$alpha <- 0.12
#q <- ggplot_gtable(q)
#plot(q)
#A1C_ps <- as.ggplot(function() plot(q))
#A1C_ps


# ----------------------------------------------- STC

#STC_pm <- STC$pm + ggtitle("STC")

#STC_ps <- STC$ps + ggtitle("STC")
#STC_ps <- STC_ps + theme_classic()
#STC_ps <- STC_ps + theme(axis.text.x = 
#                           element_text(color="white", 
#                                        size=12, angle=45, hjust = 1),
#                         axis.text.y = element_text(size = 12, color = "white")) +
#                        xlab("") + ylab("")

#STC_ps <- STC_ps + scale_y_continuous(breaks=seq(0.5,3,0.5))
#STC_ps <- STC_ps + theme(plot.title = element_text(size=22, face="bold.italic"))
#STC_ps

#q <- ggplot_build(STC_ps)
#q$data[[2]]$size <- 2
#q$data[[2]]$colour <- "#56B4E9"
#q$data[[1]]$fill <- "#56B4E9"
#q$data[[1]]$alpha <- 0.12
#q <- ggplot_gtable(q)
#plot(q)
#STC_ps <- as.ggplot(function() plot(q))
#STC_ps



# Now add the vertical code



# ----------------------------------------------- ITC

#ITC_pm <- ITC$pm + ggtitle("ITC")
#ITC_ps <- ITC$ps + ggtitle("ITC")

#ITC_ps <- ITC_ps + theme_classic()
#ITC_ps <- ITC_ps + theme(axis.text.x = 
#                           element_text(color="white", 
#                                        size=12, angle=45, hjust = 1),
#                         axis.text.y = element_text(size = 12 )) +
#                        xlab("") + ylab("")
#ITC_ps

#ITC_ps <- ITC_ps + scale_y_continuous(breaks=seq(0.5,3,0.5))
#ITC_ps <- ITC_ps + theme(plot.title = element_text(size=22, face="bold.italic"),
#                         axis.text.y = element_text(size=22))
#ITC_ps

#q <- ggplot_build(ITC_ps)
#q$data[[2]]$size <- 2
#q$data[[2]]$colour <- "#56B4E9"
#q$data[[1]]$fill <- "#56B4E9"
#q$data[[1]]$alpha <- 0.12
#q <- ggplot_gtable(q)
#plot(q)
#ITC_ps <- as.ggplot(function() plot(q))
#ITC_ps


# ----------------------------------------------- V1C

#V1C_pm <- V1C$pm + ggtitle("V1C")

#V1C_ps <- V1C$ps + ggtitle("V1C")
#V1C_ps <- V1C_ps + theme_classic()
#V1C_ps <- V1C_ps + theme(axis.text.x = 
#                           element_text(color="white", 
#                                        size=12, angle=45, hjust = 1),
#                         axis.text.y = element_text(size = 12, color = "white")) +
#                        xlab("") + ylab("")

#V1C_ps <- V1C_ps + scale_y_continuous(breaks=seq(0.5,3,0.5))
#V1C_ps <- V1C_ps + theme(plot.title = element_text(size=22, face="bold.italic"))
#V1C_ps

#q <- ggplot_build(V1C_ps)
#q$data[[2]]$size <- 2
#q$data[[2]]$colour <- "#56B4E9"
#q$data[[1]]$fill <- "#56B4E9"
#q$data[[1]]$alpha <- 0.12
#q <- ggplot_gtable(q)
#plot(q)
#V1C_ps <- as.ggplot(function() plot(q))
#V1C_ps


# ----------------------------------------------- MFC

#MFC_pm <- MFC$pm + ggtitle("MFC")

#MFC_ps <- MFC$ps + ggtitle("MFC")
#MFC_ps <- MFC_ps + theme_classic()
#MFC_ps <- MFC_ps + theme(axis.text.x = 
#                           element_text(color="white", 
#                                        size=12, angle=45, hjust = 1),
#                         axis.text.y = element_text(size = 12, color = "white")) +
#                        xlab("") + ylab("")

#MFC_ps <- MFC_ps + scale_y_continuous(breaks=seq(0.5,3,0.5))
#MFC_ps <- MFC_ps + theme(plot.title = element_text(size=22, face="bold.italic"))
#MFC_ps

#q <- ggplot_build(MFC_ps)
#q$data[[2]]$size <- 2
#q$data[[2]]$colour <- "#56B4E9"
#q$data[[1]]$fill <- "#56B4E9"
#q$data[[1]]$alpha <- 0.12
#q <- ggplot_gtable(q)
#plot(q)
#MFC_ps <- as.ggplot(function() plot(q))
#MFC_ps


# ----------------------------------------------- HIP

#HIP_pm <- HIP$pm + ggtitle("HIP")

#HIP_ps <- HIP$ps + ggtitle("HIP")
#HIP_ps <- HIP_ps + theme_classic()
#HIP_ps <- HIP_ps + theme(axis.text.x = 
#                           element_text(color="white", 
#                                        size=12, angle=45, hjust = 1),
#                         axis.text.y = element_text(size = 12, color = "white")) +
#                         xlab("") + ylab("")

#HIP_ps <- HIP_ps + scale_y_continuous(breaks=seq(0.5,3,0.5))
#HIP_ps <- HIP_ps + theme(plot.title = element_text(size=22, face="bold.italic"))
#HIP_ps

#q <- ggplot_build(HIP_ps)
#q$data[[2]]$size <- 2
#q$data[[2]]$colour <- "#56B4E9"
#q$data[[1]]$fill <- "#56B4E9"
#q$data[[1]]$alpha <- 0.12
#q <- ggplot_gtable(q)
#plot(q)
#HIP_ps <- as.ggplot(function() plot(q))
#HIP_ps


# ----------------------------------------------- STR

#STR_pm <- STR$pm + ggtitle("STR")

#STR_ps <- STR$ps + ggtitle("STR")
#STR_ps <- STR_ps + theme_classic()
#STR_ps <- STR_ps + theme(axis.text.x = 
#                           element_text(color="black", 
#                                        size=22, angle=45, hjust = 1),
#                         axis.text.y = element_text(size = 12 )) +
#                        xlab("") + ylab("")
#STR_ps

#STR_ps <- STR_ps + scale_y_continuous(breaks=seq(0.5,3,0.5))
#STR_ps <- STR_ps + theme(plot.title = element_text(size=22, face="bold.italic"),
#                         axis.text.y = element_text(size=22))
#STR_ps

#q <- ggplot_build(STR_ps)
#q$data[[2]]$size <- 2
#q$data[[2]]$colour <- "#56B4E9"
#q$data[[1]]$fill <- "#56B4E9"
#q$data[[1]]$alpha <- 0.12
#q <- ggplot_gtable(q)
#plot(q)
#STR_ps <- as.ggplot(function() plot(q))
#STR_ps


# ----------------------------------------------- AMY

#AMY_pm <- AMY$pm + ggtitle("AMY")

#AMY_ps <- AMY$ps + ggtitle("AMY")
#AMY_ps <- AMY_ps + theme_classic()
#AMY_ps <- AMY_ps + theme(axis.text.x = 
#                           element_text(color="black", 
#                                        size=22, angle=45, hjust = 1),
#                         axis.text.y = element_text(size = 12, color = "white")) +
#                        xlab("") + ylab("")

#AMY_ps <- AMY_ps + scale_y_continuous(breaks=seq(0.5,3,0.5))
#AMY_ps <- AMY_ps + theme(plot.title = element_text(size=22, face="bold.italic"))
#AMY_ps

#q <- ggplot_build(AMY_ps)
#q$data[[2]]$size <- 2
#q$data[[2]]$colour <- "#56B4E9"
#q$data[[1]]$fill <- "#56B4E9"
#q$data[[1]]$alpha <- 0.12
#q <- ggplot_gtable(q)
#plot(q)
#AMY_ps <- as.ggplot(function() plot(q))
#AMY_ps

# ----------------------------------------------- MD

#MD_pm <- MD$pm + ggtitle("MD")

#MD_ps <- MD$ps + ggtitle("MD")
#MD_ps <- MD_ps + theme_classic()
#MD_ps <- MD_ps + theme(axis.text.x = 
#                         element_text(color="black", 
#                                      size=22, angle=45, hjust = 1),
#                         axis.text.y = element_text(size = 12, color = "white")) +
#                         xlab("") + ylab("")

#MD_ps <- MD_ps + scale_y_continuous(breaks=seq(0.5,3,0.5))
#MD_ps <- MD_ps + theme(plot.title = element_text(size=22, face="bold.italic"))
#MD_ps

#q <- ggplot_build(MD_ps)
#q$data[[2]]$size <- 2
#q$data[[2]]$colour <- "#56B4E9"
#q$data[[1]]$fill <- "#56B4E9"
#q$data[[1]]$alpha <- 0.12
#q <- ggplot_gtable(q)
#plot(q)
#MD_ps <- as.ggplot(function() plot(q))
#MD_ps


# ----------------------------------------------- CBC


#CBC_pm <- CBC$pm + ggtitle("CBC")

#CBC_ps <- CBC$ps + ggtitle("CBC")
#CBC_ps <- CBC_ps + theme_classic()
#CBC_ps <- CBC_ps + theme(axis.text.x = 
# element_text(color="black", 
#  size=22, angle=45, hjust = 1),
# axis.text.y = element_text(size = 12, color = "white")) +
#xlab("") + ylab("")

#CBC_ps <- CBC_ps + scale_y_continuous(breaks=seq(0.5,3,0.5))
#CBC_ps <- CBC_ps + theme(plot.title = element_text(size=22, face="bold.italic"))
#CBC_ps

#q <- ggplot_build(CBC_ps)
#q$data[[2]]$size <- 2
#q$data[[2]]$colour <- "#56B4E9"
#q$data[[1]]$fill <- "#56B4E9"
#q$data[[1]]$alpha <- 0.12
#q <- ggplot_gtable(q)
#plot(q)
#CBC_ps <- as.ggplot(function() plot(q))
#CBC_ps



#grid_ps1 <- plot_grid(M1C_ps, DFC_ps, VFC_ps, OFC_ps,
#                     S1C_ps, IPC_ps, A1C_ps, STC_ps,
#                     ITC_ps, V1C_ps, MFC_ps, HIP_ps,
#                     nrow = 3)

#grid_ps2 <- plot_grid(STR_ps, AMY_ps, MD_ps, CBC_ps,
#                      nrow = 1)

#grid_ps <- plot_grid(grid_ps1, grid_ps2,
#                     nrow = 2, rel_heights = c(1.6, 0.6))

#tai_p_com <- plot_grid(grid_ps, pplot,
#                       ncol = 2, rel_widths = c(1, 0.6),
#                       labels = c('A', 'B'),
#                      label_size = 14) # 12 * 20

















