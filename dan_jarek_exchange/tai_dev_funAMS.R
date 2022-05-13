

tai_dev <-
  function(region) 
  {
    exp <- get_expression(structure_ids=c(region),
                          gene_ids=t1_gene, 
                          dataset='5_stages') # Get expression values
    
    exp_df <- as.data.frame(do.call(rbind, exp))                                 # new
    rownames(exp_df) <- NULL                                                     # new
    rownames(exp_df) <- c("prenatal", "infant", "child", "adolescent", "adult")  # new
    exp_t <- as.data.frame(t(exp_df))                                            # new
    exp_t <- rownames_to_column(exp_t, var = "gene")                             # new
    
    ens_gene_dev <- ensembldb::select(EnsDb.Hsapiens.v79, 
                                      keys= exp_t$gene,                          # adjusted "exp_l_t" to "exp_t"
                                      keytype = "SYMBOL", 
                                      columns = c("SYMBOL","GENEID"))
    
    ens_gene_dev <- ens_gene_dev %>% distinct(SYMBOL, .keep_all = TRUE)
    
    names(ens_gene_dev)[names(ens_gene_dev) == "SYMBOL"] <- "gene"
    
    ex_t_f <- dplyr::inner_join(ens_gene_dev, exp_t, by = "gene")                # adjusted "exp_l_t" to "exp_t"
    
    ex_t_f <- dplyr::select(ex_t_f, -gene)
    
    names(ex_t_f)[names(ex_t_f) == "GENEID"] <- "Gene_ID"
    
    HomoSapiens.PhyloMap <- read_excel(paste0(BASE, "MBE_2008_Homo_Sapiens_PhyloMap.xls"), 
                                       sheet = 1, skip = 1)
    
    HomoSapiens.PhyloMap <- HomoSapiens.PhyloMap[ , c("Phylostratum","Gene_ID")]
    
    mm <- MatchMap(HomoSapiens.PhyloMap, ex_t_f)
    
    ps <- PlotSignature( ExpressionSet = mm,
                         measure       = "TAI", 
                         TestStatistic = "FlatLineTest",
                         p.value = FALSE,
                         xlab          = "Developmental stage", 
                         ylab          = "TAI" )
    
    ps_dat <- ps$data
    
    pc <- PlotContribution(mm, 
                           legendName = "PS")
    
    value <- list(
      ps = ps,
      ps_dat = ps_dat,
      pc = pc
    ) # Create a list of output objects
    attr(value, "class") <- "tai_dev"
    value
  }

