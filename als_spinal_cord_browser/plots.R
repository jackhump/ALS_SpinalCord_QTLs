
theme_jh <- function () { 
  theme_bw(base_size=12) %+replace% 
    theme(
      panel.grid = element_blank(),
      strip.background = element_blank(),
      #panel.border = element_blank(),
      axis.line = element_line(),
      #text = element_text(color = "black"), 
      strip.text = element_text(color = "black"),
      axis.text = element_text(colour = "black"),
      panel.background  = element_blank(),
      plot.background = element_rect(fill="white", colour=NA), 
      legend.background = element_rect(fill="transparent", colour=NA),
      legend.key = element_rect(fill="transparent", colour=NA), 
      plot.title = element_text(face = "bold.italic", size = 14),
      axis.line.y = element_line(), strip.text.x = element_text(face = "bold", margin = margin(t = 2,r = 0,b = 2,l=0))
    )
}

# match abbreviated labels to long form labels
choice_df <- tibble::tibble(
  short =  c("TPM","disease","C9orf72 status","mutations","site_of_motor_onset","age_rounded","disease_duration","sex","% Astrocytes", "% Endothelial", "% Neurons", "% Microglia", "% Oligodendrocytes", "% Pericytes","RIN","seq_platform", "tissue"),
  long =  c("TPM","Disease", "C9orf72 status", "Known ALS mutations", "Site of motor symptom onset","Age at death (decade)","Disease duration (months)","Sex","% Astrocytes", "% Endothelial", "% Neurons", "% Microglia", "% Oligodendrocytes", "% Pericytes","RIN","Sequencing platform", "Region" )
)

choice_to_short <- function(x){
  if(x %in% choice_df$short){return(x) }
  if(x %in% choice_df$long){
    choice_df$short[match(x, choice_df$long)]
  }else{
    return(NA)
  }
}

choice_to_long <- function(x){
  if(x %in% choice_df$long){return(x) }
  if(x %in% choice_df$short){
    choice_df$long[match(x, choice_df$short)]
  }else{
    return(NA)
  }
}

gene_plot <- function(mygene, counts, meta,
                      x = "disease", 
                      y = "TPM", 
                      colourby = "disease", 
                      facet = "tissue", 
                      tissues = c("Cervical", "Lumbar", "Thoracic"),
                      boxplot = FALSE,
                      corline = FALSE,
                      stats = FALSE,
                      log = FALSE,
                      return_data = FALSE){
  require(dplyr)
  require(ggplot2)
  require(stringr)
  # convert to short form
  x <- choice_to_short(x)
  y <- choice_to_short(y)
  colourby <- choice_to_short(colourby)
  
  # set labels to long form
  xlab <- choice_to_long(x)
  ylab <- choice_to_long(y)
  collab <- choice_to_long(colourby)
  if( y == "TPM" & log ==TRUE){
    ylab <- "log2(TPM + 0.01)"
  }
  
  ## allow multiple genes
  if(grepl(",", mygene)){
    mygene <- unlist(str_split(mygene, ","))
    facet <- paste0("gene + ", facet )
  }
  
  
  if( all(mygene != "") ){
  # return closest match to gene
    all_genes <- row.names(counts)
    if( ! all(mygene %in% all_genes) ){
      mygene <- all_genes[ grepl(mygene, all_genes, ignore.case = TRUE ) ]
      if(length(mygene) != 1){
        message( "Gene not found!")
        return(NULL)
      }
      message(mygene)
    }
  }else{
    mygene <- NA
  }

  
  
  df <- as.data.frame(t(counts[ mygene, ])) %>%
    tibble::rownames_to_column(var = "rna_id") %>%
    tidyr::pivot_longer(names_to = "gene", values_to = "TPM", cols = !c(rna_id) )
  df <- df %>%
    dplyr::mutate( TPM = as.numeric(TPM)) %>%
    dplyr::left_join(meta, by = "rna_id") %>%
    filter(tissue %in% tissues)
  
  if(return_data == TRUE){ return(df)}
  if( y == "TPM" & log ==TRUE){
    df$TPM <- log2(df$TPM + 0.01)  
  }
  
  ## MAKE PLOT
  plot <- df %>%
    ggplot(aes_(x = as.name(x), y = as.name(y) )) + 
    facet_wrap( as.formula(paste("~", facet) ) ) + 
    theme_jh() +
    labs(x = xlab, y = ylab, colour = collab)
  
  
  
  # if categorical variable on X axis 
  #specify whether to jitter points - only if x is categorical and y is continuous
  # give option of boxplots
  if( x %in% c("tissue","seq_platform","mutations","disease","disease_c9","site_of_motor_onset","sex", "C9orf72 status")  ){
    plot <- plot + geom_jitter(width = 0.25, aes_(colour = as.name(colourby) ), height = 0)
    if(boxplot == TRUE){
      plot <- plot + geom_boxplot(aes_(group = as.name(x) ), fill = NA, outlier.color = NA)
    }
    if(stats == TRUE){
      plot <- plot + ggpubr::stat_compare_means()
    }
  }else{
    # IF CONTINUOUS
    plot <- plot + geom_point(aes_(colour = as.name(colourby) ))
    
    if(corline == TRUE){
      plot <- plot + geom_smooth(method = "lm", formula = 'y ~ x')
    }
    if(stats == TRUE){
      plot <- plot + ggpubr::stat_cor()
    }
  }
  # percent scales
  if( grepl("%", x) ){ 
    plot <- plot + scale_x_continuous(labels = scales::percent_format(accuracy = 1) ) 
  }
  if( grepl("%", y) ){ 
    plot <- plot + scale_y_continuous(labels = scales::percent_format(accuracy = 1) ) 
  }
  
  if( x == "TPM" | y == "TPM"){
    plot <- plot + labs(title = mygene) 
  }
  
  return(plot)
}

set_df <- tibble(
  test = c("CHIT1", "APOE", "APOC1")
)

set_plot <- function(sets = kelley,
         tissues = c("Cervical","Lumbar"), type = "de"){
  # get log_fc values from DE or duration results
  if(type == "de"){
    res <- de_res
  }
  if(type == "dur"){
    res <- dur_res
  }
  
  set_df <- map_df(sets, ~{tibble(gene = .x) }, .id = "set")
  
  df <- 
    map_df(res, ~{.x %>% filter(genename %in% set_df$gene)}, .id = "tissue") %>%
    filter(tissue %in% tissues) %>%
    left_join(set_df, by = c("genename" = "gene"))
  
  df %>%
    ggplot(aes(x = log_fc, y = set)) + 
    geom_jitter(size = 0.8, aes(colour = -log10(p_value)) ) + 
    facet_wrap(~tissue, ncol = 1) +
    xlim(-2,4) +
    geom_vline(xintercept = 0, linetype = 3) +
    theme_jh()
  
}

set_plot()

#gene_plot("MOBP", x = "% Oligodendrocytes", colourby = "disease_duration", tissues = "Cervical", boxplot = TRUE, corline = TRUE, stats = TRUE)


#setdiff(choice_df$short, names(metadata))
