# !/usr/bin/env Rscript
#                                 _       _
#     /\                         (_)     | |
#    /  \   _ __ ___   ___  _ __  _  ___ | |_
#   / /\ \ | '_ ` _ \ / _ \| '_ \| |/ _ \| __|
#  / ____ \| | | | | | (_) | | | | | (_) | |_
# /_/    \_\_| |_| |_|\___/|_| |_| |\___/ \__|
#                               _/ |
#                              |__/
# 30/05/2022
#
# Script Rarecurve
set.seed(1)
# Set input and output -----------------------------------------------------------
# Detect R or Rstudio
  se <- Sys.getenv()
  if (length(se[grepl("rstudio",se,ignore.case=TRUE)]) > 0 ) { inputmode <- TRUE }
  if (length(se[grepl("rstudio",se,ignore.case=TRUE)]) == 0 ) { inputmode <- FALSE }
#
# Input argument if using Rstudio
  if (inputmode == TRUE) {
    result <- "V4-unified-correct-paired-out-compo"
    region <- "V4"
  }
#
# Input argument if using R
  args = commandArgs(trailingOnly=TRUE)
  if ( inputmode == FALSE ) {
    result <- args[1]
    region <- args[2]
  }
#
# Import package and palette -----------------------------------------------------------
  pkg <- c("ggplot2","dplyr","tidyr","cowplot","FactoMineR","factoextra","reshape2","varhandle","ggrepel","ggpubr","ggsci","scales","hrbrthemes","GUniFrac","svglite","stringr","paletteer","elementalist","SciViews","vegan","parallel")
  lapply(pkg, require, character.only = TRUE)
  #palette <- c("#AD002ACC","#EEA236CC","#00468BCC","#0099B4CC","#D62768CC","#ADB6B6CC","#42B540CC","#357EBDCC","#1B1919CC","#ED0000CC","#FF7F0ECC","#925E9FCC","#D43F3ACC","#9632B8CC","#46B8DACC","#5CB85CCC","#EFC008CC")
  palet_monimixo <- paletteer_d("ggthemes::Nuriel_Stone",n=2,direction=-1)
# Set directory, create file result -----------------------------------------------------------
  if (inputmode == TRUE) {
    current <- dirname(rstudioapi::getSourceEditorContext()$path)
    setwd(current)
    result <- paste(current,"../dataDADA2/result",result, sep = "/")
  }
  if (inputmode == FALSE) {
    result <- paste("../dataDADA2/result",result, sep = "/")
  }
## Set working directory
  if (dir.exists(result) == FALSE) { dir.create(result,recursive = TRUE) }
  if (dir.exists(result) == TRUE) { setwd(result) }
## Create result files
  system("mkdir Rarecurve")
  system("mkdir Diversity")
#
# Theme unique Dark perso -------------------------------------------------------
  theme_unique_dark <- function (base_size = 12, base_family = "") {
    ret <- (theme_bw(base_size = base_size, base_family = base_family) +
              theme(text = element_text(colour = "black"),
                    title = element_text(color = "black", face = "bold"),
                    axis.ticks = element_blank(),
                    line = element_line(color = "black"),
                    rect = element_rect(fill = "white", color = "black"),
                    axis.title = element_text(color = "black", face = "bold"),
                    axis.text.y = element_blank(),
                    axis.text.x = element_blank(),
                    #axis.text.x = element_text(color = "black", size = 8, vjust = 2),
                    axis.line = element_line(color = "#969696", linetype = 1),
                    legend.background = element_rect(fill = NULL, color = NULL),
                    legend.position = "bottom",
                    legend.key = element_rect(fill = NA, color = NA, linetype = 0),
                    legend.text = element_text(size = 8),
                    legend.title = element_text(size = 8, face="bold"),
                    strip.background = element_rect(fill=NA,colour=NA,size=NA,linetype = NULL),
                    strip.text = element_text(color="black",face="bold",vjust=.5,hjust=.5),
                    panel.background = element_rect(fill = "white", color = NULL),
                    panel.border = element_blank(),
                    panel.grid = element_line(color = "#252525"),
                    panel.grid.major = element_line(color = "white"),
                    panel.grid.minor = element_line(color = "white"),
                    plot.background = element_rect(fill = "white", colour = "white", linetype = 0)))
    ret
  } 
  theme_unique_darkbis <- function (base_size = 12, base_family = "") {
    ret <- (theme_bw(base_size = base_size, base_family = base_family) +
              theme(text = element_text(colour = "black"),
                    title = element_text(color = "black", face = "bold", vjust = 0.5, hjust = 0.5),
                    axis.ticks = element_blank(),
                    line = element_line(color = "black"),
                    rect = element_rect(fill = "white", color = "black"),
                    axis.title = element_text(color = "black", face = "bold"),
                    axis.text.y = element_blank(),
                    axis.text.x = element_blank(),
                    axis.line = element_line(color = "#969696", linetype = 1),
                    legend.background = element_rect(fill = NULL, color = NULL),
                    legend.position = "none",
                    legend.key = element_rect(fill = NA, color = NA, linetype = 0),
                    legend.text = element_text(size = 8),
                    legend.title = element_text(size = 8, face="bold"),
                    strip.background = element_rect(fill=NA,colour=NA,size=NA,linetype = NULL),
                    strip.text = element_text(color="black",face="bold",vjust=.5,hjust=.5),
                    panel.background = element_rect(fill = "white", color = NULL),
                    panel.border = element_blank(),
                    panel.grid = element_line(color = "#252525"),
                    panel.grid.major = element_line(color = "white"),
                    panel.grid.minor = element_line(color = "white"),
                    plot.background = element_rect(fill = "white", colour = "white", linetype = 0)))
    ret
  } 
  theme_unique_art <- function (base_size = 12, base_family = "") {
    ret <- (theme_minimal(base_size = base_size, base_family = base_family) +
              theme(text = element_text(colour = "black"),
                    title = element_text(color = "black", face = "bold", vjust = 0.5, hjust = 0.5),
                    axis.ticks = element_blank(),
                    line = element_line(color = "black"),
                    rect = element_rect(fill = "white", color = "black"),
                    axis.title = element_text(color = "black", face = "bold"),
                    #axis.text.y = element_blank(),
                    #axis.text.x = element_blank(),
                    #axis.line = element_line(color = "#969696", linetype = 1),
                    legend.background = element_blank(),
                    legend.key = element_rect(fill = NA, color = NA, linetype = 0),
                    legend.text = element_text(size = 8),
                    legend.title = element_text(size = 8, face="bold"),
                    legend.position = "none",
                    panel.background = element_rect_round(radius = unit(0.5, "cm"),fill = "white", color = NULL,size = 1),
                    #strip.background = element_rect(fill=NA,colour=NA,size=NA,linetype = NULL),
                    strip.text = element_text(color="black",face="bold",vjust=.5,hjust=.5),
                    #panel.background = element_rect(fill = "white", color = NULL),
                    panel.border = element_blank(),
                    panel.grid = element_line(color = "#252525"),
                    panel.grid.major = element_line(color = "grey86"),
                    panel.grid.minor = element_line(color = "grey86"),
                    plot.background = element_rect(fill = "white", colour = "white", linetype = 0)))
    ret
  } 
#
# quickRareCurve fonction -------------------------------------------------
  quickRareCurve <- function (x, step = 1, sample, xlab = "Sequences",
                              ylab = "ASVs", label = TRUE, col, lty, max.cores = T, nCores = 1, ...)
  {
    require(parallel)
    x <- as.matrix(x)
    if (!identical(all.equal(x, round(x)), TRUE))
      stop("function accepts only integers (counts)")
    if (missing(col))
      col <- par("col")
    if (missing(lty))
      lty <- par("lty")
    tot <- rowSums(x) # calculates library sizes
    S <- specnumber(x) # calculates n species for each sample
    if (any(S <= 0)) {
      message("empty rows removed")
      x <- x[S > 0, , drop = FALSE]
      tot <- tot[S > 0]
      S <- S[S > 0]
    } # removes any empty rows
    nr <- nrow(x) # number of samples
    col <- rep(col, length.out = nr)
    lty <- rep(lty, length.out = nr)
    # parallel mclapply
    # set number of cores
    mc <- getOption("mc.cores", ifelse(max.cores, detectCores()/2, nCores))
    message(paste("Using ", mc, " cores"))
    out <- mclapply(seq_len(nr), mc.cores = mc, function(i) {
      n <- seq(1, tot[i], by = step)
      if (n[length(n)] != tot[i])
        n <- c(n, tot[i])
      drop(rarefy(x[i, ], n))
    })
    Nmax <- sapply(out, function(x) max(attr(x, "Subsample")))
    Smax <- sapply(out, max)
    plot(c(1, max(Nmax)), c(1, max(Smax)), xlab = xlab, ylab = ylab,
         type = "n", ...)
    if (!missing(sample)) {
      abline(v = sample)
      rare <- sapply(out, function(z) approx(x = attr(z, "Subsample"),
                                             y = z, xout = sample, rule = 1)$y)
      abline(h = rare, lwd = 0.5)
    }
    for (ln in seq_along(out)) {
      N <- attr(out[[ln]], "Subsample")
      lines(N, out[[ln]], col = col[ln], lty = lty[ln], ...)
    }
    if (label) {
      ordilabel(cbind(tot, S), labels = rownames(x), ...)
    }
    invisible(out)
  }
#
# Input  and prepare Tables ---------------------------------------------------------
## Raw
  seq_mat <- read.csv(file = paste("ASV-Table/Seq-mat.csv",sep = "/"), sep = "\t")
  amplicon <- c(grep(pattern = "OSTA", colnames(seq_mat), value = TRUE))
  seq_mat <- seq_mat %>% select(all_of(amplicon))
  seq_mat[,all_of(amplicon)] <- lapply(seq_mat[,all_of(amplicon)], as.numeric)
## Rare
  seq_mat_rare <- read.csv(file = paste("ASV-Table/Seq-mat-rare.csv",sep = "/"), sep = "\t")
  amplicon <- c(grep(pattern = "OSTA", colnames(seq_mat_rare), value = TRUE))
  seq_mat_rare <- seq_mat_rare %>% select(all_of(amplicon))
  seq_mat_rare[,all_of(amplicon)] <- lapply(seq_mat_rare[,all_of(amplicon)], as.numeric)
#
  # Prepare data inf --------------------------------------------------------
  infdataini <- read.table(file = "../../../rawdata/metadata_metaB", sep = "_", header = FALSE, col.names = c("Id","Condition","Period","Region","Replicat","empty"))
  infdataini$OSTA <- "OSTA"
  infdataini$Technologie  <- "Hiseq"
  infdataini$Amplicon <- paste(infdataini$Id,infdataini$OSTA, sep = "")
  infdataini <- infdataini %>% select(-"Id",-"OSTA")
  ## Prepare sample_df
  infdataini <- separate(infdataini, Condition, c("empty","ADN","Cycle","Zone","Fraction"),sep = "")
  infdataini <- infdataini %>% select(-"empty")
  infdataini<-as.data.frame(infdataini)
  col <- c("Amplicon","Technologie","Region","ADN","Cycle","Zone","Fraction","Period","Replicat")
  samples_df <- infdataini[,col]
  samples_df$Condition <- paste(samples_df$ADN,samples_df$Cycle,samples_df$Zone,samples_df$Fraction, sep = "")
  for (i in row.names(samples_df)) {
    if (samples_df[i,"Cycle"] == "J") { samples_df[i,"Cycle"] <- "Day"}
    if (samples_df[i,"Cycle"] == "N") { samples_df[i,"Cycle"] <- "Night"}
    if (samples_df[i,"Zone"] == "O") { samples_df[i,"Zone"] <- "Mixolimnion"}
    if (samples_df[i,"Zone"] == "A") { samples_df[i,"Zone"] <- "Monimolimnion"}
    if (samples_df[i,"Fraction"] == "G") { samples_df[i,"Fraction"] <- "Large"}
    if (samples_df[i,"Fraction"] == "P") { samples_df[i,"Fraction"] <- "Small"}
    if (samples_df[i,"Period"] == "4") { samples_df[i,"Period"] <- "04"}
    if (samples_df[i,"Period"] == "6") { samples_df[i,"Period"] <- "06"}
    if (samples_df[i,"Period"] == "9") { samples_df[i,"Period"] <- "09"}
    if (samples_df[i,"Period"] == "11") { samples_df[i,"Period"] <- "11"}
  }
  #
# Select V4 or V9 ---------------------------------------------------------
  samples_df<-filter(samples_df, Region == region)
# Create sort Condition pattern -------------------------------------------
## Cycle
  CycleDay <- samples_df %>% filter(Cycle == "Day") %>% filter(Replicat == 1)
  CycleDay <- CycleDay$Amplicon
  CycleNight <- samples_df %>% filter(Cycle == "Night") %>% filter(Replicat == 1)
  CycleNight <- CycleNight$Amplicon
## Zone
  ZoneMixolimnion <- samples_df %>% filter(`Zone` == "Mixolimnion") %>% filter(Replicat == 1)
  ZoneMixolimnion <- ZoneMixolimnion$Amplicon
  ZoneMonimolimnion <- samples_df %>% filter(`Zone` == "Monimolimnion") %>% filter(Replicat == 1)
  ZoneMonimolimnion <- ZoneMonimolimnion$Amplicon
## Fraction
  FractionSmall <- samples_df %>% filter(`Fraction` == "Small") %>% filter(Replicat == 1)
  FractionSmall <- FractionSmall$Amplicon
  FractionLarge <- samples_df %>% filter(`Fraction` == "Large") %>% filter(Replicat == 1)
  FractionLarge <- FractionLarge$Amplicon
## FractionXZone
  FractionSmallXMixolimnion <- samples_df %>% filter(`Fraction` == "Small") %>% filter(`Zone` == "Mixolimnion") %>% filter(Replicat == 1)
  FractionSmallXMixolimnion <- FractionSmallXMixolimnion$Amplicon
  FractionLargeXMixolimnion <- samples_df %>% filter(`Fraction` == "Large") %>% filter(`Zone` == "Mixolimnion") %>% filter(Replicat == 1)
  FractionLargeXMixolimnion <- FractionLargeXMixolimnion$Amplicon
  FractionSmallXMonimolimnion <- samples_df %>% filter(`Fraction` == "Small") %>% filter(`Zone` == "Monimolimnion") %>% filter(Replicat == 1)
  FractionSmallXMonimolimnion <- FractionSmallXMonimolimnion$Amplicon
  FractionLargeXMonimolimnion <- samples_df %>% filter(`Fraction` == "Large") %>% filter(`Zone` == "Monimolimnion") %>% filter(Replicat == 1)
  FractionLargeXMonimolimnion <- FractionLargeXMonimolimnion$Amplicon
## Period
  Period04 <- samples_df %>% filter(Period == "04") %>% filter(Replicat == 1)
  Period04 <- Period04$Amplicon
  Period06 <- samples_df %>% filter(Period == "06") %>% filter(Replicat == 1)
  Period06 <- Period06$Amplicon
  Period09 <- samples_df %>% filter(Period == "09") %>% filter(Replicat == 1)
  Period09 <- Period09$Amplicon
  Period11 <- samples_df %>% filter(Period == "11") %>% filter(Replicat == 1)
  Period11 <- Period11$Amplicon
#
# Rarecurve ---------------------------------------------------------------
  # Raw --------------------------------------------------------------------
  ## rarecurve-raw
    rawcurve <- t(seq_mat)
    pdf("Rarecurve/Rarecurve-Raw.pdf",width = 5.00,height = 5.00)
    curve <- quickRareCurve(rawcurve, col = "black", cex = 0.6, label = FALSE)
    dev.off()
  ## Total FractionXZone Raw
  ### SmallXMixo
    total_SmallXMixolimnion_raw <- as.data.frame(rawcurve) %>% filter(row.names(rawcurve) %in% all_of(FractionSmallXMixolimnion))
    total_SmallXMixolimnion_raw <- t(colSums(total_SmallXMixolimnion_raw)) ; row.names(total_SmallXMixolimnion_raw) <- "Small"
  ### LargeXMixo
    total_LargeXMixolimnion_raw <- as.data.frame(rawcurve) %>% filter(row.names(rawcurve) %in% all_of(FractionLargeXMixolimnion))
    total_LargeXMixolimnion_raw <- t(colSums(total_LargeXMixolimnion_raw)) ; row.names(total_LargeXMixolimnion_raw) <- "Large"
  ### SmallXMoni
    total_SmallXMonimolimnion_raw <- as.data.frame(rawcurve) %>% filter %>% filter(row.names(rawcurve) %in% all_of(FractionSmallXMonimolimnion))
    total_SmallXMonimolimnion_raw <- t(colSums(total_SmallXMonimolimnion_raw)) ; row.names(total_SmallXMonimolimnion_raw) <- "Small"
  ### LargeXMoni
    total_LargeXMonimolimnion_raw <- as.data.frame(rawcurve) %>% filter(row.names(rawcurve) %in% all_of(FractionLargeXMonimolimnion))
    total_LargeXMonimolimnion_raw <- t(colSums(total_LargeXMonimolimnion_raw)) ; row.names(total_LargeXMonimolimnion_raw) <- "Large"
  ### plot
    total_FractionXZone_raw <- rbind(total_SmallXMixolimnion_raw,total_LargeXMixolimnion_raw,total_SmallXMonimolimnion_raw,total_LargeXMonimolimnion_raw)
    Color <- c("#6FB899FF","#6FB899FF","#8175AAFF","#8175AAFF")
    pdf("Rarecurve/Rarecurve-FractionXZoneSum-Rawdata.pdf", width = 5.00, height = 5.00)
    #svglite("Rarecurve/Rarecurve-FractionXZoneSum-Rawdata.svg", width = 5.00, height = 5.00)
    curve_d <- quickRareCurve(total_FractionXZone_raw, col = Color, cex = 0.6, ylab = "ASVs", xlab = "Illumina reads")
    legend("topleft", inset=.02, title="Fraction",c("Monimolimnion","Mixolimnion"), fill = c("#8175AAFF","#6FB899FF"), horiz=TRUE, cex=0.6, bg="white")
    dev.off()
#
# Diversity Raw ---------------------------------------------------------------
  rawcurve <- t(seq_mat)
  # Shannon -----------------------------------------------------------------
    Div_Table <- samples_df
    rownames(Div_Table) <- Div_Table[,"Amplicon"]
    Divraw <- as.data.frame(diversity(rawcurve, index = "shannon"))
    colnames(Divraw) <- "Shannon"
    for (i in rownames(Divraw)) { Div_Table[i,"Shannon"] <- Divraw[i,"Shannon"] }
    Div_Table <- Div_Table %>% filter(Replicat!=2)
  ## Period
    my_comp <- list(c("04","06"),c("04","09"),c("04","11"),c("06","09"),c("06","11"),c("09","11"))
    adiv <- ggplot(Div_Table, aes(x = Period, y = Shannon)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
      stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
      theme(axis.title = element_text(face="bold", size=12), 
            axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
            title = element_text(face="bold", size=14),
            legend.title = element_text(face="bold"),
            legend.position = "right",
            legend.text = element_text(size=10)) +
      labs(x="Period",y="Diversity") #+ guides(color = FALSE)
    print(adiv)
    legend <- get_legend(adiv)
    adiv<- adiv + theme(legend.position="none")
  ## Cycle
    my_comp <- list(c("Day","Night"))
    bdiv <- ggplot(Div_Table, aes(x = Cycle, y = Shannon)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
      stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
      theme(axis.title = element_text(face="bold", size=12), 
            axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
            title = element_text(face="bold", size=14),
            legend.title = element_text(face="bold"),
            legend.position = "right",
            legend.text = element_text(size=10)) +
      labs(x="Cycle",y="Diversity") #+ guides(color = FALSE)
    print(bdiv)
    bdiv<- bdiv + theme(legend.position="none")
  ## Fraction
    my_comp <- list(c("Small","Large"))
    cdiv <- ggplot(Div_Table, aes(x = Fraction, y = Shannon)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
      stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
      theme(axis.title = element_text(face="bold", size=12), 
            axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
            title = element_text(face="bold", size=14),
            legend.title = element_text(face="bold"),
            legend.position = "right",
            legend.text = element_text(size=10)) +
      labs(x="Fraction",y="Diversity") #+ guides(color = FALSE)
    print(cdiv)
    cdiv<- cdiv + theme(legend.position="none")
  ## Zone
    my_comp <- list(c("Mixolimnion","Monimolimnion"))
    ddiv <- ggplot(Div_Table, aes(x = Zone, y = Shannon)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
      stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
      theme(axis.title = element_text(face="bold", size=12), 
            axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
            title = element_text(face="bold", size=14),
            legend.title = element_text(face="bold"),
            legend.position = "right",
            legend.text = element_text(size=10)) +
      labs(x="Zone",y="Diversity") #+ guides(color = FALSE)
    print(ddiv)
    ddiv<- ddiv + theme(legend.position="none")
  ## Coplot
    fdiv<-ggplot() + theme(panel.background = element_rect(fill="white",colour = "white", size = 0.5, linetype = "solid"))
    pdf("Diversity/Analyse-Shannon-Raw.pdf",width = 16.00,height = 9.00)
    Diversity_plot <- plot_grid(adiv,bdiv, fdiv, fdiv,fdiv,legend,cdiv,ddiv, labels=c("A","B", "", "", "", "", "C","D", ""), ncol = 3, nrow = 3, rel_widths = c(3,3,1.5), rel_heights = c(5,0.1,5))
    print(Diversity_plot)
    dev.off()
  ## ZoneXFraction
    my_comp <- list(c("Mixolimnion","Monimolimnion"))
    xidiv <- ggplot(Div_Table, aes(x = Zone, y = Shannon)) + geom_violin(aes(fill = Zone),trim = FALSE,width=0.5) + geom_boxplot(width=0.075) + geom_jitter(position=position_jitter(0.05)) +
      stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
      facet_grid(~Fraction) + scale_fill_manual(values = palet_monimixo) + theme_unique_art() +
      labs(x="Zone",y="Shannon index") #+ guides(color = FALSE)
    print(xidiv)
    ggsave("Diversity/ZoneXfraction_Shannon_Raw.svg", device = "svg", width = 5, height = 5)
  # Simpson -----------------------------------------------------------------
    Divraw <- as.data.frame(diversity(rawcurve, index = "simpson"))
    colnames(Divraw) <- "Simpson"
    for (i in rownames(Divraw)) { Div_Table[i,"Simpson"] <- Divraw[i,"Simpson"] }
  ## Period
    my_comp <- list(c("04","06"),c("04","09"),c("04","11"),c("06","09"),c("06","11"),c("09","11"))
    adiv <- ggplot(Div_Table, aes(x = Period, y = Simpson)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
      stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
      theme(axis.title = element_text(face="bold", size=12), 
            axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
            title = element_text(face="bold", size=14),
            legend.title = element_text(face="bold"),
            legend.position = "right",
            legend.text = element_text(size=10)) +
      labs(x="Period",y="Diversity") #+ guides(color = FALSE)
    print(adiv)
    legend <- get_legend(adiv)
    adiv<- adiv + theme(legend.position="none")
  ## Cycle
    my_comp <- list(c("Day","Night"))
    bdiv <- ggplot(Div_Table, aes(x = Cycle, y = Simpson)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
      stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
      theme(axis.title = element_text(face="bold", size=12), 
            axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
            title = element_text(face="bold", size=14),
            legend.title = element_text(face="bold"),
            legend.position = "right",
            legend.text = element_text(size=10)) +
      labs(x="Cycle",y="Diversity") #+ guides(color = FALSE)
    print(bdiv)
    bdiv<- bdiv + theme(legend.position="none")
  ## Fraction
    my_comp <- list(c("Small","Large"))
    cdiv <- ggplot(Div_Table, aes(x = Fraction, y = Simpson)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
      stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
      theme(axis.title = element_text(face="bold", size=12), 
            axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
            title = element_text(face="bold", size=14),
            legend.title = element_text(face="bold"),
            legend.position = "right",
            legend.text = element_text(size=10)) +
      labs(x="Fraction",y="Diversity") #+ guides(color = FALSE)
    print(cdiv)
    cdiv<- cdiv + theme(legend.position="none")
  ## Zone
    my_comp <- list(c("Mixolimnion","Monimolimnion"))
    ddiv <- ggplot(Div_Table, aes(x = Zone, y = Simpson)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
      stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
      theme(axis.title = element_text(face="bold", size=12), 
            axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
            title = element_text(face="bold", size=14),
            legend.title = element_text(face="bold"),
            legend.position = "right",
            legend.text = element_text(size=10)) +
      labs(x="Zone",y="Diversity") #+ guides(color = FALSE)
    print(ddiv)
    ddiv<- ddiv + theme(legend.position="none")
  ## Coplot
    fdiv<-ggplot() + theme(panel.background = element_rect(fill="white",colour = "white", size = 0.5, linetype = "solid"))
    pdf("Diversity/Analyse-Simpson-Raw.pdf",width = 16.00,height = 9.00)
    Diversity_plot <- plot_grid(adiv,bdiv, fdiv, fdiv,fdiv,legend,cdiv,ddiv, labels=c("A","B", "", "", "", "", "C","D", ""), ncol = 3, nrow = 3, rel_widths = c(3,3,1.5), rel_heights = c(5,0.1,5))
    print(Diversity_plot)
    dev.off()
  ## ZoneXFraction
    my_comp <- list(c("Mixolimnion","Monimolimnion"))
    xjdiv <- ggplot(Div_Table, aes(x = Zone, y = Simpson)) + geom_violin(aes(fill = Zone),trim = FALSE,width=0.5) + geom_boxplot(width=0.075) + geom_jitter(position=position_jitter(0.05)) +
      stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
      facet_grid(~Fraction) + scale_fill_manual(values = palet_monimixo) + theme_unique_art() +
      labs(x="Zones",y="Simpson index") #+ guides(color = FALSE)
    print(xjdiv)
    ggsave("Diversity/ZoneXfraction_Simpson_Raw.svg", device = "svg", width = 5, height = 5)
#
  # Chao1 -------------------------------------------------------------------
    chao_pool <- as.data.frame(estimateR(rawcurve))
    for (i in colnames(chao_pool)) { Div_Table[i,"chao"] <- chao_pool["S.chao1",i] }
    for (i in colnames(chao_pool)) { Div_Table[i,"chao_se"] <- chao_pool["se.chao1",i] }
  ## Period
    my_comp <- list(c("04","06"),c("04","09"),c("04","11"),c("06","09"),c("06","11"),c("09","11"))
    adiv <- ggplot(Div_Table, aes(x = Period, y = chao)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
      geom_errorbar(aes(ymin=chao-chao_se,ymax=chao+chao_se), width=.1,position=position_dodge(0.05)) +
      stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
      theme(axis.title = element_text(face="bold", size=12), 
            axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
            title = element_text(face="bold", size=14),
            legend.title = element_text(face="bold"),
            legend.position = "right",
            legend.text = element_text(size=10)) +
      labs(x="Period",y="Diversity") #+ guides(color = FALSE)
    print(adiv)
    legend <- get_legend(adiv)
    adiv<- adiv + theme(legend.position="none")
  ## Cycle
    my_comp <- list(c("Day","Night"))
    bdiv <- ggplot(Div_Table, aes(x = Cycle, y = chao)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
      geom_errorbar(aes(ymin=chao-chao_se,ymax=chao+chao_se), width=.1,position=position_dodge(0.05)) +
      stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
      theme(axis.title = element_text(face="bold", size=12), 
            axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
            title = element_text(face="bold", size=14),
            legend.title = element_text(face="bold"),
            legend.position = "right",
            legend.text = element_text(size=10)) +
      labs(x="Cycle",y="Diversity") #+ guides(color = FALSE)
    print(bdiv)
    bdiv<- bdiv + theme(legend.position="none")
  ## Fraction
    my_comp <- list(c("Small","Large"))
    cdiv <- ggplot(Div_Table, aes(x = Fraction, y = chao)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
      geom_errorbar(aes(ymin=chao-chao_se,ymax=chao+chao_se), width=.1,position=position_dodge(0.05)) +
      stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
      theme(axis.title = element_text(face="bold", size=12), 
            axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
            title = element_text(face="bold", size=14),
            legend.title = element_text(face="bold"),
            legend.position = "right",
            legend.text = element_text(size=10)) +
      labs(x="Fraction",y="Diversity") #+ guides(color = FALSE)
    print(cdiv)
    cdiv<- cdiv + theme(legend.position="none")
  ## Zone
    my_comp <- list(c("Mixolimnion","Monimolimnion"))
    ddiv <- ggplot(Div_Table, aes(x = Zone, y = chao)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
      geom_errorbar(aes(ymin=chao-chao_se,ymax=chao+chao_se), width=.1,position=position_dodge(0.05)) + 
      stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
      theme(axis.title = element_text(face="bold", size=12), 
            axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
            title = element_text(face="bold", size=14),
            legend.title = element_text(face="bold"),
            legend.position = "right",
            legend.text = element_text(size=10)) +
      labs(x="Zone",y="Diversity") #+ guides(color = FALSE)
    print(ddiv)
    ddiv<- ddiv + theme(legend.position="none")
  ## Coplot
    fdiv<-ggplot() + theme(panel.background = element_rect(fill="white",colour = "white", size = 0.5, linetype = "solid"))
    pdf("Diversity/Analyse-Chao-Raw.pdf",width = 16.00,height = 9.00)
    Diversity_plot <- plot_grid(adiv,bdiv, fdiv, fdiv,fdiv,legend,cdiv,ddiv, labels=c("A","B", "", "", "", "", "C","D", ""), ncol = 3, nrow = 3, rel_widths = c(3,3,1.5), rel_heights = c(5,0.1,5))
    print(Diversity_plot)
    dev.off()
  ## ZoneXFraction
    my_comp <- list(c("Mixolimnion","Monimolimnion"))
    xkdiv <- ggplot(Div_Table, aes(x = Zone, y = chao)) + geom_violin(aes(fill = Zone),trim = FALSE,width=0.5) + geom_boxplot(width=0.075) + geom_jitter(position=position_jitter(0.05)) +
      stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
      facet_grid(~Fraction) + scale_fill_manual(values = palet_monimixo) + theme_unique_art() +
      labs(x="Zone",y="Chao index") #+ guides(color = FALSE)
    print(xkdiv)
    ggsave("Diversity/ZoneXfraction_chao_Raw.svg", device = "svg", width = 5, height = 5)
#
  # ACE -------------------------------------------------------------------
    for (i in colnames(chao_pool)) { Div_Table[i,"ACE"] <- chao_pool["S.ACE",i] }
    for (i in colnames(chao_pool)) { Div_Table[i,"ACE_se"] <- chao_pool["se.ACE",i] }
  ## Period
    my_comp <- list(c("04","06"),c("04","09"),c("04","11"),c("06","09"),c("06","11"),c("09","11"))
    adiv <- ggplot(Div_Table, aes(x = Period, y = ACE)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
      geom_errorbar(aes(ymin=ACE-ACE_se,ymax=ACE+ACE_se), width=.1,position=position_dodge(0.05)) +
      stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
      theme(axis.title = element_text(face="bold", size=12),
            axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5),
            title = element_text(face="bold", size=14),
            legend.title = element_text(face="bold"),
            legend.position = "right",
            legend.text = element_text(size=10)) +
      labs(x="Period",y="Diversity") #+ guides(color = FALSE)
    print(adiv)
    legend <- get_legend(adiv)
    adiv<- adiv + theme(legend.position="none")
  ## Cycle
    my_comp <- list(c("Day","Night"))
    bdiv <- ggplot(Div_Table, aes(x = Cycle, y = ACE)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
      geom_errorbar(aes(ymin=ACE-ACE_se,ymax=ACE+ACE_se), width=.1,position=position_dodge(0.05)) +
      stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
      theme(axis.title = element_text(face="bold", size=12), 
            axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
            title = element_text(face="bold", size=14),
            legend.title = element_text(face="bold"),
            legend.position = "right",
            legend.text = element_text(size=10)) +
      labs(x="Cycle",y="Diversity") #+ guides(color = FALSE)
    print(bdiv)
    bdiv<- bdiv + theme(legend.position="none")
  ## Fraction
    my_comp <- list(c("Small","Large"))
    cdiv <- ggplot(Div_Table, aes(x = Fraction, y = ACE)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
      geom_errorbar(aes(ymin=ACE-ACE_se,ymax=ACE+ACE_se), width=.1,position=position_dodge(0.05)) +
      stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
      theme(axis.title = element_text(face="bold", size=12), 
            axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
            title = element_text(face="bold", size=14),
            legend.title = element_text(face="bold"),
            legend.position = "right",
            legend.text = element_text(size=10)) +
      labs(x="Fraction",y="Diversity") #+ guides(color = FALSE)
    print(cdiv)
    cdiv<- cdiv + theme(legend.position="none")
  ## Zone
    my_comp <- list(c("Mixolimnion","Monimolimnion"))
    ddiv <- ggplot(Div_Table, aes(x = Zone, y = ACE)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
      geom_errorbar(aes(ymin=ACE-ACE_se,ymax=ACE+ACE_se), width=.1,position=position_dodge(0.05)) + 
      stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
      theme(axis.title = element_text(face="bold", size=12), 
            axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
            title = element_text(face="bold", size=14),
            legend.title = element_text(face="bold"),
            legend.position = "right",
            legend.text = element_text(size=10)) +
      labs(x="Zone",y="Diversity") #+ guides(color = FALSE)
    print(ddiv)
    ddiv<- ddiv + theme(legend.position="none")
  ## Coplot
    fdiv<-ggplot() + theme(panel.background = element_rect(fill="white",colour = "white", size = 0.5, linetype = "solid"))
    pdf("Diversity/Analyse-ACE-Raw.pdf",width = 16.00,height = 9.00)
    Diversity_plot <- plot_grid(adiv,bdiv, fdiv, fdiv,fdiv,legend,cdiv,ddiv, labels=c("A","B", "", "", "", "", "C","D", ""), ncol = 3, nrow = 3, rel_widths = c(3,3,1.5), rel_heights = c(5,0.1,5))
    print(Diversity_plot)
    dev.off()
  ## ZoneXFraction
    my_comp <- list(c("Mixolimnion","Monimolimnion"))
    xldiv <- ggplot(Div_Table, aes(x = Zone, y = ACE)) + geom_violin(aes(fill = Zone),trim = FALSE,width=0.5) + geom_boxplot(width=0.075) + geom_jitter(position=position_jitter(0.05)) +
      stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
      facet_grid(~Fraction) + scale_fill_manual(values = palet_monimixo) + theme_unique_art() +
      labs(x="Zone",y="ACE index") #+ guides(color = FALSE)
    print(xldiv)
    ggsave("Diversity/ZoneXfraction_ACE_Raw.svg", device = "svg", width = 5, height = 5)
#
  # Obs -------------------------------------------------------------------
    for (i in colnames(chao_pool)) { Div_Table[i,"Obs"] <- chao_pool["S.obs",i] }
  ## Period
    my_comp <- list(c("04","06"),c("04","09"),c("04","11"),c("06","09"),c("06","11"),c("09","11"))
    adiv <- ggplot(Div_Table, aes(x = Period, y = Obs)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
      stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
      theme(axis.title = element_text(face="bold", size=12), 
            axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
            title = element_text(face="bold", size=14),
            legend.title = element_text(face="bold"),
            legend.position = "right",
            legend.text = element_text(size=10)) +
      labs(x="Period",y="Diversity") #+ guides(color = FALSE)
    print(adiv)
    legend <- get_legend(adiv)
    adiv<- adiv + theme(legend.position="none")
  ## Cycle
    my_comp <- list(c("Day","Night"))
    bdiv <- ggplot(Div_Table, aes(x = Cycle, y = Obs)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
      stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
      theme(axis.title = element_text(face="bold", size=12), 
            axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
            title = element_text(face="bold", size=14),
            legend.title = element_text(face="bold"),
            legend.position = "right",
            legend.text = element_text(size=10)) +
      labs(x="Cycle",y="Diversity") #+ guides(color = FALSE)
    print(bdiv)
    bdiv<- bdiv + theme(legend.position="none")
  ## Fraction(Taille)
    my_comp <- list(c("Small","Large"))
    cdiv <- ggplot(Div_Table, aes(x = Fraction, y = Obs)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
      stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
      theme(axis.title = element_text(face="bold", size=12), 
            axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
            title = element_text(face="bold", size=14),
            legend.title = element_text(face="bold"),
            legend.position = "right",
            legend.text = element_text(size=10)) +
      labs(x="Fraction",y="Diversity") #+ guides(color = FALSE)
    print(cdiv)
    cdiv<- cdiv + theme(legend.position="none")
  ## Zone
    my_comp <- list(c("Mixolimnion","Monimolimnion"))
    ddiv <- ggplot(Div_Table, aes(x = Zone, y = Obs)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
      stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
      theme(axis.title = element_text(face="bold", size=12), 
            axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
            title = element_text(face="bold", size=14),
            legend.title = element_text(face="bold"),
            legend.position = "right",
            legend.text = element_text(size=10)) +
      labs(x="Zone",y="Diversity") #+ guides(color = FALSE)
    print(ddiv)
    ddiv<- ddiv + theme(legend.position="none")
  ## Coplot
    fdiv<-ggplot() + theme(panel.background = element_rect(fill="white",colour = "white", size = 0.5, linetype = "solid"))
    pdf("Diversity/Analyse-Obs-Raw.pdf",width = 16.00,height = 9.00)
    Diversity_plot <- plot_grid(adiv,bdiv, fdiv, fdiv,fdiv,legend,cdiv,ddiv, labels=c("A","B", "", "", "", "", "C","D", ""), ncol = 3, nrow = 3, rel_widths = c(3,3,1.5), rel_heights = c(5,0.1,5))
    print(Diversity_plot)
    dev.off()
  ## ZoneXFraction
    my_comp <- list(c("Mixolimnion","Monimolimnion"))
    xmdiv <- ggplot(Div_Table, aes(x = Zone, y = Obs)) + geom_violin(aes(fill = Zone),trim = FALSE,width=0.5) + geom_boxplot(width=0.075) + geom_jitter(position=position_jitter(0.05)) +
      stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
      facet_grid(~Fraction) + scale_fill_manual(values = palet_monimixo) + theme_unique_art() +
      labs(x="Zone",y="Observed richness") #+ guides(color = FALSE)
    print(xmdiv)
    ggsave("Diversity/ZoneXfraction_Obs_Raw.svg", device = "svg", width = 5, height = 5)
#
# Diversity Rare ---------------------------------------------------------------
  rarecurve <- t(seq_mat_rare)
#
  # Shannon -----------------------------------------------------------------
    Div_Table <- samples_df
    rownames(Div_Table) <- Div_Table[,"Amplicon"]
    Divrare <- as.data.frame(diversity(rarecurve, index = "shannon"))
    colnames(Divrare) <- "Shannon"
    for (i in rownames(Divrare)) { Div_Table[i,"Shannon"] <- Divrare[i,"Shannon"] }
    Div_Table <- Div_Table %>% filter(Replicat!=2)
    Div_Table$Pielou <- c(Div_Table$Shannon / ln(specnumber(rarecurve)))
  ## Period
    my_comp <- list(c("04","06"),c("04","09"),c("04","11"),c("06","09"),c("06","11"),c("09","11"))
    adiv <- ggplot(Div_Table, aes(x = Period, y = Shannon)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
      stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
      theme(axis.title = element_text(face="bold", size=12), 
            axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
            title = element_text(face="bold", size=14),
            legend.title = element_text(face="bold"),
            legend.position = "right",
            legend.text = element_text(size=10)) +
      labs(x="Period",y="Diversity") #+ guides(color = FALSE)
    print(adiv)
    legend <- get_legend(adiv)
    adiv<- adiv + theme(legend.position="none")
  ## Cycle
    my_comp <- list(c("Day","Night"))
    bdiv <- ggplot(Div_Table, aes(x = Cycle, y = Shannon)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
      stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
      theme(axis.title = element_text(face="bold", size=12), 
            axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
            title = element_text(face="bold", size=14),
            legend.title = element_text(face="bold"),
            legend.position = "right",
            legend.text = element_text(size=10)) +
      labs(x="Cycle",y="Diversity") #+ guides(color = FALSE)
    print(bdiv)
    bdiv<- bdiv + theme(legend.position="none")
  ## Fraction
    my_comp <- list(c("Small","Large"))
    cdiv <- ggplot(Div_Table, aes(x = Fraction, y = Shannon)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
      stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
      theme(axis.title = element_text(face="bold", size=12), 
            axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
            title = element_text(face="bold", size=14),
            legend.title = element_text(face="bold"),
            legend.position = "right",
            legend.text = element_text(size=10)) +
      labs(x="Fraction",y="Diversity") #+ guides(color = FALSE)
    print(cdiv)
    cdiv<- cdiv + theme(legend.position="none")
  ## Zone
    my_comp <- list(c("Mixolimnion","Monimolimnion"))
    ddiv <- ggplot(Div_Table, aes(x = Zone, y = Shannon)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
      stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
      theme(axis.title = element_text(face="bold", size=12), 
            axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
            title = element_text(face="bold", size=14),
            legend.title = element_text(face="bold"),
            legend.position = "right",
            legend.text = element_text(size=10)) +
      labs(x="Zone",y="Diversity") #+ guides(color = FALSE)
    print(ddiv)
    ddiv<- ddiv + theme(legend.position="none")
  ## Coplot
    fdiv<-ggplot() + theme(panel.background = element_rect(fill="white",colour = "white", size = 0.5, linetype = "solid"))
    pdf("Diversity/Analyse-Shannon-Rare.pdf",width = 16.00,height = 9.00)
    Diversity_plot <- plot_grid(adiv,bdiv, fdiv, fdiv,fdiv,legend,cdiv,ddiv, labels=c("A","B", "", "", "", "", "C","D", ""), ncol = 3, nrow = 3, rel_widths = c(3,3,1.5), rel_heights = c(5,0.1,5))
    print(Diversity_plot)
    dev.off()
  ## ZoneXFraction
    my_comp <- list(c("Mixolimnion","Monimolimnion"))
    yidiv <- ggplot(Div_Table, aes(x = Zone, y = Shannon)) + geom_violin(aes(fill = Zone),trim = FALSE,width=0.5) + geom_boxplot(width=0.075) + geom_jitter(position=position_jitter(0.05)) +
      stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
      stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
      facet_grid(~Fraction) + scale_fill_manual(values = palet_monimixo) + theme_unique_art() +
      labs(x="Zones",y="Shannon index") #+ guides(color = FALSE)
    print(yidiv)
    ggsave("Diversity/ZoneXfraction_Shannon_Rare.svg", device = "svg", width = 5, height = 5)
  ## ZoneXFraction
    my_comp <- list(c("Mixolimnion","Monimolimnion"))
    zidiv <- ggplot(Div_Table, aes(x = Zone, y = Pielou)) + geom_violin(aes(fill = Zone),trim = FALSE,width=0.5) + geom_boxplot(width=0.075) + geom_jitter(position=position_jitter(0.05)) +
      stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
      stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
      facet_grid(~Fraction) + scale_fill_manual(values = palet_monimixo) + theme_unique_art() +
      labs(x="Zones",y="Pielou index") #+ guides(color = FALSE)
    print(zidiv)
    ggsave("Diversity/ZoneXfraction_Pielou_Rare.svg", device = "svg", width = 5, height = 5)
  ## ZoneXPeriods
    my_comp <- list(c("04","06"),c("04","09"),c("04","11"),c("06","09"),c("06","11"),c("09","11"))
    palet_periodcomp <- c("#7375b5FF","#8ca252FF","#f6ae34FF","#ad494aFF")
    aidiv <- ggplot(Div_Table, aes(x = Period, y = Pielou)) + geom_violin(aes(fill = Period),trim = FALSE,width=0.5,alpha=0.9) + geom_boxplot(width=0.075) + geom_jitter(position=position_jitter(0.05)) +
      stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
      stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
      facet_grid(~Zone) + scale_fill_manual(values = palet_periodcomp) + theme_unique_art() +
      labs(x="Periods",y="Pielou index") #+ guides(color = FALSE)
    print(aidiv)
    ggsave("Diversity/ZoneXPeriod_Pielou_Rare.svg", device = "svg", width = 5, height = 5)
  ## ZoneXPeriods
    my_comp <- list(c("04","06"),c("04","09"),c("04","11"),c("06","09"),c("06","11"),c("09","11"))
    palet_periodcomp <- c("#7375b5FF","#8ca252FF","#f6ae34FF","#ad494aFF")
    bidiv <- ggplot(Div_Table, aes(x = Period, y = Shannon)) + geom_violin(aes(fill = Period),trim = FALSE,width=0.5,alpha=0.9) + geom_boxplot(width=0.075) + geom_jitter(position=position_jitter(0.05)) +
      stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
      stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
      facet_grid(~Zone) + scale_fill_manual(values = palet_periodcomp) + theme_unique_art() +
      labs(x="Periods",y="Shannon index") #+ guides(color = FALSE)
    print(bidiv)
    ggsave("Diversity/ZoneXPeriod_Shannon_Rare.svg", device = "svg", width = 5, height = 5)
#  
  # Simpson -----------------------------------------------------------------
    Divrare <- as.data.frame(diversity(rarecurve, index = "simpson"))
    colnames(Divrare) <- "Simpson"
    for (i in rownames(Divrare)) { Div_Table[i,"Simpson"] <- Divrare[i,"Simpson"] }
  ## Period
    my_comp <- list(c("04","06"),c("04","09"),c("04","11"),c("06","09"),c("06","11"),c("09","11"))
    adiv <- ggplot(Div_Table, aes(x = Period, y = Simpson)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
      stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
      theme(axis.title = element_text(face="bold", size=12), 
            axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
            title = element_text(face="bold", size=14),
            legend.title = element_text(face="bold"),
            legend.position = "right",
            legend.text = element_text(size=10)) +
      labs(x="Period",y="Diversity") #+ guides(color = FALSE)
    print(adiv)
    legend <- get_legend(adiv)
    adiv<- adiv + theme(legend.position="none")
  ## Cycle
    my_comp <- list(c("Day","Night"))
    bdiv <- ggplot(Div_Table, aes(x = Cycle, y = Simpson)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
      stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
      theme(axis.title = element_text(face="bold", size=12), 
            axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
            title = element_text(face="bold", size=14),
            legend.title = element_text(face="bold"),
            legend.position = "right",
            legend.text = element_text(size=10)) +
      labs(x="Cycle",y="Diversity") #+ guides(color = FALSE)
    print(bdiv)
    bdiv<- bdiv + theme(legend.position="none")
  ## Fraction
    my_comp <- list(c("Small","Large"))
    cdiv <- ggplot(Div_Table, aes(x = Fraction, y = Simpson)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
      stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
      theme(axis.title = element_text(face="bold", size=12), 
            axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
            title = element_text(face="bold", size=14),
            legend.title = element_text(face="bold"),
            legend.position = "right",
            legend.text = element_text(size=10)) +
      labs(x="Fraction",y="Diversity") #+ guides(color = FALSE)
    print(cdiv)
    cdiv<- cdiv + theme(legend.position="none")
  ## Zone
    my_comp <- list(c("Mixolimnion","Monimolimnion"))
    ddiv <- ggplot(Div_Table, aes(x = Zone, y = Simpson)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
      stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
      theme(axis.title = element_text(face="bold", size=12), 
            axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
            title = element_text(face="bold", size=14),
            legend.title = element_text(face="bold"),
            legend.position = "right",
            legend.text = element_text(size=10)) +
      labs(x="Zone",y="Diversity") #+ guides(color = FALSE)
    print(ddiv)
    ddiv<- ddiv + theme(legend.position="none")
  ## Coplot
    fdiv<-ggplot() + theme(panel.background = element_rect(fill="white",colour = "white", size = 0.5, linetype = "solid"))
    pdf("Diversity/Analyse-Simpson-Rare.pdf",width = 16.00,height = 9.00)
    Diversity_plot <- plot_grid(adiv,bdiv, fdiv, fdiv,fdiv,legend,cdiv,ddiv, labels=c("A","B", "", "", "", "", "C","D", ""), ncol = 3, nrow = 3, rel_widths = c(3,3,1.5), rel_heights = c(5,0.1,5))
    print(Diversity_plot)
    dev.off()
  ## ZoneXFraction
    my_comp <- list(c("Mixolimnion","Monimolimnion"))
    yjdiv <- ggplot(Div_Table, aes(x = Zone, y = Simpson)) + geom_violin(aes(fill = Zone),trim = FALSE,width=0.5) + geom_boxplot(width=0.075)+ geom_jitter(position=position_jitter(0.05)) +
      stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
      facet_grid(~Fraction) + scale_fill_manual(values = palet_monimixo) + theme_unique_art() +
      labs(x="Zone",y="Simpson index") #+ guides(color = FALSE)
    print(yjdiv)
    ggsave("Diversity/ZoneXfraction_Simpson_Rare.svg", device = "svg", width = 5, height = 5)
#
  # Chao1 -------------------------------------------------------------------
    chao_pool <- as.data.frame(estimateR(rarecurve))
    for (i in colnames(chao_pool)) { Div_Table[i,"chao"] <- chao_pool["S.chao1",i] }
    for (i in colnames(chao_pool)) { Div_Table[i,"chao_se"] <- chao_pool["se.chao1",i] }
  ## Period
    my_comp <- list(c("04","06"),c("04","09"),c("04","11"),c("06","09"),c("06","11"),c("09","11"))
    adiv <- ggplot(Div_Table, aes(x = Period, y = chao)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
      geom_errorbar(aes(ymin=chao-chao_se,ymax=chao+chao_se), width=.1,position=position_dodge(0.05)) +
      stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
      theme(axis.title = element_text(face="bold", size=12), 
            axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
            title = element_text(face="bold", size=14),
            legend.title = element_text(face="bold"),
            legend.position = "right",
            legend.text = element_text(size=10)) +
      labs(x="Period",y="Diversity") #+ guides(color = FALSE)
    print(adiv)
    legend <- get_legend(adiv)
    adiv<- adiv + theme(legend.position="none")
  ## Cycle
    my_comp <- list(c("Day","Night"))
    bdiv <- ggplot(Div_Table, aes(x = Cycle, y = chao)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
      geom_errorbar(aes(ymin=chao-chao_se,ymax=chao+chao_se), width=.1,position=position_dodge(0.05)) +
      stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
      theme(axis.title = element_text(face="bold", size=12), 
            axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
            title = element_text(face="bold", size=14),
            legend.title = element_text(face="bold"),
            legend.position = "right",
            legend.text = element_text(size=10)) +
      labs(x="Cycle",y="Diversity") #+ guides(color = FALSE)
    print(bdiv)
    bdiv<- bdiv + theme(legend.position="none")
  ## Fraction
    my_comp <- list(c("Small","Large"))
    cdiv <- ggplot(Div_Table, aes(x = Fraction, y = chao)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
      geom_errorbar(aes(ymin=chao-chao_se,ymax=chao+chao_se), width=.1,position=position_dodge(0.05)) +
      stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
      theme(axis.title = element_text(face="bold", size=12), 
            axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
            title = element_text(face="bold", size=14),
            legend.title = element_text(face="bold"),
            legend.position = "right",
            legend.text = element_text(size=10)) +
      labs(x="Fraction",y="Diversity") #+ guides(color = FALSE)
    print(cdiv)
    cdiv<- cdiv + theme(legend.position="none")
  ## Zone
    my_comp <- list(c("Mixolimnion","Monimolimnion"))
    ddiv <- ggplot(Div_Table, aes(x = Zone, y = chao)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
      geom_errorbar(aes(ymin=chao-chao_se,ymax=chao+chao_se), width=.1,position=position_dodge(0.05)) + 
      stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
      theme(axis.title = element_text(face="bold", size=12), 
            axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
            title = element_text(face="bold", size=14),
            legend.title = element_text(face="bold"),
            legend.position = "right",
            legend.text = element_text(size=10)) +
      labs(x="Zone",y="Diversity") #+ guides(color = FALSE)
    print(ddiv)
    ddiv<- ddiv + theme(legend.position="none")
  ## Coplot
    fdiv<-ggplot() + theme(panel.background = element_rect(fill="white",colour = "white", size = 0.5, linetype = "solid"))
    pdf("Diversity/Analyse-Chao-Rare.pdf",width = 16.00,height = 9.00)
    Diversity_plot <- plot_grid(adiv,bdiv, fdiv, fdiv,fdiv,legend,cdiv,ddiv, labels=c("A","B", "", "", "", "", "C","D", ""), ncol = 3, nrow = 3, rel_widths = c(3,3,1.5), rel_heights = c(5,0.1,5))
    print(Diversity_plot)
    dev.off()
  ## ZoneXFraction
    my_comp <- list(c("Mixolimnion","Monimolimnion"))
    ykdiv <- ggplot(Div_Table, aes(x = Zone, y = chao)) + geom_violin(aes(fill = Zone),trim = FALSE,width=0.5) + geom_boxplot(width=0.075)+ geom_jitter(position=position_jitter(0.05)) +
      stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
      facet_grid(~Fraction) + scale_fill_manual(values = palet_monimixo) + theme_unique_art() +
      labs(x="Zone",y="Chao index") #+ guides(color = FALSE)
    print(ykdiv)
    ggsave("Diversity/ZoneXfraction_chao_Rare.svg", device = "svg", width = 5, height = 5)
#
  # ACE -------------------------------------------------------------------
    for (i in colnames(chao_pool)) { Div_Table[i,"ACE"] <- chao_pool["S.ACE",i] }
    for (i in colnames(chao_pool)) { Div_Table[i,"ACE_se"] <- chao_pool["se.ACE",i] }
  ## Period
    my_comp <- list(c("04","06"),c("04","09"),c("04","11"),c("06","09"),c("06","11"),c("09","11"))
    adiv <- ggplot(Div_Table, aes(x = Period, y = ACE)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
      geom_errorbar(aes(ymin=ACE-ACE_se,ymax=ACE+ACE_se), width=.1,position=position_dodge(0.05)) +
      stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
      theme(axis.title = element_text(face="bold", size=12), 
            axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
            title = element_text(face="bold", size=14),
            legend.title = element_text(face="bold"),
            legend.position = "right",
            legend.text = element_text(size=10)) +
      labs(x="Period",y="Diversity") #+ guides(color = FALSE)
    print(adiv)
    legend <- get_legend(adiv)
    adiv<- adiv + theme(legend.position="none")
  ## Cycle
    my_comp <- list(c("Day","Night"))
    bdiv <- ggplot(Div_Table, aes(x = Cycle, y = ACE)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
      geom_errorbar(aes(ymin=ACE-ACE_se,ymax=ACE+ACE_se), width=.1,position=position_dodge(0.05)) +
      stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
      theme(axis.title = element_text(face="bold", size=12), 
            axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
            title = element_text(face="bold", size=14),
            legend.title = element_text(face="bold"),
            legend.position = "right",
            legend.text = element_text(size=10)) +
      labs(x="Cycle",y="Diversity") #+ guides(color = FALSE)
    print(bdiv)
    bdiv<- bdiv + theme(legend.position="none")
  ## Fraction
    my_comp <- list(c("Small","Large"))
    cdiv <- ggplot(Div_Table, aes(x = Fraction, y = ACE)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
      geom_errorbar(aes(ymin=ACE-ACE_se,ymax=ACE+ACE_se), width=.1,position=position_dodge(0.05)) +
      stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
      theme(axis.title = element_text(face="bold", size=12), 
            axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
            title = element_text(face="bold", size=14),
            legend.title = element_text(face="bold"),
            legend.position = "right",
            legend.text = element_text(size=10)) +
      labs(x="Fraction",y="Diversity") #+ guides(color = FALSE)
    print(cdiv)
    cdiv<- cdiv + theme(legend.position="none")
  ## Zone
    my_comp <- list(c("Mixolimnion","Monimolimnion"))
    ddiv <- ggplot(Div_Table, aes(x = Zone, y = ACE)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
      geom_errorbar(aes(ymin=ACE-ACE_se,ymax=ACE+ACE_se), width=.1,position=position_dodge(0.05)) + 
      stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
      theme(axis.title = element_text(face="bold", size=12), 
            axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
            title = element_text(face="bold", size=14),
            legend.title = element_text(face="bold"),
            legend.position = "right",
            legend.text = element_text(size=10)) +
      labs(x="Zone",y="Diversity") #+ guides(color = FALSE)
    print(ddiv)
    ddiv<- ddiv + theme(legend.position="none")
  ## Coplot
    fdiv<-ggplot() + theme(panel.background = element_rect(fill="white",colour = "white", size = 0.5, linetype = "solid"))
    pdf("Diversity/Analyse-ACE-Rare.pdf",width = 16.00,height = 9.00)
    Diversity_plot <- plot_grid(adiv,bdiv, fdiv, fdiv,fdiv,legend,cdiv,ddiv, labels=c("A","B", "", "", "", "", "C","D", ""), ncol = 3, nrow = 3, rel_widths = c(3,3,1.5), rel_heights = c(5,0.1,5))
    print(Diversity_plot)
    dev.off()
  ## ZoneXFraction
    my_comp <- list(c("Mixolimnion","Monimolimnion"))
    yldiv <- ggplot(Div_Table, aes(x = Zone, y = ACE)) + geom_violin(aes(fill = Zone),trim = FALSE,width=0.5) + geom_boxplot(width=0.075)+ geom_jitter(position=position_jitter(0.05)) +
      stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
      facet_grid(~Fraction) + scale_fill_manual(values = palet_monimixo) + theme_unique_art() +
      labs(x="Zone",y="ACE index") #+ guides(color = FALSE)
    print(yldiv)
    ggsave("Diversity/ZoneXfraction_ACE_Rare.svg", device = "svg", width = 5, height = 5)
#
  # Obs -------------------------------------------------------------------
    for (i in colnames(chao_pool)) { Div_Table[i,"Obs"] <- chao_pool["S.obs",i] }
  ## Period
    my_comp <- list(c("04","06"),c("04","09"),c("04","11"),c("06","09"),c("06","11"),c("09","11"))
    adiv <- ggplot(Div_Table, aes(x = Period, y = Obs)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
      stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
      theme(axis.title = element_text(face="bold", size=12), 
            axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
            title = element_text(face="bold", size=14),
            legend.title = element_text(face="bold"),
            legend.position = "right",
            legend.text = element_text(size=10)) +
      labs(x="Period",y="Diversity") #+ guides(color = FALSE)
    print(adiv)
    legend <- get_legend(adiv)
    adiv<- adiv + theme(legend.position="none")
  ## Cycle
    my_comp <- list(c("Day","Night"))
    bdiv <- ggplot(Div_Table, aes(x = Cycle, y = Obs)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
      stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
      theme(axis.title = element_text(face="bold", size=12), 
            axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
            title = element_text(face="bold", size=14),
            legend.title = element_text(face="bold"),
            legend.position = "right",
            legend.text = element_text(size=10)) +
      labs(x="Cycle",y="Diversity") #+ guides(color = FALSE)
    print(bdiv)
    bdiv<- bdiv + theme(legend.position="none")
  ## Fraction
    my_comp <- list(c("Small","Large"))
    cdiv <- ggplot(Div_Table, aes(x = Fraction, y = Obs)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
      stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
      theme(axis.title = element_text(face="bold", size=12), 
            axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
            title = element_text(face="bold", size=14),
            legend.title = element_text(face="bold"),
            legend.position = "right",
            legend.text = element_text(size=10)) +
      labs(x="Fraction",y="Diversity") #+ guides(color = FALSE)
    print(cdiv)
    cdiv<- cdiv + theme(legend.position="none")
  ## Zone
    my_comp <- list(c("Mixolimnion","Monimolimnion"))
    ddiv <- ggplot(Div_Table, aes(x = Zone, y = Obs)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
      stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
      theme(axis.title = element_text(face="bold", size=12), 
            axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
            title = element_text(face="bold", size=14),
            legend.title = element_text(face="bold"),
            legend.position = "right",
            legend.text = element_text(size=10)) +
      labs(x="Zone",y="Diversity") #+ guides(color = FALSE)
    print(ddiv)
    ddiv<- ddiv + theme(legend.position="none")
  ## Coplot
    fdiv<-ggplot() + theme(panel.background = element_rect(fill="white",colour = "white", size = 0.5, linetype = "solid"))
    pdf("Diversity/Analyse-Obs-Rare.pdf",width = 16.00,height = 9.00)
    Diversity_plot <- plot_grid(adiv,bdiv, fdiv, fdiv,fdiv,legend,cdiv,ddiv, labels=c("A","B", "", "", "", "", "C","D", ""), ncol = 3, nrow = 3, rel_widths = c(3,3,1.5), rel_heights = c(5,0.1,5))
    print(Diversity_plot)
    dev.off()
  ## ZoneXFraction
    my_comp <- list(c("Mixolimnion","Monimolimnion"))
    ymdiv <- ggplot(Div_Table, aes(x = Zone, y = Obs)) + geom_violin(aes(fill = Zone),trim = FALSE,width=0.5) + geom_boxplot(width=0.075)  + geom_jitter(position=position_jitter(0.05)) +
      stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
      facet_grid(~Fraction) + scale_fill_manual(values = palet_monimixo) + theme_unique_art() +
      labs(x="Zones",y="Observed richness") #+ guides(color = FALSE)
    print(ymdiv)
    ggsave("Diversity/ZoneXfraction_Obs_Rare.svg", device = "svg", width = 5, height = 5)
  ## ZoneXPeriods
    my_comp <- list(c("04","06"),c("04","09"),c("04","11"),c("06","09"),c("06","11"),c("09","11"))
    palet_periodcomp <- c("#7375b5FF","#8ca252FF","#f6ae34FF","#ad494aFF")
    zmdiv <- ggplot(Div_Table, aes(x = Period, y = Obs)) + geom_violin(aes(fill = Period),trim = FALSE,width=0.5,alpha=0.9) + geom_boxplot(width=0.075) + geom_jitter(position=position_jitter(0.05)) +
      stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
      #stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
      facet_grid(~Zone) + scale_fill_manual(values = palet_periodcomp) + theme_unique_art() +
      labs(x="Periods",y="Shannon index") #+ guides(color = FALSE)
    print(zmdiv)
    ggsave("Diversity/ZoneXPeriod_Obs_Rare.svg", device = "svg", width = 5, height = 5)
#
  # Coplot -----------------------------------------------------------------
    yidivx <- yidiv + theme(axis.title.x=element_blank(), axis.text.x=element_blank()) + lims(y=c(0,7)) + theme(plot.margin = margin(r=1))
    bidivx <- bidiv + theme(axis.title=element_blank(), axis.text=element_blank()) + lims(y=c(0,7)) + theme(plot.margin = margin(l = 1))
    zidivx <- zidiv + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), strip.text=element_blank()) + lims(y=c(0,1.3)) + theme(plot.margin = margin(r = 1))
    aidivx <- aidiv + theme(axis.title=element_blank(), axis.text=element_blank(), strip.text=element_blank()) + lims(y=c(0,1.3)) + theme(plot.margin = margin(l = 1))
    ymdivx <- ymdiv + theme(strip.text=element_blank()) + lims(y=c(0,950)) + theme(plot.margin = margin(r = 1))
    zmdivx <- zmdiv + theme(strip.text=element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) + lims(y=c(0,950)) + theme(plot.margin = margin(l = 1))
    svglite("Diversity/All_article_Diversity.svg",width = 10,height = 9)
    Total_fig <- plot_grid(yidivx,bidivx,zidivx,aidivx,ymdivx,zmdivx, ncol = 2, nrow = 3, labels =c("A","","B","","C",""), rel_widths = c(3,3.5),rel_heights = c(3,3,3),align="v")
    print(Total_fig)
    dev.off()
#