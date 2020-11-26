#                                 _       _
#     /\                         (_)     | |
#    /  \   _ __ ___   ___  _ __  _  ___ | |_
#   / /\ \ | '_ ` _ \ / _ \| '_ \| |/ _ \| __|
#  / ____ \| | | | | | (_) | | | | | (_) | |_
# /_/    \_\_| |_| |_|\___/|_| |_| |\___/ \__|
#                               _/ |
#                              |__/
# !/usr/bin/env Rscript
# XX/XX/XXXX
#
# Script Rarecurve
# Set directory, input, output and import packages -----------------------------------------------------------
#
#output <- "Analyse-Composition-Rarefy-V4-095-199-NOfilter-Eukaryota"
#input <- "../dataPANAM/PANAM2/V4-result-095-199/OTU_distribution_tax.txt"
#region <- "V4"
#sortop <- "no"
rarecurvecolor <- no
#
args = commandArgs(trailingOnly=TRUE)

if (length(args)==2) {
cat("Enter ribosomal region (V4 or V9) : ");
region <- readLines("stdin",n=1);
cat("You entered")
str(region);
cat( "\n" )}
if (length(args)>2) {
region <- args[3]
}

output <- args[2]
input <- args[1]

# Import packages ---------------------------------------------------------
#if (!require("pacman")) {install.packages("pacman", repos="http://cran.rstudio.com/")}
#pacman::p_load("ggplot2", "readxl","dplyr","tidyr","cowplot","FactoMineR","factoextra","reshape2","varhandle","gplots","ggrepel","ggpubr","ggsci","scales","hrbrthemes","GUniFrac","svglite")

pkg <- c("parallel","ggplot2", "readxl","dplyr","tidyr","cowplot","FactoMineR","factoextra","reshape2","varhandle","ggrepel","ggpubr","ggsci","scales","hrbrthemes","GUniFrac","svglite")
lapply(pkg, require, character.only = TRUE)
input <- paste("..",input, sep = "/")
output <- paste("../result",output, sep = "/")
if (dir.exists("../result") == FALSE) { dir.create("../result") }
if (dir.exists(output) == FALSE) { dir.create(output) }
if (dir.exists(output) == TRUE) { setwd(output) }

# quickRareCurve fonction -------------------------------------------------
quickRareCurve <- function (x, step = 1, sample, xlab = "Sequences",
                            ylab = "OTUs", label = TRUE, col, lty, max.cores = T, nCores = 1, ...)
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

# Input OTU Table ---------------------------------------------------------
tableVinput <- read.csv(file = input, sep = "\t", row.names = "OTU_Id")
##0.0005% filter
if (length(args)==2) {
  cat("Should I filter 0.0005% of total OTUs (yes or no) ? : ");
  sortop <- readLines("stdin",n=1);
  cat("You entered")
  str(sortop);
  cat( "\n" )}
if (length(args)>2) {
  sortop <- args[4]
}

if (sortop == "yes") {
  amplicon <- grep(pattern = "OSTA", colnames(tableVinput), value = TRUE)
  tableVinput$SUM <- rowSums(tableVinput %>% select(all_of(amplicon)))
  tableVinput <- tableVinput %>% filter(SUM > 0.0005*sum(tableVinput$SUM)/100)
  tableVinput <- tableVinput %>% select(-"SUM")
}

# Prepare data inf --------------------------------------------------------
infdataini <- read.table(file = "../../rawdata/data-inf.txt", sep = "\t", header = FALSE,row.names = "row", col.names = c("row","Id","Condition","Technologie","Region"))
infdataini$Rep <- rep("_1", each = nrow(infdataini))
infdataini$Variable <- paste(infdataini$Condition,infdataini$Rep, sep = "")
infdataini <- infdataini %>% select(-"Rep",-"Condition")
infdataini <- separate(infdataini, Variable, c("Condition","Date","Replicat"), sep = "_")
for (i in row.names(infdataini)) { 
  if (infdataini[i,"Replicat"] == "01") infdataini[i,"Replicat"] <- "1"
  if (infdataini[i,"Replicat"] == "02") infdataini[i,"Replicat"] <- "2"
}
infdataini$OSTA <- rep("OSTA", each = nrow(infdataini))
infdataini$Variable <- paste(infdataini$Id,infdataini$OSTA, sep = "")
infdataini <- infdataini %>% select(-"Id",-"OSTA")
infdataini <- separate(infdataini, Variable, c("Cin","Amplicon"), sep = "_")
infdataini <- infdataini %>% select(-"Cin")

## FR sample (V9 in reality and not V4)
for (i in row.names(infdataini)) { if (infdataini[i,"Amplicon"] == "FROSTA") { infdataini[i,"Region"] <- "V9"}}
for (i in row.names(infdataini)) { if (infdataini[i,"Amplicon"] == "FROSTA") { infdataini[i,"Technologie"] <- "Miseq"}}

## Prepare sample_df
infdataini <- separate(infdataini, Condition, c("empty","ADN","Cycle","Zone","Fraction"),sep = "")
infdataini <- infdataini %>% select(-"empty")
infdataini<-as.data.frame(infdataini)
col <- c("Amplicon","Technologie","Region","ADN","Cycle","Zone","Fraction","Date","Replicat")
samples_df <- infdataini[,all_of(col)]
samples_df$Condition <- paste(samples_df$ADN,samples_df$Cycle,samples_df$Zone,samples_df$Fraction, sep = "")
for (i in row.names(samples_df)) {
  if (samples_df[i,"Cycle"] == "J") { samples_df[i,"Cycle"] <- "Jour"}
  if (samples_df[i,"Cycle"] == "N") { samples_df[i,"Cycle"] <- "Nuit"}
  if (samples_df[i,"Zone"] == "O") { samples_df[i,"Zone"] <- "Oxique"}
  if (samples_df[i,"Zone"] == "A") { samples_df[i,"Zone"] <- "Anoxique"}
  if (samples_df[i,"Fraction"] == "G") { samples_df[i,"Fraction"] <- "Grande"}
  if (samples_df[i,"Fraction"] == "P") { samples_df[i,"Fraction"] <- "Petite"}}

# Select V4 or V9 ---------------------------------------------------------
samples_df<-filter(samples_df, Region == region)
# Prepare Object Nopool-------------------------------------------------
  ## Prepare seq_mat ---------------------------------------------------------------------
amplicon <- c(grep(pattern = "OSTA", colnames(tableVinput), value = TRUE))
seq_mat <- tableVinput %>% select(all_of(amplicon))
seq_mat[,all_of(amplicon)] <- lapply(seq_mat[,all_of(amplicon)], as.numeric)
  ## Rarefy (seq_mat_rare) -----------------------------------------------------------------------
seq_matt <- seq_mat
seq_matt <- t(seq_matt)
seq_mat_rare <- as.data.frame(t(Rarefy(seq_matt)$otu.tab.rff)) ## Ref : Jun Chen et al. (2012). Associating microbiome composition with environmental covariates using generalized UniFrac distances. 28(16): 2106–2113.

  ## Rarecurve ---------------------------------------------------------------
    # Raw --------------------------------------------------------------------
#rarecurve-raw
rawcurve <- seq_matt
pdf("Rarecurve/Nopool/Rarecurve-Raw.pdf",width = 5.00,height = 5.00)
curve <- quickRareCurve(rawcurve, col = "black", cex = 0.6)
dev.off()
#rarecurve all
#rawcurve_all <- t(colSums(rawcurve))
#rownames(rawcurve_all) <- "sum"
#pdf("Rarecurve/Nopool/Rarecurve-all-Raw.pdf",width = 5.00,height = 5.00)
#curve <- quickRareCurve(rawcurve_all, col = "black", cex = 0.6)
#dev.off()
if (rarecurvecolor == "yes") {
#Color
colraw <- as.data.frame(rownames(rawcurve))
colnames(colraw) <- "sample"
w <- 0
for (x in colraw[,"sample"]) { w <- w+1
colraw[w,1:10] <- subset(samples_df, Amplicon == x)}
#Color Cycle
colraw$color_cycle <- rep("black", time = nrow(colraw))
for (w in rownames(colraw)) {
  if (colraw[w,"Cycle"] == "Jour") { colraw[w,"color_cycle"] <- "#D43F3AFF"}}
pdf("Rarecurve/Nopool/Rarecurve-Cycle-Rawdata.pdf", width = 5.00, height = 5.00)
curve_a <- quickRareCurve(rawcurve, col = colraw$color_cycle, cex = 0.6)
legend("topleft", inset=.02, title="Cycle",
       c("Day","Night"), fill = c("#D43F3AFF","black"), horiz=TRUE, cex=0.6)
print(curve_a)
dev.off()
#Color Oxygène
colraw$color_Ox <- rep("black", time = nrow(colraw))
for (w in rownames(colraw)) {
  if (colraw[w,"Zone"] == "Oxique") { colraw[w,"color_Ox"] <- "#46B8DAFF"}}
pdf("Rarecurve/Nopool/Rarecurve-Zone-Rawdata.pdf", width = 5.00, height = 5.00)
curve_b <- quickRareCurve(rawcurve, col = colraw$color_Ox, cex = 0.6)
legend("topleft", inset=.02, title="Zone",
       c("Oxique","Anoxique"), fill = c("#46B8DAFF","black"), horiz=TRUE, cex=0.6)
print(curve_b)
dev.off()
#Color Taille
colraw$color_Size <- rep("black", time = nrow(colraw))
for (w in rownames(colraw)) {
  if (colraw[w,"Fraction"] == "Petite") { colraw[w,"color_Size"] <- "#5CB85CFF"}}
pdf("Rarecurve/Nopool/Rarecurve-Fraction-Rawdata.pdf", width = 5.00, height = 5.00)
curve_c <- quickRareCurve(rawcurve, col = colraw$color_Size, cex = 0.6)
legend("topleft", inset=.02, title="Fraction",
       c("Petite","Grande"), fill = c("#5CB85CFF","black"), horiz=TRUE, cex=0.6)
print(curve_c)
dev.off()
#Color Replicat
colraw$color_Rep <- rep("black", time = nrow(colraw))
for (w in rownames(colraw)) {
  if (colraw[w,"Replicat"] == 1) { colraw[w,"color_Rep"] <- "#EEA236CC"}}
pdf("Rarecurve/Nopool/Rarecurve-Rep-Rawdata.pdf",width = 5.00,height = 5.00)
curve_d <- quickRareCurve(rawcurve, col = colraw$color_Rep, cex = 0.6)
legend("topleft", inset=.02, title="Replicat",
       c("Rep1","Rep2"), fill = c("#EEA236CC","black"), horiz=TRUE, cex=0.6)
print(curve_d)
dev.off()
}

    # Rare --------------------------------------------------------------------
  rarecurve <- t(seq_mat_rare)
  pdf("Rarecurve/Nopool/Rarecurve-Rarefy.pdf",width = 5.00,height = 5.00)
  curve <- quickRareCurve(rarecurve, col = "black", cex = 0.6)
  dev.off()
  #rarecurve all
  #rarecurve_all <- t(colSums(rarecurve))
  #rownames(rarecurve_all) <- "sum"
  #pdf("Rarecurve/Nopool/Rarecurve-all-Rarefy.pdf",width = 5.00,height = 5.00)
  #curve <- quickRareCurve(rarecurve_all, col = "black", cex = 0.6)
  #dev.off()
  if (rarecurvecolor == "yes") {
  #Color
  colraw <- as.data.frame(rownames(rarecurve))
  colnames(colraw) <- "sample"
  row.names(colraw) <- colraw[,"sample"]
  for (x in row.names(colraw)) {
    colraw[x,1:10] <- samples_df %>% filter(Amplicon == x)}
  #Color Cycle
  colraw$color_cycle <- rep("black", time = nrow(colraw))
  for (w in rownames(colraw)) {
    if (colraw[w,"Cycle"] == "Jour") { colraw[w,"color_cycle"] <- "#D43F3AFF"}}
  pdf("Rarecurve/Nopool/Rarecurve-Cycle-Rarefy.pdf",width = 5.00,height = 5.00)
  curve_a <- quickRareCurve(rarecurve, col = colraw$color_cycle, cex = 0.6)
  legend("topleft", inset=.02, title="Cycle",
         c("Day","Night"), fill = c("#D43F3AFF","black"), horiz=TRUE, cex=0.6)
  print(curve_a)
  dev.off()
  #Color Oxygène
  colraw$color_Ox <- rep("black", time = nrow(colraw))
  for (w in rownames(colraw)) {
    if (colraw[w,"Zone"] == "Oxique") { colraw[w,"color_Ox"] <- "#46B8DAFF"}}
  pdf("Rarecurve/Nopool/Rarecurve-Zone-Rarefy.pdf",width = 5.00,height = 5.00)
  curve_b <- quickRareCurve(rarecurve, col = colraw$color_Ox, cex = 0.6)
  legend("topleft", inset=.02, title="Zone",
         c("Oxique","Anoxique"), fill = c("#46B8DAFF","black"), horiz=TRUE, cex=0.6)
  print(curve_b)
  dev.off()
  #Color Taille
  colraw$color_Size <- rep("black", time = nrow(colraw))
  for (w in rownames(colraw)) {
    if (colraw[w,"Fraction"] == "Petite") { colraw[w,"color_Size"] <- "#5CB85CFF"}}
  pdf("Rarecurve/Nopool/Rarecurve-Fraction-Rarefy.pdf",width = 5.00,height = 5.00)
  curve_c <- quickRareCurve(rarecurve, col = colraw$color_Size, cex = 0.6)
  legend("topleft", inset=.02, title="Fraction",
         c("Petite","Grande"), fill = c("#5CB85CFF","black"), horiz=TRUE, cex=0.6)
  print(curve_c)
  dev.off()
  #Color Replicat
  colraw$color_Rep <- rep("black", time = nrow(colraw))
  for (w in rownames(colraw)) {
    if (colraw[w,"Replicat"] == 1) { colraw[w,"color_Rep"] <- "#EEA236CC"}}
  pdf("Rarecurve/Nopool/Rarecurve-Rep-Rarefy.pdf",width = 5.00,height = 5.00)
  curve_d <- quickRareCurve(rarecurve, col = colraw$color_Rep, cex = 0.6)
  legend("topleft", inset=.02, title="Replicat",
         c("Rep1","Rep2"), fill = c("#EEA236CC","black"), horiz=TRUE, cex=0.6)
  print(curve_d)
  dev.off()
  }
  ## Diversity Raw ---------------------------------------------------------------
    # Shannon -----------------------------------------------------------------
  Div_Table <- samples_df
  rownames(Div_Table) <- Div_Table[,"Amplicon"]
  Divraw <- as.data.frame(diversity(seq_matt, index = "shannon"))
  colnames(Divraw) <- "Shannon"
  for (i in rownames(Divraw)) { Div_Table[i,"Shannon"] <- Divraw[i,"Shannon"] }
  # Dates
  my_comp <- list(c("04","06"),c("04","09"),c("04","11"),c("06","09"),c("06","11"),c("09","11"))
  adiv <- ggplot(Div_Table, aes(x = Date, y = Shannon)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
    stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x="Dates",y="Diversity") #+ guides(color = FALSE)
  print(adiv)
  legend <- get_legend(adiv)
  adiv<- adiv + theme(legend.position="none")
  
  #Cycle
  my_comp <- list(c("Jour","Nuit"))
  bdiv <- ggplot(Div_Table, aes(x = Cycle, y = Shannon)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
    stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x="Cycles",y="Diversity") #+ guides(color = FALSE)
  print(bdiv)
  bdiv<- bdiv + theme(legend.position="none")
  #Fraction
  my_comp <- list(c("Petite","Grande"))
  cdiv <- ggplot(Div_Table, aes(x = Fraction, y = Shannon)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
    stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x="Fractions",y="Diversity") #+ guides(color = FALSE)
  print(cdiv)
  cdiv<- cdiv + theme(legend.position="none")
  #Zone
  my_comp <- list(c("Oxique","Anoxique"))
  ddiv <- ggplot(Div_Table, aes(x = Zone, y = Shannon)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
    stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x="Zones",y="Diversity") #+ guides(color = FALSE)
  print(ddiv)
  ddiv<- ddiv + theme(legend.position="none")
  #Coplot
  fdiv<-ggplot() + theme(panel.background = element_rect(fill="white",colour = "white", size = 0.5, linetype = "solid"))
  pdf("Diversity/Nopool/Analyse-Shannon-Raw.pdf",width = 16.00,height = 9.00)
  Diversity_plot <- plot_grid(adiv,bdiv, fdiv, fdiv,fdiv,legend,cdiv,ddiv, labels=c("A","B", "", "", "", "", "C","D", ""), ncol = 3, nrow = 3, rel_widths = c(3,3,1.5), rel_heights = c(5,0.1,5))
  print(Diversity_plot)
  dev.off()
    # Simpson -----------------------------------------------------------------
  Divraw <- as.data.frame(diversity(seq_matt, index = "simpson"))
  colnames(Divraw) <- "Simpson"
  for (i in rownames(Divraw)) { Div_Table[i,"Simpson"] <- Divraw[i,"Simpson"] }
  
  # Date
  my_comp <- list(c("04","06"),c("04","09"),c("04","11"),c("06","09"),c("06","11"),c("09","11"))
  adiv <- ggplot(Div_Table, aes(x = Date, y = Simpson)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
    stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x="Dates",y="Diversity") #+ guides(color = FALSE)
  print(adiv)
  legend <- get_legend(adiv)
  adiv<- adiv + theme(legend.position="none")
  
  #Cycle
  my_comp <- list(c("Jour","Nuit"))
  bdiv <- ggplot(Div_Table, aes(x = Cycle, y = Simpson)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
    stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x="Cycles",y="Diversity") #+ guides(color = FALSE)
  print(bdiv)
  bdiv<- bdiv + theme(legend.position="none")
  #Fraction
  my_comp <- list(c("Petite","Grande"))
  cdiv <- ggplot(Div_Table, aes(x = Fraction, y = Simpson)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
    stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x="Fractions",y="Diversity") #+ guides(color = FALSE)
  print(cdiv)
  cdiv<- cdiv + theme(legend.position="none")
  #Zone
  my_comp <- list(c("Oxique","Anoxique"))
  ddiv <- ggplot(Div_Table, aes(x = Zone, y = Simpson)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
    stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x="Zones",y="Diversity") #+ guides(color = FALSE)
  print(ddiv)
  ddiv<- ddiv + theme(legend.position="none")
  #Coplot
  fdiv<-ggplot() + theme(panel.background = element_rect(fill="white",colour = "white", size = 0.5, linetype = "solid"))
  pdf("Diversity/Nopool/Analyse-Simpson-Raw.pdf",width = 16.00,height = 9.00)
  Diversity_plot <- plot_grid(adiv,bdiv, fdiv, fdiv,fdiv,legend,cdiv,ddiv, labels=c("A","B", "", "", "", "", "C","D", ""), ncol = 3, nrow = 3, rel_widths = c(3,3,1.5), rel_heights = c(5,0.1,5))
  print(Diversity_plot)
  dev.off()
  
    # Chao1 -------------------------------------------------------------------
  chao_pool <- as.data.frame(estimateR(seq_matt))
  for (i in colnames(chao_pool)) { Div_Table[i,"chao"] <- chao_pool["S.chao1",i] }
  for (i in colnames(chao_pool)) { Div_Table[i,"chao_se"] <- chao_pool["se.chao1",i] }
  #Date
  my_comp <- list(c("04","06"),c("04","09"),c("04","11"),c("06","09"),c("06","11"),c("09","11"))
  adiv <- ggplot(Div_Table, aes(x = Date, y = chao)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
    geom_errorbar(aes(ymin=chao-chao_se,ymax=chao+chao_se), width=.1,position=position_dodge(0.05)) +
    stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x="Dates",y="Diversity") #+ guides(color = FALSE)
  print(adiv)
  legend <- get_legend(adiv)
  adiv<- adiv + theme(legend.position="none")
  
  #Cycle
  my_comp <- list(c("Jour","Nuit"))
  bdiv <- ggplot(Div_Table, aes(x = Cycle, y = chao)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
    geom_errorbar(aes(ymin=chao-chao_se,ymax=chao+chao_se), width=.1,position=position_dodge(0.05)) +
    stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x="Cycles",y="Diversity") #+ guides(color = FALSE)
  print(bdiv)
  bdiv<- bdiv + theme(legend.position="none")
  #Fraction
  my_comp <- list(c("Petite","Grande"))
  cdiv <- ggplot(Div_Table, aes(x = Fraction, y = chao)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
    geom_errorbar(aes(ymin=chao-chao_se,ymax=chao+chao_se), width=.1,position=position_dodge(0.05)) +
    stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x="Fractions",y="Diversity") #+ guides(color = FALSE)
  print(cdiv)
  cdiv<- cdiv + theme(legend.position="none")
  #Zone
  my_comp <- list(c("Oxique","Anoxique"))
  ddiv <- ggplot(Div_Table, aes(x = Zone, y = chao)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
    geom_errorbar(aes(ymin=chao-chao_se,ymax=chao+chao_se), width=.1,position=position_dodge(0.05)) + 
    stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x="Zones",y="Diversity") #+ guides(color = FALSE)
  print(ddiv)
  ddiv<- ddiv + theme(legend.position="none")
  #Coplot
  fdiv<-ggplot() + theme(panel.background = element_rect(fill="white",colour = "white", size = 0.5, linetype = "solid"))
  pdf("Diversity/Nopool/Analyse-Chao-Raw.pdf",width = 16.00,height = 9.00)
  Diversity_plot <- plot_grid(adiv,bdiv, fdiv, fdiv,fdiv,legend,cdiv,ddiv, labels=c("A","B", "", "", "", "", "C","D", ""), ncol = 3, nrow = 3, rel_widths = c(3,3,1.5), rel_heights = c(5,0.1,5))
  print(Diversity_plot)
  dev.off()
  
    # ACE -------------------------------------------------------------------
  for (i in colnames(chao_pool)) { Div_Table[i,"ACE"] <- chao_pool["S.ACE",i] }
  for (i in colnames(chao_pool)) { Div_Table[i,"ACE_se"] <- chao_pool["se.ACE",i] }
  #Dates
  my_comp <- list(c("04","06"),c("04","09"),c("04","11"),c("06","09"),c("06","11"),c("09","11"))
  adiv <- ggplot(Div_Table, aes(x = Date, y = ACE)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
    geom_errorbar(aes(ymin=ACE-ACE_se,ymax=ACE+ACE_se), width=.1,position=position_dodge(0.05)) +
    stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
    theme(axis.title = element_text(face="bold", size=12),
          axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5),
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x="Dates",y="Diversity") #+ guides(color = FALSE)
  print(adiv)
  legend <- get_legend(adiv)
  adiv<- adiv + theme(legend.position="none")
  
  #Cycle
  my_comp <- list(c("Jour","Nuit"))
  bdiv <- ggplot(Div_Table, aes(x = Cycle, y = ACE)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
    geom_errorbar(aes(ymin=ACE-ACE_se,ymax=ACE+ACE_se), width=.1,position=position_dodge(0.05)) +
    stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x="Cycles",y="Diversity") #+ guides(color = FALSE)
  print(bdiv)
  bdiv<- bdiv + theme(legend.position="none")
  #Fraction
  my_comp <- list(c("Petite","Grande"))
  cdiv <- ggplot(Div_Table, aes(x = Fraction, y = ACE)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
    geom_errorbar(aes(ymin=ACE-ACE_se,ymax=ACE+ACE_se), width=.1,position=position_dodge(0.05)) +
    stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x="Fractions",y="Diversity") #+ guides(color = FALSE)
  print(cdiv)
  cdiv<- cdiv + theme(legend.position="none")
  #Zone
  my_comp <- list(c("Oxique","Anoxique"))
  ddiv <- ggplot(Div_Table, aes(x = Zone, y = ACE)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
    geom_errorbar(aes(ymin=ACE-ACE_se,ymax=ACE+ACE_se), width=.1,position=position_dodge(0.05)) + 
    stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=90, size=10, hjust = 0.5, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x="Zones",y="Diversity") #+ guides(color = FALSE)
  print(ddiv)
  ddiv<- ddiv + theme(legend.position="none")
  #Coplot
  fdiv<-ggplot() + theme(panel.background = element_rect(fill="white",colour = "white", size = 0.5, linetype = "solid"))
  pdf("Diversity/Nopool/Analyse-ACE-Raw.pdf",width = 16.00,height = 9.00)
  Diversity_plot <- plot_grid(adiv,bdiv, fdiv, fdiv,fdiv,legend,cdiv,ddiv, labels=c("A","B", "", "", "", "", "C","D", ""), ncol = 3, nrow = 3, rel_widths = c(3,3,1.5), rel_heights = c(5,0.1,5))
  print(Diversity_plot)
  dev.off()
  
  
    # Obs -------------------------------------------------------------------
  for (i in colnames(chao_pool)) { Div_Table[i,"Obs"] <- chao_pool["S.obs",i] }
  #Date
  my_comp <- list(c("04","06"),c("04","09"),c("04","11"),c("06","09"),c("06","11"),c("09","11"))
  adiv <- ggplot(Div_Table, aes(x = Date, y = Obs)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
    stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x="Dates",y="Diversity") #+ guides(color = FALSE)
  print(adiv)
  legend <- get_legend(adiv)
  adiv<- adiv + theme(legend.position="none")
  
  #Cycle
  my_comp <- list(c("Jour","Nuit"))
  bdiv <- ggplot(Div_Table, aes(x = Cycle, y = Obs)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
    stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x="Cycles",y="Diversity") #+ guides(color = FALSE)
  print(bdiv)
  bdiv<- bdiv + theme(legend.position="none")
  #Fraction(Taille)
  my_comp <- list(c("Petite","Grande"))
  cdiv <- ggplot(Div_Table, aes(x = Fraction, y = Obs)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
    stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x="Fractions",y="Diversity") #+ guides(color = FALSE)
  print(cdiv)
  cdiv<- cdiv + theme(legend.position="none")
  #Zone
  my_comp <- list(c("Oxique","Anoxique"))
  ddiv <- ggplot(Div_Table, aes(x = Zone, y = Obs)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
    stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x="Zones",y="Diversity") #+ guides(color = FALSE)
  print(ddiv)
  ddiv<- ddiv + theme(legend.position="none")
  #Coplot
  fdiv<-ggplot() + theme(panel.background = element_rect(fill="white",colour = "white", size = 0.5, linetype = "solid"))
  pdf("Diversity/Nopool/Analyse-Obs-Raw.pdf",width = 16.00,height = 9.00)
  Diversity_plot <- plot_grid(adiv,bdiv, fdiv, fdiv,fdiv,legend,cdiv,ddiv, labels=c("A","B", "", "", "", "", "C","D", ""), ncol = 3, nrow = 3, rel_widths = c(3,3,1.5), rel_heights = c(5,0.1,5))
  print(Diversity_plot)
  dev.off()
  
  
  
  ## Diversity Rare ---------------------------------------------------------------
    # Shannon -----------------------------------------------------------------
  Div_Table <- samples_df
  rownames(Div_Table) <- Div_Table[,"Amplicon"]
  Divrare <- as.data.frame(diversity(t(seq_mat_rare), index = "shannon"))
  colnames(Divrare) <- "Shannon"
  for (i in rownames(Divrare)) { Div_Table[i,"Shannon"] <- Divrare[i,"Shannon"] }
  Div_Table <- Div_Table %>% filter(Replicat == 1)
  
  # Date
  my_comp <- list(c("04","06"),c("04","09"),c("04","11"),c("06","09"),c("06","11"),c("09","11"))
  adiv <- ggplot(Div_Table, aes(x = Date, y = Shannon)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
    stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x="Dates",y="Diversity") #+ guides(color = FALSE)
  print(adiv)
  legend <- get_legend(adiv)
  adiv<- adiv + theme(legend.position="none")
  
  #Cycle
  my_comp <- list(c("Jour","Nuit"))
  bdiv <- ggplot(Div_Table, aes(x = Cycle, y = Shannon)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
    stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x="Cycles",y="Diversity") #+ guides(color = FALSE)
  print(bdiv)
  bdiv<- bdiv + theme(legend.position="none")
  #Fraction
  my_comp <- list(c("Petite","Grande"))
  cdiv <- ggplot(Div_Table, aes(x = Fraction, y = Shannon)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
    stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x="Fractions",y="Diversity") #+ guides(color = FALSE)
  print(cdiv)
  cdiv<- cdiv + theme(legend.position="none")
  #Zone
  my_comp <- list(c("Oxique","Anoxique"))
  ddiv <- ggplot(Div_Table, aes(x = Zone, y = Shannon)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
    stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x="Zones",y="Diversity") #+ guides(color = FALSE)
  print(ddiv)
  ddiv<- ddiv + theme(legend.position="none")
  #Coplot
  fdiv<-ggplot() + theme(panel.background = element_rect(fill="white",colour = "white", size = 0.5, linetype = "solid"))
  pdf("Diversity/Nopool/Analyse-Shannon-Rare.pdf",width = 16.00,height = 9.00)
  Diversity_plot <- plot_grid(adiv,bdiv, fdiv, fdiv,fdiv,legend,cdiv,ddiv, labels=c("A","B", "", "", "", "", "C","D", ""), ncol = 3, nrow = 3, rel_widths = c(3,3,1.5), rel_heights = c(5,0.1,5))
  print(Diversity_plot)
  dev.off()
    # Simpson -----------------------------------------------------------------
  Divrare <- as.data.frame(diversity(t(seq_mat_rare), index = "simpson"))
  colnames(Divrare) <- "Simpson"
  for (i in rownames(Divrare)) { Div_Table[i,"Simpson"] <- Divrare[i,"Simpson"] }
  # Date
  my_comp <- list(c("04","06"),c("04","09"),c("04","11"),c("06","09"),c("06","11"),c("09","11"))
  adiv <- ggplot(Div_Table, aes(x = Date, y = Simpson)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
    stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x="Dates",y="Diversity") #+ guides(color = FALSE)
  print(adiv)
  legend <- get_legend(adiv)
  adiv<- adiv + theme(legend.position="none")
  
  #Cycle
  my_comp <- list(c("Jour","Nuit"))
  bdiv <- ggplot(Div_Table, aes(x = Cycle, y = Simpson)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
    stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x="Cycles",y="Diversity") #+ guides(color = FALSE)
  print(bdiv)
  bdiv<- bdiv + theme(legend.position="none")
  #Fraction
  my_comp <- list(c("Petite","Grande"))
  cdiv <- ggplot(Div_Table, aes(x = Fraction, y = Simpson)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
    stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x="Fractions",y="Diversity") #+ guides(color = FALSE)
  print(cdiv)
  cdiv<- cdiv + theme(legend.position="none")
  #Zone
  my_comp <- list(c("Oxique","Anoxique"))
  ddiv <- ggplot(Div_Table, aes(x = Zone, y = Simpson)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
    stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x="Zones",y="Diversity") #+ guides(color = FALSE)
  print(ddiv)
  ddiv<- ddiv + theme(legend.position="none")
  #Coplot
  fdiv<-ggplot() + theme(panel.background = element_rect(fill="white",colour = "white", size = 0.5, linetype = "solid"))
  pdf("Diversity/Nopool/Analyse-Simpson-Rare.pdf",width = 16.00,height = 9.00)
  Diversity_plot <- plot_grid(adiv,bdiv, fdiv, fdiv,fdiv,legend,cdiv,ddiv, labels=c("A","B", "", "", "", "", "C","D", ""), ncol = 3, nrow = 3, rel_widths = c(3,3,1.5), rel_heights = c(5,0.1,5))
  print(Diversity_plot)
  dev.off()
  
    # Chao1 -------------------------------------------------------------------
  chao_pool <- as.data.frame(estimateR(t(seq_mat_rare)))
  for (i in colnames(chao_pool)) { Div_Table[i,"chao"] <- chao_pool["S.chao1",i] }
  for (i in colnames(chao_pool)) { Div_Table[i,"chao_se"] <- chao_pool["se.chao1",i] }
  #Date
  my_comp <- list(c("04","06"),c("04","09"),c("04","11"),c("06","09"),c("06","11"),c("09","11"))
  adiv <- ggplot(Div_Table, aes(x = Date, y = chao)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
    geom_errorbar(aes(ymin=chao-chao_se,ymax=chao+chao_se), width=.1,position=position_dodge(0.05)) +
    stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x="Dates",y="Diversity") #+ guides(color = FALSE)
  print(adiv)
  legend <- get_legend(adiv)
  adiv<- adiv + theme(legend.position="none")
  
  #Cycle
  my_comp <- list(c("Jour","Nuit"))
  bdiv <- ggplot(Div_Table, aes(x = Cycle, y = chao)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
    geom_errorbar(aes(ymin=chao-chao_se,ymax=chao+chao_se), width=.1,position=position_dodge(0.05)) +
    stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x="Cycles",y="Diversity") #+ guides(color = FALSE)
  print(bdiv)
  bdiv<- bdiv + theme(legend.position="none")
  #Fraction
  my_comp <- list(c("Petite","Grande"))
  cdiv <- ggplot(Div_Table, aes(x = Fraction, y = chao)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
    geom_errorbar(aes(ymin=chao-chao_se,ymax=chao+chao_se), width=.1,position=position_dodge(0.05)) +
    stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x="Fractions",y="Diversity") #+ guides(color = FALSE)
  print(cdiv)
  cdiv<- cdiv + theme(legend.position="none")
  #Zone
  my_comp <- list(c("Oxique","Anoxique"))
  ddiv <- ggplot(Div_Table, aes(x = Zone, y = chao)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
    geom_errorbar(aes(ymin=chao-chao_se,ymax=chao+chao_se), width=.1,position=position_dodge(0.05)) + 
    stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x="Zones",y="Diversity") #+ guides(color = FALSE)
  print(ddiv)
  ddiv<- ddiv + theme(legend.position="none")
  #Coplot
  fdiv<-ggplot() + theme(panel.background = element_rect(fill="white",colour = "white", size = 0.5, linetype = "solid"))
  pdf("Diversity/Nopool/Analyse-Chao-Rare.pdf",width = 16.00,height = 9.00)
  Diversity_plot <- plot_grid(adiv,bdiv, fdiv, fdiv,fdiv,legend,cdiv,ddiv, labels=c("A","B", "", "", "", "", "C","D", ""), ncol = 3, nrow = 3, rel_widths = c(3,3,1.5), rel_heights = c(5,0.1,5))
  print(Diversity_plot)
  dev.off()
  
    # ACE -------------------------------------------------------------------
  for (i in colnames(chao_pool)) { Div_Table[i,"ACE"] <- chao_pool["S.ACE",i] }
  for (i in colnames(chao_pool)) { Div_Table[i,"ACE_se"] <- chao_pool["se.ACE",i] }
  #Date
  my_comp <- list(c("04","06"),c("04","09"),c("04","11"),c("06","09"),c("06","11"),c("09","11"))
  adiv <- ggplot(Div_Table, aes(x = Date, y = ACE)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
    geom_errorbar(aes(ymin=ACE-ACE_se,ymax=ACE+ACE_se), width=.1,position=position_dodge(0.05)) +
    stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x="Dates",y="Diversity") #+ guides(color = FALSE)
  print(adiv)
  legend <- get_legend(adiv)
  adiv<- adiv + theme(legend.position="none")
  
  #Cycle
  my_comp <- list(c("Jour","Nuit"))
  bdiv <- ggplot(Div_Table, aes(x = Cycle, y = ACE)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
    geom_errorbar(aes(ymin=ACE-ACE_se,ymax=ACE+ACE_se), width=.1,position=position_dodge(0.05)) +
    stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x="Cycles",y="Diversity") #+ guides(color = FALSE)
  print(bdiv)
  bdiv<- bdiv + theme(legend.position="none")
  #Fraction
  my_comp <- list(c("Petite","Grande"))
  cdiv <- ggplot(Div_Table, aes(x = Fraction, y = ACE)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
    geom_errorbar(aes(ymin=ACE-ACE_se,ymax=ACE+ACE_se), width=.1,position=position_dodge(0.05)) +
    stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x="Fractions",y="Diversity") #+ guides(color = FALSE)
  print(cdiv)
  cdiv<- cdiv + theme(legend.position="none")
  #Zone
  my_comp <- list(c("Oxique","Anoxique"))
  ddiv <- ggplot(Div_Table, aes(x = Zone, y = ACE)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
    geom_errorbar(aes(ymin=ACE-ACE_se,ymax=ACE+ACE_se), width=.1,position=position_dodge(0.05)) + 
    stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=90, size=10, hjust = 0.5, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x="Zones",y="Diversity") #+ guides(color = FALSE)
  print(ddiv)
  ddiv<- ddiv + theme(legend.position="none")
  #Coplot
  fdiv<-ggplot() + theme(panel.background = element_rect(fill="white",colour = "white", size = 0.5, linetype = "solid"))
  pdf("Diversity/Nopool/Analyse-ACE-Rare.pdf",width = 16.00,height = 9.00)
  Diversity_plot <- plot_grid(adiv,bdiv, fdiv, fdiv,fdiv,legend,cdiv,ddiv, labels=c("A","B", "", "", "", "", "C","D", ""), ncol = 3, nrow = 3, rel_widths = c(3,3,1.5), rel_heights = c(5,0.1,5))
  print(Diversity_plot)
  dev.off()
  
  
    # Obs -------------------------------------------------------------------
  for (i in colnames(chao_pool)) { Div_Table[i,"Obs"] <- chao_pool["S.obs",i] }
  #Date
  my_comp <- list(c("04","06"),c("04","09"),c("04","11"),c("06","09"),c("06","11"),c("09","11"))
  adiv <- ggplot(Div_Table, aes(x = Date, y = Obs)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
    stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x="Dates",y="Diversity") #+ guides(color = FALSE)
  print(adiv)
  legend <- get_legend(adiv)
  adiv<- adiv + theme(legend.position="none")
  
  #Cycle
  my_comp <- list(c("Jour","Nuit"))
  bdiv <- ggplot(Div_Table, aes(x = Cycle, y = Obs)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
    stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x="Cycles",y="Diversity") #+ guides(color = FALSE)
  print(bdiv)
  bdiv<- bdiv + theme(legend.position="none")
  #Fraction
  my_comp <- list(c("Petite","Grande"))
  cdiv <- ggplot(Div_Table, aes(x = Fraction, y = Obs)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
    stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x="Fractions",y="Diversity") #+ guides(color = FALSE)
  print(cdiv)
  cdiv<- cdiv + theme(legend.position="none")
  #Zone
  my_comp <- list(c("Oxique","Anoxique"))
  ddiv <- ggplot(Div_Table, aes(x = Zone, y = Obs)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
    stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x="Zones",y="Diversity") #+ guides(color = FALSE)
  print(ddiv)
  ddiv<- ddiv + theme(legend.position="none")
  #Coplot
  fdiv<-ggplot() + theme(panel.background = element_rect(fill="white",colour = "white", size = 0.5, linetype = "solid"))
  pdf("Diversity/Nopool/Analyse-Obs-Rare.pdf",width = 16.00,height = 9.00)
  Diversity_plot <- plot_grid(adiv,bdiv, fdiv, fdiv,fdiv,legend,cdiv,ddiv, labels=c("A","B", "", "", "", "", "C","D", ""), ncol = 3, nrow = 3, rel_widths = c(3,3,1.5), rel_heights = c(5,0.1,5))
  print(Diversity_plot)
  dev.off()
  
  
  
# Prepare Object Pool-------------------------------------------------
  ## Unify
  if (length(args)==2) {
    cat("Should I unify duplicat (yes or no) ? : ");
    UnifyYoN <- readLines("stdin",n=1);
    cat("You entered")
    str(UnifyYoN);
    cat( "\n" )}
  if (length(args)>2) {
    UnifyYoN <- args[5]
  }
  if (UnifyYoN == "yes") {
  ## Prepare seq_mat and Pool (seq_mat and seq_mat_pool) ---------------------------------------------------------------------
  ## Prepare otu_mat
  amplicon <- c(grep(pattern = "OSTA", colnames(tableVinput), value = TRUE))
  seq_mat <- tableVinput %>% select(all_of(amplicon))
  seq_mat[,all_of(amplicon)] <- lapply(seq_mat[,all_of(amplicon)], as.numeric)
  ## Prepare correspondance table
  tblc <- samples_df %>% select("Amplicon", "Date", "Replicat", "Condition")
  tblc$Condition_Date <- paste(tblc$Condition,tblc$Date, sep = "_")
  tblc <- tblc %>% select(-"Date",-"Condition")
  tblc1 <- tblc %>% filter(tblc$Replicat == 1)
  tblc2 <- tblc %>% filter(tblc$Replicat == 2)
  tblcx <- merge(x = tblc1,y = tblc2, by = "Condition_Date", all = TRUE)
  Amplicon_without_duplicat <- c(tblcx[,"Amplicon.x"][which(is.na(tblcx[,"Amplicon.y"]) == TRUE)])
  for (i in row.names(tblcx)) { if (is.na(tblcx[i,"Amplicon.y"]) == TRUE) { tblcx[i,"Amplicon.y"] <- tblcx[i,"Amplicon.x"]}}
  for (i in row.names(tblcx)) { if (is.na(tblcx[i,"Replicat.y"]) == TRUE) { tblcx[i,"Replicat.y"] <- tblcx[i,"Replicat.x"]}}
  
  ## Pool
  seq_mat_pool <- as.data.frame(x=row.names(seq_mat),row.names = row.names(seq_mat))
  colnames(seq_mat_pool) <- "OTU_Id"
  for (h in tblcx[,"Amplicon.x"] ) {
    f <- (tblcx %>% filter(Amplicon.x == h) %>% select(Amplicon.y))[[1]]
    print (h)
    print (f)
    for (r in row.names(seq_mat)) {
      if (seq_mat[r,h] * seq_mat[r,f] == 0) { seq_mat_pool[r,h] <- round((seq_mat[r,h]+seq_mat[r,f]))}
      if (seq_mat[r,h] * seq_mat[r,f] != 0) { seq_mat_pool[r,h] <- round((seq_mat[r,h]+seq_mat[r,f])/2)}
    }}
  seq_mat_pool <- seq_mat_pool %>% select(-OTU_Id) 
  ## Rarefy (seq_mat_pool_rare) ---------------------------------------------------------------------
  seq_mat_poolt <- seq_mat_pool
  seq_mat_poolt <- t(seq_mat_poolt)
  seq_mat_pool_rare <- as.data.frame(t(Rarefy(seq_mat_poolt)$otu.tab.rff)) ## Ref : Jun Chen et al. (2012). Associating microbiome composition with environmental covariates using generalized UniFrac distances. 28(16): 2106–2113.
  ## Rarecurve ---------------------------------------------------------------
    # Raw ---------------------------------------------------------------------
    #rarecurve
    pdf("Rarecurve/Pool/Rarecurve-Raw.pdf",width = 5.00,height = 5.00)
    curve <- quickRareCurve(seq_mat_poolt, col = "black", cex = 0.6)
    dev.off()
    
    #rarecurve all
    #rawcurve_all <- t(colSums(data_poolx))
    #rownames(rawcurve_all) <- "sum"
    #pdf("Rarecurve/Pool/Rarecurve-all-Raw.pdf",width = 5.00,height = 5.00)
    #curve <- quickRareCurve(rawcurve_all, col = "black", cex = 0.6)
    #dev.off()
    if (rarecurvecolor == "yes") {
    #Color
    colraw <- as.data.frame(rownames(seq_mat_poolt))
    colnames(colraw) <- "sample"
    w <- 0
    for (x in colraw[,"sample"]) { w <- w+1
    colraw[w,1:10] <- subset(samples_df, Amplicon == x)}
    #Color Cycle
    colraw$color_cycle <- rep("black", time = nrow(colraw))
    for (w in rownames(colraw)) {
      if (colraw[w,"Cycle"] == "Jour") { colraw[w,"color_cycle"] <- "#D43F3AFF"}}
    pdf("Rarecurve/Pool/Rarecurve-Cycle-Rawdata.pdf", width = 5.00, height = 5.00)
    curve_a <- quickRareCurve(seq_mat_poolt, col = colraw$color_cycle, cex = 0.6)
    legend("topleft", inset=.02, title="Cycle",
           c("Day","Night"), fill = c("#D43F3AFF","black"), horiz=TRUE, cex=0.6)
    print(curve_a)
    dev.off()
    #Color Oxygène
    colraw$color_Ox <- rep("black", time = nrow(colraw))
    for (w in rownames(colraw)) {
      if (colraw[w,"Zone"] == "Oxique") { colraw[w,"color_Ox"] <- "#46B8DAFF"}}
    pdf("Rarecurve/Pool/Rarecurve-Zone-Rawdata.pdf", width = 5.00, height = 5.00)
    curve_b <- quickRareCurve(seq_mat_poolt, col = colraw$color_Ox, cex = 0.6)
    legend("topleft", inset=.02, title="Zone",
           c("Oxique","Anoxique"), fill = c("#46B8DAFF","black"), horiz=TRUE, cex=0.6)
    print(curve_b)
    dev.off()
    #Color Taille
    colraw$color_Size <- rep("black", time = nrow(colraw))
    for (w in rownames(colraw)) {
      if (colraw[w,"Fraction"] == "Petite") { colraw[w,"color_Size"] <- "#5CB85CFF"}}
    pdf("Rarecurve/Pool/Rarecurve-Fraction-Rawdata.pdf", width = 5.00, height = 5.00)
    curve_c <- quickRareCurve(seq_mat_poolt, col = colraw$color_Size, cex = 0.6)
    legend("topleft", inset=.02, title="Fraction",
           c("Petite","Grande"), fill = c("#5CB85CFF","black"), horiz=TRUE, cex=0.6)
    print(curve_c)
    dev.off()
    }
    # Rare --------------------------------------------------------------------
    #rarecurve
    rarecurve <- t(seq_mat_pool_rare)
    pdf("Rarecurve/Pool/Rarecurve-Rarefy.pdf",width = 5.00,height = 5.00)
    curve <- quickRareCurve(rarecurve, col = "black", cex = 0.6)
    dev.off()
    
    #rarecurve all
    #rawcurve_all <- t(colSums(rarecurve))
    #rownames(rawcurve_all) <- "sum"
    #pdf("Rarecurve/Pool/Rarecurve-all-Rarefy.pdf",width = 5.00,height = 5.00)
    #curve <- quickRareCurve(rawcurve_all, col = "black", cex = 0.6)
    #dev.off()
    if (rarecurvecolor == "yes") {
    #Color
    colraw <- as.data.frame(rownames(rarecurve))
    colnames(colraw) <- "sample"
    w <- 0
    for (x in colraw[,"sample"]) { w <- w+1
    colraw[w,1:10] <- subset(samples_df, Amplicon == x)}
    #Color Cycle
    colraw$color_cycle <- rep("black", time = nrow(colraw))
    for (w in rownames(colraw)) {
      if (colraw[w,"Cycle"] == "Jour") { colraw[w,"color_cycle"] <- "#D43F3AFF"}}
    pdf("Rarecurve/Pool/Rarecurve-Cycle-Rarefy.pdf", width = 5.00, height = 5.00)
    curve_a <- quickRareCurve(rarecurve, col = colraw$color_cycle, cex = 0.6)
    legend("topleft", inset=.02, title="Cycle",
           c("Day","Night"), fill = c("#D43F3AFF","black"), horiz=TRUE, cex=0.6)
    print(curve_a)
    dev.off()
    #Color Oxygène
    colraw$color_Ox <- rep("black", time = nrow(colraw))
    for (w in rownames(colraw)) {
      if (colraw[w,"Zone"] == "Oxique") { colraw[w,"color_Ox"] <- "#46B8DAFF"}}
    pdf("Rarecurve/Pool/Rarecurve-Zone-Rarefy.pdf", width = 5.00, height = 5.00)
    curve_b <- quickRareCurve(rarecurve, col = colraw$color_Ox, cex = 0.6)
    legend("topleft", inset=.02, title="Zone",
           c("Oxique","Anoxique"), fill = c("#46B8DAFF","black"), horiz=TRUE, cex=0.6)
    print(curve_b)
    dev.off()
    #Color Taille
    colraw$color_Size <- rep("black", time = nrow(colraw))
    for (w in rownames(colraw)) {
      if (colraw[w,"Fraction"] == "Petite") { colraw[w,"color_Size"] <- "#5CB85CFF"}}
    pdf("Rarecurve/Pool/Rarecurve-Fraction-Rarefy.pdf", width = 5.00, height = 5.00)
    curve_c <- quickRareCurve(rarecurve, col = colraw$color_Size, cex = 0.6)
    legend("topleft", inset=.02, title="Fraction",
           c("Petite","Grande"), fill = c("#5CB85CFF","black"), horiz=TRUE, cex=0.6)
    print(curve_c)
    dev.off()
    }
  ## Diversity Raw ---------------------------------------------------------------
    # Shannon -----------------------------------------------------------------
  Div_Table <- samples_df
  rownames(Div_Table) <- Div_Table[,"Amplicon"]
  Divraw <- as.data.frame(diversity(seq_mat_poolt, index = "shannon"))
  colnames(Divraw) <- "Shannon"
  for (i in rownames(Divraw)) { Div_Table[i,"Shannon"] <- Divraw[i,"Shannon"] }
  Div_Table <- Div_Table %>% filter(Replicat == 1)
  # Dates
  my_comp <- list(c("04","06"),c("04","09"),c("04","11"),c("06","09"),c("06","11"),c("09","11"))
  adiv <- ggplot(Div_Table, aes(x = Date, y = Shannon)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
    stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x="Dates",y="Diversity") #+ guides(color = FALSE)
  print(adiv)
  legend <- get_legend(adiv)
  adiv<- adiv + theme(legend.position="none")
  
  #Cycle
  my_comp <- list(c("Jour","Nuit"))
  bdiv <- ggplot(Div_Table, aes(x = Cycle, y = Shannon)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
    stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x="Cycles",y="Diversity") #+ guides(color = FALSE)
  print(bdiv)
  bdiv<- bdiv + theme(legend.position="none")
  #Fraction
  my_comp <- list(c("Petite","Grande"))
  cdiv <- ggplot(Div_Table, aes(x = Fraction, y = Shannon)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
    stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x="Fractions",y="Diversity") #+ guides(color = FALSE)
  print(cdiv)
  cdiv<- cdiv + theme(legend.position="none")
  #Zone
  my_comp <- list(c("Oxique","Anoxique"))
  ddiv <- ggplot(Div_Table, aes(x = Zone, y = Shannon)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
    stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x="Zones",y="Diversity") #+ guides(color = FALSE)
  print(ddiv)
  ddiv<- ddiv + theme(legend.position="none")
  #Coplot
  fdiv<-ggplot() + theme(panel.background = element_rect(fill="white",colour = "white", size = 0.5, linetype = "solid"))
  pdf("Diversity/Pool/Analyse-Shannon-Raw.pdf",width = 16.00,height = 9.00)
  Diversity_plot <- plot_grid(adiv,bdiv, fdiv, fdiv,fdiv,legend,cdiv,ddiv, labels=c("A","B", "", "", "", "", "C","D", ""), ncol = 3, nrow = 3, rel_widths = c(3,3,1.5), rel_heights = c(5,0.1,5))
  print(Diversity_plot)
  dev.off()
    # Simpson -----------------------------------------------------------------
  Divraw <- as.data.frame(diversity(seq_mat_poolt, index = "simpson"))
  colnames(Divraw) <- "Simpson"
  for (i in rownames(Divraw)) { Div_Table[i,"Simpson"] <- Divraw[i,"Simpson"] }

  # Date
  my_comp <- list(c("04","06"),c("04","09"),c("04","11"),c("06","09"),c("06","11"),c("09","11"))
  adiv <- ggplot(Div_Table, aes(x = Date, y = Simpson)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
    stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x="Dates",y="Diversity") #+ guides(color = FALSE)
  print(adiv)
  legend <- get_legend(adiv)
  adiv<- adiv + theme(legend.position="none")
  
  #Cycle
  my_comp <- list(c("Jour","Nuit"))
  bdiv <- ggplot(Div_Table, aes(x = Cycle, y = Simpson)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
    stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x="Cycles",y="Diversity") #+ guides(color = FALSE)
  print(bdiv)
  bdiv<- bdiv + theme(legend.position="none")
  #Fraction
  my_comp <- list(c("Petite","Grande"))
  cdiv <- ggplot(Div_Table, aes(x = Fraction, y = Simpson)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
    stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x="Fractions",y="Diversity") #+ guides(color = FALSE)
  print(cdiv)
  cdiv<- cdiv + theme(legend.position="none")
  #Zone
  my_comp <- list(c("Oxique","Anoxique"))
  ddiv <- ggplot(Div_Table, aes(x = Zone, y = Simpson)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
    stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x="Zones",y="Diversity") #+ guides(color = FALSE)
  print(ddiv)
  ddiv<- ddiv + theme(legend.position="none")
  #Coplot
  fdiv<-ggplot() + theme(panel.background = element_rect(fill="white",colour = "white", size = 0.5, linetype = "solid"))
  pdf("Diversity/Pool/Analyse-Simpson-Raw.pdf",width = 16.00,height = 9.00)
  Diversity_plot <- plot_grid(adiv,bdiv, fdiv, fdiv,fdiv,legend,cdiv,ddiv, labels=c("A","B", "", "", "", "", "C","D", ""), ncol = 3, nrow = 3, rel_widths = c(3,3,1.5), rel_heights = c(5,0.1,5))
  print(Diversity_plot)
  dev.off()
  
    # Chao1 -------------------------------------------------------------------
  chao_pool <- as.data.frame(estimateR(seq_mat_poolt))
  for (i in colnames(chao_pool)) { Div_Table[i,"chao"] <- chao_pool["S.chao1",i] }
  for (i in colnames(chao_pool)) { Div_Table[i,"chao_se"] <- chao_pool["se.chao1",i] }
  #Date
  my_comp <- list(c("04","06"),c("04","09"),c("04","11"),c("06","09"),c("06","11"),c("09","11"))
  adiv <- ggplot(Div_Table, aes(x = Date, y = chao)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
    geom_errorbar(aes(ymin=chao-chao_se,ymax=chao+chao_se), width=.1,position=position_dodge(0.05)) +
    stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x="Dates",y="Diversity") #+ guides(color = FALSE)
  print(adiv)
  legend <- get_legend(adiv)
  adiv<- adiv + theme(legend.position="none")
  
  #Cycle
  my_comp <- list(c("Jour","Nuit"))
  bdiv <- ggplot(Div_Table, aes(x = Cycle, y = chao)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
    geom_errorbar(aes(ymin=chao-chao_se,ymax=chao+chao_se), width=.1,position=position_dodge(0.05)) +
    stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x="Cycles",y="Diversity") #+ guides(color = FALSE)
  print(bdiv)
  bdiv<- bdiv + theme(legend.position="none")
  #Fraction
  my_comp <- list(c("Petite","Grande"))
  cdiv <- ggplot(Div_Table, aes(x = Fraction, y = chao)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
    geom_errorbar(aes(ymin=chao-chao_se,ymax=chao+chao_se), width=.1,position=position_dodge(0.05)) +
    stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x="Fractions",y="Diversity") #+ guides(color = FALSE)
  print(cdiv)
  cdiv<- cdiv + theme(legend.position="none")
  #Zone
  my_comp <- list(c("Oxique","Anoxique"))
  ddiv <- ggplot(Div_Table, aes(x = Zone, y = chao)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
    geom_errorbar(aes(ymin=chao-chao_se,ymax=chao+chao_se), width=.1,position=position_dodge(0.05)) + 
    stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x="Zones",y="Diversity") #+ guides(color = FALSE)
  print(ddiv)
  ddiv<- ddiv + theme(legend.position="none")
  #Coplot
  fdiv<-ggplot() + theme(panel.background = element_rect(fill="white",colour = "white", size = 0.5, linetype = "solid"))
  pdf("Diversity/Pool/Analyse-Chao-Raw.pdf",width = 16.00,height = 9.00)
  Diversity_plot <- plot_grid(adiv,bdiv, fdiv, fdiv,fdiv,legend,cdiv,ddiv, labels=c("A","B", "", "", "", "", "C","D", ""), ncol = 3, nrow = 3, rel_widths = c(3,3,1.5), rel_heights = c(5,0.1,5))
  print(Diversity_plot)
  dev.off()
  
    # ACE -------------------------------------------------------------------
  for (i in colnames(chao_pool)) { Div_Table[i,"ACE"] <- chao_pool["S.ACE",i] }
  for (i in colnames(chao_pool)) { Div_Table[i,"ACE_se"] <- chao_pool["se.ACE",i] }
  #Dates
  my_comp <- list(c("04","06"),c("04","09"),c("04","11"),c("06","09"),c("06","11"),c("09","11"))
  adiv <- ggplot(Div_Table, aes(x = Date, y = ACE)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
    geom_errorbar(aes(ymin=ACE-ACE_se,ymax=ACE+ACE_se), width=.1,position=position_dodge(0.05)) +
    stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
    theme(axis.title = element_text(face="bold", size=12),
          axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5),
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x="Dates",y="Diversity") #+ guides(color = FALSE)
  print(adiv)
  legend <- get_legend(adiv)
  adiv<- adiv + theme(legend.position="none")
  
  #Cycle
  my_comp <- list(c("Jour","Nuit"))
  bdiv <- ggplot(Div_Table, aes(x = Cycle, y = ACE)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
    geom_errorbar(aes(ymin=ACE-ACE_se,ymax=ACE+ACE_se), width=.1,position=position_dodge(0.05)) +
    stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x="Cycles",y="Diversity") #+ guides(color = FALSE)
  print(bdiv)
  bdiv<- bdiv + theme(legend.position="none")
  #Fraction
  my_comp <- list(c("Petite","Grande"))
  cdiv <- ggplot(Div_Table, aes(x = Fraction, y = ACE)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
    geom_errorbar(aes(ymin=ACE-ACE_se,ymax=ACE+ACE_se), width=.1,position=position_dodge(0.05)) +
    stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x="Fractions",y="Diversity") #+ guides(color = FALSE)
  print(cdiv)
  cdiv<- cdiv + theme(legend.position="none")
  #Zone
  my_comp <- list(c("Oxique","Anoxique"))
  ddiv <- ggplot(Div_Table, aes(x = Zone, y = ACE)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
    geom_errorbar(aes(ymin=ACE-ACE_se,ymax=ACE+ACE_se), width=.1,position=position_dodge(0.05)) + 
    stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=90, size=10, hjust = 0.5, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x="Zones",y="Diversity") #+ guides(color = FALSE)
  print(ddiv)
  ddiv<- ddiv + theme(legend.position="none")
  #Coplot
  fdiv<-ggplot() + theme(panel.background = element_rect(fill="white",colour = "white", size = 0.5, linetype = "solid"))
  pdf("Diversity/Pool/Analyse-ACE-Raw.pdf",width = 16.00,height = 9.00)
  Diversity_plot <- plot_grid(adiv,bdiv, fdiv, fdiv,fdiv,legend,cdiv,ddiv, labels=c("A","B", "", "", "", "", "C","D", ""), ncol = 3, nrow = 3, rel_widths = c(3,3,1.5), rel_heights = c(5,0.1,5))
  print(Diversity_plot)
  dev.off()
  
  
    # Obs -------------------------------------------------------------------
  for (i in colnames(chao_pool)) { Div_Table[i,"Obs"] <- chao_pool["S.obs",i] }
  #Date
  my_comp <- list(c("04","06"),c("04","09"),c("04","11"),c("06","09"),c("06","11"),c("09","11"))
  adiv <- ggplot(Div_Table, aes(x = Date, y = Obs)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
    stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x="Dates",y="Diversity") #+ guides(color = FALSE)
  print(adiv)
  legend <- get_legend(adiv)
  adiv<- adiv + theme(legend.position="none")
  
  #Cycle
  my_comp <- list(c("Jour","Nuit"))
  bdiv <- ggplot(Div_Table, aes(x = Cycle, y = Obs)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
    stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x="Cycles",y="Diversity") #+ guides(color = FALSE)
  print(bdiv)
  bdiv<- bdiv + theme(legend.position="none")
  #Fraction(Taille)
  my_comp <- list(c("Petite","Grande"))
  cdiv <- ggplot(Div_Table, aes(x = Fraction, y = Obs)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
    stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x="Fractions",y="Diversity") #+ guides(color = FALSE)
  print(cdiv)
  cdiv<- cdiv + theme(legend.position="none")
  #Zone
  my_comp <- list(c("Oxique","Anoxique"))
  ddiv <- ggplot(Div_Table, aes(x = Zone, y = Obs)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
    stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x="Zones",y="Diversity") #+ guides(color = FALSE)
  print(ddiv)
  ddiv<- ddiv + theme(legend.position="none")
  #Coplot
  fdiv<-ggplot() + theme(panel.background = element_rect(fill="white",colour = "white", size = 0.5, linetype = "solid"))
  pdf("Diversity/Pool/Analyse-Obs-Raw.pdf",width = 16.00,height = 9.00)
  Diversity_plot <- plot_grid(adiv,bdiv, fdiv, fdiv,fdiv,legend,cdiv,ddiv, labels=c("A","B", "", "", "", "", "C","D", ""), ncol = 3, nrow = 3, rel_widths = c(3,3,1.5), rel_heights = c(5,0.1,5))
  print(Diversity_plot)
  dev.off()
  
  
  
  ## Diversity Rare ---------------------------------------------------------------
    # Shannon -----------------------------------------------------------------
  Div_Table <- samples_df
  rownames(Div_Table) <- Div_Table[,"Amplicon"]
  Divrare <- as.data.frame(diversity(t(seq_mat_pool_rare), index = "shannon"))
  colnames(Divrare) <- "Shannon"
  for (i in rownames(Divrare)) { Div_Table[i,"Shannon"] <- Divrare[i,"Shannon"] }
  Div_Table <- Div_Table %>% filter(Replicat == 1)

  # Date
  my_comp <- list(c("04","06"),c("04","09"),c("04","11"),c("06","09"),c("06","11"),c("09","11"))
  adiv <- ggplot(Div_Table, aes(x = Date, y = Shannon)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
    stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x="Dates",y="Diversity") #+ guides(color = FALSE)
  print(adiv)
  legend <- get_legend(adiv)
  adiv<- adiv + theme(legend.position="none")
  
  #Cycle
  my_comp <- list(c("Jour","Nuit"))
  bdiv <- ggplot(Div_Table, aes(x = Cycle, y = Shannon)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
    stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x="Cycles",y="Diversity") #+ guides(color = FALSE)
  print(bdiv)
  bdiv<- bdiv + theme(legend.position="none")
  #Fraction
  my_comp <- list(c("Petite","Grande"))
  cdiv <- ggplot(Div_Table, aes(x = Fraction, y = Shannon)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
    stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x="Fractions",y="Diversity") #+ guides(color = FALSE)
  print(cdiv)
  cdiv<- cdiv + theme(legend.position="none")
  #Zone
  my_comp <- list(c("Oxique","Anoxique"))
  ddiv <- ggplot(Div_Table, aes(x = Zone, y = Shannon)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
    stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x="Zones",y="Diversity") #+ guides(color = FALSE)
  print(ddiv)
  ddiv<- ddiv + theme(legend.position="none")
  #Coplot
  fdiv<-ggplot() + theme(panel.background = element_rect(fill="white",colour = "white", size = 0.5, linetype = "solid"))
  pdf("Diversity/Pool/Analyse-Shannon-Rare.pdf",width = 16.00,height = 9.00)
  Diversity_plot <- plot_grid(adiv,bdiv, fdiv, fdiv,fdiv,legend,cdiv,ddiv, labels=c("A","B", "", "", "", "", "C","D", ""), ncol = 3, nrow = 3, rel_widths = c(3,3,1.5), rel_heights = c(5,0.1,5))
  print(Diversity_plot)
  dev.off()
    # Simpson -----------------------------------------------------------------
  Divrare <- as.data.frame(diversity(t(seq_mat_pool_rare), index = "simpson"))
  colnames(Divrare) <- "Simpson"
  for (i in rownames(Divrare)) { Div_Table[i,"Simpson"] <- Divrare[i,"Simpson"] }
  # Date
  my_comp <- list(c("04","06"),c("04","09"),c("04","11"),c("06","09"),c("06","11"),c("09","11"))
  adiv <- ggplot(Div_Table, aes(x = Date, y = Simpson)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
    stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x="Dates",y="Diversity") #+ guides(color = FALSE)
  print(adiv)
  legend <- get_legend(adiv)
  adiv<- adiv + theme(legend.position="none")
  
  #Cycle
  my_comp <- list(c("Jour","Nuit"))
  bdiv <- ggplot(Div_Table, aes(x = Cycle, y = Simpson)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
    stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x="Cycles",y="Diversity") #+ guides(color = FALSE)
  print(bdiv)
  bdiv<- bdiv + theme(legend.position="none")
  #Fraction
  my_comp <- list(c("Petite","Grande"))
  cdiv <- ggplot(Div_Table, aes(x = Fraction, y = Simpson)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
    stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x="Fractions",y="Diversity") #+ guides(color = FALSE)
  print(cdiv)
  cdiv<- cdiv + theme(legend.position="none")
  #Zone
  my_comp <- list(c("Oxique","Anoxique"))
  ddiv <- ggplot(Div_Table, aes(x = Zone, y = Simpson)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
    stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x="Zones",y="Diversity") #+ guides(color = FALSE)
  print(ddiv)
  ddiv<- ddiv + theme(legend.position="none")
  #Coplot
  fdiv<-ggplot() + theme(panel.background = element_rect(fill="white",colour = "white", size = 0.5, linetype = "solid"))
  pdf("Diversity/Pool/Analyse-Simpson-Rare.pdf",width = 16.00,height = 9.00)
  Diversity_plot <- plot_grid(adiv,bdiv, fdiv, fdiv,fdiv,legend,cdiv,ddiv, labels=c("A","B", "", "", "", "", "C","D", ""), ncol = 3, nrow = 3, rel_widths = c(3,3,1.5), rel_heights = c(5,0.1,5))
  print(Diversity_plot)
  dev.off()
  
    # Chao1 -------------------------------------------------------------------
  chao_pool <- as.data.frame(estimateR(t(seq_mat_pool_rare)))
  for (i in colnames(chao_pool)) { Div_Table[i,"chao"] <- chao_pool["S.chao1",i] }
  for (i in colnames(chao_pool)) { Div_Table[i,"chao_se"] <- chao_pool["se.chao1",i] }
  #Date
  my_comp <- list(c("04","06"),c("04","09"),c("04","11"),c("06","09"),c("06","11"),c("09","11"))
  adiv <- ggplot(Div_Table, aes(x = Date, y = chao)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
    geom_errorbar(aes(ymin=chao-chao_se,ymax=chao+chao_se), width=.1,position=position_dodge(0.05)) +
    stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x="Dates",y="Diversity") #+ guides(color = FALSE)
  print(adiv)
  legend <- get_legend(adiv)
  adiv<- adiv + theme(legend.position="none")
  
  #Cycle
  my_comp <- list(c("Jour","Nuit"))
  bdiv <- ggplot(Div_Table, aes(x = Cycle, y = chao)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
    geom_errorbar(aes(ymin=chao-chao_se,ymax=chao+chao_se), width=.1,position=position_dodge(0.05)) +
    stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x="Cycles",y="Diversity") #+ guides(color = FALSE)
  print(bdiv)
  bdiv<- bdiv + theme(legend.position="none")
  #Fraction
  my_comp <- list(c("Petite","Grande"))
  cdiv <- ggplot(Div_Table, aes(x = Fraction, y = chao)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
    geom_errorbar(aes(ymin=chao-chao_se,ymax=chao+chao_se), width=.1,position=position_dodge(0.05)) +
    stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x="Fractions",y="Diversity") #+ guides(color = FALSE)
  print(cdiv)
  cdiv<- cdiv + theme(legend.position="none")
  #Zone
  my_comp <- list(c("Oxique","Anoxique"))
  ddiv <- ggplot(Div_Table, aes(x = Zone, y = chao)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
    geom_errorbar(aes(ymin=chao-chao_se,ymax=chao+chao_se), width=.1,position=position_dodge(0.05)) + 
    stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x="Zones",y="Diversity") #+ guides(color = FALSE)
  print(ddiv)
  ddiv<- ddiv + theme(legend.position="none")
  #Coplot
  fdiv<-ggplot() + theme(panel.background = element_rect(fill="white",colour = "white", size = 0.5, linetype = "solid"))
  pdf("Diversity/Pool/Analyse-Chao-Rare.pdf",width = 16.00,height = 9.00)
  Diversity_plot <- plot_grid(adiv,bdiv, fdiv, fdiv,fdiv,legend,cdiv,ddiv, labels=c("A","B", "", "", "", "", "C","D", ""), ncol = 3, nrow = 3, rel_widths = c(3,3,1.5), rel_heights = c(5,0.1,5))
  print(Diversity_plot)
  dev.off()
  
    # ACE -------------------------------------------------------------------
  for (i in colnames(chao_pool)) { Div_Table[i,"ACE"] <- chao_pool["S.ACE",i] }
  for (i in colnames(chao_pool)) { Div_Table[i,"ACE_se"] <- chao_pool["se.ACE",i] }
  #Date
  my_comp <- list(c("04","06"),c("04","09"),c("04","11"),c("06","09"),c("06","11"),c("09","11"))
  adiv <- ggplot(Div_Table, aes(x = Date, y = ACE)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
    geom_errorbar(aes(ymin=ACE-ACE_se,ymax=ACE+ACE_se), width=.1,position=position_dodge(0.05)) +
    stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x="Dates",y="Diversity") #+ guides(color = FALSE)
  print(adiv)
  legend <- get_legend(adiv)
  adiv<- adiv + theme(legend.position="none")
  
  #Cycle
  my_comp <- list(c("Jour","Nuit"))
  bdiv <- ggplot(Div_Table, aes(x = Cycle, y = ACE)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
    geom_errorbar(aes(ymin=ACE-ACE_se,ymax=ACE+ACE_se), width=.1,position=position_dodge(0.05)) +
    stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x="Cycles",y="Diversity") #+ guides(color = FALSE)
  print(bdiv)
  bdiv<- bdiv + theme(legend.position="none")
  #Fraction
  my_comp <- list(c("Petite","Grande"))
  cdiv <- ggplot(Div_Table, aes(x = Fraction, y = ACE)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
    geom_errorbar(aes(ymin=ACE-ACE_se,ymax=ACE+ACE_se), width=.1,position=position_dodge(0.05)) +
    stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x="Fractions",y="Diversity") #+ guides(color = FALSE)
  print(cdiv)
  cdiv<- cdiv + theme(legend.position="none")
  #Zone
  my_comp <- list(c("Oxique","Anoxique"))
  ddiv <- ggplot(Div_Table, aes(x = Zone, y = ACE)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
    geom_errorbar(aes(ymin=ACE-ACE_se,ymax=ACE+ACE_se), width=.1,position=position_dodge(0.05)) + 
    stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=90, size=10, hjust = 0.5, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x="Zones",y="Diversity") #+ guides(color = FALSE)
  print(ddiv)
  ddiv<- ddiv + theme(legend.position="none")
  #Coplot
  fdiv<-ggplot() + theme(panel.background = element_rect(fill="white",colour = "white", size = 0.5, linetype = "solid"))
  pdf("Diversity/Pool/Analyse-ACE-Rare.pdf",width = 16.00,height = 9.00)
  Diversity_plot <- plot_grid(adiv,bdiv, fdiv, fdiv,fdiv,legend,cdiv,ddiv, labels=c("A","B", "", "", "", "", "C","D", ""), ncol = 3, nrow = 3, rel_widths = c(3,3,1.5), rel_heights = c(5,0.1,5))
  print(Diversity_plot)
  dev.off()
  
  
    # Obs -------------------------------------------------------------------
  for (i in colnames(chao_pool)) { Div_Table[i,"Obs"] <- chao_pool["S.obs",i] }
  #Date
  my_comp <- list(c("04","06"),c("04","09"),c("04","11"),c("06","09"),c("06","11"),c("09","11"))
  adiv <- ggplot(Div_Table, aes(x = Date, y = Obs)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
    stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x="Dates",y="Diversity") #+ guides(color = FALSE)
  print(adiv)
  legend <- get_legend(adiv)
  adiv<- adiv + theme(legend.position="none")
  
  #Cycle
  my_comp <- list(c("Jour","Nuit"))
  bdiv <- ggplot(Div_Table, aes(x = Cycle, y = Obs)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
    stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x="Cycles",y="Diversity") #+ guides(color = FALSE)
  print(bdiv)
  bdiv<- bdiv + theme(legend.position="none")
  #Fraction
  my_comp <- list(c("Petite","Grande"))
  cdiv <- ggplot(Div_Table, aes(x = Fraction, y = Obs)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
    stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x="Fractions",y="Diversity") #+ guides(color = FALSE)
  print(cdiv)
  cdiv<- cdiv + theme(legend.position="none")
  #Zone
  my_comp <- list(c("Oxique","Anoxique"))
  ddiv <- ggplot(Div_Table, aes(x = Zone, y = Obs)) + geom_boxplot() + geom_point(aes(color = Amplicon), size = 2.5) + 
    stat_compare_means(method="wilcox.test", paired = FALSE, label = "p.signif", comparisons = my_comp) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=0, size=10, hjust = 0.5, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x="Zones",y="Diversity") #+ guides(color = FALSE)
  print(ddiv)
  ddiv<- ddiv + theme(legend.position="none")
  #Coplot
  fdiv<-ggplot() + theme(panel.background = element_rect(fill="white",colour = "white", size = 0.5, linetype = "solid"))
  pdf("Diversity/Pool/Analyse-Obs-Rare.pdf",width = 16.00,height = 9.00)
  Diversity_plot <- plot_grid(adiv,bdiv, fdiv, fdiv,fdiv,legend,cdiv,ddiv, labels=c("A","B", "", "", "", "", "C","D", ""), ncol = 3, nrow = 3, rel_widths = c(3,3,1.5), rel_heights = c(5,0.1,5))
  print(Diversity_plot)
  dev.off()
  
  
  
  }