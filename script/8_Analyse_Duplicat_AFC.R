#                                 _       _
#     /\                         (_)     | |
#    /  \   _ __ ___   ___  _ __  _  ___ | |_
#   / /\ \ | '_ ` _ \ / _ \| '_ \| |/ _ \| __|
#  / ____ \| | | | | | (_) | | | | | (_) | |_
# /_/    \_\_| |_| |_|\___/|_| |_| |\___/ \__|
#                               _/ |
#                              |__/
# XX/XX/XXXX
#
# Script Composition

# Set directory, input, output and import packages -----------------------------------------------------------
#!/usr/bin/env Rscript
# Rscript 8_Analyse_Duplicat_AFC.R V4 ../dataPANAM/PANAM2/V4-result-095/OTU_distribution_tax.txt Analyse-Composition-Rarefy-V4-095-Vsearch
#output <- "Analyse-Composition-Rarefy-V4"
#input <- "../dataPANAM/PANAM2/V4-result-095/OTU_distribution_tax.txt"
#region <- "V4"
args = commandArgs(trailingOnly=TRUE)

if (length(args)==2) {
cat("Enter ribosomal region (V4 or V9) : ");
region <- readLines("stdin",n=1);
cat("You entered")
str(region);
cat( "\n" )}
if (length(args)==4) {
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

system("mkdir Figure-Sum")
system("mkdir AFC-Duplicat")
system("mkdir Composition")
system("mkdir HistOnly")
system("mkdir AFC-Parasite")
system("mkdir TableOnly")
system("mkdir Biplot")
system("mkdir Rarecurve")
system("mkdir Rarecurve/Nopool")
system("mkdir Rarecurve/Pool")
system("mkdir Diversity")
system("mkdir Diversity/Nopool")
system("mkdir Diversity/Pool")


# Input OTU Table ---------------------------------------------------------
tableVinput <- read.csv(file = input, sep = "\t")

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
  mc <- getOption("mc.cores", ifelse(max.cores, detectCores(), nCores))
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

# Prepare data inf --------------------------------------------------------
infdataini <- read.table(file = "../../rawdata/data-inf.txt", sep = "\t", header = FALSE, col.names = c("row","Id","Conditions","Technologies","Regions"))
infdataini <- infdataini %>% select(-"row")
infdataini$Rep <- rep("_1", each = nrow(infdataini))
infdataini$Variable <- paste(infdataini$Condition,infdataini$Rep, sep = "")
infdataini <- infdataini %>% select(-"Rep",-"Conditions")
infdataini <- separate(infdataini, Variable, c("Conditions","Dates","Replicats"), sep = "_")
x <- 0
for (i in infdataini$Replicats) { 
  x <- x+1
  if (i == "01") infdataini[x,6] <- "1"
  if (i == "02") infdataini[x,6] <- "2"
}
infdataini$OSTA <- rep("OSTA", each = nrow(infdataini))
infdataini$Variable <- paste(infdataini$Id,infdataini$OSTA, sep = "")
infdataini <- infdataini %>% select(-"Id",-"OSTA")
infdataini <- separate(infdataini, Variable, c("Cin","Variable"), sep = "_")
infdataini <- infdataini %>% select(-"Cin")

# FR sample (V9 in reality and not V4)
for (i in row.names(infdataini)) { if (infdataini[i,"Variable"] == "FROSTA") { infdataini[i,"Regions"] <- "V9"}}
for (i in row.names(infdataini)) { if (infdataini[i,"Variable"] == "FROSTA") { infdataini[i,"Technologies"] <- "Miseq"}}

# Prepare sample_df
infdataini[infdataini == "DJOG"] <- "D;J;O;G"
infdataini[infdataini == "DJOP"] <- "D;J;O;P"
infdataini[infdataini == "DJAG"] <- "D;J;A;G"
infdataini[infdataini == "DJAP"] <- "D;J;A;P"
infdataini[infdataini == "DNOG"] <- "D;N;O;G"
infdataini[infdataini == "DNOP"] <- "D;N;O;P"
infdataini[infdataini == "DNAG"] <- "D;N;A;G"
infdataini[infdataini == "DNAP"] <- "D;N;A;P"
infdataini<-as.data.frame(infdataini)
infdataini <- separate(infdataini, Conditions, c("ADN","Jour_Nuit","Oxique_Anoxique","Grande_Petite"),sep = ";")
pattern <- c("Variable","Technologies","Regions","ADN","Jour_Nuit","Oxique_Anoxique","Grande_Petite","Dates","Replicats")
samples_df <- infdataini[,all_of(pattern)]
samples_df$Conditions <- paste(samples_df$ADN,samples_df$Jour_Nuit,samples_df$Oxique_Anoxique,samples_df$Grande_Petite, sep = "")
colnames(samples_df) <- c("sample","Technologies","Regions","ADN","Cycle","Fraction-Oxygène","Fraction-Taille","Dates","Replicats","Conditions")
w <- 0
for (i in samples_df$Cycle) { w <- w +1
if (i == "J") { samples_df[w,"Cycle"] <- "Jour"}
if (i == "N") { samples_df[w,"Cycle"] <- "Nuit"}}
w <- 0
for (i in samples_df$`Fraction-Oxygène`) { w <- w +1
if (i == "O") { samples_df[w,"Fraction-Oxygène"] <- "Oxique"}
if (i == "A") { samples_df[w,"Fraction-Oxygène"] <- "Anoxique"}}
w <- 0
for (i in samples_df$`Fraction-Taille`) { w <- w +1
if (i == "G") { samples_df[w,"Fraction-Taille"] <- "Grande"}
if (i == "P") { samples_df[w,"Fraction-Taille"] <- "Petite"}}

# Select V4 or V9 ---------------------------------------------------------
samples_df<-filter(samples_df, Regions == region)

# Prepare  Object -------------------------------------------------
  # Rarefy ---------------------------------------------------------------------
pattern <- c(grep(pattern = "OSTA", colnames(tableVinput), value = FALSE, fixed = FALSE))
otu_mat <- tableVinput %>% select(OTU_Id,all_of(pattern))
row.names(otu_mat) <- otu_mat$OTU_Id
otu_mat <- otu_mat %>% select(-OTU_Id)
otu_mat[,] <- lapply(otu_mat[,], as.numeric)
otu_matx <- t(otu_mat)
sequence_mat_rare <- t(Rarefy(otu_matx)$otu.tab.rff)


  # PA-AB ------------------------------------------------------------------
otu_mat_rare<- sequence_mat_rare
pattern <- c(grep(pattern = "OSTA", colnames(otu_mat_rare), value = FALSE, fixed = FALSE))
otu_mat_rare[,all_of(pattern)][otu_mat_rare[,all_of(pattern)] != 0] <- 1

#Stat Rarefy
  #Sequence
  avRarefyS <- as.data.frame(colSums(otu_mat))
  apRarefyS <- as.data.frame(colSums(sequence_mat_rare))
  statRarefy <- cbind(avRarefyS,apRarefyS)
  #OTU
  otu_mat_uniq <- otu_mat
  otu_mat_uniq[otu_mat_uniq != 0] <- 1
  avRarefyO <- as.data.frame(colSums(otu_mat_uniq))
  apRarefyO <- as.data.frame(colSums(otu_mat_rare))
  statRarefy <- cbind(statRarefy,avRarefyO,apRarefyO)
  statRarefy["Total",]<-colSums(statRarefy)
  colnames(statRarefy) <- c("avRarefy-Sequence","apRarefy-Sequence","avRarefy-OTU","apRarefy-OTU")
  write.table(statRarefy, file = "TableOnly/StatRarefy_withDuplicat.txt", sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)
  
  
  # Rarecurve ---------------------------------------------------------------
  #rarecurveop <- readline(prompt="Should I calculate the rarefaction curves (yes or no) ? : ")
  if (length(args)==2) {
  cat("Should I calculate the rarefaction curves (yes or no) (It takes a lot of time and resources) ? : ");
  rarecurveop <- readLines("stdin",n=1);
  cat("You entered")
  str(rarecurveop);
  cat( "\n" )}
  if (length(args)==4) {
    rarecurveop <- args[4]
  }
  
  
  if (rarecurveop == "yes") {
  #rarecurve Rarefy
  rarecurve <- t(sequence_mat_rare)
  pdf("Rarecurve/Nopool/Rarecurve-Rarefy.pdf",width = 5.00,height = 5.00)
  curve <- quickRareCurve(rarecurve, col = "black", cex = 0.6)
  dev.off()
  #rarecurve all
  #rarecurve_all <- t(colSums(rarecurve))
  #rownames(rarecurve_all) <- "sum"
  #pdf("Rarecurve/Nopool/Rarecurve-all-Rarefy.pdf",width = 5.00,height = 5.00)
  #curve <- quickRareCurve(rarecurve_all, col = "black", cex = 0.6)
  #dev.off()
  #Color
  colraw <- as.data.frame(rownames(rarecurve))
  colnames(colraw) <- "sample"
  row.names(colraw) <- colraw[,"sample"]
  for (x in row.names(colraw)) {
    colraw[x,1:10] <- samples_df %>% filter(sample == x)}
  #Color Cycle
  colraw$color_cycle <- rep("black", time = nrow(colraw))
  for (w in rownames(colraw)) {
    if (colraw[w,"Cycle"] == "Jour") { colraw[w,"color_cycle"] <- "#D43F3AFF"}}
  pdf("Rarecurve/Nopool/Rarecurve-Cycle-Rarefy.pdf",width = 5.00,height = 5.00)
  curve_a <- quickRareCurve(rarecurve, col = colraw$color_cycle, cex = 0.6)
  legend("topleft", inset=.02, title="Day/Night cycle",
         c("Day","Night"), fill = c("#D43F3AFF","black"), horiz=TRUE, cex=0.6)
  print(curve_a)
  dev.off()
  #Color Oxygène
  colraw$color_Ox <- rep("black", time = nrow(colraw))
  for (w in rownames(colraw)) {
    if (colraw[w,"Fraction-Oxygène"] == "Oxique") { colraw[w,"color_Ox"] <- "#46B8DAFF"}}
  pdf("Rarecurve/Nopool/Rarecurve-Ox-Rarefy.pdf",width = 5.00,height = 5.00)
  curve_b <- quickRareCurve(rarecurve, col = colraw$color_Ox, cex = 0.6)
  legend("topleft", inset=.02, title="Fraction Ox/Anox",
         c("Oxique","Anoxique"), fill = c("#46B8DAFF","black"), horiz=TRUE, cex=0.6)
  print(curve_b)
  dev.off()
  #Color Taille
  colraw$color_Size <- rep("black", time = nrow(colraw))
  for (w in rownames(colraw)) {
    if (colraw[w,"Fraction-Taille"] == "Petite") { colraw[w,"color_Size"] <- "#5CB85CFF"}}
  pdf("Rarecurve/Nopool/Rarecurve-Size-Rarefy.pdf",width = 5.00,height = 5.00)
  curve_c <- quickRareCurve(rarecurve, col = colraw$color_Size, cex = 0.6)
  legend("topleft", inset=.02, title="Fraction Taille",
         c("Petite","Grande"), fill = c("#5CB85CFF","black"), horiz=TRUE, cex=0.6)
  print(curve_c)
  dev.off()
  #Color Replicat
  colraw$color_Rep <- rep("black", time = nrow(colraw))
  for (w in rownames(colraw)) {
    if (colraw[w,"Replicats"] == 1) { colraw[w,"color_Rep"] <- "#EEA236CC"}}
  pdf("Nopool/Rarecurve-Rep-Rarefy.pdf",width = 5.00,height = 5.00)
  curve_d <- quickRareCurve(rarecurve, col = colraw$color_Rep, cex = 0.6)
  legend("topleft", inset=.02, title="Replicat",
         c("Rep1","Rep2"), fill = c("#EEA236CC","black"), horiz=TRUE, cex=0.6)
  print(curve_d)
  dev.off()
  
  #rarecurve no-rarefy
  rawcurve <- otu_matx
  pdf("Rarecurve/Nopool/Rarecurve-Raw.pdf",width = 5.00,height = 5.00)
  curve <- quickRareCurve(rawcurve, col = "black", cex = 0.6)
  dev.off()
  #rarecurve all
  #rawcurve_all <- t(colSums(rawcurve))
  #rownames(rawcurve_all) <- "sum"
  #pdf("Rarecurve/Nopool/Rarecurve-all-Raw.pdf",width = 5.00,height = 5.00)
  #curve <- quickRareCurve(rawcurve_all, col = "black", cex = 0.6)
  #dev.off()
  #Color
  colraw <- as.data.frame(rownames(rawcurve))
  colnames(colraw) <- "sample"
  w <- 0
  for (x in colraw[,"sample"]) { w <- w+1
  colraw[w,1:10] <- subset(samples_df, sample == x)}
  #Color Cycle
  colraw$color_cycle <- rep("black", time = nrow(colraw))
  for (w in rownames(colraw)) {
    if (colraw[w,"Cycle"] == "Jour") { colraw[w,"color_cycle"] <- "#D43F3AFF"}}
  pdf("Rarecurve/Nopool/Rarecurve-Cycle-Rawdata.pdf", width = 5.00, height = 5.00)
  curve_a <- quickRareCurve(rawcurve, col = colraw$color_cycle, cex = 0.6)
  legend("topleft", inset=.02, title="Day/Night cycle",
         c("Day","Night"), fill = c("#D43F3AFF","black"), horiz=TRUE, cex=0.6)
  print(curve_a)
  dev.off()
  #Color Oxygène
  colraw$color_Ox <- rep("black", time = nrow(colraw))
  for (w in rownames(colraw)) {
    if (colraw[w,"Fraction-Oxygène"] == "Oxique") { colraw[w,"color_Ox"] <- "#46B8DAFF"}}
  pdf("Rarecurve/Nopool/Rarecurve-Ox-Rawdata.pdf", width = 5.00, height = 5.00)
  curve_b <- quickRareCurve(rawcurve, col = colraw$color_Ox, cex = 0.6)
  legend("topleft", inset=.02, title="Fraction Ox/Anox",
         c("Oxique","Anoxique"), fill = c("#46B8DAFF","black"), horiz=TRUE, cex=0.6)
  print(curve_b)
  dev.off()
  #Color Taille
  colraw$color_Size <- rep("black", time = nrow(colraw))
  for (w in rownames(colraw)) {
    if (colraw[w,"Fraction-Taille"] == "Petite") { colraw[w,"color_Size"] <- "#5CB85CFF"}}
  pdf("Rarecurve/Nopool/Rarecurve-Size-Rawdata.pdf", width = 5.00, height = 5.00)
  curve_c <- quickRareCurve(rawcurve, col = colraw$color_Size, cex = 0.6)
  legend("topleft", inset=.02, title="Fraction Taille",
         c("Petite","Grande"), fill = c("#5CB85CFF","black"), horiz=TRUE, cex=0.6)
  print(curve_c)
  dev.off()
  #Color Replicat
  colraw$color_Rep <- rep("black", time = nrow(colraw))
  for (w in rownames(colraw)) {
    if (colraw[w,"Replicats"] == 1) { colraw[w,"color_Rep"] <- "#EEA236CC"}}
  pdf("Rarecurve/Nopool/Rarecurve-Rep-Rawdata.pdf",width = 5.00,height = 5.00)
  curve_d <- quickRareCurve(rawcurve, col = colraw$color_Rep, cex = 0.6)
  legend("topleft", inset=.02, title="Replicat",
         c("Rep1","Rep2"), fill = c("#EEA236CC","black"), horiz=TRUE, cex=0.6)
  print(curve_d)
  dev.off()
  }
  
# OTU ---------------------------------------------------------------------
  # AFC Plot ----------------------------------------------------------------
  dt <- as.data.frame(otu_mat_rare)
  res.ca <- CA (dt, graph = FALSE,ncp = 2 )
  pdf("AFC-Duplicat/AFC_OTU.pdf",width = 16.00,height = 9.00)
  fviz_ca_col(res.ca, repel = TRUE, col.col = "lightgreen") + theme_unique_dark()
  dev.off()
  p <- get_ca_col(res.ca)
  coord <- p$coord
  coord <- as.data.frame(coord)
  coord$sample <- row.names(coord)
  data <- merge(x = coord, y = samples_df, by= "sample")
  fviz_screeplot (res.ca, addlabels = TRUE, ylim = c(0, 50))
  
  Xseq <- fviz_screeplot (res.ca, addlabels = TRUE, ylim = c(0, 50))
  Xseq <- Xseq$data
  Dim1Seq <- paste("Dim 1 [",round(Xseq %>% filter(dim == 1) %>% select(eig),1),"%]", sep = "")
  Dim2Seq <- paste("Dim 2 [",round(Xseq %>% filter(dim == 2) %>% select(eig),1),"%]", sep = "")
  
  # Plot AFC ----------------------------------------------------------------
    # Overview* -------------------------------------------------------------------
  pdf("AFC-Duplicat/Analyse-OTU.pdf",width = 16.00,height = 9.00)
  a <- ggplot(data, aes(y = `Dim 1`, x = `Dim 2`)) + geom_point(aes(color = Conditions, shape = Dates), size = 2) + 
    geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
    geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
    facet_grid(~Dates, switch = "x") + theme_bw() +
    geom_label_repel(aes(label = sample,fill = Conditions), color = 'white',size = 2.5, segment.size = 0,segment.color = "black", alpha = 0.8) +
    geom_polygon(aes(color = Conditions)) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=90, size=10, hjust = 1, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(title="Analyse Factorielle des Correspondances",x=Dim1Seq,y=Dim2Seq) +
    guides(color = FALSE)
  print(a)
  dev.off()
  
  svglite("AFC-Duplicat/Analyse-OTU-04.svg",width = 8.00,height = 6.00)
  a <- ggplot(data %>% filter(Dates == "04"), aes(y = `Dim 1`, x = `Dim 2`)) + geom_point(aes(color = Conditions, shape = Dates), size = 3) + 
    geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
    geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
    geom_label_repel(aes(label = sample,fill = Conditions), color = 'white',size = 4, segment.size = 0.5,segment.color = "black", alpha = 0.8) +
    #geom_polygon(aes(color = Conditions)) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=90, size=10, hjust = 1, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x=Dim1Seq,y=Dim2Seq) +
    guides(shape = FALSE, color = FALSE) + labs(color = "Conditions")
  print(a)
  dev.off()
  
  svglite("AFC-Duplicat/Analyse-OTU-06.svg",width = 8.00,height = 6.00)
  a <- ggplot(data %>% filter(Dates == "06"), aes(y = `Dim 1`, x = `Dim 2`)) + geom_point(aes(color = Conditions, shape = Dates), size = 3) + 
    geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
    geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
    geom_label_repel(aes(label = sample,fill = Conditions), color = 'white',size = 4, segment.size = 0.5,segment.color = "black", alpha = 0.8) +
    #geom_polygon(aes(color = Conditions)) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=90, size=10, hjust = 1, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x=Dim1Seq,y=Dim2Seq) +
    guides(shape = FALSE, color = FALSE) + labs(color = "Conditions")
  print(a)
  dev.off()
  
  svglite("AFC-Duplicat/Analyse-OTU-09.svg",width = 8.00,height = 6.00)
  a <- ggplot(data %>% filter(Dates == "09"), aes(y = `Dim 1`, x = `Dim 2`)) + geom_point(aes(color = Conditions, shape = Dates), size = 3) + 
    geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
    geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
    geom_label_repel(aes(label = sample,fill = Conditions), color = 'white',size = 4, segment.size = 0.5,segment.color = "black", alpha = 0.8) +
    #geom_polygon(aes(color = Conditions)) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=90, size=10, hjust = 1, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x=Dim1Seq,y=Dim2Seq) +
    guides(shape = FALSE, color = FALSE) + labs(color = "Conditions")
  print(a)
  dev.off()
  
  svglite("AFC-Duplicat/Analyse-OTU-11.svg",width = 8.00,height = 6.00)
  a <- ggplot(data %>% filter(Dates == "11"), aes(y = `Dim 1`, x = `Dim 2`)) + geom_point(aes(color = Conditions, shape = Dates), size = 3) + 
    geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
    geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
    geom_label_repel(aes(label = sample,fill = Conditions), color = 'white',size = 4, segment.size = 0.5,segment.color = "black", alpha = 0.8) +
    #geom_polygon(aes(color = Conditions)) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=90, size=10, hjust = 1, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x=Dim1Seq,y=Dim2Seq) +
    guides(shape = FALSE, color = FALSE) + labs(color = "Conditions")
  print(a)
  dev.off()
  

# Sequence ---------------------------------------------------------------------
  # AFC Plot ----------------------------------------------------------------
  dt <- as.data.frame(sequence_mat_rare)
  res.ca <- CA (dt, graph = FALSE,ncp = 2 )
  pdf("AFC-Duplicat/AFC_Seq.pdf",width = 16.00,height = 9.00)
  fviz_ca_col(res.ca, repel = TRUE, col.col = "lightgreen") + theme_unique_dark()
  dev.off()
  p <- get_ca_col(res.ca)
  coord <- p$coord
  coord <- as.data.frame(coord)
  coord$sample <- row.names(coord)
  data <- merge(x = coord, y = samples_df, by= "sample")
  fviz_screeplot (res.ca, addlabels = TRUE, ylim = c(0, 50))
  
  Xseq <- fviz_screeplot (res.ca, addlabels = TRUE, ylim = c(0, 50))
  Xseq <- Xseq$data
  Dim1Seq <- paste("Dim 1 [",round(Xseq %>% filter(dim == 1) %>% select(eig),1),"%]", sep = "")
  Dim2Seq <- paste("Dim 2 [",round(Xseq %>% filter(dim == 2) %>% select(eig),1),"%]", sep = "")
  
  # Plot AFC ----------------------------------------------------------------
    # Overview* -------------------------------------------------------------------
  pdf("AFC-Duplicat/Analyse-Seq.pdf",width = 16.00,height = 9.00)
  a <- ggplot(data, aes(y = `Dim 1`, x = `Dim 2`)) + geom_point(aes(color = Conditions, shape = Dates), size = 2) + 
    geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
    geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
    facet_grid(~Dates, switch = "x") + theme_bw() +
    geom_label_repel(aes(label = sample,fill = Conditions), color = 'white',size = 2.5, segment.size = 0,segment.color = "black", alpha = 0.8) +
    geom_polygon(aes(color = Conditions)) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=90, size=10, hjust = 1, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(title="Analyse Factorielle des Correspondances",x=Dim1Seq,y=Dim2Seq) +
    guides(color = FALSE)
  print(a)
  dev.off()
  
  svglite("AFC-Duplicat/Analyse-Seq-04.svg",width = 8.00,height = 6.00)
  a <- ggplot(data %>% filter(Dates == "04"), aes(y = `Dim 1`, x = `Dim 2`)) + geom_point(aes(color = Conditions, shape = Dates), size = 3) + 
    geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
    geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
    geom_label_repel(aes(label = sample,fill = Conditions), color = 'white',size = 4, segment.size = 0.5,segment.color = "black", alpha = 0.8) +
    #geom_polygon(aes(color = Conditions)) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=90, size=10, hjust = 1, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x=Dim1Seq,y=Dim2Seq) +
    guides(shape = FALSE, color = FALSE) + labs(color = "Conditions")
  print(a)
  dev.off()
  
  svglite("AFC-Duplicat/Analyse-Seq-06.svg",width = 8.00,height = 6.00)
  a <- ggplot(data %>% filter(Dates == "06"), aes(y = `Dim 1`, x = `Dim 2`)) + geom_point(aes(color = Conditions, shape = Dates), size = 3) + 
    geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
    geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
    geom_label_repel(aes(label = sample,fill = Conditions), color = 'white',size = 4, segment.size = 0.5,segment.color = "black", alpha = 0.8) +
    #geom_polygon(aes(color = Conditions)) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=90, size=10, hjust = 1, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x=Dim1Seq,y=Dim2Seq) +
    guides(shape = FALSE, color = FALSE) + labs(color = "Conditions")
  print(a)
  dev.off()
  
  svglite("AFC-Duplicat/Analyse-Seq-09.svg",width = 8.00,height = 6.00)
  a <- ggplot(data %>% filter(Dates == "09"), aes(y = `Dim 1`, x = `Dim 2`)) + geom_point(aes(color = Conditions, shape = Dates), size = 3) + 
    geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
    geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
    geom_label_repel(aes(label = sample,fill = Conditions), color = 'white',size = 4, segment.size = 0.5,segment.color = "black", alpha = 0.8) +
    #geom_polygon(aes(color = Conditions)) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=90, size=10, hjust = 1, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x=Dim1Seq,y=Dim2Seq) +
    guides(shape = FALSE, color = FALSE) + labs(color = "Conditions")
  print(a)
  dev.off()
  
  svglite("AFC-Duplicat/Analyse-Seq-11.svg",width = 8.00,height = 6.00)
  a <- ggplot(data %>% filter(Dates == "11"), aes(y = `Dim 1`, x = `Dim 2`)) + geom_point(aes(color = Conditions, shape = Dates), size = 3) + 
    geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
    geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
    geom_label_repel(aes(label = sample,fill = Conditions), color = 'white',size = 4, segment.size = 0.5,segment.color = "black", alpha = 0.8) +
    #geom_polygon(aes(color = Conditions)) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=90, size=10, hjust = 1, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10)) +
    labs(x=Dim1Seq,y=Dim2Seq) +
    guides(shape = FALSE, color = FALSE) + labs(color = "Conditions")
  print(a)
  dev.off()
  
  
