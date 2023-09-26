# !/usr/bin/env Rscript
#                                 _       _
#     /\                         (_)     | |
#    /  \   _ __ ___   ___  _ __  _  ___ | |_
#   / /\ \ | '_ ` _ \ / _ \| '_ \| |/ _ \| __|
#  / ____ \| | | | | | (_) | | | | | (_) | |_
# /_/    \_\_| |_| |_|\___/|_| |_| |\___/ \__|
#                               _/ |
#                              |__/
# 09/06/2022
#
# Script Metatrans_Table analyze
Cstack_info()["size"]
# Set input and output -----------------------------------------------------------
# Detect R or Rstudio
se <- Sys.getenv()
if (length(se[grepl("rstudio",se,ignore.case=TRUE)]) > 0 ) { inputmode <- TRUE }
if (length(se[grepl("rstudio",se,ignore.case=TRUE)]) == 0 ) { inputmode <- FALSE }
#
# Input argument if using Rstudio
if (inputmode == TRUE) {
input <- "main_table.mapping.unique.raw.noHuman.noConta.noMetazoa.annot.tsv"
output <- "V4-unified-correct-paired-out-compo"
rarefaction <- FALSE
}
#
# Input argument if using R
args = commandArgs(trailingOnly=TRUE)
if ( inputmode == FALSE ) {
  input <- args[1]
  output <- args[2]
}
print(input)
print(output)

#
# Import package -----------------------------------------------------------
pkg <- c("ggplot2","dplyr","tidyr","cowplot","FactoMineR","factoextra","reshape2","varhandle","ggrepel","ggpubr","ggsci","scales","hrbrthemes","data.table","svglite","treemap", "VennDiagram","stringr","paletteer","elementalist","gtools","vegan")
lapply(pkg, require, character.only = TRUE)

# Set directory, create file result -----------------------------------------------------------
if (inputmode == TRUE) {
  current <- dirname(rstudioapi::getSourceEditorContext()$path)
  setwd(current)
  output <- paste(current,"../dataDADA2/result",output,"Metatrans", sep = "/")
}
if (inputmode == FALSE) {
  output <- paste("../dataDADA2/result",output,"Metatrans", sep = "/")
}
# Set working directory
if (dir.exists(output) == FALSE) {dir.create(output)}
if (dir.exists(output) == TRUE) { setwd(output) }
# Create result files
dir.create("Duplicat")
dir.create("Betadisp")

# Theme unique Dark perso -------------------------------------------------------
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
                  legend.title = element_text(size = 12, face="bold"),
                  panel.background = element_rect_round(radius = unit(0.5, "cm"),fill = "white", color = NULL,size = 1),
                  #legend.position = "none",
                  #strip.background = element_rect(fill=NA,colour=NA,size=NA,linetype = NULL),
                  strip.text = element_text(size = 12,color="black",face="bold",vjust=.5,hjust=.5),
                  #panel.background = element_rect(fill = "white", color = NULL),
                  panel.border = element_blank(),
                  panel.grid = element_line(color = "#252525"),
                  panel.grid.major = element_line(color = "grey86"),
                  panel.grid.minor = element_line(color = "grey86"),
                  plot.background = element_rect(fill = "white", colour = "white", linetype = 0)))
  ret
}
# Input Tables ---------------------------------------------------------
tableVinput <- fread(file = paste("../../../../rawdata",input,sep = "/"), sep = "\t")
# Prepare data inf --------------------------------------------------------
infdataini <- read.table(file = "../../../../rawdata/metadata_metaT", sep = ";", header = TRUE)
infdataini$Echantillon <- sapply(strsplit(as.character(infdataini$Echantillon), "_"), function(x) x[[2]])
infdataini$Rep <- "_1"

infdataini$Variable <- paste0(infdataini$Ref..collab,infdataini$Rep)
infdataini <- infdataini %>% select(-"Rep")
infdataini <- separate(infdataini, Variable, c("Condition","Periods","Replicat"), sep = "_")
infdataini[,"Replicat"][infdataini[,"Replicat"] == "01"] <- "1"
infdataini[,"Replicat"][infdataini[,"Replicat"] == "02"] <- "2"
infdataini$Ref
## Prepare sample_df
infdataini <- separate(infdataini, Condition, c("empty","Cycle","Zone","Fraction"),sep = "")
infdataini <- infdataini %>% select(-"empty")
infdataini<-as.data.frame(infdataini)
col <- c("Echantillon","Ref..collab","Cycle","Zone","Fraction","Periods","Replicat")
samples_df <- infdataini[,col]
samples_df$Condition <- paste(samples_df$Cycle,samples_df$Zone,samples_df$Fraction, sep = "")
for (i in row.names(samples_df)) {
if (samples_df[i,"Cycle"] == "J") { samples_df[i,"Cycle"] <- "Day"}
if (samples_df[i,"Cycle"] == "N") { samples_df[i,"Cycle"] <- "Night"}
if (samples_df[i,"Zone"] == "O") { samples_df[i,"Zone"] <- "Mixolimnion"}
if (samples_df[i,"Zone"] == "A") { samples_df[i,"Zone"] <- "Monimolimnion"}
if (samples_df[i,"Fraction"] == "G") { samples_df[i,"Fraction"] <- "Large"}
if (samples_df[i,"Fraction"] == "P") { samples_df[i,"Fraction"] <- "Small"}}
rm(infdataini)

# MetaT table -------------------------------------------------------------
  tableVinput_abundance <- tableVinput
  dt <- setDF(tableVinput_abundance,rownames = tableVinput_abundance$unigene) %>% select(-unigene, -PfamDom, -ko, -ko_BestHitNotSignif, -GOterms)
  dt_PA <- dt
  dt_PA[dt_PA > 0 ] <- 1
  MetaT_Summary <- as.data.frame(colSums(dt_PA))
  rm(dt_PA)
  colnames(MetaT_Summary) <- "unigene_raw_count"
  MetaT_Summary$Ref..collab <- row.names(MetaT_Summary)
  MetaT_Summaryx <- merge(MetaT_Summary,samples_df,by="Ref..collab")
  rm(MetaT_Summary)
# Normalization TPM -------------------------------------------------------
  if (exists("rarefaction") == TRUE) { if ( rarefaction == TRUE ) {
  tableVinput_abundance_normalize <- setDF(tableVinput_abundance,rownames = tableVinput_abundance$unigene) %>% select(-unigene, -PfamDom, -ko, -ko_BestHitNotSignif, -GOterms)
  tableVinput_abundance_normalize["scalingFactor",] <- colSums(tableVinput_abundance_normalize) / 1000000
  for (i in colnames(tableVinput_abundance_normalize)) {
    tableVinput_abundance_normalize[1:nrow(tableVinput_abundance_normalize)-1,i] <- tableVinput_abundance_normalize[1:nrow(tableVinput_abundance_normalize)-1,i] / tableVinput_abundance_normalize[nrow(tableVinput_abundance_normalize),i] }
  tableVinput_abundance_normalize <- tableVinput_abundance_normalize[1:nrow(tableVinput_abundance_normalize)-1,]
  Unigene_table_c <- tableVinput_abundance_normalize
  }}else{
  Unigene_table_c <- setDF(tableVinput_abundance,rownames = tableVinput_abundance$unigene) %>% select(-unigene, -PfamDom, -ko, -ko_BestHitNotSignif, -GOterms)
  }

#
# Analyse Duplicat --------------------------------------------------------
  # Abundance ---------------------------------------------------------------
    # NMDS --------------------------------------------------------------------
    ## Prepare table
  dt <- Unigene_table_c
  dist <- vegdist(t(dt), method ="bray")
  res.nmds <- metaMDS(dist, parralel = 10)
  coord <- scores(res.nmds, display = "sites")
  data <- merge(y = coord, x = MetaT_Summaryx, by.y = "row.names", by.x = "Ref..collab")
    # Plot NMDS ----------------------------------------------------------------
  svglite("Duplicat/NMDS_abundance_byPeriods.svg",width = 8.00,height = 8.00)
  i <- ggplot(data, aes(y = `NMDS2`, x = `NMDS1`)) + geom_point(aes(color = Condition, shape = Periods), size = 2) + 
    geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
    geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
    facet_wrap(~Periods, strip.position="bottom") + theme_bw() +
    geom_label_repel(aes(label = Echantillon,fill = Condition), color = 'white',size = 2.5, segment.size = 0,segment.color = "black", alpha = 0.8) +
    geom_polygon(aes(color = Condition)) + theme_unique_art() +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10),
          strip.placement = "outside") +
    annotate("text", x = min(scores(res.nmds, display = "sites")), y = max(scores(res.nmds, display = "sites"))+1.2, label = paste0("stress: ", format(res.nmds$stress, digits = 4)), hjust = 0) +
    guides(color = "none") + scale_fill_paletteer_d("ggsci::default_jco") + scale_color_paletteer_d("ggsci::default_jco")
  print(i)
  dev.off()
  
  svglite("Duplicat/NMDS_abundance_byPeriodsXZoneXfraction.svg",width = 8.00,height = 8.00)
  j <- ggplot(data, aes(y = `NMDS2`, x = `NMDS1`)) + geom_point(aes(color = Zone,shape=Fraction), size = 2) + 
    geom_hline(yintercept=0, linetype=2, color = "black", size=0.2) +
    geom_vline(xintercept=0, linetype=2, color = "black", size=0.2) +
    geom_label_repel(aes(label = Echantillon,fill = Zone), color = 'white',size = 2.5, segment.size = 0,segment.color = "black", alpha = 0.8) +
    stat_ellipse(aes(alpha = Zone,linetype = Zone,color=Zone),geom = "polygon",type = "norm") +
    theme_unique_art() +
    facet_wrap(.~Periods) +
    scale_linetype_manual(values = c(2,2,2,2)) +
    scale_alpha_manual(values = c(0.1,0.1,0.1,0.1)) +
    scale_fill_manual(values = c("#7375b5FF","#8ca252FF","#f6ae34FF","#ad494aFF")) +
    scale_color_manual(values = c("#7375b5FF","#8ca252FF","#f6ae34FF","#ad494aFF")) +
    theme(axis.title = element_text(face="bold", size=12), 
          axis.text.x = element_text(angle=0, size=10, hjust = 1, vjust=0.5), 
          title = element_text(face="bold", size=14),
          legend.title = element_text(face="bold"),
          legend.position = "right",
          legend.text = element_text(size=10),
          strip.placement = "outside") +
    annotate("text", x = min(scores(res.nmds, display = "sites")), y = max(scores(res.nmds, display = "sites"))+1.2, label = paste0("stress: ", format(res.nmds$stress, digits = 4)), hjust = 0) +
    guides(color = "none",fill = guide_legend(title = "Zone",override.aes = aes(label = "")))
  print(j)
  dev.off()
#
  # Clear env ---------------------------------------------------------------
rm(dt)
rm(tableVinput)
rm(samples_df)
# Remove duplicat ---------------------------------------------------------
    tuniq <- unique(MetaT_Summaryx$Condition)
    puniq <- unique(MetaT_Summaryx$Periods)
    MetaT_Summary <- tibble()
    for (i in tuniq) {
      for (j in puniq) {
        k <- MetaT_Summaryx %>% filter(Condition == i) %>% filter(Periods == j)
        if (nrow(k) > 1) { interest <- k %>% filter(Replicat == 1) }
        else { interest <- k }
        MetaT_Summary <- rbind(MetaT_Summary, interest)
      }
    }
    oldcol <- colnames(Unigene_table_c)
    w <- 0
    for (i in oldcol) { f <- as.character(MetaT_Summary %>% filter(Ref..collab == i) %>% select(Echantillon))
      if (f != "character(0)") { 
        if (w == 0) { newcol <- f 
          selectcol <- i }
        else { newcol <- c(newcol,f)
          selectcol <- c(selectcol,i)}}
    w <- w+1 }
    Unigene_table_c <- Unigene_table_c %>% select(all_of(selectcol))
    colnames(Unigene_table_c) <- newcol
#

# Sample Ordination : beta disp -------------------------------------------------------
    dt <- t(Unigene_table_c)
    dt <- merge(dt, MetaT_Summary %>% select(Echantillon, Cycle, Periods, Fraction, Zone), by.x = "row.names", by.y = "Echantillon")
    dtx <- dt %>% select(-Cycle,-Periods, -Fraction, -Zone, -Row.names)
    dist <- vegdist(dtx, method ="bray")
    NMDS_tab <- metaMDS(dist, parralel = 10)
    # Plots
    Condition <- dt$Zone
    plot(NMDS_tab, display = "sites")
    ordiellipse(NMDS_tab,groups = Condition,label = T)
    orditorp(NMDS_tab,display="sites")
    # ggplot
    scores(NMDS_tab, display = "sites") %>%
      cbind(dt %>% select(Row.names, Cycle, Periods, Fraction, Zone)) %>%
      ggplot(aes(x = NMDS1, y = NMDS2)) +
      geom_point(aes(color = Zone)) +
      facet_grid(~Fraction) +
      stat_ellipse(geom = "polygon", aes(group = Zone, color = Zone), alpha = 0.1) +
      annotate("text", x = min(scores(NMDS_tab,display = "sites")), y = max(scores(NMDS_tab,display = "sites"))+1.2, label = paste0("stress: ", format(NMDS_tab$stress, digits = 4)), hjust = 0) +
      theme_bw() + scale_color_jco() # + geom_text(aes(label = Amplicon))
    ggsave("Betadisp/NMDS.svg", device = "svg", width = 5, height = 5)
    # Adonis
    capture.output(adonis2(formula = dist~dt$Zone+dt$Fraction+dt$Cycle+dt$Periods,distance = "bray", parallel = 10),file="Betadisp/Adonis.txt")
    # Anosim
    capture.output(anosim(dist, dt$Zone, permutations = 9999, parallel = 10),file="Betadisp/anosim_transcrit_Zone.txt")
    capture.output(anosim(dist, dt$Fraction, permutations = 9999, parallel = 10),file="Betadisp/anosim_transcrit_Fraction.txt")
    capture.output(anosim(dist, dt$Periods, permutations = 9999, parallel = 10),file="Betadisp/anosim_transcrit_Periods.txt")
    capture.output(anosim(dist, dt$Cycle, permutations = 9999, parallel = 10),file="Betadisp/anosim_transcrit_Cycle.txt")
    # Betadisper
    mod <- betadisper(dist,dt$Zone)
    labs <- paste("Dimension", 1:4, "(", 
                  round(100*mod$eig / sum(mod$eig), 2), "%)")
    svglite("Betadisp/Zone_PCoA.svg", width = 5.00, height = 5.00)
    plot(mod, cex=1, pch=15:17, cex.lab=1.25,
         xlab=labs[1], ylab=labs[2], main = "Zone",
         hull=TRUE, ellipse=FALSE, conf=0.95, lwd=2)
    dev.off()
    # Periods
    mod <- betadisper(dist,dt$Periods)
    labs <- paste("Dimension", 1:4, "(", 
                  round(100*mod$eig / sum(mod$eig), 2), "%)")
    svglite("Betadisp/Periods_PCoA.svg", width = 5.00, height = 5.00)
    plot(mod, cex=1, pch=15:17, cex.lab=1.25,
         xlab=labs[1], ylab=labs[2], main = "Periods",
         hull=TRUE, ellipse=FALSE, conf=0.95, lwd=2)
    dev.off()
    # Fraction
    mod <- betadisper(dist,dt$Fraction)
    labs <- paste("Dimension", 1:4, "(", 
                  round(100*mod$eig / sum(mod$eig), 2), "%)")
    svglite("Betadisp/Fraction_PCoA.svg", width = 5.00, height = 5.00)
    plot(mod, cex=1, pch=15:17, cex.lab=1.25,
         xlab=labs[1], ylab=labs[2], main = "Fraction",
         hull=TRUE, ellipse=FALSE, conf=0.95, lwd=2)
    dev.off()
    # Cycle
    mod <- betadisper(dist,dt$Cycle)
    labs <- paste("Dimension", 1:4, "(", 
                  round(100*mod$eig / sum(mod$eig), 2), "%)")
    svglite("Betadisp/Cycle_PCoA.svg", width = 5.00, height = 5.00)
    plot(mod, cex=1, pch=15:17, cex.lab=1.25,
         xlab=labs[1], ylab=labs[2], main = "Cycle",
         hull=TRUE, ellipse=FALSE, conf=0.95, lwd=2)
    dev.off()
#
# Save IMG ----------------------------------------------------------------
if (inputmode == FALSE) {save.image("R_Duplicat_image.RData")}
