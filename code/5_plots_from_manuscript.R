# generate plots for Gordon et al, Nature Neuroscience, 2021

# setup -----------------------------------------------------------------------------------------------------------
options(stringsAsFactors = FALSE)
library(RColorBrewer)
library(gplots)
library(WGCNA)
library(tidyverse)
library(patchwork)
library(readxl)
library(viridis)
library(ggrepel)
library(tidytext)
library(pzfx)
library(multcomp)
library(ggsignif)
library(ComplexHeatmap)


outputFolder <- file.path("output")

tables_folder <- file.path(outputFolder, "tables")
dir.create(outputFolder, showWarnings = F)
dir.create(tables_folder, showWarnings = F)

pub_theme <- theme_classic(base_size = 6) +
  theme(axis.text = element_text(size = 6, color = "black"),
        strip.background = element_rect(colour = NA),
        strip.text = element_text(size = 6),
        plot.title = element_text(size = 6),
        axis.line = element_line(colour = 'black', size = 1/2.13), 
        axis.ticks = element_line(colour = "black", size = 1/2.13),
        axis.title.x = element_text(vjust=-0.75),
        legend.text = element_text(size = 6))
theme_set(pub_theme)


#load brain span data
load(file.path("data","brainSpan_expr_subset.rdata"))

load(file.path("output/","processed_data.rdata"))
load(file.path("data","gene_annotation.Rdata"))

#functions ---------------------------------------------------------
OR <- function(q,k,m,t) {
  q #<-  ## Intersection of test list and reference list, aka number of white balls drawn
  m #<-  ## All genes in reference list, aka number of draws
  k #<-  ## All genes in test list, aka number white balls
  t #<-  ## Total number of genes assessed, aka black plus white balls
  
  fisher.out <- fisher.test(matrix(c(q, k-q, m-q, t-m-k+q), 2, 2),conf.int=TRUE)
  OR <- fisher.out$estimate
  pval <- fisher.out$p.value
  upCI <- fisher.out$conf.int[1]
  downCI <- fisher.out$conf.int[2]
  
  output <- c(OR,pval,upCI,downCI)
  names(output) <- c("OR","Fisher p","-95%CI","+95%CI")
  return(output)
}

ORA <- function(testpath,refpath,testbackground,refbackground) {
  q <- length(intersect(testpath,refpath)) ## overlapped pathway size
  k <- length(intersect(refpath,testbackground))  ## input gene set
  m <- length(intersect(testpath,refbackground)) ## input module
  t <- length(intersect(testbackground,refbackground)) ## Total assessed background (intersect reference and test backgrounds)
  
  empvals <- OR(q,k,m,t)
  
  tmpnames <- names(empvals)
  empvals <- as.character(c(empvals,q,k,m,t,100*signif(q/k,3)))
  names(empvals) <- c(tmpnames,"Overlap","Reference List","Input List","Background","% List Overlap")
  return(empvals)
}

geneSetEnrichment <- function(test_matrix, ref_matrix){
  ORA.P = ORA.OR = matrix(NA,nrow = dim(test_matrix)[2], ncol = dim(ref_matrix)[2]);
  
  testbackground = refbackground = rownames(test_matrix)
  
  for (i in 1:dim(ref_matrix)[2]) {
    for (j in 1:dim(test_matrix)[2]) {
      testGenes <- rownames(test_matrix)[test_matrix[, j] > 0]
      refgnenes <-  rownames(ref_matrix)[ref_matrix[,i] > 0]
      result = ORA(testGenes,refgnenes,testbackground,refbackground);
      ORA.OR[j,i] = result[1];
      ORA.P[j,i] = result[2];
    }
  }
  
  colnames(ORA.P) = colnames(ORA.OR) = colnames(ref_matrix);
  rownames(ORA.P) = rownames(ORA.OR) = colnames(test_matrix);
  FDRmat <- matrix(p.adjust(as.numeric(data.matrix(ORA.P)), method = "BH"), nrow = nrow(ORA.P), ncol = ncol(ORA.P))
  rownames(FDRmat) = rownames(ORA.P)
  colnames(FDRmat) = colnames(ORA.P)
  ORA.P = matrix(as.numeric(ORA.P), nrow = nrow( ORA.P), ncol = ncol(ORA.P))
  ORA.OR = matrix(as.numeric(ORA.OR), nrow = nrow(ORA.OR), ncol = ncol(ORA.OR))
  rownames(ORA.P) <- rownames(ORA.OR) <- rownames(FDRmat)
  colnames(ORA.P) <- colnames(ORA.OR) <- colnames(FDRmat)
  
  Pmat <- ORA.P
  ORmat <- ORA.OR
  
  dispMat <- log2(ORmat)
  dispMat[FDRmat > 0.05] <- 0
  rownames(dispMat) <- rownames(ORmat)
  tMat <- matrix(signif(ORmat,2), nrow = nrow(dispMat), ncol = ncol(dispMat))
  tMat[FDRmat > 0.05] <- ""
  FDRmat[is.na(FDRmat)] <- 1
  tMat[is.na(FDRmat)] <- ""
  tMat[ORmat < 1] <- ""
  list(dispMat = dispMat ,tMat = tMat, FDRmat = FDRmat)
}

plot_switch_hCS <- function(genes){
  datExpr_switch <- datExpr_reg_batch %>%
    as.data.frame() %>%
    rownames_to_column("ensembl_gene_id") %>%
    left_join(geneAnno, by = "ensembl_gene_id") %>%
    filter(hgnc_symbol %in% genes) %>%
    gather("SampleID", "value", matches("MEF|Pool")) %>%
    left_join(datMeta_hcs, by = "SampleID") %>%
    group_by(hgnc_symbol) %>%
    mutate(scaled_expression = scale(value))
  
  ggplot(datExpr_switch, aes(x = Differentiation.day, y = scaled_expression)) +
    geom_rect(aes(xmin=250, xmax=300, ymin=-Inf, ymax=Inf), color = "grey90", fill = "grey90") +
    geom_jitter(width = 0.1, aes(color = hgnc_symbol, shape = hgnc_symbol), size = 0.6, show.legend = F) +
    geom_smooth(aes(color = hgnc_symbol, linetype = hgnc_symbol), method = 'loess', formula = y~x, size =0.6) +
    scale_y_continuous(breaks = seq(-10,10,2)) +
    scale_x_continuous(breaks = seq(0,700,100)) +
    scale_color_brewer(palette = "Set1", name = "") +
    scale_linetype(name = "", guide = FALSE) +
    scale_shape(name = "", guide = FALSE ) +
    labs(x = "Differentiation Day", y = "Scaled normalized expression") +
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5, size = 6),
          axis.text.y = element_text(size = 6))
}

plot_switch_bs <- function(genes){
  datExpr_switch <- datExpr_reg_cortical %>%
    as.data.frame() %>%
    rownames_to_column("ensembl_gene_id") %>%
    left_join(geneAnno, by = "ensembl_gene_id") %>%
    filter(hgnc_symbol %in% genes) %>%
    gather("SampleID", "value", matches("HSB")) %>%
    left_join(datMetaCortical, by = "SampleID") %>%
    group_by(hgnc_symbol) %>%
    mutate(scaled_expression = scale(value))
  
  ggplot(datExpr_switch, aes(x = Period, y = scaled_expression)) +
    labs(x = "BrainSpan Stage", y = "Scaled normalized expression") +
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5, size = 6),
          axis.text.y = element_text(size = 6))  +
    geom_vline(xintercept = 7.5 , linetype = 1, color = "grey 70") +
    geom_jitter(width = 0.1, aes(color = hgnc_symbol, shape = hgnc_symbol), size = 0.6, show.legend = F) +
    geom_smooth(aes(color = hgnc_symbol, linetype = hgnc_symbol), method = 'loess', formula = y~x, size =0.6) +
    scale_y_continuous(breaks = seq(-10,10,2)) +
    scale_x_continuous(breaks = seq(2,15,1)) +
    scale_color_brewer(palette = "Set1", name = "") +
    scale_linetype(name = "", guide = FALSE) +
    scale_shape(name = "", guide = FALSE)
}

plot_comparison <- function(genes){
  
  
  p1 <- plot_switch_hCS(genes) +
    guides(color = guide_legend(override.aes = list(fill = NA))) +
    theme(legend.key.size =  unit(.1, "in"), legend.position = "top")
  p2 <-   plot_switch_bs(genes) +
    guides(color = guide_legend(override.aes = list(fill = NA))) +
    theme(legend.key.size =  unit(.1, "in"), legend.position = "top")
  
  list(p1,p2)
}

singleGenePlot <- function(g, dataset = "hcs"){
  genes <- toupper(g)
  single_genes_anno <- geneAnno %>%
    filter(hgnc_symbol %in% genes)
  if (dataset == "hcs") {
    current_expr <- datExpr_reg_batch
    current_meta <- datMeta
    sample_regex <- "MEF|Pool"
    x_axis <- "Differentiation.day.original"
    xlab_pretty <- "Differentiation Day"
    
  } else if (dataset == "bs") {
    current_expr <- datExpr_reg_cortical
    current_meta <- datMetaCortical
    sample_regex <- "HSB"
    x_axis <- "Period"
    xlab_pretty <- "BrainSpan Stage"
  }
  
  plot_data <- current_expr %>%
    as.data.frame() %>%
    rownames_to_column("ensembl_gene_id") %>%
    inner_join(single_genes_anno, by = "ensembl_gene_id") %>%
    mutate(hgnc_symbol = factor(hgnc_symbol, levels = genes)) %>%
    gather("SampleID", "value", matches(sample_regex)) %>%
    left_join(current_meta, by = "SampleID")
  gene_display <- set_names(unique(as.character(plot_data$hgnc_symbol))) %>%
    ifelse(. == "POU3F2",paste(.,"\n(BRN2)") ,.) %>%
    ifelse(. == "BCL11B",paste(.,"\n(CTIP2)") ,.)
  marker_plot <- ggplot(plot_data,aes_string(x = x_axis, y = "value")) 
  if (dataset == "hcs") {
    marker_plot <- marker_plot +
      geom_rect(aes(xmin = 250, xmax = 300, ymin = -Inf, ymax = Inf), color = "grey90", fill = "grey90") 
  } else if (dataset == "bs") {
    marker_plot <- marker_plot +
      geom_vline(xintercept = 7.5 , linetype = 1, color = "grey 70")
  }
  marker_plot <- marker_plot +
    geom_jitter(width = 0.1, color = "grey50", size = 0.3) +
    facet_wrap(~hgnc_symbol,scales = "free_y", labeller = as_labeller(gene_display)) +
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 6, margin = margin(t = 0, r = 0, b = 0, l = 0)),
          strip.background = element_rect(colour = NA)) +
    geom_smooth(method = 'loess', formula = y~x, color = "grey10", size = 0.5) +
    labs(x = xlab_pretty, y = "Normalized expression")
}

plot_genes <- function(datExpr,datMeta, genes, scale_expr = T){
  if ("Period" %in% names(datMeta)) {
    xlab <-  "Period"
    sample_regex <- "HSB"
    xlab_pretty <- "BrainSpan Stage"
    x_breaks <- seq(2,15,1)
  } else {
    xlab <-  "Differentiation.day"
    sample_regex <- "MEF|Pool"
    xlab_pretty <- "Differentiation Day"
    x_breaks <- seq(0,700,100)
  }
  datExpr_tidy <- datExpr %>%
    as.data.frame() %>%
    rownames_to_column("ensembl_gene_id") %>%
    left_join(geneAnno, by = "ensembl_gene_id") %>%
    filter(hgnc_symbol %in% genes)  %>%
    gather("SampleID", "value", matches(sample_regex)) %>%
    left_join(datMeta, by = "SampleID") %>%
    group_by(hgnc_symbol)
  
  if(scale_expr){
    datExpr_tidy <- datExpr_tidy %>%
      mutate(scaled_expression =  scale(value))
  } else {
    datExpr_tidy <- datExpr_tidy %>%
      mutate(scaled_expression =  value)
  }
  
  
  p1 <- ggplot(datExpr_tidy, aes_string(x = xlab, y = "scaled_expression"))
  
  if ("Period" %in% names(datMeta)) {
    p1 <- p1 +
      geom_vline(xintercept = 7.5 , linetype = 1, color = "grey 70")
  } else {
    p1 <- p1 +
      geom_rect(aes(xmin = 250, xmax = 300, ymin = -Inf, ymax = Inf),
                color = "grey90", fill = "grey90")
  }
  p1 +
    geom_smooth(method = 'loess', formula = y~x, size =0.5, color = "grey10") +
    geom_jitter(width = 0.1, color = "grey50", size = 0.3) +
    scale_y_continuous(breaks = seq(-10,10,2)) +
    scale_x_continuous(breaks = x_breaks) +
    scale_color_brewer(palette = "Set1", name = "") +
    scale_linetype(name = "") +
    scale_shape(name = "") +
    labs(x = xlab_pretty, y = "Scaled normalized expression") +
    facet_wrap(~hgnc_symbol, nrow = 1) +
    theme(axis.text.y = element_text(size = 6),
          axis.text.x = element_text(size = 6, angle = -90, hjust = 0, vjust = 0.5),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 6, margin = margin(t = 0, r = 0, b = 0, l = 0)),
          strip.background = element_rect(colour = NA))
}


# Transition mapping ------------------------------------------------------------------------------------------------------------
# run transition mapping script first to obtain the input for this plot
hypmat_file <- load(file.path(outputFolder,"data", "hypermatAll.rdata"))


nr <- dim(hypermat.all)[1]
nc <- dim(hypermat.all)[2]

jet.colors = colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"));
colormap = jet.colors(100)


minhypermat = min(hypermat.all,na.rm=TRUE);
maxhypermat = max(hypermat.all,na.rm=TRUE)
png(file.path(outputFolder,"iPSC_brainSpan_day25_RRHOMap.png"),width=2.25,height=2.25,units="in",res=300);
par(mfrow = c(nr, nc), mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0))
for (k in 1:nr){
  for (i in 1:nc) {
    ## Scale the output to the max across all hypermats so that the colorbar is the same for all hypermats
    image(hypermat.all[k, i, , ],
          xlab = '',
          ylab = '',
          axes = FALSE,
          col = colormap,
          zlim = c(minhypermat,maxhypermat))
  }
}
dev.off()
## Function to plot color bar
## Modified from http://www.colbyimaging.com/wiki/statistics/color-bars
color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title = '') {
  scale = (length(lut)-1)/(max-min);
  plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='');
  mtext(title,2,2.7,cex=1.33333);
  axis(2, round(ticks,0), las=1,cex.lab=2, tick = FALSE);
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min;
    rect(0,y,10,y+1/scale, col=lut[i], border=NA);
  }
}
png(file.path(outputFolder,"iPSC_brainSpanAll_RRHOMap.png",sep = ""),width = 1.5,height = 6,units = "in", res = 300)
color.bar(lut = jet.colors(100),min = minhypermat,max = maxhypermat,nticks = 6, title = "-log10(Nominal P-value)")
dev.off()

# DNAmAge ---------------------------------------------------------------------------------------------------------
#DNA methylation age was calculated using https://horvath.genetics.ucla.edu/html/dnamage/
#raw methylion data can be found in GEO at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE150122

load(file.path("data", "DNAmAge_data.rdata"))
replicated_samples <- datSampleRaw$External.Sample.ID[duplicated(datSampleRaw$External.Sample.ID)]
datSample <- datSampleRaw %>%
  group_by(External.Sample.ID, Differentiation.day,Line, Sex ) %>%
  summarise(DNAmAge = mean(DNAmAge)) %>% 
  mutate(replicated = ifelse(External.Sample.ID %in% replicated_samples, "Replicated\nsample", "single"))

r <- cor.test(datSample$DNAmAge,datSample$Differentiation.day)

r_p_val <- formatC(r$p.value, format = "e", digits = 2) %>%
  str_split("e") %>%
  unlist

p1 <- ggplot(datSample, aes(x = Differentiation.day, y = DNAmAge)) +
  geom_smooth(method = "lm", color = "grey40") +
  annotate(geom = "text", x = 400, y = 8.5, size = 2, hjust = 0,label = paste("r =", round(r$estimate,2))) +
  annotate(geom = "text", x = 400, y = 7.8, size = 2, hjust = 0,
           label = substitute(
             paste("p = ", part1,"e"^part2),
             list(
               part1 = r_p_val[1],
               part2 = r_p_val[2])
           ),
  ) +
  geom_point(aes(color = Line), size = 1) +
  scale_color_discrete(name = "hiPS cell line") +
  scale_x_continuous(limits = c(0, 700), breaks = seq(0,700,100)) +
  labs(x = "Differentiation Day", title = "DNA methylation age") +
  guides(color = guide_legend(keywidth = 0.1, keyheight = 0.05,default.unit = "inch")) +
  theme(legend.justification = "top",
    legend.margin = margin(0,20,0,0),
    legend.background = element_blank(),
    legend.position = c(0.25, 1),
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(size = 6, color = "black"),
    axis.text.y = element_text(size = 6, color = "black"),
    axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)))

pdf(file = file.path(outputFolder,"DNAmAge_scatter.pdf"), w = 2.2, h = 2.2, useDingbats = FALSE)
p1 
dev.off()

# number of Samples -----------------------------------------------------------------------------------------------

stages <- sort(unique(as.numeric(datMeta$Day.grouped)))
breaks <- sapply(1:(length(stages) - 1),  function(i){
  stages[i] + ((stages[i + 1] - stages[i])/2)
})
rects <- data.frame(xstart = c(0,breaks),
                    xend = c(breaks, Inf),
                    col = rep_len(c("A","B"), length.out = length(stages)))

p_rna_samples <- ggplot() +
  geom_rect(data = rects, aes(xmin = xstart, xmax = xend, fill = col), ymin = -Inf, ymax = Inf, alpha = 0.1) +
  geom_jitter(dat = datMeta, aes(x = Differentiation.day.original, y = CellLine, shape = Sex),
              width = 0, height = 0.15, size = 1) +
  labs(x = "Differentiation day", y = "hiPSC line") +
  scale_x_continuous(breaks = seq(0, 700, 100)) +
  scale_shape_manual(values = c(1,16)) +
  scale_fill_manual(values = c("white", "black"),guide = FALSE) +
  theme_minimal(base_size = 6) +
  theme(axis.text = element_text(size = 6, color = "black"),
        axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5),
        axis.title.x = element_text(margin = margin(t = 5)),
        axis.title.y = element_text(margin = margin(r = 5)),
        legend.position = "none",
        plot.title = element_text(size = 6)) +
  labs(title = "RNA Expression")

p_meth_samples <- ggplot(datSample, aes(x = Differentiation.day, y = Line, shape = Sex)) +
  geom_jitter(width = 0, height = 0.15, size = 1, aes(color = replicated)) +
  scale_color_manual(values = c("single"="black", "Replicated\nsample" = "steelblue4"), 
                     name = "") +
  labs(x = "Differentiation day", y = "Individual") +
  scale_x_continuous(breaks = seq(0, 700, 100)) +
  scale_shape_manual(values = c(1,16)) +
  theme_minimal(base_size = 6) +
  theme(axis.text = element_text(size = 6, color = "black"),
        axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5),
        axis.title.x = element_text(margin = margin(t = 5)),
        axis.title.y = element_blank(),
        legend.text = element_text(size = 6),
        legend.margin = margin(0,0,0,0),
        plot.title = element_text(size = 6)) +
  guides(shape = guide_legend( keywidth = 0.1, keyheight = 0.1,
                               default.unit = "inch", order = 1), 
         color = guide_legend( keywidth = 0.1, keyheight = 0.1,
                               default.unit = "inch")) +
  labs(title = "Methylation")
pdf(file.path(outputFolder, "Samples.pdf"), w = 4, h = 2)
p_rna_samples | p_meth_samples
dev.off()

datMeta_rnaseq <- datMeta %>%
  mutate(Induction = as.numeric(as.factor(Induction))) %>%
  mutate(batch = as.numeric(as.factor(batch))) %>%
  mutate(Individual = gsub("-*", "", CellLine)) %>% 
  rowid_to_column(var = "Sample") %>%
  dplyr::select(Sample, CellLine, Differentiation = Induction, Individual,
                Differentiation.day = Differentiation.day.original,
                Day.grouped, Sex, Batch = batch, ancestryPC1, ancestryPC2,
                SeqPC1, SeqPC2, SeqPC3,SeqPC4,SeqPC5)

datMeta_methylation <- datSampleRaw %>%
  rowid_to_column(var = "Sample") %>%
  mutate(Sample = as.numeric(as.factor(External.Sample.ID))) %>%
  mutate("Technical replication" = ifelse(is.na(TechnicalReplicate)," " , T)) %>%
  arrange(Sample) %>%
  dplyr::select(Sample, Line, Sex,
                Differentiation.day, chip.ID,DNAmAge,  "Technical replication")


supp_table1 <- list("meta data RNAseq" = datMeta_rnaseq,
                    "meta data methylation" = datMeta_methylation)

writexl::write_xlsx(supp_table1, path = file.path(tables_folder, "Supp Table 1 - Sample data.xlsx"))
# PCA -------------------------------------------------------------------------------------------------------------
mdsG = cmdscale(dist(t(datExpr_reg_batch)), eig = TRUE)
mdsPlots <- cbind(datMeta, as.data.frame(mdsG$points)[match(datMeta$SampleID, rownames(mdsG$points)), ])
colnames(mdsPlots)[c((ncol(mdsPlots) - 1):ncol(mdsPlots))] <- c("MDS1", "MDS2")

pc_n <- 10
thisdat.expr <- t(scale(t(datExpr_reg_batch),scale=F))
PC.expr <- prcomp(thisdat.expr)
topPC.expr <- PC.expr$rotation[,1:pc_n]
varexp <- (PC.expr$sdev)^2 / sum(PC.expr$sdev^2)
topvar <- varexp[1:pc_n]
xy_labs <- paste(colnames(topPC.expr)," (",signif(100*topvar[1:pc_n],2),"%)",sep = "")
pairs_data <- datMeta %>%
  dplyr::select(CellLine,  Differentiation.day, Sex, batch) %>%
  mutate_at(vars(-Differentiation.day),list(~as.numeric(as.factor(.)))) %>%
  mutate(Differentiation.day = as.numeric(Differentiation.day))
all(rownames(topPC.expr) == datMeta$SampleID)
expr_pc_cors <- cor(topPC.expr,pairs_data)
expr_pc_rsqr <- sapply(pairs_data, function(cvr){
  apply(topPC.expr,2,function(pc){
    summary(lm(pc~as.factor(cvr)))$adj.r.squared
  })
})

rownames(expr_pc_rsqr) <- rownames(expr_pc_cors) <- paste("Expression", xy_labs)
colnames(expr_pc_rsqr) <- colnames(expr_pc_cors) <- c("Individual", "Differentiation day", "Sex", "Batch")

pdf(file.path(outputFolder, "PC_cor.pdf"), h = 2.5, w = 3 )

corrplot::corrplot(
  t(expr_pc_rsqr),
  method = "ellipse",
  tl.col	= "black",
  mar = c(0, 1.5, 0, 0),
  cl.align.text = "l",
  tl.cex = 0.5,
  cl.cex = 0.5,
  cl.length = 3
)
dev.off()

PC_plot <- topPC.expr %>%
  as.data.frame() %>%
  rownames_to_column("SampleID") %>%
  left_join(datMeta, by = "SampleID")
pdf(file.path(outputFolder, "LT_PC_plot_1.pdf"), width = 2.3, h = 2.2, useDingbats = FALSE)
ggplot(PC_plot,aes(x=PC1,y=PC2,color=Differentiation.day)) +
  geom_point(size=1) +
  theme(legend.position=c(0.94,0.8),
        legend.title.align = 0.5,
        legend.key.size = unit(0.075, "in"),
        legend.background = element_blank(),
        legend.text = element_text(size = 6),
        plot.title = element_text(hjust = 0.5),
        legend.title = element_blank(),
        axis.text = element_text(color = "black", size = 6)) +
  scale_color_viridis_c("Differentiation\nday") +
  labs(x = xy_labs[1], y= xy_labs[2], title = " ") +
  coord_cartesian(clip = "off")
dev.off()

# Dendogram -------------------------------------------------------------------------------------------------------
datMeta1 <- datMeta %>%
  dplyr::select(Line = CellLine, Differentiation.day, Sex, Batch = batch) %>%
  mutate_all(~(as.numeric(as.factor(.))))  %>%
  mutate(Line  = numbers2colors(Line, signed = FALSE, colors = viridis_pal()(6)),
         Sex = numbers2colors(Sex, signed = FALSE, colors = brewer.pal(6,"Paired")[5:6]),
         Batch = numbers2colors(Batch, signed = FALSE, colors = brewer.pal(4,"Paired")[3:4]),
         Differentiation.day = numbers2colors(Differentiation.day, signed = FALSE, colors =  colorRampPalette(c("white", "dodgerblue4"))( 20 )))

hc <- hclust(dist(t(datExpr_reg_batch)), method = "average")
pdf(file.path(outputFolder, "sample_dendogram.pdf"), width = 3.5, height = 2)
plotDendroAndColors( hc,
                     datMeta1,
                     hang = -1,
                     dendroLabels = FALSE, 
                     groupLabels = gsub("\\."," ",names(datMeta1)),
                     main = NULL,
                     marAll = c(0, 5,  2, 0),
                     cex.dendroLabels =  0.5,
                     autoColorHeight = F,
                     colorHeight = 0.48,
                     cex.colorLabels = 0.5,
                     cex.rowText = 0.5,
                     cex.axis = 0.5,
                     cex.lab = 0.5,
                     ylab = NULL,
                     axes = F
)
dev.off()

datMeta2 <- datMeta %>%
  dplyr::select(Individual = CellLine, Differentiation.day, Sex, batch) %>%
  mutate(Individual2 = Individual, Differentiation.day2 = Differentiation.day, Sex2 = Sex, batch2 = batch)  %>% 
  mutate_at(c("Individual" , "Differentiation.day", "Sex", "batch"), ~(as.numeric(as.factor(.))))  %>%
  mutate(Individual  = numbers2colors(Individual, signed = FALSE, colors = viridis_pal()(6)),
         Sex = numbers2colors(Sex, signed = FALSE, colors = brewer.pal(6,"Paired")[5:6]),
         batch = numbers2colors(batch, signed = FALSE, colors = brewer.pal(4,"Paired")[3:4]),
         Differentiation.day = numbers2colors(Differentiation.day, signed = FALSE, colors =  colorRampPalette(c("white", "dodgerblue4"))( 20 )))

lg1 <- Legend(at = unique(datMeta2$Individual2), 
              legend_gp = gpar(fill = unique(datMeta2$Individual)), 
              labels_gp = gpar(cex = 0.5),
              title_gp = gpar(cex = 0.5, fontface = 1),
              title = "Line", 
              grid_height = unit(2, "mm"), grid_width = unit(3, "mm"),
              nrow = 1)
lg2 <- Legend(at = c(25,200,400,600),
              col_fun = circlize::colorRamp2(c(25, 650), c("white", "dodgerblue4")), 
              labels_gp = gpar(cex = 0.5),
              title_gp = gpar(cex = 0.5, fontface = 1),
              title = "Differentiation day", 
              grid_height = unit(2, "mm"), grid_width = unit(3, "mm"),
              direction = "horizontal")
lg3 <- Legend(at = unique(datMeta2$Sex2), 
              legend_gp = gpar(fill = unique(datMeta2$Sex)), 
              labels_gp = gpar(cex = 0.5),
              title_gp = gpar(cex = 0.5, fontface = 1),
              title = "Sex", 
              grid_height = unit(2, "mm"), grid_width = unit(3, "mm"),
              nrow = 1)
lg4 <- Legend(at = unique(datMeta2$batch2), 
              legend_gp = gpar(fill = unique(datMeta2$batch)), 
              labels_gp = gpar(cex = 0.5),
              title_gp = gpar(cex = 0.5, fontface = 1),
              title = "Batch",
              grid_height = unit(2, "mm"), grid_width = unit(3, "mm"),
              nrow = 1)

dendo_legends <- packLegend(list = list(lg1, lg2, lg3, lg4), direction = "h", max_width = unit(3, "in"))
pdf(file.path(outputFolder, "sample_dendogram_legned.pdf"), width = 3.5, height = 1)
draw(dendo_legends)
dev.off()

# Reproducibility -------------------------------------------
# Correlation between samples from same/different individuals 

# get all samples from a specifc Cell-line/day
induction.replicates <- datMeta %>% 
  unite(ID, CellLine, Day.grouped) %>% 
  group_by(ID) %>% 
  filter(n() > 1) 

if(nrow(induction.replicates) > 0) {
  induction.replicates <- induction.replicates %>% 
    mutate(dummy = paste0("V",1:n())) %>% 
    ungroup() %>% 
    dplyr::select(ID, SampleID, dummy) %>% 
    spread(ID, SampleID) %>% 
    dplyr::select(-dummy)
  
  induction.combos <- lapply(induction.replicates, function(cell_line){
    combos <- combn(na.omit(cell_line),2)
  }) %>% as.data.frame()
  
  induction.cor <- apply(induction.combos, 2, function(samples){
    cor(datExpr_reg_batch[,samples[1]], datExpr_reg_batch[,samples[2]], method = "s")
  })
} else{
  induction.cor <- set_names(rep(NA, length(unique(datMeta$Day.grouped))), paste0("X_",unique(datMeta$Day.grouped)))
}

# reproducibility between cell lines (within a timepoint)
samples_by_day <- datMeta %>% 
  group_by(Day.grouped) %>% 
  mutate(dummy = paste0("V",1:n())) %>% 
  ungroup() %>% 
  dplyr::select(Day.grouped, SampleID, dummy) %>% 
  spread(Day.grouped, SampleID) %>% 
  dplyr::select(-dummy)

all_combos_day  <- lapply(samples_by_day, function(day){
  combos <- combn(na.omit(day),2)
}) %>% as.data.frame()

unique_samples_by_day <- lapply(1:ncol(all_combos_day),function(i){
  current <- all_combos_day[,i, drop = F]
  rev_current <- current[2:1,, drop = F]
  if ((current %in% induction.combos) | (rev_current %in% induction.combos)) {
    return(NULL)
  } else {
    return(current)
  }
})
stopifnot(sum(sapply(unique_samples_by_day,is.null)) == ncol(induction.combos))
unique_samples_by_day <- unique_samples_by_day[!sapply(unique_samples_by_day,is.null)] %>% 
  as.data.frame()

unique_samples_cor <- apply(unique_samples_by_day, 2, function(samples){
  cor(datExpr_reg_batch[,samples[1]], datExpr_reg_batch[,samples[2]], method = "s")
})

all_cor <- data.frame(type = "Within Individual", cor = induction.cor, day = gsub(".*_([0-9]{3}).*","\\1",names(induction.cor))) %>% 
  bind_rows(data.frame(type = "Between Individuals", cor = unique_samples_cor, day = gsub("X([0-9]{3}).*","\\1",names(unique_samples_cor)))) 


pdf(file.path(outputFolder,"reproducibility_cor.pdf"), h = 1.5, w = 4.5)
all_cor %>% 
  bind_rows() %>%
  mutate(cor = ifelse(is.na(cor),0, cor)) %>%
  mutate(day = gsub("^0", "", day)) %>% 
  mutate(day = factor(day, levels = str_sort(unique(day), numeric = T))) %>% 
  ggplot(aes(x = day, y = cor)) +
  geom_boxplot(aes(color = type), outlier.colour = NA, fill = NA,
               position = position_dodge2(preserve = 'single'),
               show.legend = FALSE
  ) +

  geom_jitter(aes(color = type),size = 0.3,
              position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.1)) +
  labs(x = "Differentiation Day", y = "Spearman's correlation", title = "Reproducibility")  +
  scale_color_brewer("", palette = "Set1", drop = FALSE) +
  scale_y_continuous(limits = c(0.65,1)) +
  theme(
    legend.position = c(0.85,0.23), 
    legend.margin = margin(0,0,0,0), 
    legend.box.margin = margin(0,0,0,0),
    plot.title = element_text(hjust = 0.5),
    legend.spacing.x = unit(-1, "mm"), 
    legend.key.height = unit(4,"mm"),
    legend.direction = "vertical", 
    legend.spacing.y = unit(2,"mm"),
    legend.background = element_blank()
  )
dev.off()

# Partition Variance ---------------------------------------------------------------------------------------------

library(variancePartition)

form <- as.formula(paste("~ (1|Day.grouped) + (1|CellLine) + (1|batch) + (1|Sex) + (1|Induction) +",
                         paste0("SeqPC",1:5, collapse = " + "),
                         " + ancestryPC1 + ancestryPC2") )

varPart <- fitExtractVarPartModel( datExpr_reg_batch, form, datMeta )


median_var <- map_df(as.data.frame(varPart), function(x){
  median(x*100)
}) %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("variable") %>% 
  arrange(desc(V1)) %>% 
  mutate(V1 = case_when(V1 < 0.01 ~"<0.01", 
                        V1 < 0.5 ~ as.character(round(V1, 2)),
                        T ~ as.character((round(V1, 1))))) %>% 
  mutate(V1 = paste0(V1, "%")) %>% 
  mutate(variable = plyr::mapvalues(variable,  from = c("Day.grouped", "CellLine", "Induction", 
                                                        "ancestryPC1","ancestryPC2", "batch" , 
                                                        paste0("SeqPC", 1:5)),
                                    to = c("Differentiation day", "Cell line", "Differentiation", 
                                           "Genetic ancestry PC1", "Genetic ancestry PC2", "Batch",
                                           paste0("Technical PC", 1:5))
  ))
fct_lvls <- c(median_var$variable[median_var$variable != "Residuals"] , "Residuals")

p11 <- as.data.frame(varPart) %>% 
  rownames_to_column("ensembl_gene_id") %>% 
  gather("variable", "value", -ensembl_gene_id) %>% 
  mutate(variable = plyr::mapvalues(variable,  from = c("Day.grouped", "CellLine", "Induction", 
                                                        "ancestryPC1","ancestryPC2", "batch" , 
                                                        paste0("SeqPC", 1:5)),
                                    to = c("Differentiation day", "Cell line", "Differentiation", 
                                           "Genetic ancestry PC1", "Genetic ancestry PC2", "Batch",
                                           paste0("Technical PC", 1:5))))%>% 
  mutate(value  = value * 100, 
         variable = factor(variable, levels = fct_lvls)) %>% 
  ggplot(aes(x = variable, y = value)) +
  geom_violin(scale = "width", fill = "grey50") +
  theme(axis.text.x = element_text(angle = -90, hjust  = 0, vjust = 0.5),
        axis.title.x = element_blank(),
        legend.position = "none") +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.5, outlier.shape = NA) +
  annotate("text", x = median_var$variable , y = 120, 
           label = median_var$V1, size = 2, angle = -90, 
           hjust = 0) +
  theme(plot.margin = margin(t = 10)) +
  scale_y_continuous(limits = c(0,120), breaks = c(0,25,50,75,100)) +
  coord_cartesian(clip = 'on') +
  labs(y = "Variance explained (%)")


pdf(file.path(outputFolder,"Variance_partioned.pdf"), h = 2.5, w = 2.5)
p11
dev.off()

# Single Gene Expression ------------------------------------------------------------------------------------------
cell_type_genes2 <- list("Radial glia" = c("PAX6", "GLI3"),
                         "Intermediate progenitors" = c("EOMES","NEUROD1"),
                         #"Neurons" = c("MAP2","NCAM1"),
                         "Deep layer cortical neurons" = c("TBR1","BCL11B"),
                         "Superficial layer cortical neurons" = c("SATB2", "POU3F2"),
                         
                         "Early stage astrocytes" = c("TMSB15A", "NNAT"),
                         "Late stage astrocytes" = c("AQP4","GFAP")
)



marker_plot <- map2(cell_type_genes2,names(cell_type_genes2),function(g,n){
  singleGenePlot(g) +
    theme(strip.text = element_text(size = 6, face = "italic"), 
          axis.title.y = element_blank()) +
    ggtitle(n)
}) 
for(i in c(1,3,5)){
  marker_plot[[i]] <- marker_plot[[i]] +
    labs(x = "")
 
}

pdf(file.path(outputFolder,"cell_type_markers_single_genes.pdf"),width = 5, height = 3)
wrap_plots(marker_plot, byrow = F)
dev.off()

marker_plot_bs <- map2(cell_type_genes2,names(cell_type_genes2),function(g,n){
  singleGenePlot(g, dataset = "bs") +
    ggtitle(n) +
    theme(strip.text = element_text(size = 6, face = "italic"), 
          axis.title.y = element_blank()) +
    scale_y_continuous(breaks = scales::breaks_pretty(n=3))
})

all_marker_plots_bs <- (marker_plot_bs[[1]] + theme(axis.title.x = element_blank()) +
                          marker_plot_bs[[2]] + theme(axis.title.x = element_blank())) /
  (marker_plot_bs[[3]] + theme(axis.title.x = element_blank() ) +
     marker_plot_bs[[4]] + theme(axis.title.x = element_blank())) /
  (marker_plot_bs[[5]] +  marker_plot_bs[[6]] ) 

pdf(file.path(outputFolder,"cell_type_markers_single_genes_bs.pdf"),width = 3.5, height = 3)
all_marker_plots_bs
dev.off()


# non preserved cell types
cell_type_genes_np <- list("Inhibitory Neurons" = c("GAD1", "GABBR1" ),
                           "OPCs" =c("OLIG2", "PLP1")
)

marker_plot_np<- map2(cell_type_genes_np,names(cell_type_genes_np),function(g,n){
  singleGenePlot(g) +
    ggtitle(n) +
    scale_y_continuous(breaks = scales::breaks_pretty(n=3)) +
    theme(strip.text = element_text(size = 6, face = "italic"))
  
})

marker_plot_bs_np <- map2(cell_type_genes_np,names(cell_type_genes_np),function(g,n){
  singleGenePlot(g, dataset = "bs") +
    theme(strip.text = element_text(size = 6, face = "italic"))
  
})

all_marker_plots_np <- marker_plot_np[[1]]   +marker_plot_np[[2]]+
  marker_plot_bs_np[[1]]  + marker_plot_bs_np[[2]]+
  plot_layout( nrow = 2) & theme(axis.title.y = element_blank())

pdf(file.path(outputFolder,"cell_type_markers_single_genes_np.pdf"),width = 3.5, height = 2.25)
all_marker_plots_np
dev.off()

# activity genes
activity_genes <- c("NPAS4", "ARC","FOS","EGR1")

activity_marker_plot_hcs <- singleGenePlot(activity_genes, ) +
  facet_wrap(~hgnc_symbol, nrow = 1, scale = "free_y") +
  scale_y_continuous(breaks = scales::breaks_pretty(n = 3))


activity_marker_plot_bs <- singleGenePlot(activity_genes, dataset = "bs") +
  facet_wrap(~hgnc_symbol, nrow = 1, scale = "free_y") + 
  scale_y_continuous(breaks = seq(-12,12,4))



activity_marker_plot <- activity_marker_plot_hcs/
  activity_marker_plot_bs & theme(axis.title.y = element_blank(), 
                                  strip.text = element_text(size = 6, face = "italic")
  )

pdf(file.path(outputFolder,"activity_markers_single_genes.pdf"),width = 3.5, height = 2.25)
activity_marker_plot
dev.off()


# astrocyte transition --------------------------------------------------------------------------------------------
sloan_data_raw <- read_excel(file.path("data","Sloan2018_astrocyte_genes.xlsx"))

sloan_data <- sloan_data_raw %>%
  gather("stage","hgnc_symbol") %>%
  left_join(geneAnno, by = "hgnc_symbol") %>%
  filter(ensembl_gene_id %in% rownames(datExpr_reg_batch)) %>% 
  mutate(stage = ifelse(stage == "Mature genes", "Postnatal genes",stage))

datMeta_by_day <- datMeta %>%
  arrange(Differentiation.day.original)

sample_by_day <- datMeta_by_day %>%
  dplyr::select(SampleID) %>%
  unlist()

astro_expr <- datExpr_reg_batch[sloan_data$ensembl_gene_id,sample_by_day]

all(rownames(astro_expr) == sloan_data$ensembl_gene_id)

astro_median_expr <- map(set_names(sort(unique(datMeta$Day.grouped))), function(day){
  current_samples <- datMeta %>%
    filter(Day.grouped == day) %>%
    dplyr::select(SampleID) %>%
    unlist()
  apply(astro_expr[,current_samples],1,median)
}) %>% 
  bind_cols() %>%
  as.matrix() %>%
  `rownames<-`(rownames(astro_expr))

row_anno <- rowAnnotation(
  "Astrocyte genes" = sloan_data$stage, 
  col = list("Astrocyte genes" = c("Fetal genes" = RColorBrewer::brewer.pal(3,"Set1")[1],
                                   "Postnatal genes" = RColorBrewer::brewer.pal(3,"Set1")[2])), 
  show_annotation_name = F, 
  
  annotation_legend_param = list(
    "Astrocyte genes" = list(nrow = 2,
                             labels_gp = gpar(fontsize = 6),
                             title_gp = gpar(fontsize = 6, fontface = "plain")))
)

col_anno <- HeatmapAnnotation(
  Stage = colnames(astro_median_expr), 
  col = list("Stage" = setNames(paste0("grey", seq(93,1,-8)),
                                colnames(astro_median_expr))), 
  show_annotation_name = F, 
  annotation_legend_param = list(
    "Stage" = list(nrow = 3,
                   labels_gp = gpar(fontsize = 6),
                   title_gp = gpar(fontsize = 6, fontface = "plain"))),
  which = "col") 

pdf(file.path(outputFolder, "astrocyte_transition.pdf"), w = 2.5, h = 3.5)

astro_hm <- Heatmap(t(apply(astro_median_expr,1,scale)), 
                    cluster_rows = F, 
                    cluster_columns = F, 
                    show_row_names = F, 
                    col = colorRampPalette( c("green","black", "magenta"))(200),
                    left_annotation = row_anno, 
                    top_annotation = col_anno,
                    heatmap_legend_param = list(title = "Scale Expression", 
                                                legend_direction = "horizontal",  
                                                legend_width = unit(1,"in"),
                                                at = c(-4, 0, 4),
                                                labels_gp = gpar(fontsize = 6),
                                                title_gp = gpar(fontsize = 6, fontface = "plain")
                    )
                    
)
draw(astro_hm, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")

dev.off()


# metabolic stress genes ------------------------------------------------------------------------------------------
stress_genes <- c("PGK1","ALDOA","BNIP3","YIP5","ARCN1","GORASP2")

stress_hcs <- plot_genes(datExpr_reg_batch, datMeta, stress_genes, scale_expr = F) %+%
  facet_wrap(~hgnc_symbol, nrow = 1) +
  theme(axis.title.y = element_blank(),  
        strip.text = element_text(size = 6, face = "italic"))+
  scale_y_continuous(limits =c(3,12), breaks = seq(4,10,2)) 

stress_bs <- plot_genes(datExpr_reg_cortical, datMetaCortical, stress_genes, scale_expr = F) %+%
  facet_wrap(~hgnc_symbol, nrow = 1) +
  theme(axis.title.y = element_blank(),  
        strip.text = element_text(size = 6, face = "italic"))+
  scale_y_continuous(limits =c(3,12), breaks = seq(4,10,2)) 


pdf(file.path(outputFolder,"metabolic_stress_single_genes.pdf"),width = 4, height = 2)
stress_hcs / plot_spacer()/stress_bs +
  plot_layout(heights = c(1,0.01,1))

dev.off()


stree_hcs_data <- stress_hcs$data %>%
  mutate(data_set = "hCS")
stress_bs_data <- stress_bs$data %>%
  mutate(data_set = "bs")
stress_data <- bind_rows(stree_hcs_data,stress_bs_data) %>%
  dplyr::select(hgnc_symbol, value, data_set)

#modules from pollen 2019
pollen_wgcna_genes <- read_excel(file.path("data","Pollen_2019_table_s3.xlsx"), skip = 753)

stress_modules <- c("ER_stress" = "organoid.human.ME.darkred",
                    "Glycolysis" = "organoid.Sloan.human.ME.paleturquoise")

stress_module_genes <-   map(stress_modules , function(mod){
  pollen_wgcna_genes %>%
    rename(mod_name = "Name of Module in Other Figures and Tables") %>%
    filter(mod_name == mod) %>%
    left_join(geneAnno, by = c("Gene" = "hgnc_symbol")) %>%
    filter(!is.na(ensembl_gene_id))
}) %>%
  bind_rows(.id = "mod")

hcs_genes <- intersect(stress_module_genes$ensembl_gene_id, rownames(datExpr_reg_batch))
mods_hcs <- stress_module_genes[match(hcs_genes, stress_module_genes$ensembl_gene_id),]

stress_expr_hcs <- t(datExpr_reg_batch[hcs_genes,])
stopifnot(all(hcs_genes == mods_hcs$ensembl_gene_id))
res_hcs <- moduleEigengenes(stress_expr_hcs,mods_hcs$mod)
eigen_values_hcs <- res_hcs$eigengenes %>%
  mutate(SampleID = rownames(stress_expr_hcs)) %>%
  gather("ME", "value", -SampleID) %>%
  left_join(datMeta, by = "SampleID") %>%
  mutate(ME = gsub("ME", "", ME)) %>%
  mutate(ME = gsub("_", " ", ME))
limits <-  range(eigen_values_hcs$value) + c(-0.1,0.1)
limits <-  round(limits,digits = 2)

clrs <- RColorBrewer::brewer.pal(4, "Set1")[3:4]
p_stress_hcs <- ggplot(eigen_values_hcs, aes(x = Differentiation.day, y = value, color = ME)) +
  geom_rect(aes(xmin=250, xmax=300, ymin=-Inf, ymax=Inf), color = "grey90", fill = "grey90") +
  geom_smooth(method = 'loess', formula = y~x, size =0.6) +
  geom_point(size = 0.3)  +
  labs(x = "Differentiation Day", y = "Module Eigengene") +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5, size = 6),
        axis.text.y = element_text(size = 6)) +
  scale_y_continuous(limits = limits) +
  scale_x_continuous(breaks = seq(0,700,100)) +
  scale_color_manual(values = clrs,"Module") +
  guides(color=guide_legend(override.aes=list(fill=NA)))


genes_bs <- intersect(stress_module_genes$ensembl_gene_id, rownames(datExpr_reg_cortical))
mods_bs <- stress_module_genes[match(genes_bs, stress_module_genes$ensembl_gene_id),]

stress_expr_bs <- t(datExpr_reg_cortical[genes_bs,])
stopifnot(all(genes_bs == mods_bs$ensembl_gene_id))
res_bs <- moduleEigengenes(stress_expr_bs,mods_bs$mod)
eigen_values_bs <- res_bs$eigengenes %>%
  mutate(SampleID = rownames(stress_expr_bs)) %>%
  gather("ME", "value", -SampleID) %>%
  left_join(datMetaCortical, by = "SampleID") %>%
  mutate(ME = gsub("ME", "", ME)) %>%
  mutate(ME = gsub("_", " ", ME))

p_stress_bs <- ggplot(eigen_values_bs, aes(x = Period, y = value, color = ME)) +
  geom_vline(xintercept = 7.5 , linetype = 1, color = "grey 70") +
  geom_smooth(method = 'loess', formula = y~x, size = 0.6) +
  geom_point(size = 0.3) +
  labs(x = "BrainSpan Stage", y = "") +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5, size = 6),
        axis.text.y = element_text(size = 6))  +
  scale_y_continuous(limits = limits) +
  scale_x_continuous(breaks = seq(2,15,1)) +
  scale_color_manual(values = clrs,"Module") +
  guides(color=guide_legend(override.aes=list(fill=NA)))


pdf(file.path(outputFolder,"metabolic_stress_ME.pdf"),width =3.5, height = 1.25)
p_stress_hcs + p_stress_bs + plot_layout(guides = "collect") &
  theme(legend.key.size =  unit(.15, "in"),
        legend.box.margin = margin(0,0,0,0),
        legend.position = "right",
        legend.box.spacing = unit(0,"pt"),
        plot.margin = margin(5.5,3.5,5.5,5.5))

dev.off()


# Context modules trajectories-------------------------------------------------------------------------------------------------
isDark <- function(colr) { (sum( col2rgb(colr) * c(299, 587,114))/1000 < 123) }

Context_MM <- read.csv(file.path("data", "Stein2015_Gene_module_membership.csv"))

load(file.path("data","context_gene_anno.rdata"))


Context_MM$ensembl_gene_id <- geneAnnoRaw[match(Context_MM$ENTREZ_GENE_ID,geneAnnoRaw$entrezgene),"ensembl_gene_id"]
commmon_ensmbl <- intersect(rownames(datExpr_reg_batch),Context_MM$ensembl_gene_id)
datExpr_Cntxt <- datExpr_reg_batch[commmon_ensmbl,]
Context_MM <- Context_MM[match(commmon_ensmbl,Context_MM$ensembl_gene_id),]
all(rownames(datExpr_Cntxt) == Context_MM$ensembl_gene_id)
ME <- moduleEigengenes(t(datExpr_Cntxt),Context_MM$moduleColor)
MEs <- ME$averageExpr
MEs$SampleID <- colnames(datExpr_Cntxt)
MEs_plot <- MEs %>%
  left_join(datMeta,by = "SampleID") %>%
  gather("Module", "value", matches("AE")) %>%
  mutate(Module = gsub(".E","",Module)) %>%
  mutate(Module = gsub("^white$","wheat",Module)) %>%
  mutate(text_color = if_else(isDark(Module), "grey80", "grey20"))

module_groups <- list(c("salmon", "cyan", "green", "lightgreen", "purple","darkgreen"),
  c("pink", "tan", "yellow", "lightcyan", "black")
 
)

context_plots <- map(module_groups, function(grp){
  plot_subset <- subset(MEs_plot,Module %in% grp) %>%
    mutate(Module = factor(Module, levels = grp))
  
  
  current_plot <- ggplot(plot_subset,aes(x = Differentiation.day, y = value)) +
    theme(strip.background = element_rect(fill = "white"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "none",
          plot.title = element_text(size = 6)) +
    
    geom_rect(aes(xmin = 250, xmax = 300, ymin = -Inf, ymax = Inf), color = "grey90", fill = "grey90") +
    geom_jitter(aes(color = Module), size = 0.75) +
    labs(x = "Differentiation day", y = "Scaled mean expression") +
    geom_smooth(method = 'loess', formula  = 'y ~ x', aes(color = Module), size = 0.75) +
    scale_color_manual(values = grp)
  
  label_info <- layer_data(current_plot,3) %>%
    group_by(colour) %>%
    arrange(desc(x)) %>%
    slice(1) %>%
    mutate(text_color = if_else(isDark(colour), "grey90", "grey20"))
  current_plot <- current_plot +
    geom_label_repel(data = label_info,
                     aes(x = x, y=y, label = colour, fill = colour),
                     color =  label_info$text_color,
                    
                     nudge_x      = 350,
                     direction = "y",
                     label.size = NA,
                     hjust        = 0.5,
                     segment.size = 0.2,
                     label.padding	= 0.15,
                     size = 2) +
    theme(legend.position = "none") +
    scale_x_continuous(limits = c(0,800), breaks = seq(0,600,100)) +
    scale_fill_manual(values = label_info$colour) +
    scale_y_continuous(breaks = seq(-2,2,1), limits = c(-1.5, 1.5)) +
    coord_cartesian(clip = "off")
})


pdf(file.path(outputFolder, "Trajectories_modules_subset.pdf"), height = 2, width = 6, useDingbats = FALSE)
context_plots[[1]] + labs(title ="Neuronal modules") +
  context_plots[[2]] + labs(title = "Glial Modules")  + plot_layout(ncol = 2)
dev.off()

# context module overlap ------------------------------------------------------------------------------------------
# run WGCNA using the designated scripn before running this section
attach(file.path("data","wgcna.Rdata"))

wgcna_matrix_raw <-  sapply(unique(geneInfo$Initially.Assigned.Module.Color), function(mod){
  mod_membership <- 1*(geneInfo$Initially.Assigned.Module.Color == mod)
  names(mod_membership) <- geneInfo$GeneSymbol
  mod_membership})
mm_modules <- unique(Context_MM$moduleColor) %>%
  .[. != "grey"]
context_matrix <- sapply(mm_modules, function(mod){
  as.numeric(Context_MM$moduleColor == mod)
})
rownames(context_matrix) <- Context_MM$geneSymbol
common_genes <- intersect(rownames(context_matrix), rownames(wgcna_matrix_raw))
wgcna_matrix_gs <- wgcna_matrix_raw[common_genes,]
context_matrix_gs <- context_matrix[common_genes,]

wgcna_or <- geneSetEnrichment( context_matrix_gs, wgcna_matrix_gs)
dispMat <- wgcna_or[["dispMat"]]
tMat <- wgcna_or[["tMat"]]
dispMat2 <- dispMat
dispMat2[dispMat2 < 0] <- 0
neuronal <- c("cyan", "darkgreen", "darkred", "green", "lightgreen", "paleturquoise", "purple", "salmon")
glial <- c("black", "lightcyan", "orange", "pink","tan","yellow")
cols <- rep('black', nrow(dispMat2))
cols[row.names(dispMat2) %in% neuronal] <- brewer.pal(4,"Set1")[3]
cols[row.names(dispMat2) %in% glial] <- brewer.pal(4,"Set1")[4]
pdf(file.path(outputFolder, "overlap_context_modules.pdf"), height = 5, width = 3.5)
heatmap.2(dispMat2,
          trace = "none",
          col = colorpanel(n = 100, low = "white", high = "blue3"),
          cellnote = tMat,
          notecex = 0.5,
          notecol = "white",
          cexRow = 0.6,
          cexCol = 0.6,
          margins = c(5, 5),
          density.info = "none",
          key = F,
          offsetRow = 0,
          offsetCol  = 0,
          cex.lab = 1,
          colRow = cols,
          lhei = c(0.2,1),
          lwid = c(0.2,1)
)
dev.off()

# DEG ---------------------------------------------------------------------------------------------------
# Dun differential expression before running this section
load(file.path(outputFolder, "data","LongTerm_pariedVoom_results.rdata"))

n_DEG <- lapply(hcs_lt_voom, function(x){
  x %>%
    rownames_to_column("ensembl_gene_id") %>%
    dplyr::filter(adj.P.Val < 0.05) %>%
    dplyr::select(ensembl_gene_id) %>%
    unlist() %>%
    length()
}) %>% unlist()

comparisons <- c("Day.grouped200-Day.grouped025",
                 "Day.grouped400-Day.grouped200")

n_DEG_direction <- lapply(comparisons, function(x){
  hcs_lt_voom[[x]] %>%
    rownames_to_column("ensembl_gene_id") %>%
    dplyr::filter(adj.P.Val < 0.05) %>%
    mutate(direction = ifelse(logFC < 0, "down", "up")) %>%
    group_by(direction) %>%
    summarise(n = n()) %>%
    mutate(cond = x)
}) %>% bind_rows()

plot_data <- n_DEG_direction %>%
  mutate(n = ifelse(direction == "up", n, -n),
         cond = gsub("Day.grouped","",cond)) %>%
  mutate(cond = gsub("-", " vs. ", cond)) %>%
  mutate(cond = gsub("0(\\d{2})","\\1", cond)) %>%
  mutate(cond = paste("Differentiation day\n",cond))
limits <- c(min(plot_data$n) - 750, max(plot_data$n) + 750)

pdf(file.path(outputFolder, "n_degs.pdf"), h = 1, w = 4.25)
ggplot(plot_data, aes(x = 1, y = n, fill = direction)) +
  geom_bar(stat = "identity", position = "identity") +
  scale_fill_manual(values =  RColorBrewer::brewer.pal(4, "Set1")[3:4], name = "", labels = c("down regulated","up regulated")) +
  geom_hline(yintercept = 0, size = 1) +
  scale_y_continuous(breaks = NULL,name = "",limits = limits) +
  scale_x_discrete(name = "") +
  theme(axis.line = element_blank(),
        legend.key.size = unit(3,"mm"),
        plot.background = element_blank(),
        panel.background = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank(),
        legend.position = c(0.5,-0.1),
        panel.grid.major.y = element_blank(),
        legend.direction = "horizontal",
        strip.text = element_text(size = 6)) +
  coord_flip() +
  geom_text(data = filter(plot_data, direction == "up"),aes(label = paste(n, "\ngenes")),hjust = -0.1, size = 2) +
  geom_text(data = filter(plot_data, direction == "down"),aes(label = paste(abs(n), "\ngenes")),hjust = 1.1, size = 2) +
  facet_wrap(~cond)
dev.off()

pretty_comparisons <- comparisons %>%
  gsub("Day.grouped","",.) %>%
  gsub("-0*", " vs. ",.)

supp_table2 <- lapply(comparisons, function(x){
  hcs_lt_voom[[x]] %>%
    rownames_to_column("ensembl_gene_id") %>%
    left_join(geneAnno, by = "ensembl_gene_id") %>%
    dplyr::select(ensembl_gene_id, hgnc_symbol, chromosome_name,band,
                  logFC, AveExpr, t, P.Value, adj.P.Val, B)
}) %>%
  setNames(pretty_comparisons)

writexl::write_xlsx(supp_table2, path = file.path(tables_folder, "Supp Table 2 - DEGs.xlsx"))
# fgsea -----------------------------------------------------------------------------------------------------------
# Dun differential expression before running this section
load(file = file.path(outputFolder,"data", "fgsea_results.rdata"))
load(file.path(outputFolder,"data","genedat.rdata"))

geneDat_full <- geneDat %>%
  left_join(geneAnno,by = "ensembl_gene_id") %>%
  dplyr::select(ensembl_gene_id, hgnc_symbol, entrezgene_id)

supp_table3 <- map(names(fgsea_res),function(f){
  res <- fgsea_res[[f]] %>%
    mutate(pathway = gsub("GO_", "", pathway)) %>%
    mutate(pathway = gsub("_", " ", pathway)) %>%
    filter(padj < 0.05) %>%
    arrange(desc(NES))
  
  
  hgnc_le <- map(res$leadingEdge,function(le){
    data.frame(entrezgene_id = as.numeric(le)) %>%
      left_join(geneDat_full,by = "entrezgene_id") %>%
      dplyr::select(hgnc_symbol) %>%
      unlist() %>%
      paste(collapse = "|")
  })
  res$leadingEdge <- unlist(hgnc_le)
  return(res)
}) %>% setNames(pretty_comparisons)
writexl::write_xlsx(supp_table3, path = file.path(tables_folder,"Supp Table 3 - GSEA.xlsx"))

n_to_plot <- 3
data_gsea_plot <- fgsea_res %>%
  map(~arrange(.,NES)) %>%
  map(~slice(.,1:n_to_plot,(n()-n_to_plot+1):n())) %>%
  bind_rows(.id = "Comparison") %>%
  mutate(Comparison = gsub(".*(\\d{3}).*(\\d{3}).*", "day \\1 vs. day \\2", Comparison)) %>%
  mutate(Comparison = gsub(" 0", " ", Comparison)) %>%
  dplyr::select(-leadingEdge) %>%
  mutate(pathway = gsub("GO_", "", pathway)) %>%
  mutate(pathway = str_to_sentence(gsub("_", " ", pathway))) %>%
  mutate(pathway = gsub("dna","DNA", pathway, ignore.case = T)) %>%
  mutate(pathway = gsub("rna","RNA", pathway, ignore.case = T)) %>%
  mutate(pathway = str_wrap(pathway, width = 30)) %>%
  mutate(pathway = reorder_within(pathway, NES, Comparison))

p_fgsea <- ggplot(data_gsea_plot, aes(x = pathway, y = NES, fill = -log10(padj))) +
  coord_flip() +
  geom_col() +
  theme_classic(base_size = 6) +
  facet_wrap(~Comparison, scales = "free_y", ncol = 1)+
  geom_hline(yintercept = 0) +
  scale_y_continuous(breaks=seq(-3,3,2) )+
  theme(strip.background = element_rect(),
        legend.position = "bottom",
        strip.text.x = element_text(size = 6),
        legend.key.height = unit(0.1, "in"),
        axis.title.x = element_text(hjust = 1, color = "black"),
        axis.text =  element_text(size = 6, color = "black"), 
        legend.margin = margin(0,0,0,0)) +
  labs(y = "Normalized enrichment score", x = "") +
  scale_x_reordered() +
  scale_fill_gradient(low = "#FCBBA1",high = "#DE2D26", name = expression(paste(-log[10](FDR))), limits = c(0,NA))

pdf(file.path(outputFolder, "gsea_bargraph.pdf"), width = 2.3, h = 3.75)
p_fgsea
dev.off()

# switch-overs -----------------------------------------------------------------------------------------------------
datMeta_hcs <- datMeta
plots <- list()
plots <-  c(plots,plot_comparison(c("HDAC1", "HDAC2")))
plots <-  c(plots,plot_comparison(c("HDAC11", "HDAC2")))
plots <-  c(plots,plot_comparison(c("GRIN2A", "GRIN2B")))
plots <-  c(plots,plot_comparison(c("GRIN2C", "GRIN2D")))
for (i in 1:length(plots)) {
  plots[[i]] <- plots[[i]] + labs(y = "")
  if (!i %in% c(3,4,7,8)) {
    plots[[i]] <- plots[[i]] + labs(x = "")
  }
}

plots <- lapply(plots, function(p){
  p + theme(plot.margin = margin(0,0,0,0), 
            legend.text = element_text(size = 6, face = "italic"), 
            legend.box.margin = margin(-20,-10,-10,-10), 
            plot.background = element_blank())
  
})
l1 <- plots[[1]] + plots[[2]] +  plot_layout( guides = "collect") &theme(legend.position = "bottom", legend.box.margin = margin(-13,0,0,0))
l2 <- plots[[3]] + plots[[4]]+ plot_layout( guides = "collect") &theme(legend.position = "bottom", legend.box.margin = margin(-5,0,-5,0))
l3 <-  plots[[5]] + plots[[6]] +  plot_layout(guides = "collect") &theme(legend.position = "bottom", legend.box.margin = margin(-13,0,0,0))
l4 <- plots[[7]] + plots[[8]]  + plot_layout( guides = "collect") &theme(legend.position = "bottom", legend.box.margin = margin(-5,0,-5,0))

pdf(file.path(outputFolder,"switchover_hdacs.pdf"), width = 2.5, height = 2.5)
l1 / l2
dev.off()

pdf(file.path(outputFolder,"switchover_nmda.pdf"), width = 2.5, height = 2.5)
l3 / l4
dev.off()

# NMDA western blot quantification --------------------------------------------------------------------------------------------------------------
wb_quant <- read_excel(file.path("data","NMDAR2AB_WB.xlsx"), range = "Sheet1!H2:K15") %>% 
  rename(GRIN2A = R2A, GRIN2B = R2B) %>% 
  gather("gene", "value", GRIN2A, GRIN2B ) 
p_wb <-   ggplot(wb_quant, aes(x = Day, y = value, color = gene)) +
  stat_smooth(geom = "line", alpha = 0.5, size =1,
              method = "loess", se = FALSE, formula = 'y ~ x',
              show.legend=T) +
  geom_point(size = 1) +
  labs(y = expression("GRIN / "*beta*"-actin ratio (AU)"),
       x = "Differentiation Day") +
  scale_color_brewer(palette = "Set1", name = "") +
  theme(legend.position = "none",
        legend.background = element_blank(),
        legend.key.width =  unit(0, "lines"),
        legend.margin = margin(l = -2, r = 5)
  ) +
  facet_wrap(~gene, scale = "free_y")

pdf(file.path(outputFolder,"nmda_wb_graph.pdf"), w = 2.5, h = 1.45)
p_wb
dev.off()

# RNA editing -----------------------------------------------------------------------------------------------------
load(file.path("data","RNAediting.rdata"))
colors2labels <- function(mod_colors, new_prefix){
  color_order <- WGCNA::labels2colors(1:length(mod_colors))
  data.frame(original = mod_colors) %>%
    arrange(match(original, color_order)) %>%
    rownames_to_column("id") %>%
    mutate(id = ifelse(original == "grey", 0, id)) %>%
    mutate(mod_name = paste0(new_prefix, id)) %>%
    dplyr::select(-id)
}

bs_mod_names <- colors2labels(unique(bs_modMember$moduleColor), "BSeditM")
hcs_mod_names <- colors2labels(unique(hcs_modMember$moduleColor), "hCSeditM")
p_bs_ME <- bs_ME %>%
  rownames_to_column("SampleID") %>%
  gather("Module", "ME", -SampleID) %>%
  left_join(datMetaCortical, by = "SampleID") %>%
  mutate(original = gsub("ME", "", Module)) %>%
  left_join(bs_mod_names, by = "original") %>%
  ggplot(aes(x = Period, y = ME)) +
  geom_vline(xintercept = 7.5 , linetype = 1, color = "grey 70") +
  geom_smooth(method = "loess", span = 0.75, color = "grey10",size = 0.5) +
  geom_point(size = 0.3) +
  labs(x = "BrainSpan Stage", y = "Module eigengene") +
  scale_x_continuous(breaks = seq(2,15,1)) +
  facet_wrap(~mod_name, scales = "free_y") +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_text(hjust = 0, size = 6),
        axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5)) 
pdf(file.path(outputFolder,"RNAedit_bs_mod_traj.pdf"), w = 3, h = 1.25)
p_bs_ME
dev.off()


edit_module_pres  <- mp_edit$preservation$Z[[1]][[2]] %>%
  rownames_to_column("original") %>%
  filter(original != "gold") %>%
  left_join(bs_mod_names, by = "original") %>%
  ggplot(aes(x = mod_name, y = Zsummary.pres)) +
  geom_col(fill = "grey70", width = 0.4) +
  labs(y = "Z summary", x = "Module") +
  scale_y_continuous(limits = c(0,10)) +
  facet_wrap(~"Module Preservation")


pdf(file.path(outputFolder,"RNAedit_bs_mod_pres.pdf"), w = 1.25, h = 1.25)
edit_module_pres
dev.off()


p_hcs_ME <- hcs_ME %>%
  rownames_to_column("SampleID") %>%
  gather("Module", "ME", -SampleID) %>%
  filter(Module != "MEgrey") %>%
  left_join(datMeta, by = "SampleID") %>%
  mutate(original = gsub("ME", "", Module)) %>%
  left_join(hcs_mod_names, by = "original") %>%
  filter(!is.na(Day.grouped)) %>%
  ggplot(aes(x = Differentiation.day, y = ME)) +
  geom_rect(aes(xmin = 250, xmax = 300, ymin = -Inf, ymax = Inf), color = "grey90", fill = "grey90") +
  geom_smooth(method = "loess", span = 0.75, color = "grey10", size = 0.5) +
  geom_point(size = 0.3) +
  scale_x_continuous(breaks = seq(0,700,100)) +
  labs(x = "Differentiation day", y = "Module eigengene") +
  facet_wrap(~mod_name, scales = "free_y", nrow = 1) +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_text(hjust = 0),
        axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5))

pdf(file.path(outputFolder,"RNAedit_hcs_mod_traj.pdf"), w = 4, h = 1.25)
p_hcs_ME +
  theme(axis.title.x = element_text(vjust=-1.2))
dev.off()


bs_matrix <-  sapply(unique(bs_modMember$moduleColor), function(mod){
  mod_membership <- 1*(bs_modMember$moduleColor == mod)
  names(mod_membership) <- bs_modMember$editing_site
  mod_membership
})

hcs_matrix <- unique(hcs_modMember$moduleColor) %>%
  .[. != "grey"] %>%
  sapply(function(mod){
    mod_membership <- 1*(hcs_modMember$moduleColor == mod)
    names(mod_membership) <- hcs_modMember$editing_site
    mod_membership
  })

idx <- intersect(rownames(bs_matrix) , rownames(hcs_matrix))
module_or <- geneSetEnrichment(bs_matrix[idx,], hcs_matrix[idx,])

or_mat <- module_or$dispMat %>%
  as.data.frame() %>%
  rownames_to_column("brainSpan") %>%
  gather("hCS", "log2(OR)", -brainSpan)
tmat <- module_or$tMat %>%
  as.data.frame() %>%
  setNames(colnames(module_or$dispMat )) %>%
  mutate(brainSpan = rownames(module_or$dispMat))%>%
  gather("hCS", "text", -brainSpan)

fdr_mat <- module_or$FDRmat %>%
  as.data.frame() %>%
  setNames(colnames(module_or$dispMat )) %>%
  mutate(brainSpan = rownames(module_or$dispMat))%>%
  gather("hCS", "fdr", -brainSpan)


or_plot_data <- or_mat %>%
  left_join(tmat, by = c("brainSpan", "hCS")) %>%
  left_join(fdr_mat, by = c("brainSpan", "hCS")) %>%
  left_join(hcs_mod_names, by = c("hCS"  = "original")) %>%
  dplyr::select(-hCS) %>%
  rename(hCS = mod_name) %>%
  left_join(bs_mod_names, by = c("brainSpan"  = "original")) %>%
  dplyr::select(-brainSpan) %>%
  rename(brainSpan = mod_name) %>%
  mutate(sig = case_when(fdr < 0.005 ~ "***",
                         fdr < 0.01 ~ "**",
                         fdr< 0.05 ~ "*",
                         TRUE ~ "")) %>%
  mutate(text = ifelse(text == "", text, paste(text, sig, sep = "\n"))) %>%
  mutate(hCS = factor(hCS, levels = str_sort(unique(hCS), decreasing = T, numeric = T)))


or_plot <- ggplot(or_plot_data,aes(y = hCS, x = brainSpan, fill = `log2(OR)`, label = text)) +
  geom_tile() +
  geom_text(size = 2) +
  scale_fill_gradient(limits = c(0,NA), na.value = NA, low = "white", high = "blue2", 
                      name= expression(paste(log[2](OR)))) +
  theme(axis.text = element_text(size = 6),
        axis.title.y = element_text(vjust = 2), 
        legend.key.width = unit(0.1 ,"in")) +
  labs(x = "BrainSpan")
pdf(file.path(outputFolder, "editing_or_mods.pdf"), w = 2.25, h = 1.75)
or_plot
dev.off()

# editing genes ---------------------------------------------------------------------------------------------------
editing_genes <- c("ADAR","ADARB1", "ADARB2", "FXR1", "FMR1")
eg_hcs <- plot_genes(datExpr_reg_batch, datMeta, editing_genes) %+%
  facet_wrap(~hgnc_symbol, nrow = 1, scales = "free_y") +
  ggtitle("hCS expression")
eg_bs <- plot_genes(datExpr_reg_cortical, datMetaCortical, editing_genes) %+%
  facet_wrap(~hgnc_symbol, nrow = 1, scales = "free_y") +
  theme(strip.text = element_blank()) +
  ggtitle("BrainSpan expression")

editing_genes_plot <- (eg_hcs + labs(y = "")) +
  plot_spacer() +
  (eg_bs + labs(y =""))*
  theme(axis.title.y = element_text(size = 11),
        plot.title = element_blank()) &
  theme(strip.text = element_text(size = 6, face = "italic")) 

pdf(file.path(outputFolder, "editing_genes_expr.pdf"), w = 5, h = 2.25)
editing_genes_plot+plot_layout(heights = c(1,0.01,1))
dev.off()

current_genes_ensembl <- geneAnno %>%
  filter(hgnc_symbol %in% editing_genes) %>%
  dplyr::select(ensembl_gene_id, hgnc_symbol) %>%
  as.data.frame()
gene_expr <- datExpr_reg_batch[current_genes_ensembl$ensembl_gene_id, ] %>%
  t()
colnames(gene_expr) <- current_genes_ensembl$hgnc_symbol
me <- hcs_ME[rownames(gene_expr),!names(hcs_ME) %in% c("MEgrey")] %>%
  as.matrix

colnames(me) <- hcs_mod_names %>%
  filter(original != "grey") %>%
  arrange(match(original,gsub("ME","", colnames(me)))) %>%
  dplyr::select(mod_name) %>%
  unlist
me <- me[, str_sort(colnames(me), numeric = T)]
cor_mat <- psych::corr.test(gene_expr, me, adjust= "BH", method = "s")
pdf(file.path(outputFolder, "editing_ME_enz_cor.pdf"), w = 3, h = 3.5)
corrplot::corrplot(
  cor_mat$r,
  p.mat = cor_mat$p,
  method = "ellipse",
  tl.col	= "black",
  mar = c(0, 0, 3, 0),
  cl.align.text = "l",
  cl.cex = 1.3 ,
  cl.length = 3,
  tl.cex = 1.3,
  insig = "label_sig",
  sig.level = c(.001, .01, .05),
  pch.cex = 1, 
  pch.col = "black",
  cl.pos="n", 
  addgrid = F
)
dev.off()

pdf(file.path(outputFolder, "editing_ME_enz_cor_legend.pdf"), w = 1.5, h = 3)
color.bar(lut = colorRampPalette(c("#67001F", "#B2182B", "#D6604D", 
                                   "#F4A582", "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE", 
                                   "#4393C3", "#2166AC", "#053061"))(200),
          min = -1,max = 1,nticks = 3, title = " ")
dev.off()

# E-phys ----------------------------------------------------------------------------------------------------------
epyhs_file <- file.path("data","NMDA_ephys.pzfx")

response_amp <- read_pzfx(path = epyhs_file, table = "Response Amplitudes")
response_amplitude_baseline <- read_pzfx(path = epyhs_file, table = "Response Amplitude Baseline") %>%
  gather("measurement","Amp", -DIV, na.rm = T )

r_ephys <- cor.test(response_amplitude_baseline$Amp, response_amplitude_baseline$DIV)
r_p_val_ephys <- formatC(r_ephys$p.value, format = "e", digits = 2) %>%
  str_split("e") %>%
  unlist

pdf(file.path(outputFolder, "e-phsy_baseline.pdf"), h = 1.4, w = 1.8)
ggplot(response_amplitude_baseline, aes(x = DIV, y = Amp)) +
  geom_smooth(method = "lm", color = alpha("black", 0.4), se = F) +
  geom_point(size = 1) +
  geom_vline(xintercept = 250, color = "black", linetype = 2) +
  annotate("segment", x = 240, xend = 140, y = 300, yend = 300,
           arrow = arrow(length = unit(0.05, "in"), type = "closed")) +
  annotate("segment", x = 260, xend = 360, y = 300, yend = 300,
           arrow = arrow(length = unit(0.05, "in"), type = "closed")) +
  annotate("text", label = "Early", x = -20, y = 300, hjust = 0,size = 2,) +
  annotate("text", label = "Late", x = 390, y = 300, hjust = 0,size = 2,) +
  annotate(geom = "text", x = -20, y = 255, size = 2, hjust = 0,label = paste("r =", round(r_ephys$estimate,2))) +
  annotate(geom = "text", x = -20, y = 220, size = 2, hjust = 0,
           label = substitute(
             paste("p = ", part1,"e"^part2),
             list(
               part1 = r_p_val_ephys[1],
               part2 = as.numeric(r_p_val_ephys[2]))
           ),
  ) +
  labs(x = "Differentiation Day", y = "NMDA Amp (pA)")
dev.off()

ifn_all <- readxl::read_excel(file.path("data" ,"NMDA_IFN Responses.xlsx")) %>%
  setNames(.[2,]) %>%
  mutate(stage = ifelse(DIV %in% c("Young", "Old"), DIV, NA)) %>%
  fill(stage) %>%
  mutate(DIV = as.numeric(DIV)) %>%
  filter(!is.na(DIV)) %>%
  mutate_at(vars(matches("pA|%")), as.numeric) %>%
  rowid_to_column("SampleID") %>%
  mutate(SampleID = letters[SampleID]) %>%
  gather("condition","pA", matches("pA")) %>% 
  mutate(log_pa = log(pA)) %>% 
  mutate(`Change` = ifelse(condition == "Baseline (pA)", NA_real_, `% Change`))

ifn_for_beta <- ifn_all %>% 
  filter(!is.na(Change)) %>% 
  mutate(Change = Change/100)
logit_model <- betareg::betareg(Change ~ DIV, data = ifn_for_beta, link = "logit") 
summary(logit_model)
logit_pval <- summary(logit_model)$coefficients$mean[2,4] %>% 
  formatC(format = "e", digits = 2) %>%
  str_split("e") %>%
  unlist



pdf(file.path(outputFolder, "e-phsy_IFN_reduction.pdf"), h = 1.4, w = 1.8)
color_for_ifn <- brewer.pal(3,"Set1")[3]
ggplot(ifn_for_beta, aes(x = DIV, y = Change*100)) +
  geom_line(aes(y = predict(logit_model, ifn_for_beta)*100), size = 1, color = alpha(color_for_ifn, 0.5)) + 
  geom_point(size = 1, color = color_for_ifn) +
  geom_vline(xintercept = 250, color = "black", linetype = 2) +
  annotate("segment", x = 240, xend = 140, y = 100, yend = 100,
           arrow = arrow(length = unit(0.05, "in"), type = "closed")) +
  annotate("segment", x = 260, xend = 360, y = 100, yend = 100,
           arrow = arrow(length = unit(0.05, "in"), type = "closed")) +
  annotate("text", label = "Early", x = -20, y = 100, hjust = 0,size = 2,) +
  annotate("text", label = "Late", x = 450, y = 100, hjust = 0,size = 2,) +
  annotate(geom = "text", x = 570, y = 88, size = 2, hjust = 1, label = substitute(
    paste("p = ", part1,"e"^part2),
    list(
      part1 = logit_pval[1],
      part2 = as.numeric(logit_pval[2]))
  )
  ) +
  labs(x = "Differentiation Day", y = "% Change in Amp") +
  scale_y_continuous(breaks = seq(0,100, 25), labels = seq(0,100,25)) +
  expand_limits(y = 0)
dev.off()

# disease genes ---------------------------------------------------------------------------------------------------
samples_by_day <- datMeta %>% 
  arrange(Differentiation.day.original) %>% 
  as.data.frame() %>% 
  `rownames<-`(.$SampleID)
lgd = Legend(labels = unique(samples_by_day$Day.grouped), 
             legend_gp = gpar(fill = paste0("grey", seq(93,1,-8))),
             labels_gp = gpar(fontsize = 6),
             title_gp = gpar(fontsize = 6), 
             title_gap = unit(2, "pt"),
             gap = unit(1, "pt"),
             by_row = F,
             title = "Differentiation day", nrow = 2, title_position = "topleft")
pdf(file.path(outputFolder,  "day_legend.pdf"), w = 3,h = 1)
draw(lgd)
dev.off()

plot_disorder_genes <- function(expr_disorder, k, condition, wrap = 30, highlight_genes = NA){
  
  expr_heatmap <- t(apply(expr_disorder[,samples_by_day$SampleID], 1, scale))
  
  ha <- HeatmapAnnotation(
    "Differentiation\nday" = samples_by_day$Day.grouped, 
    col = list("Differentiation\nday" = setNames(paste0("grey", seq(93,1,-8)),unique(samples_by_day$Day.grouped))), 
    annotation_name_gp = gpar(fontsize = 6), 
    show_legend = F, 
    annotation_name_side = "left"
  )
  
  expr_clustered <- expr_heatmap %>% 
    dist() %>% 
    hclust 
  
  gene_cluster_id <- dendextend::cutree(expr_clustered, k = k, order_clusters_as_data = FALSE) %>% 
    as.data.frame() %>% 
    setNames("cluster") %>% 
    rownames_to_column("ensembl_gene_id") %>% 
    left_join(geneAnno, by = "ensembl_gene_id")
  assign(paste0(condition, "gene_clusters"), gene_cluster_id[,c(2,1,3)], pos = .GlobalEnv)
  cluster_eigens <- WGCNA::moduleEigengenes(t(expr_disorder[gene_cluster_id$ensembl_gene_id,]),
                                            gene_cluster_id$cluster)$eigengenes
  assign(paste0(condition, "cluster_eigens"), cluster_eigens, pos = .GlobalEnv)
  
  hub_genes <- map_df(set_names(1:k), function(i){
    cluster_genes <- gene_cluster_id %>%
      filter(cluster == i)
    cluster_expr <- expr_disorder[cluster_genes$ensembl_gene_id,]
    
    cor(t(cluster_expr), cluster_eigens[,i, drop = F]) %>%
      as.data.frame() %>%
      setNames("ME_cor") %>%
      rownames_to_column("ensembl_gene_id") %>%
      left_join(geneAnno, by = "ensembl_gene_id") %>%
      arrange(desc(ME_cor)) %>%
      slice(1:5)
  })
  
  # add * to familal genes in AD and PD
  if(!anyNA(highlight_genes)){
    hl_genes <-  geneAnno %>% 
      filter(hgnc_symbol %in% highlight_genes) %>% 
      mutate(hgnc_symbol = paste0(hgnc_symbol,"*")) %>% 
      bind_rows(hub_genes) %>% 
      mutate(hgnc_symbol = ifelse(duplicated(ensembl_gene_id, fromLast = T), paste0(hgnc_symbol, "^"),hgnc_symbol ))%>% 
      filter(!duplicated(ensembl_gene_id))
  } else {
    hl_genes <- hub_genes
  }
  
  rown_anno <- rowAnnotation(link = anno_mark(at = match(hl_genes$ensembl_gene_id, rownames(expr_heatmap)),
                                              labels = hl_genes$hgnc_symbol, 
                                              labels_gp = gpar(fontsize = 6, fontface = "italic"), 
                                              padding = unit(1, "pt"))) 
  hm_h <- ifelse(nrow(expr_disorder) < 25, 1, 2)
  
  ht <- Heatmap(expr_heatmap, 
                height = unit(hm_h, "in"),
                width  = unit(0.9, "in"),
                cluster_columns = F, 
                col = c("green","black", "magenta"), 
                row_gap = unit(1, "pt"),
                column_gap = unit(1, "pt"),
                right_annotation = rown_anno,
                row_split = k, 
                column_split = factor(c(rep( "Prenatal", 34), rep( "Postnatal",28)), levels = c("Prenatal", "Postnatal")),
                show_row_names = F,  
                top_annotation = ha, 
                column_title_gp = gpar(fontsize = 6),
                row_title_gp = gpar(fontsize = 6), 
                heatmap_legend_param = list(
                  legend_width = unit(1,"in"), 
                  title = "Scaled expression",
                  title_gp = gpar(fontsize = 6) , 
                  labels_gp = gpar(fontsize = 6), 
                  title_position = "leftcenter", 
                  grid_height = unit(0.075, "in"), 
                  direction = "horizontal"
                ),
  )
  pdf(file.path(outputFolder, paste0(condition, "_hi_conf_genes.pdf")), w = 3,h = 3.5)
  draw(ht, heatmap_legend_side = "bottom")
  dev.off()
  
  cluster_eigens_tidy <- cluster_eigens %>% 
    as.data.frame() %>% 
    mutate(SampleID = colnames(expr_disorder)) %>%  
    gather("cluster", "eigengene", -SampleID) %>% 
    left_join(datMeta, by = "SampleID") %>% 
    mutate(cluster = gsub("ME", "cluster ", cluster))
  
  cluster_anno <- data.frame(Differentiation.day.original = 30, 
                             eigengene = max(cluster_eigens_tidy$eigengene), 
                             label = paste0("C",  1:k), 
                             cluster = str_sort(unique(cluster_eigens_tidy$cluster, numeric = T)))
  
  p1 <- ggplot(cluster_eigens_tidy, aes(x = Differentiation.day.original, y = eigengene)) +
    geom_jitter(width = 0.1, color = "grey50", size = 0.3) +
    geom_smooth(method = 'loess', formula = y~x, color = "grey10", size = 0.5) +
    scale_y_continuous(breaks=seq(-4, 4, by = 0.4)) +
    facet_wrap(~cluster, ncol = 1, strip.position = "top") +
    geom_text(data = cluster_anno, aes(label = label), size = 2, vjust = 1, nudge_x = 80, nudge_y = 0.1) +
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5),
          legend.position = "none",
          strip.text = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 6, margin = margin(t = 0, r = 0, b = 0, l = 0)),
          strip.background = element_rect(colour = NA), 
          axis.title.x = element_text(hjust = 1),
          panel.spacing = unit(8, "points")) +
    labs(x = "Differentiation day", y = "eigengene") 
  
  h1 <- ifelse(k>4, k*0.55, k*0.55+0.2)
  pdf(file.path(outputFolder, paste0(condition, "_cluster_eigen_gene.pdf")), h = h1, w = 1)
  print(p1)
  dev.off()
  
}
#SCZ associated genes from Wang, D., et al. Comprehensive functional genomic resource and integrative model for the human brain. Science 362, eaat8464 (2018).
scz_hi_conf <- read.csv(file.path("data","disorder_genes", "SCZ.csv"))
expressed_hi_conf <- intersect(scz_hi_conf$ensembl_names,rownames(datExpr_reg_batch))
scz_hi_conf_expr <- datExpr_reg_batch[expressed_hi_conf,]

plot_disorder_genes(scz_hi_conf_expr, 5, "SCZ", wrap = 30)

#ASD genes downloaded from SFARI database release 06-23-2020
asd_sfari_all <- read.csv(file.path("data","disorder_genes","ASD.csv"))
asd_sfari <- asd_sfari_all %>%
  filter(gene.score <= 2 | syndromic == 1) %>%
  dplyr::rename(hgnc_symbol = gene.symbol) %>%
  left_join(geneAnno[,1:2], by = "hgnc_symbol")

expressed_sfari <- intersect(asd_sfari$ensembl_gene_id,rownames(datExpr_reg_batch))
sfari_expr <- datExpr_reg_batch[expressed_sfari,]
plot_disorder_genes(sfari_expr, 5, "ASD", wrap = 30)

#ID genes from Parikshak, N. N. et al. Integrative functional genomic analyses implicate specific molecular pathways and circuits in autism. Cell 155, 1008–1021 (2013).
id_genes <- read_excel(file.path("data","disorder_genes", "ID.xlsx"),  skip = 2, sheet = 1, ) %>%
  dplyr::select("ID All") %>%
  rename(hgnc_symbol = "ID All") %>%
  filter(!is.na(hgnc_symbol)) %>%
  left_join(geneAnno, by = "hgnc_symbol") %>%
  filter(!is.na(ensembl_gene_id))

expressed_id <- intersect(id_genes$ensembl_gene_id,rownames(datExpr_reg_batch))
id_expr <- datExpr_reg_batch[expressed_id,]

plot_disorder_genes(id_expr, 4, "id", wrap = 30)

#AD genes from  Jansen, I. E. et al. Genome-wide meta-analysis identifies new loci and functional pathways influencing Alzheimer’s disease risk. Nat. Genet. 51, 404–413 (2019).
#and Cacace, R., Sleegers, K. & Van Broeckhoven, C. Molecular genetics of early-onset Alzheimer’s disease revisited. Alzheimers Dement. 12, 733–748 (2016).
ad_genes_magma <- readxl::read_excel(file.path("data","disorder_genes","AD.xlsx"), sheet = 18, skip = 3) %>%
  left_join(geneAnno, by = c("Gene Name" = "hgnc_symbol"))
ad_genes_fuma <- readxl::read_excel(file.path("data","disorder_genes","AD.xlsx"), sheet = 13, skip = 3) %>% 
  left_join(geneAnno, by = c("Gene" = "hgnc_symbol")) %>% 
  mutate(to_keep = (`posMap #SNPs` > 0) + (`eqtlMap #SNPs` > 0)  +(ciMap == "Yes")) %>% 
  filter(to_keep >= 2)
mendelian_ad_genes <- geneAnno %>% 
  filter(hgnc_symbol %in% c("APOE","PSEN1", "PSEN2", "APP")) 
all_ad_genes <- unique(c(ad_genes_fuma$`Ensembl ID`, ad_genes_magma$ensembl_gene_id, mendelian_ad_genes$ensembl_gene_id))

expressed_ad <- intersect(all_ad_genes,rownames(datExpr_reg_batch))
ad_expr <- datExpr_reg_batch[expressed_ad,]

plot_disorder_genes(ad_expr, 4, "AD", wrap= 30, highlight_genes = c("APOE","PSEN1", "PSEN2", "APP"))


#PD genes from Chang, D. et al. A meta-analysis of genome-wide association studies identifies 17 new Parkinson’s disease risk loci. Nat. Genet. 49, 1511–1516 (2017).
#and Farrer, M. J. Genetics of Parkinson disease: paradigm shifts and future prospects. Nat. Rev. Genet. 7, 306–318 (2006).

pd_genes_gwas <- readxl::read_excel(file.path("data","disorder_genes","PD.xlsx")) %>% 
  left_join(geneAnno, by = c("Candidate Gene" = "hgnc_symbol"))%>% 
  filter(!is.na(ensembl_gene_id) & !duplicated(ensembl_gene_id)) %>% 
  dplyr::select(ensembl_gene_id) %>% 
  unlist() %>% 
  unname()
#familial pd gene from https://www.nature.com/articles/nrg1831/tables/2
pd_familial_genes <- geneAnno %>% 
  filter(hgnc_symbol %in% c("PINK1", "SNCA", "LRRK2", "PRKN", "UCHL1", "PARK7")) %>% 
  dplyr::select(ensembl_gene_id) %>% 
  unlist() %>% 
  unname()

pd_genes <- union(pd_genes_gwas , pd_familial_genes)

expressed_pd <- intersect(pd_genes,rownames(datExpr_reg_batch))
pd_expr <- datExpr_reg_batch[expressed_pd,]

plot_disorder_genes(pd_expr, 4, "pd", wrap = 35, highlight_genes = c("PINK1", "SNCA", "LRRK2", "PRKN", "UCHL1", "PARK7"))


#FTD-PSP 
# FTD genes from Greaves, C. V. & Rohrer, J. D. An update on genetic frontotemporal dementia. J. Neurol. 266, 2075–2086 (2019).
ftd <- c("C9orf72","GRN","MAPT", "VCP", "CHMP2B","TARDBP","FUS","SQSTM1","CHCHD10","TBK1","OPTN","CCNF","TIA1") 
ftd_genes <- geneAnno %>% 
  filter( hgnc_symbol %in% ftd)
#PSP genes from Chen, J. A. et al. Joint genome-wide association study of progressive supranuclear palsy identifies novel susceptibility loci and genetic correlation to neurodegenerative diseases. Mol. Neurodegener. 13, 41 (2018).
psp <- c("MAPT","MOBP","STX6","RUNX2","SLCO1A2","DUSP10","SP1","ASAP1","WDR63")
psp_genes <- geneAnno %>% 
  filter( hgnc_symbol %in% psp)

psp.ftd_genes <- union(ftd_genes$ensembl_gene_id, psp_genes$ensembl_gene_id)

expressed_psp.ftd <- intersect(psp.ftd_genes,rownames(datExpr_reg_batch))
psp.ftd_expr <- datExpr_reg_batch[expressed_psp.ftd,]


plot_disorder_genes(psp.ftd_expr, 2, "psp.ftd", wrap = 35)

#epilepsy genes from Polioudakis, D. et al. A single-cell transcriptomic atlas of human neocortical development during mid-gestation. Neuron 103, 785–801.e788 (2019).
epilepsy_genes <- read.delim(file.path("data","disorder_genes", "Epilepsy.txt")) %>%
  left_join(geneAnno, by = c("Gene" = "hgnc_symbol"))

expressed_epilepsy <- intersect(epilepsy_genes$ensembl_gene_id,rownames(datExpr_reg_batch))
epilepsy_expr <- datExpr_reg_batch[expressed_epilepsy,]

plot_disorder_genes(epilepsy_expr, 3, "epilepsy", wrap = 20)

list("ASD gene clusters" = ASDgene_clusters, 
     "ID gene clusters" = idgene_clusters, 
     "SCZ gene clusters" = SCZgene_clusters, 
     "Alzheimer gene clusters" = ADgene_clusters, 
     "Parkinson gene clusters" = pdgene_clusters, 
     "PSP-FTD gene clusters" = psp.ftdgene_clusters, 
     "Epilepsy gene clusters" = epilepsygene_clusters) %>% 
  writexl::write_xlsx(file.path(outputFolder, "tables", "Supp table 5 - disorder gene clusters.xlsx"))
