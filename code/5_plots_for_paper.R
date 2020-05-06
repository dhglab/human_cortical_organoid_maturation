# plots for Paper

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


outputFolder <- file.path("output")

tables_folder <- file.path(outputFolder, "tables")
dir.create(outputFolder, showWarnings = F)
dir.create(tables_folder, showWarnings = F)

pub_theme <- theme_classic(base_size = 8) +
  theme(axis.text = element_text(size = 8, color = "black"),
        strip.background = element_rect(colour = NA),
        strip.text = element_text(size = 8),
        plot.title = element_text(size = 8),
        axis.title.x = element_text(vjust=-0.75),
        legend.text = element_text(size = 8))
theme_set(pub_theme)


#prdct <- readRDS(file = "inVitro_LT_methyaltion/analysis/predicted_birth_age.rds")
#load brain span data
load(file.path("data","brainSpan_all_star_cqn_noramlized.rdata"))
datMetaCortical <- datMeta[grep("C$",datMeta$Region),] %>%
  rownames_to_column("SampleID") %>%
  mutate(Period = as.numeric(levels(Period)[Period]))
#datExprRawCortical <- datExprRaw[rownames(datExprCQN$counts),datMetaCortical$SampleID]
datExpr_reg_cortical <- datExpr.reg[,datMetaCortical$SampleID]

load(file.path("data","processed_data.rdata"))
load(file.path("data","gene_annotation.Rdata"))

#deidentified_codes <- read.csv("analysis/rsem/1_longTerm/deidentified codes.csv")

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
  # tMat[FDRmat < 0.05] <- paste(tMat[FDRmat < 0.05], "*", sep = "")
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
    # geom_vline(xintercept = prdct["mean"] , linetype = 1, color = "grey 70") +
    # geom_vline(xintercept = prdct["lwr_CI"] , linetype = 3, color = "grey 70") +
    # geom_vline(xintercept = prdct["upr_CI"] , linetype = 3, color = "grey 70") +
    geom_rect(aes(xmin=250, xmax=300, ymin=-Inf, ymax=Inf), color = "grey90", fill = "grey90") +
    geom_jitter(width = 0.1, aes(color = hgnc_symbol, shape = hgnc_symbol), size = 1) +
    geom_smooth(aes(color = hgnc_symbol, linetype = hgnc_symbol), method = 'loess', formula = y~x, size =0.6) +
    scale_y_continuous(breaks = seq(-10,10,2)) +
    scale_x_continuous(breaks = seq(0,700,100)) +
    scale_color_brewer(palette = "Set1", name = "") +
    scale_linetype(name = "") +
    scale_shape(name = "")+
    theme_classic(base_size = 8) +
    labs(x = "Differentiation Day", y = "Scaled normalized expression") +
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5, size = 8),
          axis.text.y = element_text(size = 8))
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
    theme_classic(base_size = 8) +
    labs(x = "BrainSpan Stage", y = "Scaled normalized expression") +
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5, size = 8),
          axis.text.y = element_text(size = 8))  +
    geom_vline(xintercept = 7.5 , linetype = 1, color = "grey 70") +
    geom_jitter(width = 0.1, aes(color = hgnc_symbol, shape = hgnc_symbol), size = 1) +
    geom_smooth(aes(color = hgnc_symbol, linetype = hgnc_symbol), method = 'loess', formula = y~x, size =0.6) +
    scale_y_continuous(breaks = seq(-10,10,2)) +
    scale_x_continuous(breaks = seq(2,15,1)) +
    scale_color_brewer(palette = "Set1", name = "") +
    scale_linetype(name = "") +
    scale_shape(name = "")
}

plot_comparison <- function(genes){


  p1 <- plot_switch_hCS(genes) + theme(legend.position = "none")
  p2 <-   plot_switch_bs(genes) +
    guides(color = guide_legend(override.aes = list(fill = NA))) +
    theme(legend.key.size =  unit(.1, "in"), legend.position = "right")
  # wi = 8
  # hei = 3
  # pdf(paste0(outputFolder,"switchover_",genes[1],"_",genes[2],".pdf"), width = wi, height = hei)
  # print(p1+p2)
  # dev.off()
  list(p1,p2)
}

singleGenePlot <- function(g){
  genes <- toupper(g)
  single_genes_anno <- geneAnno %>%
    filter(hgnc_symbol %in% genes)
  plot_data <- datExpr_reg_batch %>%
    as.data.frame() %>%
    rownames_to_column("ensembl_gene_id") %>%
    inner_join(single_genes_anno, by = "ensembl_gene_id") %>%
    mutate(hgnc_symbol = factor(hgnc_symbol, levels = genes)) %>%
    gather("SampleID", "value", matches("MEF|Pool")) %>%
    left_join(datMeta, by = "SampleID")
  gene_display <- set_names(unique(as.character(plot_data$hgnc_symbol))) %>%
    ifelse(. == "POU3F2",paste(.,"\n(BRN2)") ,.) %>%
    ifelse(. == "BCL11B",paste(.,"\n(CTIP2)") ,.)
  print(gene_display)
  marker_plot <- ggplot(plot_data,aes(x = Differentiation.day.original, y = value)) +
    # geom_vline(xintercept = prdct["mean"] , linetype = 1, color = "grey 70") +
    # geom_vline(xintercept = prdct["lwr_CI"] , linetype = 3, color = "grey 70") +
    # geom_vline(xintercept = prdct["upr_CI"] , linetype = 3, color = "grey 70") +
    geom_rect(aes(xmin=250, xmax=300, ymin=-Inf, ymax=Inf), color = "grey90", fill = "grey90") +
    # geom_vline(xintercept = 250 , linetype = 2, color = "grey 70") +
    # geom_vline(xintercept = 300 , linetype = 2, color = "grey 70") +
    geom_jitter(width = 0.1, color = "grey50", size = 0.3) +
    facet_wrap(~hgnc_symbol,scales = "free_y", labeller = as_labeller(gene_display)) +
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 8, margin = margin(t = 0, r = 0, b = 0, l = 0)),
          strip.background = element_rect(colour = NA)) +
    geom_smooth(method = 'loess', formula = y~x, color = "grey10", size = 0.5) +
    labs(x = "Differentiation day", y = "Normalized expression")
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
    theme(axis.text.y = element_text(size = 8),
          axis.text.x = element_text(size = 8, angle = -90, hjust = 0, vjust = 0.5),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 8, margin = margin(t = 0, r = 0, b = 0, l = 0)),
          strip.background = element_rect(colour = NA))
}


# RRHO ------------------------------------------------------------------------------------------------------------


hypmat_file <- load(file.path(outputFolder, "data", "hypermatAll.rdata"))

comparisons <- names(brainSpanVoom)
nr <- dim(hypermat.all)[1]
nc <- dim(hypermat.all)[2]

jet.colors = colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"));
colormap = jet.colors(100)


minhypermat = min(hypermat.all,na.rm=TRUE);
maxhypermat = max(hypermat.all,na.rm=TRUE)
png(file.path(outputFolder,"iPSC_brainSpan_day25_RRHOMap.png"),width=2.5,height=2.5,units="in",res=300);
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
  mtext(title,2,2.7,cex=1);
  axis(2, round(ticks,0), las=1,cex.lab=1);
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min;
    rect(0,y,10,y+1/scale, col=lut[i], border=NA);
  }
}
png(file.path(outputFolder,"iPSC_brainSpanAll_RRHOMap.png",sep = ""),width = 1.5,height = 6,units = "in", res = 300)
color.bar(lut = jet.colors(100),min = minhypermat,max = maxhypermat,nticks = 6, title = "-log10(Nominal P-value)")
dev.off()

# DNAmAge ---------------------------------------------------------------------------------------------------------
load(file.path("data", "DNAmAge_data.rdata"))
datSample <- datSampleRaw %>%

  group_by(External.Sample.ID, Differentiation.day,Line, Sex ) %>%
  summarise(DNAmAge = mean(DNAmAge))


r <- cor.test(datSample$DNAmAge,datSample$Differentiation.day)

r_p_val <- formatC(r$p.value, format = "e", digits = 2) %>%
  str_split("e") %>%
  unlist


p1 <- ggplot(datSample, aes(x = Differentiation.day, y = DNAmAge)) +
  geom_smooth(method = "lm", color = "grey40") +
  annotate(geom = "text", x = 0, y = 8.5, size = 2.8, hjust = 0,label = paste("r =", round(r$estimate,2))) +
  annotate(geom = "text", x = 0, y = 7.5, size = 2.8, hjust = 0,
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
    labs(x = "Differentiation Day") +
  guides(color = guide_legend(keywidth = 0.1, keyheight = 0.05,default.unit = "inch")) +
  theme(legend.position = "bottom",
        plot.margin = margin(0,0,0,0),
        legend.margin = margin(0,20,0,0),
        axis.text.x = element_text(angle = -90, hjust  = 0, vjust = 0.5, size = 8, color = "black"),
        axis.text.y = element_text(size = 8, color = "black"),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))

pdf(file = file.path(outputFolder,"DNAmAge_scatter.pdf"), w = 2.5, h = 2.5, useDingbats = FALSE)
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
  theme_minimal(base_size = 8) +
  theme(axis.text = element_text(size = 8),
        axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5),
        axis.title.x = element_text(margin = margin(t = 5)),
        axis.title.y = element_text(margin = margin(r = 5)),
        legend.position = "none",
        plot.title = element_text(size = 8)) +
  labs(title = "RNA Expression")

p_meth_samples <- ggplot(datSample, aes(x = Differentiation.day, y = Line, shape = Sex)) +
  geom_jitter(width = 0, height = 0.15, size = 1) +
  labs(x = "Differentiation day", y = "Individual") +
  scale_x_continuous(breaks = seq(0, 700, 100)) +
  scale_shape_manual(values = c(1,16)) +
  theme_minimal(base_size = 8) +
  theme(axis.text = element_text(size = 8),
        axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5),
        axis.title.x = element_text(margin = margin(t = 5)),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.title = element_text(size = 8)) +
  guides(shape = guide_legend( keywidth = 0.1, keyheight = 0.1,
    default.unit = "inch"
  )) +
  labs(title = "Methylation")
pdf(file.path(outputFolder, "Samples.pdf"), w = 4, h = 2)
p_rna_samples | p_meth_samples
dev.off()

datMeta_rnaseq <- datMeta %>%
  mutate(Induction = as.numeric(as.factor(Induction))) %>%
  mutate(batch = as.numeric(as.factor(batch))) %>%
  rowid_to_column(var = "Sample") %>%
  dplyr::select(Sample, CellLine, Differentiation = Induction,
         Differentiation.day = Differentiation.day.original,
         Day.grouped, Sex, Batch = batch, ethnicityPC1, ethnicityPC2,
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

pdf(file.path(outputFolder, "PC_cor.pdf"), h = 2.5, w = 4 )
#par(ps = 8)
corrplot::corrplot(
  t(expr_pc_rsqr),
  method = "ellipse",
  tl.col	= "black",
  mar = c(0, 0, 0, 0),
  cl.align.text = "l",
  tl.cex = 0.66,
  cl.cex = 0.66,
  cl.length = 3
)
dev.off()

PC_plot <- topPC.expr %>%
  as.data.frame() %>%
  rownames_to_column("SampleID") %>%
  left_join(datMeta, by = "SampleID")
pdf(file.path(outputFolder, "LT_PC_plot_1.pdf"), width = 3, height = 2)
ggplot(PC_plot,aes(x=PC1,y=PC2,color=Differentiation.day)) +
  geom_point(size=1) +
  theme(legend.position="right",
        legend.key.size = unit(0.1, "in"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        axis.text = element_text(color = "black", size = 8)) +
  # guides(col = guide_legend(size = 0.5)) +
  scale_color_viridis_c("Differentiation\nday") +
  scale_shape("Individual") +
  labs(x = xy_labs[1], y= xy_labs[2])
dev.off()

# Dendogram -------------------------------------------------------------------------------------------------------
datMeta1 <- datMeta %>%
  dplyr::select(Individual = CellLine, Differentiation.day, Sex, batch) %>%
  mutate_all(~(as.numeric(as.factor(.))))  %>%
  mutate(Individual  = numbers2colors(Individual, signed = FALSE, colors = viridis_pal()(6)),
         Sex = numbers2colors(Sex, signed = FALSE, colors = brewer.pal(6,"Paired")[5:6]),
         batch = numbers2colors(batch, signed = FALSE, colors = brewer.pal(4,"Paired")[3:4]),
         Differentiation.day = numbers2colors(Differentiation.day, signed = FALSE, colors =  colorRampPalette(c("white", "dodgerblue4"))( 20 )))

hc <- hclust(dist(t(datExpr_reg_batch)), method = "average")
pdf(file.path(outputFolder, "sample_dendogram.pdf"), width = 3.5, height = 2)
plotDendroAndColors( hc,
                     datMeta1,
                     hang = -1,
                     dendroLabels = FALSE, #paste(datMeta$CellLine,datMeta$Differentiation.day.original,sep = "-"),
                     groupLabels = gsub("\\."," ",names(datMeta1)),
                     main = NULL,
                     marAll = c(0, 5,  2, 0),
                     cex.dendroLabels =  0.5,
                     autoColorHeight = F,
                     colorHeight = 0.48,
                     cex.colorLabels = 0.666,
                     cex.rowText = 0.666,
                     cex.axis = 0.666,
                     cex.lab = 0.666,
                     ylab = NULL,
                     axes = F
)

dev.off()



# Single Gene Expression ------------------------------------------------------------------------------------------
cell_type_genes2 <- list("Radial Glia" = c("PAX6", "GLI3"),
                         "Inermediate Progenitors" = c("EOMES","NEUROD1"),
                         #"Neurons" = c("MAP2","NCAM1"),

                         "Upper Layer Neurons" = c("SATB2", "POU3F2"),
                         "Deep Layer Neurons" = c("TBR1","BCL11B"),
                         "Early stage Astrocytes" = c("TMSB15A", "NNAT"),
                         "Late stage Astrocytes" = c("AQP4","GFAP")
)



marker_plot <- map2(cell_type_genes2,names(cell_type_genes2),function(g,n){
  singleGenePlot(g) +
    ggtitle(n)
})

pdf(file.path(outputFolder,"cell_type_markers_single_genes.pdf"),width = 3.5, height = 4)
(marker_plot[[1]] + theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
    marker_plot[[2]] + theme(axis.title.x = element_blank(), axis.title.y = element_blank())) /
  (marker_plot[[3]] + theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
     marker_plot[[4]] + theme(axis.title.x = element_blank(), axis.title.y = element_blank())) /
  (marker_plot[[5]] + theme(axis.title.y = element_blank()) +
     marker_plot[[6]] + theme(axis.title.y = element_blank())) &
  theme(strip.text = element_text(size = 8, face = "italic"))
dev.off()


# astrocyte transition --------------------------------------------------------------------------------------------
sloan_data_raw <- read_excel(file.path(outputFolder, "data","Sloan2018_astrocyte_genes.xlsx"))

sloan_data <- sloan_data_raw %>%
  gather("stage","hgnc_symbol") %>%
  left_join(geneAnno, by = "hgnc_symbol") %>%
  filter(ensembl_gene_id %in% rownames(datExpr_reg_batch))

datMeta_by_day <- datMeta %>%
  arrange(Differentiation.day.original)

sample_by_day <- datMeta_by_day %>%
  dplyr::select(SampleID) %>%
  unlist()

astro_expr <- datExpr_reg_batch[sloan_data$ensembl_gene_id,sample_by_day]

all(rownames(astro_expr) == sloan_data$ensembl_gene_id)

astro_median_expr <- map_df(set_names(sort(unique(datMeta$Day.grouped))), function(day){
  current_samples <- datMeta %>%
    filter(Day.grouped == day) %>%
    dplyr::select(SampleID) %>%
    unlist()
  apply(astro_expr[,current_samples],1,median)
}) %>%
  as.matrix() %>%
  `rownames<-`(rownames(astro_expr))

row_anno <- data.frame("Astrocyte_genes" = sloan_data$stage) %>%
  `row.names<-`(rownames(astro_expr))

col_anno <- set_names(colnames(astro_median_expr)) %>%
  as.data.frame() %>%
  setNames("day") %>%
  mutate(Stage = case_when(
    as.numeric(day) < 250 ~ "prenatal",
    as.numeric(day) < 300 ~ "transition",
    T ~ "postnatal"
  )) %>%
  as.data.frame() %>%
  `rownames<-`(.[,1]) %>%
  .[,-1, drop = F]

anno_colors = list(
  "Astrocyte_genes" = setNames(RColorBrewer::brewer.pal(3,"Set1")[1:2], unique(sloan_data$stage)),
  "Stage" = setNames(c("grey70", "grey50", "grey30"),c("prenatal","transition","postnatal")))

pdf(file.path(outputFolder, "astrocyte_transition.pdf"), w = 4, h = 2.5)
pheatmap::pheatmap(astro_median_expr,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         scale = "row",
         color = colorRampPalette( c("green","black", "magenta"))(200),
         annotation_row = row_anno,
         annotation_col = col_anno,
         annotation_colors = anno_colors,
         show_rownames = FALSE,
         fontsize = 8,
         annotation_names_row = FALSE,
         annotation_names_col = FALSE,

)
dev.off()
# metabolic stress genes ------------------------------------------------------------------------------------------

stress_genes <- c("PGK1","ALDOA","BNIP3","YIP5","ARCN1","GORASP2")

stress_hcs <- plot_genes(datExpr_reg_batch, datMeta, stress_genes, scale_expr = F) %+%
  facet_wrap(~hgnc_symbol, ncol = 1) +
  ggtitle(expression(italic("in vitro")))+
  scale_y_continuous(limits =c(3,12), breaks = seq(4,10,2))

stress_bs <- plot_genes(datExpr_reg_cortical, datMetaCortical, stress_genes, scale_expr = F) %+%
  facet_wrap(~hgnc_symbol, ncol = 1) +
  theme(axis.title.y = element_blank()) +
  scale_y_continuous(limits =c(3,12), breaks = seq(4,10,2)) +
  ggtitle(expression(italic("in vivo")))

pdf(file.path(outputFolder,"metabolic_stress_single_genes.pdf"),width = 3.5, height = 4.42)
((stress_hcs + labs(y ="Normalized expression")) +
  (stress_bs) ) & theme(strip.text = element_text(size = 8, face = "italic"))

dev.off()


stree_hcs_data <- stress_hcs$data %>%
  mutate(data_set = "hCS")
stress_bs_data <- stress_bs$data %>%
  mutate(data_set = "bs")
stress_data <- bind_rows(stree_hcs_data,stress_bs_data) %>%
  dplyr::select(hgnc_symbol, value, data_set)

# map_df(set_names(unique(stress_data$hgnc_symbol)), function(gene){
#   current_data <- stress_data %>%
#     filter(hgnc_symbol==gene)
#   lm1 <- lm(value~data_set, data = filter(stress_data,hgnc_symbol==gene))
#   broom::tidy(lm1)
# },.id = "hgnc_symbol") %>%
#   filter(term != "(Intercept)")

#modules from pollen 2019
pollen_wgcna_genes <- read_excel(file.path(base_folder,"metabolic_stress", "Pollen_2019_table_s3.xlsx"), skip = 753)

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
  geom_point(size = 1)  +
  labs(x = "Differentiation Day", y = "Module Eigengene") +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5, size = 8),
        axis.text.y = element_text(size = 8)) +
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
  geom_point(size = 1) +
  labs(x = "BrainSpan Stage", y = "") +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5, size = 8),
        axis.text.y = element_text(size = 8))  +
  scale_y_continuous(limits = limits) +
  scale_x_continuous(breaks = seq(2,15,1)) +
  scale_color_manual(values = clrs,"Module") +
  guides(color=guide_legend(override.aes=list(fill=NA)))


pdf(file.path(outputFolder,"metabolic_stress_ME.pdf"),width = 4.2, height = 1.5)
p_stress_hcs + p_stress_bs + plot_layout(guides = "collect") &
  theme(legend.key.size =  unit(.15, "in"),
        legend.box.margin = margin(0,0,0,0),
        legend.position = "right")

dev.off()


# Context modules trajectories-------------------------------------------------------------------------------------------------
isDark <- function(colr) { (sum( col2rgb(colr) * c(299, 587,114))/1000 < 123) }

Context_MM <- read.csv(file.path(outputFolder, "data", "Stein2015_Gene_module_membership.csv"))
module_descript <- read.table("~/Aaron/Databases/Stein_ModuleMembership/Module_description.tsv", header = F, sep = "\t") %>%
  setNames(c("module", "description"))
if (!file.exists(file.path(outputFolder, "data","context_gene_anno.rdata"))){
  library(biomaRt)
  humanMart <- useMart(biomart = "ensembl",dataset = "hsapiens_gene_ensembl") #change if working with differnet spieces
  getinfo <- c("ensembl_gene_id","hgnc_symbol","entrezgene_id","chromosome_name",
               "strand","start_position", "end_position","gene_biotype")
  geneAnnoRaw <- getBM(attributes = getinfo,filters = c("entrezgene_id"),values = Context_MM$ENTREZ_GENE_ID,mart = humanMart)
  save(geneAnnoRaw,file = file.path(outputFolder, "data","context_gene_anno.rdata"))
} else {
  load(file.path(outputFolder, "data","context_gene_anno.rdata"))
}

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

module_groups <- list(#c("magenta", "red", "greenyellow", "grey60"),
  c("salmon", "cyan", "green", "lightgreen", "purple","darkgreen"), #"darkred", "paleturquoise"
  c("pink", "tan", "yellow", "lightcyan", "black")
  #c("midnightblue", "skyblue", "turquoise"),
  #c("darkturquoise", "darkolivegreen")
)

context_plots <- map(module_groups, function(grp){
  plot_subset <- subset(MEs_plot,Module %in% grp) %>%
    mutate(Module = factor(Module, levels = grp))
  table_subset <- subset(module_descript, module %in% Hmisc::capitalize(grp))

  current_plot <- ggplot(plot_subset,aes(x = Differentiation.day, y = value)) +
    theme(strip.background = element_rect(fill = "white"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "none",
          plot.title = element_text(size = 8)) +
    # geom_vline(xintercept = prdct["mean"] , linetype = 1, color = "grey 70") +
    # geom_vline(xintercept = prdct["lwr_CI"] , linetype = 3, color = "grey 70") +
    # geom_vline(xintercept = prdct["upr_CI"] , linetype = 3, color = "grey 70") +
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
                     # fontface = "bold",
                     nudge_x      = 350,
                     direction = "y",
                     label.size = NA,
                     hjust        = 1,
                     segment.size = 0.2,
                     label.padding	= 0.15,
                     size = 2.8) +
    theme(legend.position = "none") +
    scale_x_continuous(limits = c(0,800), breaks = seq(0,600,100)) +
    scale_fill_manual(values = label_info$colour) +
    scale_y_continuous(breaks = seq(-2,2,1), limits = c(-1.5, 1.5))
})

pdf(file.path(outputFolder, "Trajectories_modules_subset.pdf"), height = 4, width = 3.5, useDingbats = FALSE)
context_plots[[1]] + labs(title ="Neuronal modules") +
  context_plots[[2]] + labs(title = "Glial Modules")  + plot_layout(ncol = 1)
dev.off()

# context module overlap ------------------------------------------------------------------------------------------
# sig_mods <- read.csv(file.path(outputFolder,"wgcna/8__ModuleXGenptype.csv")) %>%
#   filter(fdr < 0.1) %>%
#   mutate(X = gsub("ME", "", X))
attach(file.path(outputFolder,"data","wgcna.Rdata"))
#geneInfo <- geneInfo
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

# subset_modules  <- read.table("~/Aaron/Databases/Stein_ModuleMembership/Module_description.tsv", header = F, sep = "\t") %>%
#   setNames(c("module", "description")) %>%
#   dplyr::select(module) %>%
#   unlist()
wgcna_or <- geneSetEnrichment( context_matrix_gs, wgcna_matrix_gs)
dispMat <- wgcna_or[["dispMat"]]
tMat <- wgcna_or[["tMat"]]
dispMat2 <- dispMat
dispMat2[dispMat2 < 0] <- 0
neuronal <- c("cyan", "darkgreen", "darkred", "green", "lightgreen", "paleturquoise", "purple", "salmon")
glial <- c("black", "lightcyan", "orange", "pink","tan","yellow")
cols <- rep('black', nrow(dispMat2))
cols[row.names(dispMat2) %in% neuronal] <- brewer.pal(3,"Set1")[1]
cols[row.names(dispMat2) %in% glial] <- brewer.pal(3,"Set1")[2]
pdf(file.path(outputFolder, "overlap_context_modules.pdf"), height = 4, width = 3.2)
par(ps = 8)
heatmap.2(dispMat2,
          trace = "none",
          col = colorpanel(n = 100, low = "white", high = "blue3"),
          cellnote = tMat,
          notecex = 0.75,
          notecol = "white",
          cexRow = 1,
          cexCol = 1,
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

load(file.path(outputFolder, "data", "LongTerm_pariedVoom_results.rdata"))

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

pdf(file.path(outputFolder, "n_degs.pdf"), h = 1, w = 4.5)
ggplot(plot_data, aes(x = 1, y = n, fill = direction)) +
  geom_bar(stat = "identity", position = "identity") +
  scale_fill_manual(values = c("chartreuse3", "firebrick3")) +
  geom_hline(yintercept = 0, size = 1) +
  scale_y_continuous(breaks = NULL,name = "",limits = limits) +
  scale_x_discrete(name = "") +
  theme(legend.position = "none",
        axis.line = element_blank(),
        # axis.text.y = element_text(size = 14, face = "bold"),
        # axis.title = element_text(size = 18, face = "bold"),
        # strip.text = element_text(size = 18, face = "bold"),
        #strip.background = element_rect(fill = "grey90", color = "white"),
        panel.grid.major.y = element_blank(),
        strip.text = element_text(size = 8)) +
  coord_flip() +
  geom_text(data = filter(plot_data, direction == "up"),aes(label = paste(n, "\ngenes")),hjust = -0.1, size = 2.8) +
  geom_text(data = filter(plot_data, direction == "down"),aes(label = paste(abs(n), "\ngenes")),hjust = 1.1, size= 2.8) +
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
  # write_csv(res,
  #           path = file.path(gsea_dir,paste0("gsea_res_",f,".csv")))
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
  mutate(pathway = str_wrap(pathway, width = 40)) %>%
  #mutate(stage = factor(stage, levels = str_sort(unique(stage), numeric = T))) %>%
  mutate(pathway = reorder_within(pathway, NES, Comparison))

p_fgsea <- ggplot(data_gsea_plot, aes(x = pathway, y = NES, fill = -log10(padj))) +
  coord_flip() +
  geom_col() +
  theme_classic(base_size = 8) +
  facet_wrap(~Comparison, scales = "free_y", ncol = 1)+#, space = "free_y") +
  geom_hline(yintercept = 0) +
  scale_y_continuous(breaks=seq(-3,3,2) )+
  theme(strip.background = element_rect(),
        # strip.text.x = element_blank(),
        strip.text.x = element_text(size = 8),
        axis.text =  element_text(size = 8), legend.margin = margin(0,0,0,0)) +
  labs(y = "Normalized enrichment score", x = "") +
  scale_x_reordered() +
  scale_fill_gradient(low = "#FCBBA1",high = "#DE2D26", name = "-log10(fdr)", limits = c(0,NA))
pdf(file.path(outputFolder, "gsea_bargraph.pdf"), width = 4, h = 4 )
p_fgsea
dev.off()

# switchovers -----------------------------------------------------------------------------------------------------
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
  p + theme(plot.margin = margin(0.0,0.0,0.0,0.0, "cm"))
})
l1 <- plots[[1]] + plots[[2]]  + plot_spacer() + plot_layout(widths = c(1,1,0.1))
l2 <- plots[[3]] + plots[[4]]+ plot_spacer() + plot_layout(widths = c(1,1,0.1))
l3 <- plots[[5]] + plots[[6]] + plot_spacer() + plot_layout(widths = c(1,1,0.1))
l4 <- plots[[7]] + plots[[8]] + plot_spacer() + plot_layout(widths = c(1,1,0.1))
pdf(file.path(outputFolder,"switchover_all.pdf"), width = 8, height = 3.5)
((l1 / l2) | (l3/l4)) & theme(legend.text = element_text(size = 8, face = "italic"))

dev.off()
# WB --------------------------------------------------------------------------------------------------------------
library(pzfx)
library(multcomp)
library(ggsignif)
f <- file.path(outputFolder,"data",  "NR2AB analysis_final.pzfx")
nr2b_raw <- read_pzfx(f,"NR2B final_raw value")
nr2b <- nr2b_raw %>%
  gather("day", "value") %>%
  mutate(NR_subunit = "GRIN2B")

nr2a_raw <- read_pzfx(f,"NR2A final_raw value")
nr2a <- nr2a_raw %>%
  gather("day", "value") %>%
  mutate(NR_subunit = "GRIN2A")

all_NRs <- rbind(nr2a,nr2b) %>%
  unite("id", "NR_subunit", "day", remove = F) %>%
  na.omit() %>%
  mutate_all(list(~gsub("-","_",.))) %>%
  mutate(value = as.numeric(value))

# av <- aov(formula = as.formula("value~0+id"), data = all_NRs)
# tukey <- TukeyHSD(av)[[1]] %>%
  # as.data.frame() %>%
  # rownames_to_column("contrast")

lm1 <- lm("value~0+id", data = all_NRs)
contr <- t(limma::makeContrasts(GRIN2A_D200 - GRIN2A_D100,
                                GRIN2A_D300_400 - GRIN2A_D100,
                                GRIN2B_D200 - GRIN2B_D100,
                                GRIN2B_D300_400 - GRIN2B_D100,
                                GRIN2A_D100 - GRIN2B_D100,
                                GRIN2A_D200 - GRIN2B_D200,
                                GRIN2A_D300_400 - GRIN2B_D300_400, levels = unique(all_NRs$id)))
glht_res1 <- summary(glht(lm1, linfct = contr), test = adjusted("BH"))
p_vals <- glht_res1$test$pvalues %>%
  p.adjust() %>%
  as.data.frame() %>%
  setNames("p_adj") %>%
  rownames_to_column("contrast") %>%
  mutate(pretty_padj = formatC(p_adj, digits = 2, format = "e"))

# summary(glht_res1)
# summary(lm1)
sig_df <- tribble(
  ~x,    ~xend,    ~y,    ~contrast,
  0.85,   1.15,    3.4,    "GRIN2A_D100 - GRIN2B_D100",
  2.85,   3.15,    3.4,    "GRIN2A_D300_400 - GRIN2B_D300_400",
  0.85,   2.85,    3.75,    "GRIN2A_D300_400 - GRIN2A_D100",
  1.15,   3.15,     4,     "GRIN2B_D300_400 - GRIN2B_D100"
) %>%
  left_join(p_vals, by = "contrast") %>%
  rowwise() %>%
  mutate(pos = mean(c(x,xend))) %>%
  mutate(annotation = case_when(`p_adj` < 0.005 ~ "***",
                                `p_adj` < 0.01 ~ "**",
                                `p_adj` < 0.05 ~ "*",
                                TRUE ~ ""))


p_wb <- ggplot(all_NRs, aes(x = day, y = value, color = NR_subunit)) +
  geom_boxplot(outlier.shape = NA, width = 0.65,
               position = position_dodge(width = 0.65),
               fill = NA, size = 0.4,
               show.legend = FALSE) +
  geom_point(position = position_jitterdodge(dodge.width = 0.65, jitter.width = 0.1),
             size = 1, show.legend = FALSE) +
  stat_summary(fun.y = median, geom = "line", aes(group = NR_subunit),
               position = position_dodge(width = 0.65), size = 0.5) +
  scale_color_brewer(palette = "Set1", name = "") +
  theme_classic() +
  geom_segment(data = sig_df[c(3,4),], aes(x = x, xend = xend, y = y, yend = y),
               color = c("black")#,"black",RColorBrewer::brewer.pal(3,"Set1")[1:2])
  ) +
  geom_text(data = sig_df, aes(x = pos, y = y + 0.01, label = annotation),
            color = c("black")#,"black",RColorBrewer::brewer.pal(3,"Set1")[1:2])
  ) +
  labs(x = "Differentiation day", y = expression("GRIN/"*beta*"-actin ratio (AU)")) +
  scale_x_discrete(labels = c("Day 77", "Days 160-200", "Days 295-468")) +
  theme(text = element_text(size = 8),
        axis.text = element_text(size = 8),
        axis.line = element_line(size = 0.4),
        legend.key.size =  unit(.1, "in"),
        axis.text.x = element_text(angle = - 45, hjust = 0, vjust = 0.5),
        legend.margin = margin(t = 0, unit = 'cm'),
        legend.text = element_text(size = 8, face = "italic"))
pdf(file.path(outputFolder,"wb_graph.pdf"), w = 2.5, h = 2.5)
p_wb
dev.off()

# RNA editing -----------------------------------------------------------------------------------------------------
load(file.path(outputFolder, "data","RNAediting.rdata"))
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
        strip.text.x = element_text(hjust = 0, size = 8),
        axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5))
pdf(file.path(outputFolder,"RNAedit_bs_mod_traj.pdf"), w = 3, h = 1.5)
p_bs_ME
dev.off()


edit_module_pres  <- mp_edit$preservation$Z[[1]][[2]] %>%
  rownames_to_column("original") %>%
  filter(original != "gold") %>%
  left_join(bs_mod_names, by = "original") %>%
  ggplot(aes(x = mod_name, y = Zsummary.pres)) +
  geom_col(fill = "grey50", width = 0.5) +
  labs(y = "Z summary", x = "Module") +
  scale_y_continuous(limits = c(0,10))


pdf(file.path(outputFolder,"RNAedit_bs_mod_pres.pdf"), w = 1.75, h = 1.5)
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

pdf(file.path(outputFolder,"RNAedit_hcs_mod_traj.pdf"), w = 4.5, h = 1.5)
p_hcs_ME
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
  mutate(hCS = factor(hCS, labels = str_sort(unique(hCS), decreasing = T, nuemric = T)))


or_plot <- ggplot(or_plot_data,aes(y = hCS, x = brainSpan, fill = `log2(OR)`, label = text)) +
  geom_tile() +
  geom_text(size = 2) +
  scale_fill_gradient(limits = c(0,NA), na.value = NA, low = "white", high = "blue2") +
  theme(axis.text = element_text(size = 8),
        axis.title.y = element_text(vjust = 2))
pdf(file.path(outputFolder, "editing_or_mods.pdf"), w = 3, h = 2)
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

pdf(file.path(outputFolder, "editing_genes_expr.pdf"), w = 8, h = 2.5)
(eg_hcs + labs(y = expression(atop("",atop(italic("In vitro"),"expression")))))/
  (eg_bs + labs(y = expression(atop("",atop(italic("In vivo"),"expression"))))) *
  theme(axis.title.y = element_text(size = 11),
        plot.title = element_blank()) &
  theme(strip.text = element_text(size = 8, face = "italic"))
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
  pch.cex = 1, pch.col = "black",
  cl.pos="n"
)
dev.off()


pdf(file.path(outputFolder, "editing_ME_enz_cor_legend.pdf"), w = 1.5, h = 3)
plot.new()
corrplot::colorlegend(colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                                         "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                                         "#4393C3", "#2166AC", "#053061"))(200),
                      labels = c(-1,0,1))

dev.off()

# E-phys ----------------------------------------------------------------------------------------------------------
epyhs_file <- file.path(base_folder,"e-phys","NMDA_Spheroid_OldYoung_Updated.pzfx")

response_amp <- read_pzfx(path = epyhs_file, table = "Response Amplitudes")
response_amplitude_baseline <- read_pzfx(path = epyhs_file, table = "Response Amplitude Baseline") %>%
  gather("measurement","Amp", -DIV, na.rm = T )

r_ephys <- cor.test(response_amplitude_baseline$Amp, response_amplitude_baseline$DIV)
r_p_val_ephys <- formatC(r_ephys$p.value, format = "e", digits = 2) %>%
  str_split("e") %>%
  unlist

pdf(file.path(outputFolder, "e-phsy_baseline.pdf"), h = 1.95, w = 2)
ggplot(response_amplitude_baseline, aes(x = DIV, y = Amp)) +
  geom_smooth(method = "lm", color = "black", se = F) +
  geom_point(size = 1) +
  geom_vline(xintercept = 250, color = "black", linetype = 2) +
  annotate("segment", x = 240, xend = 140, y = 300, yend = 300,
           arrow = arrow(length = unit(0.05, "in"), type = "closed")) +
  annotate("segment", x = 260, xend = 360, y = 300, yend = 300,
           arrow = arrow(length = unit(0.05, "in"), type = "closed")) +
  annotate("text", label = "Early", x = -20, y = 300, hjust = 0,size = 2.8,) +
  annotate("text", label = "Late", x = 390, y = 300, hjust = 0,size = 2.8,) +
  annotate(geom = "text", x = -20, y = 255, size = 2.8, hjust = 0,label = paste("r =", round(r_ephys$estimate,2))) +
  annotate(geom = "text", x = -20, y = 220, size = 2.8, hjust = 0,
           label = substitute(
             paste("p = ", part1,"e"^part2),
             list(
               part1 = r_p_val_ephys[1],
               part2 = as.numeric(r_p_val_ephys[2]))
           ),
  ) +
  labs(x = "Differentiation day", y = "NMDA Amp (pA)")
dev.off()

ifn_reduction <- read_pzfx(path = epyhs_file, table = "% reduction by IFN") %>%
  gather("stage_orig", "percent", na.rm = T) %>%
  mutate("Stage" =  ifelse(stage_orig == "Baseline Old", "Late", "Early"))


ifn_reduction %>%
  group_by(Stage) %>%
  summarise(mean_reduction = mean(percent))
wilcox.test(percent ~ Stage, ifn_reduction)

pdf(file.path(outputFolder, "e-phsy_IFN_reduction.pdf"), h = 1.95, w = 1.5)
ggplot(ifn_reduction, aes(x = Stage, y = percent)) +
  geom_boxplot(outlier.color = NA, fill = NA) +
  geom_jitter(height = 0, width = 0.15, size = 1) +
  labs(x = "", y = "% IFN Redcution") +
  scale_y_continuous(limits = c(0,100), breaks = seq(0,100,20)) +
  annotate("segment", x = 1, xend = 2, y = 95, yend = 95) +
  annotate("text", x = 1.5, y = 97, label = "***")
dev.off()
