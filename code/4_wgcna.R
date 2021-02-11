# Run WGCNA on LT data
options(stringsAsFactors = FALSE)

packages <- c("impute","dynamicTreeCut","flashClust","Hmisc",
              "WGCNA","gridExtra","grid","gtable","tidyverse", "sva")

sapply(packages, require, character.only = TRUE)



output_dir <-  "output"

dir.create(output_dir, recursive = T,showWarnings = F)

load(file.path(data, "gene_annotation.Rdata"))
load(file.path(output_dir, "processed_data.rdata"))

datExpr <- t(datExpr_reg_batch)

#remove genes with to many missing values
gsg <- goodSamplesGenes(datExpr, verbose = 3)
print(gsg$allOK)
if (!gsg$allOK) {
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes) > 0)
    printFlush(paste("Removing genes:", paste(names(datExpr)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples) > 0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 <- datExpr
  datExpr = datExpr0[gsg$goodSamples, gsg$goodGenes]
}


powers = c(c(1:10), seq(from = 12, to = 20, by = 2))

# pick soft threshold. First value with R^2 > 0.8
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab = "Soft Threshold (power)",ylab = "Scale Free Topology Model Fit,signed R^2",type = "n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels = powers, cex = cex1, col = "red")
# this line corresponds to using an R^2 cut-off of h
abline(h = 0.8, col = "red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab = "Soft Threshold (power)",ylab = "Mean Connectivity", type = "n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, cex = cex1, col = "red")

softPower <-  12

#calculate adjacency matrix
adjacency <-  adjacency(datExpr, power = softPower, type = "signed",
                      corFnc = "bicor",corOptions = "use='pairwise.complete.obs'");


# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency)
dissTOM = 1 - TOM
save(dissTOM, file = paste0(output_dir,"dissTOM.Rdata" ))

####### create modules
ds = 4
mms = 100
dthresh = 0.1

# Call the hierarchical clustering function
geneTree = flashClust(as.dist(dissTOM), method = "average");

tree = cutreeHybrid(dendro = geneTree, pamStage = FALSE,
                    minClusterSize = mms, cutHeight = 0.9999,
                    deepSplit = ds, distM = as.matrix(dissTOM))
merge <- mergeCloseModules(exprData = datExpr,colors = tree$labels, cutHeight = dthresh)
mColorh <- cbind(labels2colors(merge$colors))
mLabelh <- c("Merged Colors")


pdf(paste0(output_dir,5,"_",condition,"_Final_ModuleDendro.pdf"),height = 10, width = 16)
plotDendroAndColors(geneTree, mColorh, groupLabels = mLabelh,addGuide = TRUE, dendroLabels = FALSE, main = paste("Signed bicor network with power = ",softPower,"mms=",mms,"ds=",ds,"dthresh=",dthresh))
dev.off()

mergedColors = labels2colors(merge$colors);

# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;
MEList = moduleEigengenes(datExpr, colors = mergedColors, softPower = softPower, nPC = 1)
MEs = MEList$eigengenes
rownames(MEs) <- rownames(datExpr)
MEs = orderMEs(MEs)
moduleColors = mergedColors
table(moduleColors)
KMEs <- signedKME(datExpr, MEs, outputColumnName = "kME", corFnc = "bicor")
geneAnnoSubset <- geneAnno[match(rownames(t(datExpr)),geneAnno$ensembl_gene_id),]

geneInfo <- geneAnnoSubset %>%
  dplyr::select(ensembl_gene_id,contains("symbol")) %>%
  bind_cols(as.data.frame(x = moduleColors)) %>%
  bind_cols(KMEs)

colnames(geneInfo)[1] = "Ensembl.Gene.ID"
colnames(geneInfo)[2] = "GeneSymbol"
colnames(geneInfo)[3] = "Initially.Assigned.Module.Color"
write.csv(geneInfo, file = paste0(output_dir, "5.1_",condition,"geneModules.csv"))
