options(stringsAsFactors=FALSE)
library(WGCNA)
library(plyr)
library(cqn)
library(sva)
library(tidyverse)


#Load Expression, Meta and QC data  and format them ----------------------------------------
load(file = "data/raw_data.rdata")

datMetaAll <- datMeta %>% 
  rename(Differentiation.day.original = Differentiation.day) %>% 
  mutate(Differentiation.day = as.numeric(Differentiation.day.original)) %>% 
  mutate(Differentiation.day = plyr::round_any(Differentiation.day, 25)) %>% 
  as.data.frame() %>% 
  `rownames<-`(.$SampleID)
datMetaAll <- droplevels(datMetaAll)

datQCAll <- datQC[match(rownames(datMetaAll),rownames(datQC)),]
datExprAll <- datExpr[,match(rownames(datMetaAll),colnames(datExpr))] 

# Filter genes with less than 10 counts in 20% of  samples per ID ------------
for (btch in c(1,2)) {
  
  datMeta <- datMetaAll[datMetaAll$batch == paste0("batch", btch),]
  datExpr <- datExprAll[,datMeta$SampleID]
  datQC <- datQCAll[datMeta$SampleID,]
  
  
  # Filter genes with less than 10 counts in 20% of  samples per ID ------------
  pres <- apply(datExpr > 10,1,sum)
  idx <-  which(pres > 0.2*dim(datExpr)[2])
  datExpr <- datExpr[idx,]
  # load gene annotation data (retirved from biomart)
  load("data/gene_annotation.Rdata")
  geneAnno <-  geneAnno %>%
    filter(ensembl_gene_id %in% rownames(datExpr)) %>%
    filter(gene_biotype != "rRNA") %>% 
    filter(!duplicated(ensembl_gene_id))
  
  datExpr <- datExpr[match(geneAnno$ensembl_gene_id,rownames(datExpr)),]
  
  stopifnot(all(geneAnno$ensembl_gene_id==rownames(datExpr)))
  
  # Normalize using CQN
  
  datExpr_preNorm <- datExpr
  datExprCQN <- cqn(datExpr,x=geneAnno$percentage_gene_gc_content, lengths= geneAnno$transcript_length,sqn=F ,verbose = T)
  RPKM.cqn <- datExprCQN$y + datExprCQN$offset
  datExpr.htg <- RPKM.cqn
  
  
  # calculate Principal Component of technical data 
  ## get PCs of the sequencing statistics
  large_qc<- which(apply(datQC,2,mean)>10^4)
  datQC[,large_qc] = log10(datQC[,large_qc])
  PC.datQC <- prcomp(na.omit(t(scale((datQC),scale=T))), center=T)
  varexp <- (PC.datQC$sdev)^2 / sum(PC.datQC$sdev^2)
  
  topPC.datQC <- PC.datQC$rotation[,1:5]
  colnames(topPC.datQC) <- c("SeqPC1","SeqPC2" ,"SeqPC3","SeqPC4","SeqPC5")
  
  
  datMeta <- cbind(datMeta,topPC.datQC)
  assign(paste0("datMeta_batch",btch),datMeta)
  assign(paste0("datExprCQN_batch",btch),datExpr.htg)
}

datMeta <- bind_rows(datMeta_batch1,datMeta_batch2) %>%  
  mutate(Day.grouped = case_when(Differentiation.day.original <= 100 ~ plyr::round_any(as.numeric(Differentiation.day.original),25),
                                 Differentiation.day.original <= 400 ~ plyr::round_any(as.numeric(Differentiation.day.original),50),
                                 Differentiation.day.original <= 600 ~ plyr::round_any(as.numeric(Differentiation.day.original),100), 
                                 TRUE ~ 600)) %>% 
  mutate(Day.grouped = str_pad(Day.grouped,3,"left",0))


genes_to_keep <- intersect(rownames(datExprCQN_batch1),rownames(datExprCQN_batch2))

datExpr.cqn <- cbind(datExprCQN_batch1[genes_to_keep,], datExprCQN_batch2[genes_to_keep,])
datExpr <- datExprAll[genes_to_keep,]


# remove outliers -------------------------------------------------------------------------------------------------
sdout <- 3
normadj <- (0.5 + 0.5 * bicor(datExpr.cqn, use='pairwise.complete.obs'))^2
netsummary <- fundamentalNetworkConcepts(normadj); 
K <- netsummary$Connectivity; Z.K <- (K-mean(K))/sqrt(var(K))
C <- netsummary$ClusterCoef; Z.C = (C - mean(C))/sqrt(var(C))
outliers <- (Z.K > mean(Z.K)+sdout*sd(Z.K))|(Z.K < mean(Z.K)-sdout*sd(Z.K))
print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep="")); print(colnames(datExpr)[outliers]); print(table(outliers))

datMeta <- datMeta[!outliers, ]
datExpr.cqn <- datExpr.cqn[,datMeta$SampleID]
datExpr <- datExpr[,datMeta$SampleID]

# regress out covariates and batch correct -----------------
X =  model.matrix(~ Differentiation.day+Sex+ethnicityPC1+ethnicityPC2+SeqPC1+SeqPC2+SeqPC3+SeqPC4+SeqPC5, data = datMeta)
Y = datExpr.cqn
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)
regress_indx <- grep("Intercept|Differentiation.day",colnames(X), invert = T)
to_regress = (as.matrix(X[,regress_indx]) %*% (as.matrix(beta[regress_indx,])))  #Sex + PCs
datExpr_reg = datExpr.cqn - t(to_regress)


datExpr_reg <- datExpr_reg[,datMeta$SampleID]
covars_model <- model.matrix( ~ Differentiation.day ,  data = datMeta)
datExpr_reg_batch <- ComBat(datExpr_reg, batch = datMeta$batch, mod = covars_model)

save(datMeta, datExpr_reg_batch, datExpr, file = file.path("output","processed_data.rdata"))

