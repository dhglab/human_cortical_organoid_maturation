options(stringsAsFactors = FALSE)
library(edgeR)
library(limma)
library(tidyverse)

output_dir <- "output"

load(file.path(output_dir,"processed_data.rdata"))
datExpr.lt <- DGEList(datExpr)
datExpr.lt <- calcNormFactors(datExpr.lt, method = "TMM")

design <- "~ 0 +  Day.grouped  + batch + ethnicityPC1 + ethnicityPC2 + SeqPC1 + SeqPC2 + SeqPC3 + SeqPC4 + SeqPC5"

mod <- model.matrix(as.formula(design),data = datMeta)
colnames(mod) <- make.names(colnames(mod))
days <- grep("Day.*[0-9]{3}",colnames(mod),value = T)
days_combi <- combn(days,2)
make_combinations <- apply(days_combi,2, function(x) paste0(x[2],"-",x[1]))

contr.matrix <- makeContrasts(
  contrasts = make_combinations,
  levels = colnames(mod))

first_voom <- voom(datExpr.lt, mod)
corfit <- duplicateCorrelation(first_voom, mod, block = datMeta$CellLine)
voomExpr <- voom(datExpr.lt, mod, block = datMeta$CellLine, correlation = corfit$consensus)
lmSA <- lmFit(voomExpr,mod,block = datMeta$CellLine, correlation = corfit$consensus)

fit2 <- contrasts.fit(lmSA, contr.matrix)
efit <- eBayes(fit2)

#save(efit, file = paste0(output_dir,"voom.Rdata"))



hcs_lt_voom <- list()
for (cont in colnames(efit$contrasts)) {
  print(cont)
  hcs_lt_voom[[cont]] <- topTable(efit, coef = cont, number = Inf,resort.by = "logFC")
}
save(hcs_lt_voom,file = file.path(output_dir,"LongTerm_pariedVoom_results.rdata"))
load(file.path(output_dir,"LongTerm_pariedVoom_results.rdata"))
#send data to hoffman2
#system(paste("scp", file.path(getwd(), outputFolder, "LongTerm_pariedVoom_results.rdata"), "aarong@hoffman2.idre.ucla.edu:/u/project/geschwind/aarong/CIRM_iPSC/analysis/rsem/1_longTerm/RRHO"))


#load(paste0(outputFolder,"LongTerm_pariedVoom_results.rdata"))



comparisons <- c("Day.grouped200-Day.grouped025", 
                 "Day.grouped400-Day.grouped200")

# fgsea -----------------------------------------------------------------------------------------------------------

library(biomaRt)
library(fgsea)
library(tidytext)
human_mart = useMart("ensembl",dataset = "hsapiens_gene_ensembl")
getinfo <- c( "ensembl_gene_id", "entrezgene_id")
geneDat <- getBM(attributes = getinfo,filters = "ensembl_gene_id",
                 values = rownames(datExpr_reg_batch),mart = human_mart)
go_pathways <- gmtPathways("/Users/aarongordon/Aaron/Databases/GSEA/v7/c5.all.v7.0.entrez.gmt")


fgsea_res <- map(hcs_lt_voom[comparisons], function(de){
  de_ranked <- de %>%
    rownames_to_column("ensembl_gene_id") %>%
    arrange(logFC) %>%
    
    left_join(geneDat, by = "ensembl_gene_id") %>%
    dplyr::select(entrezgene_id,logFC) %>%
    drop_na()
  ranked_stat <- de_ranked$logFC
  names(ranked_stat) <- de_ranked$entrezgene_id
  ranked_stat <- ranked_stat[!duplicated( names(ranked_stat))]
  
  fgseaRes <- fgsea(go_pathways, stats = ranked_stat, nperm = 1000000,minSize = 30, maxSize = 500)
})
save(fgsea_res, file = file.path(output_dir, "fgsea_results.rdata"))
