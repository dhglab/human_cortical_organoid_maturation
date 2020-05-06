# run RRHO on hCS long term culutres comparing to brainSpan data

options(stringsAsFactors=FALSE)

output_dir <-  "output"


#load RRHO function


color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
  scale = (length(lut)-1)/(max-min);
  plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='');
  mtext(title,2,2.3,cex=0.8);
  axis(2, round(ticks,0), las=1,cex.lab=0.8);
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min;
    rect(0,y,10,y+1/scale, col=lut[i], border=NA);
  }
}

#Rank Rank Hypergeometric Overlap based on Plaisier et al., Nucleic Acids Research, 2010
RRHO <- function(list1,list2,stepsize,labels,plots=FALSE,outputdir=NULL,annot=NULL) {
  ## list 1 is a data.frame from experiment 1 with two columns, column 1 is the Gene Identifier, column 2 is the signed ranking value (e.g. signed -log10(p-value) or fold change)
  ## list 2 is a data.frame from experiment 2 with two columns, column 1 is the Gene Identifier, column 2 is the signed ranking value (e.g. signed -log10(p-value) or fold change)
  ## stepsize indicates how many genes to increase by in each algorithm iteration
  cat("\n\nTransition Mapping\n");
  cat(" by Jason Stein, Luis de la Torre Ubietta, Daniel Geschwind\n");
  cat("\n\n");
  if (length(list1[,1])!=length(unique(list1[,1])))
    stop('Non-unique gene identifier found in list1');
  if (length(list2[,1])!=length(unique(list2[,1])))
    stop('Non-unique gene identifier found in list2');

  list1 = list1[order(list1[,2],decreasing=TRUE),];
  list2 = list2[order(list2[,2],decreasing=TRUE),];
  nlist1 = length(list1[,1]);
  nlist2 = length(list2[,1]);

  ## Number of genes on the array
  N = max(nlist1,nlist2);

  hypermat = matrix(data=NA,nrow=length(seq(1,nlist1,stepsize)),ncol=length(seq(1,nlist2,stepsize)));
  countx = county = 0;
  ##Loop over the experiments
  for (i in seq(1,nlist1,stepsize)) {
    countx = countx + 1;
    for (j in seq(1,nlist2,stepsize)) {
      county = county + 1;
      ## Parameters for the hypergeometric test
      k = length(intersect(list1[1:i,1],list2[1:j,1]));
      s = length(list1[1:i,1]);
      M = length(list2[1:j,1]);
      hypermat[countx,county] = -phyper(k-1,M,N-M,s,lower.tail=FALSE,log.p=TRUE) * log10(exp(1));
    }
    county=0;
  }

  overlapmat = matrix(data=NA,nrow=length(seq(1,nlist1,stepsize)),ncol=length(seq(1,nlist2,stepsize)));
  countx = county = 0;
  ##Loop over the experiments
  for (i in seq(1,nlist1,stepsize)) {
    countx = countx + 1;
    for (j in seq(1,nlist2,stepsize)) {
      county = county + 1;
      ## Parameters for the hypergeometric test
      k = length(intersect(list1[1:i,1],list2[1:j,1]));
      s = length(list1[1:i,1]);
      overlapmat[countx,county] = k/s;
    }
    county=0;
  }

  ## Convert hypermat to a vector and Benjamini Yekutieli FDR correct
  hypermatvec = matrix(hypermat,nrow=nrow(hypermat)*ncol(hypermat),ncol=1);
  hypermat.byvec = p.adjust(10^-hypermatvec,method="BY");
  hypermat.by = matrix(-log10(hypermat.byvec),nrow=nrow(hypermat),ncol=ncol(hypermat));

  if (plots) {
    ## Function to plot color bar
    ## Modified from http://www.colbyimaging.com/wiki/statistics/color-bars
    color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
      scale = (length(lut)-1)/(max-min);
      plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='');
      mtext(title,2,2.3,cex=0.8);
      axis(2, round(ticks,0), las=1,cex.lab=0.8);
      for (i in 1:(length(lut)-1)) {
        y = (i-1)/scale + min;
        rect(0,y,10,y+1/scale, col=lut[i], border=NA);
      }
    }

    pdf(paste(outputdir,paste("RRHOMap",labels[1],"__VS__",labels[2],".pdf",sep=""),sep="/"));
    jet.colors = colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"));
    layout(matrix(c(rep(1,5),2), 1, 6, byrow = TRUE));
    image(hypermat,xlab='',ylab='',col=jet.colors(100),axes=FALSE,main="Rank Rank Hypergeometric Overlap Map");
    mtext(labels[2],2,0.5);
    mtext(labels[1],1,0.5);
    color.bar(jet.colors(100),min=min(hypermat,na.rm=TRUE),max=max(hypermat,na.rm=TRUE),nticks=6,title="-log10(Nominal P-value)");
    dev.off();

    pdf(paste(outputdir,paste("OverlapMap",labels[1],"__VS__",labels[2],".pdf",sep=""),sep="/"));
    jet.colors = colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"));
    layout(matrix(c(rep(1,5),2), 1, 6, byrow = TRUE));
    image(overlapmat,xlab='',ylab='',col=jet.colors(100),axes=FALSE,main="Rank Rank Percent Overlap Map");
    mtext(labels[2],2,0.5);
    mtext(labels[1],1,0.5);
    color.bar(jet.colors(100),min=min(overlapmat,na.rm=TRUE),max=max(overlapmat,na.rm=TRUE),nticks=6,title="Overlap Percentage");
    dev.off();

    ## Make a rank scatter plot
    list2ind = match(list1[,1],list2[,1]);
    list1ind = 1:nlist1;
    corval = cor(list1ind,list2ind,method="spearman");
    pdf(paste(outputdir,paste("RankScatter",labels[1],"__VS__",labels[2],".pdf",sep=""),sep="/"));
    plot(list1ind,list2ind,xlab=paste(labels[1],"(Rank)"),ylab=paste(labels[2],"(Rank)"),pch=20,main=paste("Rank-Rank Scatter (rho = ",signif(corval,digits=3),")",sep=""),cex=0.5);
    model = lm(list2ind~list1ind);
    lines(predict(model),col="red",lwd=3);
    dev.off();

    ## Make a Venn Diagram for the most significantly associated points
    ## Upper Right Corner (Downregulated in both)
    maxind.ur = which(max(hypermat[ceiling(nrow(hypermat)/2):nrow(hypermat),ceiling(ncol(hypermat)/2):ncol(hypermat)],na.rm=TRUE)==hypermat,arr.ind=TRUE);
    indlist1.ur = seq(1,nlist1,stepsize)[maxind.ur[1]];
    indlist2.ur = seq(1,nlist2,stepsize)[maxind.ur[2]];
    genelist.ur = intersect(list1[indlist1.ur:nlist1,1],list2[indlist2.ur:nlist2,1]);
    ## Lower Right corner (Upregulated in both)
    maxind.lr = which(max(hypermat[1:(ceiling(nrow(hypermat)/2)-1),1:(ceiling(ncol(hypermat)/2)-1)],na.rm=TRUE)==hypermat,arr.ind=TRUE);
    indlist1.lr = seq(1,nlist1,stepsize)[maxind.lr[1]];
    indlist2.lr = seq(1,nlist2,stepsize)[maxind.lr[2]];
    genelist.lr = intersect(list1[1:indlist1.lr,1],list2[1:indlist2.lr,1]);

    ## Write out the gene lists of overlapping
    write.table(genelist.ur,paste(outputdir,"/RRHO_GO_MostDownregulated",labels[1],"__VS__",labels[2],".csv",sep=""),row.names=F,quote=F,col.names=F);
    write.table(genelist.lr,paste(outputdir,"/RRHO_GO_MostUpregulated",labels[1],"__VS__",labels[2],".csv",sep=""),row.names=F,quote=F,col.names=F);

  }
  return(list(hypermat=hypermat,hypermat.by=hypermat.by));
}

#load hCS data
load(file.path(output_dir,"LongTerm_pariedVoom_results.rdata"))


# load brainSPan data
load(file.path("data","brainSpan_pariedVoom_results.rdata"))

stepsize=200
brainSpanVoom_earliest <- brainSpanVoom[grep("Period2", names(brainSpanVoom))]
hcs_lt_voom_earliest <- hcs_lt_voom[grep("Day.grouped025",names(hcs_lt_voom))]


for (j in 1:length(brainSpanVoom_earliest)){
  ref <- as.data.frame(brainSpanVoom_earliest[[j]])
  ref$geneID <- rownames(ref)

  for (i in c(1:length(hcs_lt_voom_earliest))){
    name_i <- names(hcs_lt_voom_earliest)[i]
    name_j <- names(brainSpanVoom_earliest)[j]
    print(paste(name_i,"......",name_j))

    test <- as.data.frame(hcs_lt_voom_earliest[[i]])
    test$geneID <- rownames(test)

    ref <- ref[intersect(rownames(ref),rownames(test)),c("geneID","logFC")]
    test <- test[intersect(rownames(ref),rownames(test)),c("geneID","logFC")]

    if(j==1 & i==1){
      hypermat.all = array(data=NA,dim=c(length(hcs_lt_voom_earliest),length(brainSpanVoom_earliest),length(seq(1,length(ref[,1]),stepsize)),length(seq(1,length(test[,1]),stepsize))));
    }

    tmp <- RRHO(test,ref,stepsize = stepsize,labels = c(paste0("iPSC_DE_",name_i),paste0("brianSpan_DE_",name_j)),plots = F)
    hypermat.all[i,j,,] = tmp$hypermat;

  }
}
save(hypermat.all,file=file.path(output_dir,"hypermatAll.rdata"))
