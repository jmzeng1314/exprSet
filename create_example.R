suppressPackageStartupMessages(library(airway))
suppressPackageStartupMessages(library(DESeq))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(edgeR)) 
suppressPackageStartupMessages(library(pasilla))
suppressPackageStartupMessages(library(pasilla))
suppressPackageStartupMessages(library(Biobase))
suppressPackageStartupMessages(library(CLL))
data(sCLLex)
data(pasillaGenes) 
data(airway)
if(T){
  exprSet=exprs(sCLLex)
  write.table(exprs(sCLLex),file='sCLLex.expression.txt',quote = F,sep = '\t',row.names = T)
  groupInfo=pData(sCLLex)[,c(1,2)]
  colnames(groupInfo)=c('sampleID','group')
  groupInfo$sampleID=rownames(groupInfo)
  write.table(groupInfo,file='sCLLex.groupInfo.txt',quote = F,sep = '\t',row.names = F)
  library(hgu95av2.db)
  probe2id=toTable(hgu95av2ENTREZID)
  probe2symbol=toTable(hgu95av2SYMBOL)
  probe2genename=toTable(hgu95av2GENENAME)
  geneInfo=merge(probe2id,probe2symbol,by='probe_id')
  write.table(geneInfo,file='sCLLex.geneInfo.txt',quote = F,sep = '\t',row.names = F)
  save(exprSet,geneInfo,groupInfo,file = 'sCLLex.Rdata')
}

if(T){
  exprSet=assays(airway)$counts
  write.table(assays(airway)$counts,file='airway.expression.txt',quote = F,sep = '\t',row.names = T)
  groupInfo=as.data.frame(colData(airway)[,c(1,3)])
  colnames(groupInfo)=c('sampleID','group')
  groupInfo$sampleID=rownames(groupInfo)
  
  write.table(groupInfo,file='airway.groupInfo.txt',quote = F,sep = '\t',row.names = F)
  
  geneAnno=read.table('hg.geneAnno',stringsAsFactors = F)
  colnames(geneAnno)=c('ensembl','type','symbol')
  geneAnno$ensembl=unlist(lapply(geneAnno$ensembl,function(x)strsplit(x,'\\.')[[1]][1]))
  geneInfo=unique(geneAnno) 
  
  write.table(geneInfo,file='airway.geneInfo.txt',quote = F,sep = '\t',row.names = F)
  save(exprSet,geneInfo,groupInfo,file = 'airway.Rdata')
}






