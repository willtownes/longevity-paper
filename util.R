#utility functions
library(dplyr)
library(Matrix)
library(caret)
library(ROCR)

org_match<-list(yeast="Saccharomyces cerevisiae",
                worm="Caenorhabditis elegans",
                fly="Drosophila melanogaster",
                mouse="Mus musculus")

load_genage<-function(species=c("yeast","worm","fly","mouse"),fpath="./data/genage_models.csv"){
  d<-read.csv(fpath)
  species<-match.arg(species)
  org<-org_match[[species]]
  res<-d %>% subset(organism==org & longevity.influence %in% c("Pro-Longevity","Anti-Longevity")) %>% transform(symbol=toupper(symbol)) %>% mutate(pro_longevity=as.numeric(longevity.influence=="Pro-Longevity")) %>% select(symbol,pro_longevity)
  res<-res %>% group_by(symbol) %>% summarise(pro_longevity=round(mean(pro_longevity))) #removes duplicates by averaging them
  as.data.frame(res)
}

load_goterms<-function(species=c("yeast","worm","fly","mouse"),fpath="./data/goterms.txt",sparse=TRUE){
  sp<-match.arg(species)
  go<-read.table(fpath,header=TRUE)#,stringsAsFactors=FALSE)
  go<-go[go$species==sp,c("external_gene_name","go_id")]
  go<-subset(go,external_gene_name != "" & go_id != "")
  #transform to a binary design matrix (go terms in columns, genes in rows)
  #reshape2::acast(go,external_gene_name~go_id,fill=0,fun.aggregate=length)
  i<-factor(go$external_gene_name)
  j<-factor(go$go_id)
  X<-sparseMatrix(as.integer(i),as.integer(j))
  rownames(X)<-toupper(levels(i)); colnames(X)<-levels(j)
  if(sparse){ return(X) } else { return(as.matrix(X)) }
}

#for any caret model, extract the true and false positive rates that can be used to make ROC
caret_roc<-function(fit,Xtst,ytst,ytrn){
  p1<-prediction(predict(fit,type="prob")[,2], ytrn)
  p2<-prediction(predict(fit,Xtst,type="prob")[,2], ytst)
  p1<-performance(p1,"tpr","fpr")
  p2<-performance(p2,"tpr","fpr")
  d1<-data.frame(tpr=p1@y.values[[1]],fpr=p1@x.values[[1]],fold="train")
  d2<-data.frame(tpr=p2@y.values[[1]],fpr=p2@x.values[[1]],fold="test")
  rbind(d1,d2)
}

roc2auc<-function(fpr,tpr){ 
  #trapezoid rule for area under the curve
  m<-(tpr[-1]+tpr[-length(tpr)])/2
  sum(m*diff(fpr)) 
}
