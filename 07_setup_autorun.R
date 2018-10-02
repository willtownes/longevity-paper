#This script should be run after 01_data_loading.Rmd
#required data sources in the data subfolder:
#archs4worm.rds
#archs4yeast.rds
#goterms.txt
#genage_models.csv
#deleteome_responsive_mutants.rds
#worm_scrna_atlas_merged_umi.txt
library(SingleCellExperiment)
library(caret)
source("./util.R")
fp<-file.path
normalize_gxp_counts<-function(X){
  #X either a counts matrix or counts scaled by gene length
  log2(1+t(1e6*t(X)/colSums(X)))
}
preprocess<-function(X,binary=FALSE,min_nonzero=5){
  X<-as.matrix(X)
  X<-X[,colSums(X!=0)>=min_nonzero]
  sds<-apply(X,2,sd)
  X<-X[,sds>1e-12]
  if(!binary){ 
    X<-.5*scale(X) 
  }
  X
}
species<-c("worm","yeast")
predictors<-c("go","gxp_archs4","gxp_other1","go_gxp_archs4","go_gxp_other1")
#yeast other1 is deleteome measured genes
#worm other1 is merged UMIs of worm cell atlas
d<-expand.grid(outcome="genage_pro_anti", predictors=predictors, species=species)

sp<-"yeast"
sv<-load_genage(sp)
rownames(sv)<-toupper(sv$symbol)
archs4<-readRDS(paste0("./data/archs4",sp,".rds"))
archs4<-assay(archs4,"counts")/rowData(archs4)$length
archs4<-normalize_gxp_counts(archs4)
if(sp=="yeast"){
  other1<-readRDS("./data/deleteome_responsive_mutants.rds")
  other1<-assay(other1,"logratio")
} else if(sp=="worm"){
  other1<-read.table("./data/worm_scrna_atlas_merged_umi.txt")
  other1<-normalize_gxp_counts(X)
}
go<-load_goterms(sp)
gg<-intersect(rownames(sv),rownames(archs4))
gg<-intersect(gg,rownames(other1))
gg<-intersect(gg,rownames(go))
y<-sv[gg,"pro_longevity"]
names(y)<-gg
archs4<-preprocess(archs4[gg,])
other1<-preprocess(other1[gg,])
go<-preprocess(go[gg,],binary=TRUE)
# train/test split
train_idx<-createDataPartition(y,p=.75,list=FALSE)
ytrn<-y[train_idx]
ytst<-y[-train_idx]


#algs<-c("svmradial","knn","glmnet","gbm")
#filter GO terms with <5 occurrences