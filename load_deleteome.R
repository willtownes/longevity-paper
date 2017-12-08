#Read in data
d <- read.delim("data/deleteome_responsive_mutants_ex_wt_var_controls.txt")
gdata<-d[2:nrow(d),1:3]
d_a<-d[2:nrow(d),seq(from=4,to=ncol(d),by=3)]
d_m<-d[2:nrow(d),seq(from=5,to=ncol(d),by=3)]
d_p<-d[2:nrow(d),seq(from=6,to=ncol(d),by=3)]
#exclude the three control experiments at the end
d_a<-d_a[,1:(ncol(d_a)-3)]
d_m<-d_m[,1:(ncol(d_m)-3)]
d_p<-d_p[,1:(ncol(d_p)-3)]
cnames<-colnames(d_a)
#extract symbol of gene that was perturbed
perts<-sapply(strsplit(cnames,"\\."),function(t){toupper(t[1])})
table(perts %in% gdata$geneSymbol) 
#39 perturbation genes not measured
#661 perturbation genes were measured
