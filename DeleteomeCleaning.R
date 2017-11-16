#this script originally created by Sheila Gaynor
#slight modifications by Will Townes

#Read in data
deleteome_responsive_mutants_ex_wt_var_controls <- read.delim("data/deleteome_responsive_mutants_ex_wt_var_controls.txt")

#Convert to scaled, numeric
deleteomeNumeric <- deleteome_responsive_mutants_ex_wt_var_controls[2:6124,]
for (i in seq(5,2112,3)){ #only interested in ratio
  deleteomeNumeric[,i] <- as.numeric(as.character(deleteomeNumeric[,i]))
}
deleteomeNumeric[,seq(5,2112,3)] <- scale(deleteomeNumeric[,seq(5,2112,3)])
deleteomeNumeric$geneSymbol <- as.character(deleteomeNumeric$geneSymbol)

#Remove duplicated rows
deleteomeData <- deleteomeNumeric[-which(duplicated(deleteomeNumeric$geneSymbol)), ]
row.names(deleteomeData) <- deleteomeData$geneSymbol

#Take only the gene measure of interest (ratio)
#Now take the row/col indices and apply them to the square matrix
#Adjust the col names to be able to match to rows
deleteomeData <- as.data.frame(deleteomeData[,seq(5,2112,3)])
cutGeneName <- function(x) strsplit(names(deleteomeData)[x],'[.]')[[1]][1]
names(deleteomeData) <- toupper(sapply(1:length(names(deleteomeData)),cutGeneName))

#Take only matching rows/columns
nameOverlap <- intersect(names(deleteomeData),row.names(deleteomeData))
cols <- names(deleteomeData) %in% nameOverlap 
rows <- row.names(deleteomeData) %in% nameOverlap 
deleteome <- deleteomeData[rows,cols]

#Make sure order matches
match(names(deleteome),row.names(deleteome))
deleteData <- deleteome[match(names(deleteome),row.names(deleteome)),]
table(names(deleteData)==row.names(deleteData))

for (cols in 1:661){
  deleteData[,cols] <- as.numeric(as.character(deleteData[,cols]))
}