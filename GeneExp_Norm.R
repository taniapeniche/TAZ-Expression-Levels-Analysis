setwd("[path]/data")
load("[path]/AllData.RData")

#load all the gene expression
geneExp_Norm<-geneExp

#select the samples that are from tumors
for (i in 1:length(geneExp_Norm)){
  #Create a vector with all the ID;
  #Separate the ID in the points, using a function to do that
  #index the position 4
  pref_Norm<-sapply(names(geneExp_Norm[[i]]),function(x)unlist(strsplit(x,"\\."))[4])
  #Choose the ones that are tumoral (01)
  #names+grep(function that search for patterns)+the pattern + the data where to search
  sample_ids<-names(pref_Norm)[grep("1[0:1:2:3:4]",pref_Norm)]
  #In geneExp_tum we will update the list with only the the samples that that are tumoral
  geneExp_Norm[[i]]<-geneExp_Norm[[i]][sample_ids]
}

#Save the file
save(geneExp_Norm, file = "geneExp_Norm.RData")

#Boxplot
setwd("[path]/Results")
load("[path]/geneExp_Norm.RData")

plotData_Norm<- lapply(geneExp_Norm, as.numeric)
boxplot(plotData_Norm,col="#8EBC2B", outline=FALSE,
        ylab="Expression Levels in Normal Tissue Cells (log2 FPKMs)", xlab="Cancer types",
        main="Expression Levels of TAZ Gene in  Normal Cells")
