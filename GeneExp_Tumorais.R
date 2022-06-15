setwd("[path]/data")
load("[path]/AllData.RData")

#load all the gene expression
geneExp_Tum<-geneExp

#select the samples that are from tumors
for (i in 1:length(geneExp_Tum)){
  #Create a vector with all the ID;
  #Separate the ID in the points, using a function to do that
  #index the position 4
  pref_Tums<-sapply(names(geneExp_Tum[[i]]),function(x)unlist(strsplit(x,"\\."))[4])
  #Choose the ones that are tumoral (01)
  #names+grep(function that search for patterns)+the pattern + the data where to search
  sample_ids<-names(pref_Tums)[grep("0[1:2:3:4:5:6:7:8:9]",pref_Tums)]
  #In geneExp_tum we will update the list with only the the samples that that are tumoral
  geneExp_Tum[[i]]<-geneExp_Tum[[i]][sample_ids]
}

#Save the file
save(geneExp_Tum, file = "geneExp_Tum.RData")

#Boxplot
setwd("[path]/Results")
load("[path]/AllData.RData")

plotData_Tum<- lapply(geneExp_Tum, as.numeric)
boxplot(plotData_Tum,col="#2FABCE", las=2, outline=FALSE,
        ylab="Expression Levels in Tumor Cells (log2 FPKMs)", xlab="Cancer types",
        main="Expression Levels of TAZ Gene in Tumor Cells")

