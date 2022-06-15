# Get all the file names and cancer types
setwd("/Users/taniapeniche/Desktop/Project/data")
fileNames <- dir(pattern ="fpkm")
cancerTypes <- gsub("TCGA-", "", gsub("_htseq_fpkm-uq.tab", "", fileNames))
# Create a R list to save all information
geneExp <- vector("list", length=length(cancerTypes))
names(geneExp) <- cancerTypes

# Loop to read each file and save information
for(i in 1:length(cancerTypes)){
  # read file
  curData <- read.table(fileNames[i], sep="\t", header=T)
  #Get the expression levels of your gene
  curData$X1=as.factor(curData$X1)
  curExp<-subset(curData,X1=="ENSG00000102125.14")
  # Save the information in R list
  geneExp[[i]] <- log2(curExp[-1])
}
#Save your data for downstream analyses
save(geneExp, file = "geneExp.RData")


#boxplot
setwd("/Users/taniapeniche/Desktop/Project/Results")
#load("/Users/taniapeniche/Desktop/Project/data/AllData.RData")

plotData <- lapply(geneExp, as.numeric)
boxplot(plotData,col="#CA4DB5", outline=FALSE, las=2,
        ylab="Expression Levels (log2 FPKMs)", xlab="Cancer types",
        main="Expression Levels of TAZ Gene")
