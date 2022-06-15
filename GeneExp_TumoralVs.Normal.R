load("[path]/AllData.RData")
setwd("[path]/Results")


#Select the data####
indexinclude<-which(unlist(lapply(geneExp_Norm,function(x)length(x)))>9)
geneExp_Norm_II<-geneExp_Norm[indexinclude]
geneExp_Tum_II<-geneExp_Tum[indexinclude]

plotData_Norm_II<-lapply(geneExp_Norm_II, as.numeric)
plotData_Tum_II<-lapply(geneExp_Tum_II, as.numeric)

#Boxplot of all the cancers that have normal and tumoral samples#####
#Add colors to distinguish both the origin of the samples

plotData_Norm_III<-plotData_Norm_II
plotData_Tum_III<-plotData_Tum_II

#Add an N e T to the names
names(plotData_Norm_III)<-paste(names(plotData_Norm_II),"_N",sep="")
names(plotData_Tum_III)<-paste(names(plotData_Tum_II),"_T",sep="")

#Binding the lists + order by name
NormVs.Tum<-c(plotData_Norm_III,plotData_Tum_III)
NormVs.Tum<-NormVs.Tum[order(names(NormVs.Tum))]

names_NormVsTum<-c(names(plotData_Norm_II)," ")

boxplot(NormVs.Tum, las=2,
        col = c("#8EBC2B","#2FABCE"), outline=FALSE,
        ylab="Expression Levels (log2 FPKMs)", xlab="Cancer types", 
        main="TAZ Expression Levels in Tumor and Normal Cells")

legend("bottomright", legend = c("Tumor Cells", "Normal Cells"), 
       col = c("#8EBC2B","#2FABCE"),
       pch = 15, bty = "n", pt.cex = 2, cex = 0.6,  horiz = F, inset = c(0.0001, 0.001))


#Calcule p-valeu and save it in one vector####

p.values<-c(t.test(plotData_Tum_II$BLCA,plotData_Norm_II$BLCA)[3],
            t.test(plotData_Tum_II$BRCA,plotData_Norm_II$BRCA)[3],
            t.test(plotData_Tum_II$COAD,plotData_Norm_II$COAD)[3],
            t.test(plotData_Tum_II$ESCA,plotData_Norm_II$ESCA)[3],
            t.test(plotData_Tum_II$HNSC,plotData_Norm_II$HNSC)[3],
            t.test(plotData_Tum_II$KICH,plotData_Norm_II$KICH)[3],
            t.test(plotData_Tum_II$KIRC,plotData_Norm_II$KIRC)[3],
            t.test(plotData_Tum_II$KIRP,plotData_Norm_II$KIRP)[3],
            t.test(plotData_Tum_II$LIHC,plotData_Norm_II$LIHC)[3],
            t.test(plotData_Tum_II$LUAD,plotData_Norm_II$LUAD)[3],
            t.test(plotData_Tum_II$LUSC,plotData_Norm_II$LUSC)[3],
            t.test(plotData_Tum_II$PRAD,plotData_Norm_II$PRAD)[3],
            t.test(plotData_Tum_II$READ,plotData_Norm_II$READ)[3],
            t.test(plotData_Tum_II$STAD,plotData_Norm_II$STAD)[3],
            t.test(plotData_Tum_II$THCA,plotData_Norm_II$THCA)[3],
            t.test(plotData_Tum_II$UCEC,plotData_Norm_II$UCEC)[3])


#Multiplicity of p-values problem####
p.values<-unlist(p.values)
adjustedp.value<-p.adjust(p.values)


#Building a data frame
p.values_I<-as.data.frame(adjustedp.value)
p.values_<-as.data.frame(p.values)
cancer_types<-as.data.frame(names(geneExp_Norm_II))
p.values_II<-cbind.data.frame(cancer_types,p.values_,p.values_I)
colnames(p.values_II)<-c("Cancer Types","p-Value", "Adjusted p-Value")


#LogFC change calcule
mean_Norm<-unlist(lapply(plotData_Norm_II,function(x)mean(na.omit(x))))
mean_Tum<-unlist(lapply(plotData_Tum_II,function(x)mean(na.omit(x))))

LogFC<-mean_Tum-mean_Norm
View(LogFC)
barplot(LogFC)

p.Value_LogFC<-cbind.data.frame((signif(adjustedp.value, digits = 2)),signif(LogFC, digits = 2))
colnames(p.Value_LogFC)<-c("Adj. p-Value", "LogFC")

Table_1<-p.Value_LogFC
Table_1$`Adj. p-Value`<-as.character(Table_1$`Adj. p-Value`)


#barplot

count_Norm_all<- vector("list", length=length(plotData_Norm))
names(count_Norm_all) <- names(plotData_Norm)

for ( i in 1:length (count_Norm_all)){
  count_Norm<-length(plotData_Norm[[i]])
  count_Norm_all[[i]]<-count_Norm
}

count_Tum_all<- vector("list", length=length(plotData_Tum))
names(count_Tum_all) <- names(plotData_Tum)

for ( i in 1:length (count_Tum_all)){
  count_Tum<-length(plotData_Tum[[i]])
  count_Tum_all[[i]]<-count_Tum
}

all_samples<-c(count_Norm_all,count_Tum_all)
all_samples<-all_samples[order(names(all_samples))]

names_plot<-c(names(geneExp)," ")

barplot(unlist(all_samples), las=2,
        col = c("#2FABCE","#8EBC2B"))
legend("topright", legend = c("Normal Sample","Tumor Sample") , 
       col = c("#2FABCE","#8EBC2B"),
       bty = "n", pch=20 , pt.cex = 2, cex = 0.8, horiz = FALSE, inset = c(0.5, 0.5))



#Saving all the data
save(geneExp,plotData,
     geneExp_Norm,geneExp_Norm_II,plotData_Norm,plotData_Norm_II,plotData_Norm_III,
     geneExp_Tum,geneExp_Tum_II,plotData_Tum,plotData_Tum_II,plotData_Norm_III,
     indexinclude,p.values,p.values_I,p.values_II,cancer_types,NormVs.Tum,
     mean_Norm,mean_Tum,LogFC, p.Value_LogFC,
     file = "AllData.RData")
