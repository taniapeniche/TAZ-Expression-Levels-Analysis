#import + twd
load("[path]/AllData.RData")
setwd("[path]/Results")

#open w pdf
pdf("Different Expresion Values.pdf")
par(mar=c(6,5,5,1))
#All Cells
plotData <- lapply(geneExp, as.numeric)
boxplot(plotData,col="#E5A2DD", outline=FALSE, las=2, cex.axis= 0.7,
        ylab="TAZ Expression Levels (log2 FPKMs)", xlab="Cancer types",
        main="TAZ expresion Level in All The Cells")

#Normal Cells
plotData_Norm<- lapply(geneExp_Norm, as.numeric)
boxplot(plotData_Norm,col="#8EBC2B", outline=FALSE, las=2, cex.axis= 0.7,
        ylab="TAZ Expression Levels in Normal Tissue Cells (log2 FPKMs)", xlab="Cancer types",
        main="TAZ expresion Level in Normal Tissue Cells")

#Tumor Cells
plotData_Tum<- lapply(geneExp_Tum, as.numeric)
boxplot(plotData_Tum,col="#2FABCE", outline=FALSE, las=2, cex.axis=0.7,
        ylab="TAZ Expression Levels in Tumor Cells (log2 FPKMs)", xlab="Cancer types",
        main="TAZ expresion Level in Tumor Cells")


#TumorVs.Normal
boxplot(NormVs.Tum, las=2,
        col = c("#8EBC2B","#2FABCE"), outline=FALSE, las=2,
        cex.axis= 0.6,
        labels = c("",names(geneExp_Tum_II)),
        ylab="TAZ Expression Levels (log2 FPKMs)", xlab="Cancer types", 
        main="TAZ Expression Levels in Tumor and Normal Tissue Cells")

legend("topright", legend = c("Normal Cells","Tumor Cells"), 
       col = c("#8EBC2B","#2FABCE"), 
       pch = 15, bty = "n", pt.cex = 2, cex = 1,  horiz = F, inset = c(0.0001, 0.001))

#sample_size
barplot(unlist(all_samples), las=2,
        title = "Number of samples per Cancer and per type of Cell",
        col = c("#8EBC2B", "#2FABCE"))
legend("topright", legend = c("Normal Sample","Tumor Sample") , 
       col = c("#8EBC2B","#2FABCE"),
       pch = 15, bty = "n", pt.cex = 2, cex = 1,  horiz = F, inset = c(0.0001, 0.001))

#plot_LogFC
plot(LogFC)


for (i in 1:length(data_all)){
  surv_data <- Surv(data_all[[i]]$time, data_all[[i]]$status)
  surv_fit <- survfit(surv_data ~ 1, data=data_all[[i]])
  a<-ggsurvplot(surv_fit, conf.int = F,risk.table = TRUE, xlab = "Time (days)", censor = T,data=data_all[[i]])
  print(a)
}


for (i in 1:length(data_all)){
  surv_data <- Surv(data_all[[i]]$time, data_all[[i]]$status)
  surv_fit_groups <- survfit(surv_data ~ group, data=data_all[[i]])  
  b<-ggsurvplot(surv_fit_groups, conf.int = F,risk.table = TRUE, pval = TRUE, p.adjust = TRUE, p.adjust.methods="FDR", xlab = "Time (days)", title = names_SURV[i] , censor = T,data=data_all[[i]] )
  print(b)
} 



dev.off()



#tables
table1<-tableGrob(p.Value_LogFC)
table2<-tableGrob(SV_final)
table3<-tableGrob(AWS_final)

pdf("Tables_result.pdf", height = 15, width = 20)
grid.arrange(table1)
layout(matrix(c(1,2), 1, 2, byrow = TRUE))
grid.arrange(table2)
grid.arrange(table3)

dev.off()

#survivalreport

# List of ggsurvplots
splots <- list()
splots <- c(survplots["LAML"], survplots["LGG"], survplots["SKCM"], 
            survplots["THCA"], survplots["UCS"], survplots["UCEC"])

# Arrange multiple ggsurvplots and print the output
pdf("Surv_repor.pdf",height = 20 , width = 16)
par(mar=c(10,9,9,5))
arrange_ggsurvplots(splots, print = TRUE, 
                    ncol = 2, nrow = 3, risk.table.height = 0.2)
dev.off()
