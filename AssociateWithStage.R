load("/Users/taniapeniche/Desktop/Project/PanCancer_GDC_clinicalInfo.RData")
load("/Users/taniapeniche/Desktop/Project/data/AllData.RData")
load("/Users/taniapeniche/Desktop/Project/data/AllData_II.RData")

# Association of expression levels and tumor stage####

##cancer with the stage
have_stage<-which(sapply(clinicalInfo, function(x) "ajcc_pathologic_tumor_stage" %in% (colnames(x))))
clinicalInfo_Stage<-clinicalInfo[have_stage]
patients_info_stage<-patients_info[have_stage]

# Get tumor stage and simplify to early and advanced stages
tumorStage_all <- vector("list", length=length(patients_info_stage))
names(tumorStage_all) <- names(patients_info_stage)


for(i in 1:length(tumorStage_all)){
  tumorStage <- as.vector(clinicalInfo_Stage[[i]][patients_info_stage[[i]],"ajcc_pathologic_tumor_stage"])
  table(tumorStage)
  tumorStage <- ifelse(tumorStage %in% c("Stage IV","Stage III"), "advancedStage", "earlyStage")
  table(tumorStage)
  tumorStage_all[[i]]<-tumorStage
}
tumorStage_all<-within(tumorStage_all, rm("DLBC"))

# Create contigency table
contTable_all<-vector("list", length=length(tumorStage_all))
names(contTable_all) <- names(tumorStage_all)
geneExpStatus_all_stage<-geneExpStatus_all[names(tumorStage_all)]

for(i in 1:length(contTable_all)){
  contTable <- table(geneExpStatus_all_stage[[i]], tumorStage_all[[i]])
  contTable_all[[i]]<-contTable
}

# Apply fisher's test
fisher_result_all<-vector("list", length=length(contTable_all))
names(fisher_result_all) <- names(contTable_all)
fisher_result_p.Value<-vector("list")
fisher_result_estimate<-vector("list")

for(i in 1:length(fisher_result_all)){
  fisher_result <- fisher.test(contTable_all[[i]])
  fisher_result_p.Value[i]<-fisher_result$p.value
  fisher_result_p.Value<-as.numeric(fisher_result_p.Value)
  fisher_result_estimate[i]<-fisher_result$estimate
  fisher_result_estimate<-as.numeric(fisher_result_estimate)
  fisher_result_all[[i]]<-fisher_result
}

fisher_result_adjustp.Value<-p.adjust(fisher_result_p.Value, method="fdr")

AWS_final<-data.frame(round(fisher_result_p.Value, digits = 2), signif(fisher_result_estimate, digits = 2), signif(fisher_result_adjustp.Value, digits = 2))
row.names(AWS_final)<-names(have_stage_II)
colnames(AWS_final)<-c("p-Value","Odds Ratio", "Adj. p-Value")

Table_3<-AWS_final


save(have_stage,clinicalInfo_Stage,patients_info_stage,
     tumorStage_all,geneExpStatus_all_stage,
     contTable_all,fisher_result_all, fisher_result_p.Value,
     fisher_result_adjustp.Value, fisher_result_estimate,
     AWS_final,
     file = "AllData_III.RData")



clinicalInfo_THCA<-as.data.frame(clinicalInfo[["THCA"]])
