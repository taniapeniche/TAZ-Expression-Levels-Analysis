load("[path]/PanCancer_GDC_clinicalInfo.RData")
load("[path]/AllData.RData")
load("[path]/AllData_II.RData")

#Install packages
install.packages("survival")
install.packages("survminer")

# Load R packages
library(survival);library(survminer)


#Get patients IDs + expression level
patients_info <- vector("list", length=length(geneExp_Tum_levels))
names(patients_info) <- names(geneExp_Tum)

for (i in 1:length(patients_info)){
  #Create a vector with all the patients ID in gene tum;
  #Select only the 12 first characters (substring)
  #Change "\\." for "-"
  patientsIDs<-substring(names(geneExp_Tum_levels[[i]]), 1,12)
  patientsIDs<-gsub("\\.", "-", patientsIDs)
  names(geneExp_Tum_levels[[i]])<-patientsIDs
  patientsIDs<-intersect(patientsIDs,row.names(clinicalInfo[[i]]))
  geneExp_Tum_levels[[i]]<-geneExp_Tum_levels[[i]][patientsIDs]
  patients_info[[i]]<-patientsIDs
}


####Get survival information####
time_all <- vector("list", length=length(patients_info))
names(time_all) <- names(patients_info)
vitalStatus_all <- vector("list", length=length(patients_info))
names(vitalStatus_all) <- names(patients_info)

#time
for (i in 1:length(clinicalInfo)){
  time<-as.data.frame(clinicalInfo[[i]][patients_info[[i]],c("days_to_death","days_to_last_followup")])
  time[time== "[Not Applicable]" ] <- NA
  time[time== "-Inf" ] <- NA
  time[,1]<-ifelse(is.na(time[,1]), time[,2], time[,1])
  time[,2]<-NULL
  time$days_to_death<-as.numeric(time$days_to_death)/365
  time[is.na(time)] <- 0
  time_all[[i]]<-time
}

#vitalStatus
for (i in 1:length(vitalStatus_all)){
  vitalStatus <- as.data.frame(clinicalInfo[[i]][patients_info[[i]],"vital_status"])
  vitalStatus[vitalStatus== "Alive"] <- 1
  vitalStatus[vitalStatus== "Dead" ] <- 0
  vitalStatus_all[[i]]<-vitalStatus
}

# Classify patients according to expression levels of gene
#Expression levels where on«btain in previous analysis geneExp_Tum
geneExpStatus_all<-vector("list", length=length(geneExp_Tum_levels))
names(geneExpStatus_all)<-names(geneExp_Tum_levels)

for (i in 1:length(geneExp_Tum_levels)){
  geneExpStatus <- t(as.data.frame(ifelse(geneExp_Tum_levels[[i]] > median(as.numeric(geneExp_Tum_levels[[i]])), "High", "Low")))
  row.names(geneExpStatus)<-patients_info[[i]]
  geneExpStatus_all[[i]]<-geneExpStatus
}

#Create survival object
surv_data_all <- vector("list", length=length(patients_info))
names(surv_data_all) <- names(patients_info)
data_all<- vector("list", length=length(patients_info))
names(data_all) <- names(patients_info)

for (i in 1:length(surv_data_all)){
  data <- cbind(time_all[[i]], vitalStatus_all[[i]], geneExpStatus_all[[i]])
  colnames(data)<-c("time", "status", "group")
  data$time<-as.numeric(data$time)
  data$status<-as.numeric(data$status)
  data_all[[i]]<-data
  surv_data <- Surv(data$time,data$status)
  surv_data_all[[i]]<-surv_data
}

#Estimate survival curve considering only 1 group with all samples
surv_fit_all<- vector("list", length=length(surv_data_all))
names(surv_fit_all) <- names(surv_data_all)

for (i in 1:length(surv_fit_all)){
  surv_fit <- survfit(surv_data_all[[i]] ~ 1, data=data_all[[i]])
  surv_fit_all[[i]]<-surv_fit
}


# Plot survival curve
#for (i in 1:length(surv_fit_all)){
  #ggsurvplot(surv_fit_all[[i]], conf.int = F,risk.table = TRUE, xlab = "Time (days)", censor = T,data=data_all[[i]])
#}

pdf("Surv_all.pdf") 
for (i in 1:length(data_all)){
  surv_data <- Surv(data_all[[i]]$time, data_all[[i]]$status)
  surv_fit <- survfit(surv_data ~ 1, data=data_all[[i]])
  a<-ggsurvplot(surv_fit, conf.int = F,risk.table = TRUE, xlab = "Time (years)", censor = T,data=data_all[[i]])
  print(a)
}
dev.off()

#Estimate survival curves for two groups: high and low expression
surv_fit_groups_all<- vector("list", length=length(surv_fit_all))
names(surv_fit_groups_all) <- names(surv_fit_all)

for (i in 1:length(surv_fit_groups_all)){
  surv_fit_groups <- survfit(surv_data_all[[i]] ~ group, data=data_all[[i]])
  surv_fit_groups_all[[i]]<-surv_fit_groups
}

# Plot survival curves for two groups: high and low expression
#for (i in 1:length(surv_fit_groups_all)){
  #ggsurvplot(surv_fit_groups_all[[i]], conf.int = F,risk.table = TRUE, xlab = "Time (years)", censor = T,data=data_all[[i]], pval=T)
#}
survplots<-vector("list", length=length(surv_fit_groups_all))
names(survplots) <- names(surv_fit_groups_all)

pdf("Surv_groups.pdf")
for (i in 1:length(data_all)){
  surv_data <- Surv(data_all[[i]]$time, data_all[[i]]$status)
  surv_fit_groups <- survfit(surv_data ~ group, data=data_all[[i]])  
  survplots[[i]]<-ggsurvplot(surv_fit_groups, conf.int = F,risk.table = TRUE,pval = TRUE, p.adjust = TRUE, p.adjust.methods="FDR", xlab = "Time (years)", title = names_SURV[i] ,censor = T,data=data_all[[i]])
  print(survplots[[i]])
}   
dev.off()


#Compare two survival curves (get statistics)
surv_groups_all<- vector("list", length=length(surv_data_all))
names(surv_groups_all) <- names(surv_data_all)
pValue_all<-vector("list", length=length(surv_data_all))
names(pValue_all) <- names(surv_data_all)
pValue_sv<-("list")

for (i in 1:length(surv_groups_all)){
  surv_groups <- survdiff(surv_data_all[[i]] ~ group, data=data_all[[i]])
  surv_groups_all[[i]]<-surv_groups
  pValue_sv[i] <- pchisq(surv_groups$chisq, length(surv_groups$n)-1, lower.tail = FALSE)
  pValue_sv<-as.numeric(pValue_sv)
}

pValue_sv_adjust<-p.adjust(pValue_sv, method="fdr")
names_SURV<-names(data_all)

#Createa a table
SV_final<-data.frame(scientific(pValue_sv, digits = 2),scientific(pValue_sv_adjust, digits = 2))
row.names(SV_final)<-names(geneExp_Tum_levels)
colnames(SV_final)<-c("p-Value", "Adj. p-Value")

Table_2<-SV_final

save(patients_info,
     time_all,vitalStatus_all,geneExpStatus_all,
     data_all,surv_data_all,surv_fit_all,
     surv_fit_groups_all,surv_groups_all,
     pValue_sv,pValue_sv_adjust, SV_final,
     file = "AllData_II.RData")


