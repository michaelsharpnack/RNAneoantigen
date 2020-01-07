######################## analysis ########################################


#WBs and SBs
netmhc_out_all_WB <- list()
netmhc_out_all_SB <- list()
for(u in 1:length(netmhc_out_all)){
  netmhc_out_all_WB[[u]] <- netmhc_out_all[[u]][rowMins(as.matrix(netmhc_out_all[[u]][,1:(dim(netmhc_out_all[[u]])[2]-1)]),na.rm = TRUE) < 2 & rowMins(as.matrix(netmhc_out_all[[u]][,1:(dim(netmhc_out_all[[u]])[2]-1)]),na.rm=TRUE) >= 0.5,]
  netmhc_out_all_SB[[u]] <- netmhc_out_all[[u]][rowMins(as.matrix(netmhc_out_all[[u]][,1:(dim(netmhc_out_all[[u]])[2]-1)]),na.rm=TRUE) < 0.5,]
}

peps.per.patient <- vector(mode='numeric',length(netmhc_out_all)) #number of peptides per patient
names(peps.per.patient) <- patient.sample
peps.per.patient.wb <- peps.per.patient #number of weak binding peptides per patient
peps.per.patient.sb <- peps.per.patient #number of strong binding peptides per patient
for(count in 1:length(netmhc_out_all)){
  peps.per.patient[count] <- dim(netmhc_out_all[[count]])[1] 
  peps.per.patient.wb[count] <- dim(netmhc_out_all_WB[[count]])[1] 
  peps.per.patient.sb[count] <- dim(netmhc_out_all_SB[[count]])[1] 
}

#peptide frequency  
peptides.all <- rownames(netmhc_out_all[[1]]) 
peptides.all.WB <- rownames(netmhc_out_all_WB[[1]]) 
peptides.all.SB <- rownames(netmhc_out_all_SB[[1]]) 
for(u in 1:length(netmhc_out_all)-1){
  peptides.all <- c(peptides.all,rownames(netmhc_out_all[[u+1]]))
  peptides.all.WB <- c(peptides.all.WB,rownames(netmhc_out_all_WB[[u+1]]))
  peptides.all.SB <- c(peptides.all.SB,rownames(netmhc_out_all_SB[[u+1]]))
}  
peptides.all <- table(peptides.all)
peptides.all.WB <- table(peptides.all.WB)
peptides.all.SB <- table(peptides.all.SB)

bind.matrix <- matrix(0,nrow=length(union(names(peptides.all.SB),names(peptides.all.WB))),ncol=length(patient.sample))
rownames(bind.matrix) <- sort(union(names(peptides.all.SB),names(peptides.all.WB)))
colnames(bind.matrix) <- patient.sample
for(count in 1:dim(bind.matrix)[2]){
  bind.matrix[rownames(netmhc_out_all_WB[[count]]),count] <- 1
  bind.matrix[rownames(netmhc_out_all_SB[[count]]),count] <- 2
}


hla_fraction_bind <- vector(mode='numeric',0)
npeps_hla <- vector(mode='numeric',0)
npeps <- vector(mode='numeric',length(netmhc_out_all))
npeps_frac <- vector(mode='numeric',length(netmhc_out_all))
npeps_bind <- vector(mode='numeric',length(netmhc_out_all))
for(u in 1:length(netmhc_out_all)){
  hla_fraction_bind <- c(hla_fraction_bind,colSums(as.matrix(netmhc_out_all[[u]][,1:(dim(netmhc_out_all[[u]])[2]-1)]) < 2,na.rm=TRUE)/dim(netmhc_out_all[[u]])[1])
  npeps_hla <- c(npeps_hla,rep(dim(netmhc_out_all[[u]])[1],dim(netmhc_out_all[[u]])[2]-1))
  npeps[u] <- dim(netmhc_out_all[[u]])[1]
  npeps_bind[u] <- dim(netmhc_out_all_SB[[u]])[1]+dim(netmhc_out_all_WB[[u]])[1]
  npeps_frac[u] <- (dim(netmhc_out_all_SB[[u]])[1]+dim(netmhc_out_all_WB[[u]])[1])/dim(netmhc_out_all[[u]])[1]
}

#mean number of binders per hla allele
mean.hla.bind <- vector(mode='numeric',length(unique(names(hla_fraction_bind))))
for(i in 1:length(unique(names(hla_fraction_bind)))){
  mean.hla.bind[i] <- mean(hla_fraction_bind[names(hla_fraction_bind) == unique(names(hla_fraction_bind))[i]])
}
names(mean.hla.bind) <- unique(names(hla_fraction_bind))

pepbind.overlap.tcia <- list()
tcia.pep.patient <- list()
pep.overlap.tcia <- list()
for(i in 1:length(patient.sample)){
  pepbind.overlap.tcia[[i]] <- intersect(tcia[tcia$patientBarcode == substr(patient.sample[i],1,12),]$peptide,rownames(bind.matrix)[bind.matrix[,i] > 0])
  tcia.pep.patient[[i]] <- tcia[tcia$patientBarcode == substr(patient.sample[i],1,12),]$peptide
  pepmet <- fread(paste(patient.sample[i],"_",j,"mer_neoantmetadata.txt",sep=""))
  pepmet <- pepmet[-grep('_',pepmet[[4]]),]
  pep.count <- table(pepmet[[4]])
  pep.overlap.tcia[[i]] <- intersect(tcia[tcia$patientBarcode == substr(patient.sample[i],1,12),]$peptide,names(pep.count))
  
}

#normalize for hla alleles
hla.alleles.per.patient <- vector(mode='numeric',length(patient.sample))
for(i in 1:length(patient.sample)){
  hla.alleles.per.patient[i] <- dim(netmhc_out_all_SB[[i]])[2]-1
}

#hla allele specific bound peptides for each patient
pep.bind.perhla <- list()
for(i in 1:length(patient.sample)){
  pep.bind.perhla[[i]] <- colSums(netmhc_out_all[[i]][,1:(dim(netmhc_out_all[[i]])[2]-1)] < 2,na.rm=TRUE)
}





