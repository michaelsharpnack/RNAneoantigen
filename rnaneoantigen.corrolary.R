#load in corrolary data for RNA neoantigen data


setwd(wd)

#neoantigen data
neoantigens.data <- fread('corrollary_files/luad_neoantigens.csv')
neoantigens.data$Individual_ID <- substr(neoantigens.data$Individual_ID,6,17)
neoantigens.data <- as.matrix(neoantigens.data)
rownames(neoantigens.data) <- neoantigens.data[,1]
neoantigens.data <- neoantigens.data[intersect(substr(patients,1,12),rownames(neoantigens.data)),]
hla.neoantigens <- neoantigens.data[,18]
clinical.data <- neoantigens.data[,c(27:30)]
neoantigens.data <- neoantigens.data[,-c(1,18,27:30)]
class(neoantigens.data) <- 'numeric'

#load in maf file
maf <- list()
wd <- '/Users/michaelsharpnack/Dropbox/rnaeoantigen/corrollary_files/LUAD.maf/'
maf.to.pep <- function(wd,patient.sample)
  setwd(wd)
  for(u in 1:length(patient.sample)){
    maf[[u]] <- fread(paste(substr(patient.sample[u],1,12),'-01.maf.txt',sep=''))
  }
  #create bed file from maf file
  
  

hla.neoantigens.mat <- matrix(NA,nrow=length(patient.sample),ncol=12)
rownames(hla.neoantigens.mat) <- substr(patient.sample,1,12)
colnames(hla.neoantigens.mat) <- c('Exome: HLA-A1','Exome: HLA-A2','Exome: HLA-B1','Exome: HLA-B2','Exome: HLA-C1','Exome: HLA-C2',
                                   'RNAseq: HLA-A1', 'RNAseq: HLA-A2', 'RNAseq: HLA-B1', 'RNAseq: HLA-B2', 'RNAseq: HLA-C1', 'RNAseq: HLA-C2')
for(i in 1:length(patient.sample)){
  split.temp <- strsplit(hla.neoantigens[i],',')[[1]]
  hla.neoantigens.mat[i,1:length(grep('HLA-A',split.temp))] <- split.temp[grep('HLA-A',split.temp)]
  hla.neoantigens.mat[i,3:(length(grep('HLA-B',split.temp))+2)] <- split.temp[grep('HLA-B',split.temp)]
  hla.neoantigens.mat[i,5:(length(grep('HLA-C',split.temp))+4)] <- split.temp[grep('HLA-C',split.temp)]
  try(hla.neoantigens.mat[i,7:(length(grep('HLA-A',colnames(netmhc_out_all[[i]])))+6)] <- colnames(netmhc_out_all[[i]])[grep('HLA-A',colnames(netmhc_out_all[[i]]))])
  try(hla.neoantigens.mat[i,9:(length(grep('HLA-B',colnames(netmhc_out_all[[i]])))+8)] <- colnames(netmhc_out_all[[i]])[grep('HLA-B',colnames(netmhc_out_all[[i]]))])
  try(hla.neoantigens.mat[i,11:(length(grep('HLA-C',colnames(netmhc_out_all[[i]])))+10)] <- colnames(netmhc_out_all[[i]])[grep('HLA-C',colnames(netmhc_out_all[[i]]))])
}
hla.neoantigens.mat <- sub('\\*','',hla.neoantigens.mat)

#timer immune estimation data
timer <- fread('corrollary_files/immuneEstimation.txt')
timer <- as.matrix(timer)
rownames(timer) <- timer[,1]
timer <- timer[,-1]
class(timer) <- 'numeric'
timer <- timer[substr(rownames(timer),14,15) == '01',]
rownames(timer) <- substr(rownames(timer),1,12)
timer <- timer[intersect(substr(patients,1,12),rownames(timer)),]

#tumor purity data
tumor.purity <- as.matrix(fread('corrollary_files/tumorpurity.csv'))
rownames(tumor.purity) <- substr(tumor.purity[,1],1,12)
tumor.purity <- tumor.purity[,-c(1:2)]
class(tumor.purity) <- 'numeric'
tumor.purity <- tumor.purity[intersect(substr(patients,1,12),rownames(tumor.purity)),]

neoantigens.data <- cbind(neoantigens.data,tumor.purity[,5],timer)

npeps_cor <- matrix(nrow=dim(neoantigens.data)[2],ncol=3)
for(l in 1:dim(neoantigens.data)[2]){
  npeps_cor[l,1] <- cor(neoantigens.data[,l],npeps[,1])
  npeps_cor[l,2] <- cor(neoantigens.data[,l],npeps[,2])
  npeps_cor[l,3] <- cor(neoantigens.data[,l],npeps[,1]/npeps[,2])
}
rownames(npeps_cor) <- colnames(neoantigens.data)
colnames(npeps_cor) <- c('npeps','npeps_bind','npeps_frac')

#TCIA neoantigen data
tcia <- fread('/Users/michaelsharpnack/Dropbox/rnaeoantigen/corrollary_files/TCIA-NeoantigensData copy.txt')
tcia <- tcia[,-5]
tcia <- tcia[nchar(tcia[[4]]) == 8,]
  
#input to protein blast:
tcia.pblast <- data.frame(matrix(nrow=nrow(tcia)*2,ncol=1))
x = 0
for(i in 1:dim(tcia)[1]){
  if(i != 1 && (tcia$gene[i] == tcia$gene[i-1] & tcia$patientBarcode[i] == tcia$patientBarcode[i-1])){x = x+1} else {x = 0}
  tcia.pblast[[1]][2*i-1] <- paste('>',tcia$patientBarcode[i],tcia$gene[i],x)
  tcia.pblast[[1]][2*i] <- tcia$peptide[i]
}
write.csv(tcia.pblast,file='/Users/michaelsharpnack/Dropbox/rnaeoantigen/corrollary_files/tcia.pblast.csv')




