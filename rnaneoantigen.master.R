######################## rnaneoantigen.master.R ################################
library(data.table)
library(matrixStats)
library(reshape2)
library(biomaRt)

wd <- '/Volumes/Michael-backup3/rnaneoantigen/'
setwd(wd)
source('RNAneoantigen/rnaneoantigen.loader.R')
source('RNAneoantigen/rnaneoantigen.annotate.R')
patients <- fread('tcga.luad.patients.csv',header=FALSE)[[1]]
gtf = gtf.load(wd)
ccds = framed.load(wd)
peplengths <- c(8,9)

netmhc.out <- list()
blastout.out <- list()
for(i in 1:length(peplengths)){
  peplength <- peplengths[i]
  blastout.out.temp <- list()
  netmhc.out.temp <- list()
  for(j in 1:length(patients)){
    print(c(i,j))
    patient <- patients[j]
    netmhc <- netmhc.loader(peplength,patient)
    blastout <- blast(patient,peplength,wd,max.mutations = 2,max.gaps = 1,onlytophits = TRUE)
    blastout.fasta <- fastafromblast(patient,wd,blastout,netmhc[[2]],netmhc[[3]],peplength)
    blastout <- junctioncaller(blastout.fasta,wd,patient,netmhc[[2]],netmhc[[3]],peplength)
    blastout.fasta <- fastafromblast(patient,wd,blastout,netmhc[[2]],netmhc[[3]],peplength)
    blastout.fasta.variant <- variantcaller(patient,wd,blastout.fasta,peplength,max.mutations = 2, max.gaps = 1)
    blastout.fasta.variant.gtf <- gtf.blast(blastout.fasta.variant,gtf[[1]],gtf[[2]],gtf[[3]])
    blastout.fasta.variant.gtf.framed <- framed(ccds[[2]],blastout.fasta.variant.gtf)
    blastout.out.temp[[j]] <- blastout.fasta.variant.gtf.framed
    netmhc.out.temp[[j]] <- netmhc
    rm(blastout,blastout.fasta,blastout.fasta.variant,blastout.fasta.variant.gtf,blastout.fasta.variant.gtf.framed)
    save.image(paste('RData/rnaneoantigen.master.RData',sep=''))
  }
  blastout.out[[i]] <- blastout.out.temp
  netmhc.out[[i]] <- netmhc.out.temp
}



