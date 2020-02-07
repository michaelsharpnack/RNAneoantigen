################################ rnaneoantigen.pepfinder.analysis.R ################################
library(data.table)
library(matrixStats)
neoantburden = matrix(0,nrow=)
for(i in 1:length(peps.found.tumor.kmer[[2]])){
  peps.found.tumor.kmer
  
  col.names <- c('Pos','HLA','Peptide','Core','Of','Gp','Gl','Ip','ll','lcore','Identity','Score','Aff(nM)','%Rank','BindLevel')
  #load netmhcpan output
  netmhc <- read.table(paste('Desktop/Desktop/netmhcpan_out/',patients.netmhc[i],'_clean',sep=""),header=FALSE,fill=TRUE,col.names = col.names,sep=',')
  netmhc <- data.frame(netmhc[[3]],netmhc[[2]],netmhc[[13]])
  colnames(netmhc) <- c('V1','V2','V3')
  netmhc <- dcast(netmhc,V1 ~ V2,fun.aggregate = mean,value.var = 'V3')
  rownames(netmhc) <- netmhc[[1]]
  netmhc <- netmhc[,-1]
  temp <- peps.found.normal.kmer[[2]][[i]][[2]][[1]]
  netmhc <- cbind(netmhc,matrix(0,nrow=dim(netmhc)[1],ncol=2))
  netmhc[temp[[2]],7] <- temp[[1]]
  temp <- peps.found.tumor.kmer[[2]][[i]][[2]][[1]]
  netmhc[temp[[2]],8] <- temp[[1]]
}

neoantigens <- matrix(0,nrow=39,ncol=4)
for(i in 1:39){
  if(is.null(dim(netmhc.out[[i]])) == FALSE){
    neoantigens[i,1] <- sum(colMins(netmhc.out[[i]][,1:(dim(netmhc.out[[i]])[2]-2)]) < 2 & netmhc.out[[i]][,dim(netmhc.out[[i]])[2]-1] > 0)
    neoantigens[i,2] <- sum(colMins(netmhc.out[[i]][,1:(dim(netmhc.out[[i]])[2]-2)]) < 2 & netmhc.out[[i]][,dim(netmhc.out[[i]])[2]] > 0)
    neoantigens[i,3] <- sum(colMins(netmhc.out[[i]][,1:(dim(netmhc.out[[i]])[2]-2)]) < 2 & netmhc.out[[i]][,dim(netmhc.out[[i]])[2]-1] > 0 & netmhc.out[[i]][,dim(netmhc.out[[i]])[2]] > 0)
    neoantigens[i,4] <- neoantigens[i,3]/neoantigens[i,1]
  }
}

filesize = matrix(0,nrow=39,ncol=2)
for(i in 1:39){
  print(i)
  temp <- dir(paste(your_file_dir,'TCGA_neoantigens_',9,'mer/','tumor','/',sep=''))
  your_file_location <- paste(your_file_dir,'TCGA_neoantigens_',9,'mer/','tumor','/',temp[i],sep='')
  filesize[i,1] = as.integer(system(paste('wc -l < ',your_file_location,sep = ''),intern = TRUE)[[1]])
  
  temp <- dir(paste(your_file_dir,'TCGA_neoantigens_',9,'mer/','normal','/',sep=''))
  your_file_location <- paste(your_file_dir,'TCGA_neoantigens_',9,'mer/','normal','/',temp[i],sep='')
  filesize[i,2] = as.integer(system(paste('wc -l < ',your_file_location,sep = ''),intern = TRUE)[[1]])
}
neoants.Nsearched <- vector(mode='numeric',39)
for(i in 1:39){
  if(is.null(dim(netmhc.out[[i]])) == FALSE){
    neoants.Nsearched[i] <- dim(netmhc.out[[i]])[1]
  }
}
total.expr <- matrix(0,nrow=39,ncol=3)
for(i in 1:39){
  total.expr[i,1] <- peps.found.tumor.kmer[[2]][[i]][[2]][[2]]
  total.expr[i,2] <- peps.found.normal.kmer[[2]][[i]][[2]][[2]]
  total.expr[i,3] <- sum(peps.found.tumor.kmer[[2]][[i]][[2]][[1]][,1])
}
for(i in 1:39){
  temp.pep <- union(peps.found.tumor.kmer[[2]][[i]][[2]][[1]][,2],peps.found.normal.kmer[[2]][[i]][[2]][[1]][,2])
  if(length(temp.pep) > 0){
    pep.out <- vector(mode='character',length(temp.pep)*2)
    for(j in 1:length(temp.pep)){
      pep.out[2*(j-1)+1] = paste('>',temp.pep[j])
      pep.out[2*j] = temp.pep[j]
    }
  }
  write.table(pep.out,sep=',',quote = FALSE,row.names = FALSE,col.names = FALSE,file=paste(your_file_dir,'peps.tumornormal/',patients.netmhc[i],sep=''))
}




