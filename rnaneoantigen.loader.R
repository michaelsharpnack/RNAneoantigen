################################ rnaneoantigen.loader.R   Load in files  ########################################################

peptide.to.seq <- function(peplength,netmhc_out,pepmet,patient){
  #this function creates a list of contigs for each peptide. Steps:
  #1. find sub-sequence that codes for peptide within rnaseq read
  #2. Center reads on peptide coding sequence (maximum possible contig length)
  #2. Use results from #1 as seed to add on reads
  
  #create pep to codon table
  codon.table <- list()
  codon.table[["F"]] <- c('TTC','TTT')
  codon.table[['L']] <- c('CTA','CTC','CTG','CTT','TTA','TTG')
  codon.table[['I']] <- c('ATA','ATC','ATT')
  codon.table[['M']] <- 'ATG'
  codon.table[['V']] <- c('GTA','GTC','GTG','GTT')
  codon.table[['S']] <- c('TCA','TCC','TCG','TCT','AGT','AGC')
  codon.table[['P']] <- c('CCA','CCC','CCG','CCT')
  codon.table[['T']] <- c('ACA','ACC','ACG','ACT')
  codon.table[['A']] <- c('GCA','GCC','GCG','GCT')
  codon.table[['Y']] <- c('TAC','TAT')
  codon.table[['H']] <- c('CAT','CAC')
  codon.table[['Q']] <- c('CAA','CAG')
  codon.table[['N']] <- c('AAT','AAC')
  codon.table[['K']] <- c('AAA','AAG')
  codon.table[['D']] <- c('GAT','GAC')
  codon.table[['E']] <- c('GAA','GAG')
  codon.table[['C']] <- c('TGT','TGC')
  codon.table[['W']] <- 'TGG'
  codon.table[['R']] <- c('CGA','CGC','CGG','CGT','AGA','AGG')
  codon.table[['G']] <- c('GGA','GGC','GGG','GGT')
  
  bind.names <- rownames(netmhc_out[rowMins(as.matrix(netmhc_out),na.rm=TRUE) < 2,])
  pepmet <- pepmet[(which(pepmet[[4]] %in% bind.names)),]
  pepmet <- pepmet[order(pepmet[[4]]),]
  pepmet <- pepmet[,c(3:4)]
  
  pep.contig <- vector(mode='character',length(bind.names))
  names(pep.contig) <- bind.names
  peponly.contig <- pep.contig
  contig.pepstart <- pep.contig
  #loop through peptides
  for(y in 1:length(bind.names)){
    #print(paste('analyzing peptide',y,bind.names[y]))
    
    #unique rna reads corresponding for peptide y
    reads.temp <- pepmet[[1]][which(pepmet[[2]] == bind.names[y])]
    unique.reads.freq <- table(reads.temp)
    unique.reads <- unique(reads.temp)
    
    pep.rna.seq <- matrix(NA,nrow=length(unique.reads),ncol=4)  #RNAseq locations and sequences for each peptide
    #figure out read corresponding to each peptide
    for(o in 1:length(unique.reads)){
      n = 1
      for(z in 1:nchar(unique.reads[o])){
        #print(paste('analyzing RNA read position',z))
        if(substr(unique.reads[o],z,z+2) %in% codon.table[[substr(bind.names[y],n,n)]]){
          #print('peptide matched, stop on read and check subsequent peptides')
          n.temp <- n+1
          for(m in n:(peplength-1)){
            if(substr(unique.reads[o],z+m*3,z+m*3+2) %in% codon.table[[substr(bind.names[y],n.temp,n.temp)]]){
              if(m < peplength-1){n.temp <- n.temp+1}
            } else {
              #print('subsequent peptides do no match, move down RNA read')
              break
            }
          }
          if(n.temp == peplength){
            #print('match found!')
            #print(paste(as.numeric(unique.reads.entries[[2]][o]),as.numeric(unique.reads.entries[[2]][o])+z-1,z))
            pep.rna.seq[o,] <- c(z,substr(unique.reads[o],z,z+peplength*3-1),unique.reads[o],unique.reads.freq[o])
            break
          }
        }
      }
    }
    pep.rna.seq <- pep.rna.seq[order(as.numeric(pep.rna.seq[,1])),]
    read.length.max <- max(nchar(pep.rna.seq[,3]))
    contig.length <- peplength*3+(read.length.max-peplength*3)*2
    pep.rna.seq.full <- matrix(NA,nrow=dim(pep.rna.seq)[1],ncol=contig.length)
    for(counter in 1:dim(pep.rna.seq.full)[1]){
      read.length <- nchar(pep.rna.seq[counter,3])
      pep.rna.seq.full[counter,(read.length.max-peplength*3-as.numeric(pep.rna.seq[counter,1])+2):(read.length.max-peplength*3-as.numeric(pep.rna.seq[counter,1])+1+read.length)] <- strsplit(pep.rna.seq[counter,3],'')[[1]]
    }
    pep.rna.seq.unique <- list()
    for(counter in 1:dim(pep.rna.seq.full)[2]){
      pep.temp <- cbind(pep.rna.seq[,4],pep.rna.seq.full[,counter])[is.na(pep.rna.seq.full[,counter]) == FALSE,]
      if(is.null(dim(pep.temp))){
        pep.rna.seq.unique[[counter]] <- as.numeric(pep.temp[1])
        names(pep.rna.seq.unique[[counter]]) <- pep.temp[2]
      } else {
        pep.rna.seq.unique[[counter]] <- vector(mode='numeric',length(unique(pep.temp[,2])))
        names(pep.rna.seq.unique[[counter]]) <- unique(pep.temp[,2])
        for(counter2 in 1:length(pep.rna.seq.unique[[counter]])){
          pep.rna.seq.unique[[counter]][counter2] <- sum(as.numeric(pep.temp[pep.temp[,2] == unique(pep.temp[,2])[counter2],1]))
        }
      }
    }
    #create the contig with winner takes all approach
    contig  <- vector(mode='character',length(pep.rna.seq.unique))
    for(counter in 1:length(pep.rna.seq.unique)){
      contig[counter] <- names(which.max(pep.rna.seq.unique[[counter]]))
    }
    peponly.contig[y] <- paste(contig[(read.length.max-peplength*3+1):(read.length.max)],collapse='')
    contig.pepstart[y] <- read.length.max-peplength*3+1-sum(nchar(contig[1:read.length.max]) == 0)
    pep.contig[y] <- paste(contig,collapse = '')
    
  }
  contig.write <- vector(mode='character',length(pep.contig)*2)
  for(i in 1:length(pep.contig)){
    contig.write[2*i-1] <- paste(">",names(pep.contig)[i],sep=' ')
    contig.write[2*i] <- pep.contig[i]
  }
  write.table(contig.write,file=paste('contig/',patient,peplength,sep=''),quote=FALSE,row.names = FALSE,col.names=FALSE)
  pep.contig.save <- cbind(pep.contig,peponly.contig,contig.pepstart)
  rownames(pep.contig.save) <- names(pep.contig)
  
  return(list(contig.write,pep.contig.save))
}


netmhc.loader <- function(peplength,patient){
  col.names <- c('Pos','HLA','Peptide','Core','Of','Gp','Gl','Ip','ll','lcore','Identity','Score','Aff(nM)','%Rank','BindLevel')

  #load netmhcpan output
  netmhc <- read.table(paste('TCGA_neoantigens_',peplength,'mer/',patient,"_",peplength,"mer_neoantigens.txt",sep=""),header=FALSE,fill=TRUE,col.names = col.names,sep=',')
  netmhc <- data.frame(netmhc[[3]],netmhc[[2]],netmhc[[13]])
  colnames(netmhc) <- c('V1','V2','V3')
  netmhc <- dcast(netmhc,V1 ~ V2)
  rownames(netmhc) <- netmhc[[1]]
  netmhc <- netmhc[,-1]
  netmhc <- netmhc[-grep('X',rownames(netmhc)),]
      
  #load in peptide counts
  pepmet <- fread(paste('TCGA_neoantigens_',peplength,'mer/',patient,"_",peplength,"mer_neoantmetadata.txt",sep=""))
  pepmet <- pepmet[-grep('_',pepmet[[4]]),]
  pep.count <- table(pepmet[[4]])
  if(length(pep.count) != nrow(netmhc)){print('warning: number of peptides in metadata and neoantigens run through netmhcpan are not equal')}
  pep.count <- pep.count[intersect(names(pep.count),rownames(netmhc))]
  contig <- peptide.to.seq(peplength,netmhc,pepmet,patient)
  netmhc <- cbind(netmhc,as.vector(pep.count))
  
  return(list(netmhc,contig[[1]],contig[[2]]))
}



