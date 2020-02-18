#process your pdb input fasta into a vector of unique sequences of proteins

#pdb.unique <- vector(mode='character',0)
#for(i in 1:length(pdb.seq)){
#  if(i %% 1000 == TRUE){print(i)}
#  temp <- try(grepl(pdb.seq[i],pdb.seq[i+1],perl=TRUE))
#  if(temp == FALSE){
#    pdb.unique <- c(pdb.unique,pdb.seq[i])
#  }
#}


setwd('/users/PAS1203/osu1042/rnaneoantigen/unique_kmers')
args = commandArgs(trailingOnly=TRUE)

unique.kmer <- function(proteinfasta,kmer_length,aminoacid){ #generate unique peptides from a single protein fasta sequence
  unique.kmers <- vector(mode='character',0)
  protein.temp <- strsplit(proteinfasta,'')[[1]]
  #aa.location <- str_locate_all(proteinfasta,'A')[[1]]
  #protein.aa.temp <- matrix(NA,nrow=dim(aa.location)[1],ncol=kmer_length)
  #for(i in 1:kmer_length){
  #  protein.aa.temp[,i] <- protein.temp[aa.location[,1]+i-1]
  #} 
  #protein.aa.unique <- unique(protein.aa.temp)
  #return(matrix(do.call(paste0, as.data.frame(protein.aa.unique))))
  for(i in 1:(length(protein.temp)-kmer_length+1)){
    if(protein.temp[i] == aminoacid){
      unique.kmers <- union(unique.kmers,substr(proteinfasta,i,i+kmer_length-1))
    }
  }
  return(unique.kmers)
}


fasta.to.kmer <- function(pdb,kmer_length,aminoacid){
  ptm <- proc.time()
  unique.kmers <- vector(mode='character',0)
  for(i in 1:length(pdb)){
    if(nchar(pdb[i]) >= kmer_length){
      unique.kmers.temp <- unique.kmer(pdb[i],kmer_length,aminoacid)
    }
    unique.kmers <- c(unique.kmers,unique.kmers.temp)
  }
  print(proc.time() - ptm)
  return(unique(unique.kmers))
}

kmer_length = 8
unique.peptides <- fasta.to.kmer(pdb.unique,kmer_length,args[1])
save.image(paste('unique.kmers',kmer_length,'.',args[1],'.RData',sep=''))









