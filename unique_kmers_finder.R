########################### unique_kmers_finder.R ####################################

#this loop converts the RData file into the list of peptides
#setwd('/Volumes/Michael-backup3/rnaneoantigen/references/unique_kmers/RData/')
#files <- dir()
#for(i in 1:length(files)){
#  print(i)
#  load(files[i])
#  write.table(unique.peptides,file=paste('/Volumes/Michael-backup3/rnaneoantigen/references/unique_kmers/',substr(files[i],1,nchar(files[i])-6),sep=''),col.names = FALSE,row.names = FALSE,quote=FALSE)
#}

setwd('/Volumes/Michael-backup3/rnaneoantigen/references/unique_kmers/')

reference.peptide.search <- function(peptide){
  setwd('reference_peptides/')
  temp <- try(system(paste('grep ',peptide,' unique.kmers',nchar(peptide),'.',substr(peptide,1,1),sep=''),intern=TRUE)[[1]],silent = TRUE)
  setwd('..')
  return(temp)
}

load('pepfinder_01212020.RData') 
patients <- seq(from=1,to=39,by=1)
for(j in 5:length(patients)){
  print(j)
  peptides <- union(peps.found.tumor.kmer[[1]][[patients[j]]][[1]],peps.found.tumor.kmer[[2]][[patients[j]]][[1]])
  peptides <- peptides[nchar(peptides) >= 8]
  peptide.reference.found <- vector(mode='character',length(peptides))
  for(i in 1:length(peptides)){
    #if(i %% 100 == 0){print(i)}
    peptide.reference.found[i] <- reference.peptide.search(peptides[i])[[1]]
  }
  save.image(paste('peptides_found/',patients[j],'.RData',sep=''))
}
