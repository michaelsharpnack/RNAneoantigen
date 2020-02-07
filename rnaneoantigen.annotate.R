##################################### annotation ####################################################
#annotate your blast output with genes/exons
#first, load and process gtf file
gtf.load <- function(wd){
  gtf <- read.table(paste(wd,'/references/gencode.v22.annotation.gtf',sep=''),sep=";",fill=TRUE,header = F,col.names=c(1:24),colClasses = c("character",rep('NULL',23)))
  gtf.mat <- matrix(NA,nrow=length(gtf[[1]]),ncol=6)
  gtf.temp <- strsplit(gtf[[1]],'\t')
  for(i in 1:length(gtf[[1]])){
    gtf.mat[i,] <- gtf.temp[[i]][c(1,3,4,5,7,9)]
  }
  gtf.gene <- gtf.mat[grep('gene',gtf.mat[,2]),]
  gtf.exon <- gtf.mat[grep('exon',gtf.mat[,2]),]
  rm(gtf,gtf.mat,gtf.temp)
  #split gtf objects into chromosomes
  gtf.gene.list <- list()
  gtf.exon.list <- list()
  for(i in 1:length(unique(gtf.gene[,1]))){
    gtf.gene.list[[i]] <- gtf.gene[gtf.gene[,1] == unique(gtf.gene[,1])[i],]
    gtf.exon.list[[i]] <- gtf.exon[gtf.exon[,1] == unique(gtf.exon[,1])[i],]
  }
  names(gtf.gene.list) <- unique(gtf.gene[,1])
  names(gtf.exon.list) <- unique(gtf.exon[,1])
  
  mart <- useMart("ENSEMBL_MART_ENSEMBL")
  mart <- useDataset("hsapiens_gene_ensembl", mart)
  
  return(list(gtf.gene.list,gtf.exon.list,mart))
  rm(gtf.exon,gtf.gene)
}

#blast output of contigs
blast <- function(patient,peplength,wd,max.mutations = 2,max.gaps = 1,onlytophits = TRUE){
  setwd(paste('references/',sep=''))
  #blast your files against reference sequence
  blastout<- system(paste('./blastn -strand both -reward 1 -penalty -1 -gapextend 1 -gapopen 2 -db ','human_grch38.fa',' -query ',wd,'/contig/',patient,peplength,' -outfmt 6',sep=''),intern=TRUE)
  temp <- strsplit(blastout,'\t')
  blastout.mat <- matrix(NA,nrow=length(temp),ncol=length(temp[[1]]))
  for(tempcount in 1:length(temp)){
    blastout.mat[tempcount,] <- temp[[tempcount]]
  }
  col.names.blast <- c('qseqid','sseqid','% identical matches','alignment length','number of mismatches','number of gap openings','query start','query end','subject start','subject end','evalue','bitscore')
  colnames(blastout.mat) <- col.names.blast
  blastout <- blastout.mat
  rm(blastout.mat,temp)
  setwd('..')
  #filter out blast results with multiple mutations
  blastout <- blastout[as.numeric(blastout[,5]) <= max.mutations,]
  #max number of gaps in the alignment
  blastout <- blastout[as.numeric(blastout[,6]) <= max.gaps,]
  #only take hits that have top score (keeps alignments that are ties)
  blastout.temp <- matrix(NA,nrow=0,ncol=dim(blastout)[2])
  for(i in 1:length(unique(blastout[,1]))){
    blastout.temp <- rbind(blastout.temp,blastout[blastout[,1] == unique(blastout[,1])[i] & as.numeric(blastout[,12]) == max(as.numeric(blastout[blastout[,1] == unique(blastout[,1])[i],12])),])
  }
  blastout <- blastout.temp
  
  return(blastout)
}

junctioncaller <- function(blastout,wd,patient,contig,pep.contig,peplength){
  #reblast subsequences with long seqments that do not align and insert them into blastout.mat
  reblast.front <- blastout[as.numeric(blastout[,7]) > 10,]
  reblast.back <- blastout[(max(as.numeric(blastout[,8]))-as.numeric(blastout[,8])) >= 10,]
  
  contig.write <- vector(mode='character',dim(reblast.front)[1]*2)
  for(i in 1:dim(reblast.front)[1]){
    contig.write[2*i-1] <- paste(">",reblast.front[i,1],sep=' ')
    contig.write[2*i] <- substr(reblast.front[i,17],1,as.numeric(reblast.front[i,7])-1)
  }
  write.table(contig.write,file=paste('contig.front/',patient,peplength,sep=''),quote=FALSE,row.names = FALSE,col.names=FALSE)
  
  contig.write <- vector(mode='character',dim(reblast.back)[1]*2)
  for(i in 1:dim(reblast.back)[1]){
    contig.write[2*i-1] <- paste(">",reblast.back[i,1],sep=' ')
    contig.write[2*i] <- substr(reblast.back[i,17],as.numeric(reblast.back[i,8])+1,nchar(reblast.back[i,17]))
  }
  write.table(contig.write,file=paste('contig.back/',patient,peplength,sep=''),quote=FALSE,row.names = FALSE,col.names=FALSE)
  
  setwd(paste('references/',sep=''))
  #blast your files against reference sequence
  blastout.temp <- system(paste('./blastn -strand both -reward 1 -penalty -1 -gapextend 1 -gapopen 2 -db ','human_grch38.fa',' -query ',wd,'/contig.front/',patient,peplength,' -outfmt 6',sep=''),intern=TRUE)
  temp <- strsplit(blastout.temp,'\t')
  blastout.mat <- matrix(NA,nrow=length(temp),ncol=length(temp[[1]]))
  for(tempcount in 1:length(temp)){
    blastout.mat[tempcount,] <- temp[[tempcount]]
  }
  col.names.blast <- c('qseqid','sseqid','% identical matches','alignment length','number of mismatches','number of gap openings','query start','query end','subject start','subject end','evalue','bitscore')
  colnames(blastout.mat) <- col.names.blast
  blastout.front <- blastout.mat
  rm(blastout.mat,temp)
  setwd('..')
  
  setwd(paste('references/',sep=''))
  #blast your files against reference sequence
  blastout.temp <- system(paste('./blastn -strand both -reward 1 -penalty -1 -gapextend 1 -gapopen 2 -db ','human_grch38.fa',' -query ',wd,'/contig.back/',patient,peplength,' -outfmt 6',sep=''),intern=TRUE)
  temp <- strsplit(blastout.temp,'\t')
  blastout.mat <- matrix(NA,nrow=length(temp),ncol=length(temp[[1]]))
  for(tempcount in 1:length(temp)){
    blastout.mat[tempcount,] <- temp[[tempcount]]
  }
  col.names.blast <- c('qseqid','sseqid','% identical matches','alignment length','number of mismatches','number of gap openings','query start','query end','subject start','subject end','evalue','bitscore')
  colnames(blastout.mat) <- col.names.blast
  blastout.back <- blastout.mat
  rm(blastout.mat,temp)
  setwd('..')
  
  
  
}

#getfasta of your new blasted segments and compare them with your contigs
fastafromblast <- function(patient,wd,blastout,contig,pep.contig,peplength){
  #write bed file
  strand <- rep('-',dim(blastout)[1])
  #add strand information
  strand[as.numeric(blastout[,9]) < as.numeric(blastout[,10])] <- '+'
  blastout.temp <- blastout
  #reverse strand to find correct sequence in reference genome
  blastout.temp[,c(9,10)][as.numeric(blastout[,9]) > as.numeric(blastout[,10]),] <- blastout[,c(10,9)][as.numeric(blastout[,9]) > as.numeric(blastout[,10]),]
  blastout.temp[,9] <- as.numeric(blastout.temp[,9])-1
  write.table(cbind(blastout.temp[,c(2,9,10,1,4)],strand),file=paste(wd,'bedfiles/',patient,sep=''),row.names=FALSE,quote = FALSE,sep='\t',col.names=FALSE)
  setwd('references/')
  fasta <- system(paste('./bedtools getfasta -fi ',wd,'references/human_grch38.fa -bed ',wd,'bedfiles/',patient,sep=''),intern=TRUE)
  setwd('..')
  
  #add in the contigs column
  contig.named <- contig[seq(from=2,to=length(contig),by=2)]
  names(contig.named) <- substr(contig[seq(from=1,to=length(contig),by=2)],3,2+peplength)
  
  #add in the pep.contigs column
  fasta.temp <- cbind(strand,fasta[seq(from=2,to=length(fasta),by=2)],pep.contig[blastout[,1],c(2,3)],contig.named[blastout[,1]])
  colnames(fasta.temp) <- c('Strand','GRC38 reference sequence','RNAseq peptide contig','Peptide start on which base of contig','Contig sequence')
  
  #reverse - stranded sequence and complement
  complementer <- list()
  complementer[['A']] <- 'T'
  complementer[['C']] <- 'G'
  complementer[['G']] <- 'C'
  complementer[['T']] <- 'A'
  
  for(i in 1:dim(fasta.temp)[1]){
    if(fasta.temp[i,1] == '-'){
      temp.bases <- rev(strsplit(toupper(fasta.temp[i,2]),'')[[1]])
      for(j in 1:length(temp.bases)){
        temp.bases[j] <- complementer[[temp.bases[j]]]
      }
      fasta.temp[i,2] <- paste(temp.bases,collapse='')
    }
  }
  
  return(cbind(blastout,fasta.temp))
  
}



#second, match gtf features to blast locations
gtf.blast <- function(blastout,gtf.gene.list,gtf.exon.list,mart){
  blastout.temp <- matrix(NA,nrow=dim(blastout)[1],ncol=8)
  for(i in 1:dim(blastout)[1]){
    if(as.numeric(blastout[i,10]) > as.numeric(blastout[i,9])){
      temp <- gtf.gene.list[[blastout[i,2]]]
      temp <- temp[temp[,5] == '+',]
      temp1 <- temp[(as.numeric(blastout[i,9]) > as.numeric(temp[,3]) & as.numeric(blastout[i,9]) < as.numeric(temp[,4])) | (as.numeric(blastout[i,10]) > as.numeric(temp[,3]) & as.numeric(blastout[i,10]) < as.numeric(temp[,4])),6]
      if(length(temp1) == 1){blastout.temp[i,1] <- gsub("\\.[0-9]*$", "",substr(temp1,9,nchar(temp1)[1]))
      } else if(length(temp1) == 2){blastout.temp[i,c(1:2)] <- gsub("\\.[0-9]*$", "",substr(temp1,9,nchar(temp1)[1]))}
      
      temp <- gtf.exon.list[[blastout[i,2]]]
      temp <- temp[temp[,5] == '+',]
      temp1 <- temp[(as.numeric(blastout[i,9]) > as.numeric(temp[,3]) & as.numeric(blastout[i,9]) < as.numeric(temp[,4])) | (as.numeric(blastout[i,10]) > as.numeric(temp[,3]) & as.numeric(blastout[i,10]) < as.numeric(temp[,4])),6]
      if(length(temp1) == 1){blastout.temp[i,3] <- substr(temp1,9,nchar(temp1)[1])
      } else if(length(temp1) == 2){blastout.temp[i,c(3:4)] <- substr(temp1,9,nchar(temp1)[1])}
      
    } else {
      temp <- gtf.gene.list[[blastout[i,2]]]
      temp <- temp[temp[,5] == '-',]
      temp1 <- temp[(as.numeric(blastout[i,9]) > as.numeric(temp[,3]) & as.numeric(blastout[i,9]) < as.numeric(temp[,4])) | (as.numeric(blastout[i,10]) > as.numeric(temp[,3]) & as.numeric(blastout[i,10]) < as.numeric(temp[,4])),6]
      if(length(temp1) == 1){blastout.temp[i,1] <- gsub("\\.[0-9]*$", "",substr(temp1,9,nchar(temp1)[1]))
      } else if(length(temp1) == 2){blastout.temp[i,c(1:2)] <- gsub("\\.[0-9]*$", "",substr(temp1,9,nchar(temp1)[1]))}
       
      temp <- gtf.exon.list[[blastout[i,2]]]
      temp <- temp[temp[,5] == '-',]
      temp1 <- temp[(as.numeric(blastout[i,9]) > as.numeric(temp[,3]) & as.numeric(blastout[i,9]) < as.numeric(temp[,4])) | (as.numeric(blastout[i,10]) > as.numeric(temp[,3]) & as.numeric(blastout[i,10]) < as.numeric(temp[,4])),6]
      if(length(temp1) == 1){blastout.temp[i,3] <- substr(temp1,9,nchar(temp1)[1])
      } else if(length(temp1) == 2){blastout.temp[i,c(3:4)] <- substr(temp1,9,nchar(temp1)[1])}
      
    }
  }
  
  #use biomart to add in gene symbols
  ens <- c(blastout.temp[,1][is.na(blastout.temp[,1]) == FALSE],blastout.temp[,2][is.na(blastout.temp[,2]) == FALSE])
  ensLookup <- unique(ens)
  annotLookup <- getBM(
    mart=mart,
    attributes=c("ensembl_transcript_id", "ensembl_gene_id", "gene_biotype", "external_gene_name"),
    filter="ensembl_gene_id",
    values=ensLookup,
    uniqueRows=TRUE)
  annotLookup <- unique(annotLookup[,c(2,3,4)])
  for(i in 1:dim(blastout.temp)[1]){
    j = sum(is.na(blastout.temp[i,c(1:2)]) == FALSE)
    if(j == 1){
      
      temp2 <- annotLookup$gene_biotype[annotLookup$ensembl_gene_id == blastout.temp[i,1]]
      if(length(temp2) > 0){blastout.temp[i,5] <- temp2}
      
      temp2 <- annotLookup$external_gene_name[annotLookup$ensembl_gene_id == blastout.temp[i,1]]
      if(length(temp2) > 0){blastout.temp[i,6] <- temp2}
      
    } else if(j == 2){
      
      temp2 <- annotLookup$gene_biotype[annotLookup$ensembl_gene_id == blastout.temp[i,1]]
      if(length(temp2) > 0){blastout.temp[i,5]}
      
      temp2 <- annotLookup$external_gene_name[annotLookup$ensembl_gene_id == blastout.temp[i,1]]
      if(length(temp2) > 0){blastout.temp[i,6]}
      
      temp2 <- annotLookup$gene_biotype[annotLookup$ensembl_gene_id == blastout.temp[i,2]]
      if(length(temp2) > 0){blastout.temp[i,7]}
      
      temp2 <- annotLookup$external_gene_name[annotLookup$ensembl_gene_id == blastout.temp[i,2]]
      if(length(temp2) > 0){blastout.temp[i,8]}
    }
  }
  colnames(blastout.temp) <- c('Gencode gene 1','Gencode gene 2','Gencode exon 1','Gencode exon 2','Gene type 1','Hugo Symbol 1','Gene type 2','Hugo symbol 2')
  blastout <- cbind(blastout,blastout.temp)
  
  
  return(blastout)
}



variantcaller <- function(patient,wd,blastout,peplength,max.mutations = 2, max.gaps = 1){
  
  #get start and stop of peptide coding sequence on the reference mapping
  start.coding <- as.numeric(blastout[,16])-as.numeric(blastout[,7])+1
  coding.bases <- as.numeric(blastout[,4])-start.coding+1

  variants <- list()
  #filter out any mappings that do not cover the peptide coding sequence -- will have to modify this step after we account for gapped alignment
  blastout.temp <- blastout[start.coding >= 1 & coding.bases >= peplength*3,]
  #add annotation columns
  annotation <- matrix(NA,nrow=dim(blastout.temp)[1],ncol=5*(max.mutations+max.gaps))
  annotation.labels <- c('Start_position','End_position','Variant_type','Reference Allele','Tumor Allele')
  colnames(annotation) <- rep(annotation.labels,max.mutations+max.gaps)
  for(i in 1:dim(annotation)[2]){
    colnames(annotation)[i] <- paste(colnames(annotation)[i],ceiling(i/5),sep=' ')
  }
  dim.temp <- dim(blastout.temp)[2]+1
  blastout.temp <- cbind(blastout.temp,vector(mode='character',dim(blastout.temp)[1]),annotation)
  colnames(blastout.temp)[dim.temp] <- 'Reference peptide'
  
  #first pull out any contigs that matched exactly into list[[1]]
  variants[[1]] <- blastout.temp[as.numeric(blastout.temp[,5]) == 0 & as.numeric(blastout.temp[,6]) == 0,]
  
  #second pull out any point mutations - sequences with 0 gap openings into list[[2]]
  variants[[2]] <- blastout.temp[as.numeric(blastout.temp[,5]) > 0 & as.numeric(blastout.temp[,6]) == 0,]
  #check if mutation is within peptide coding region, if not -- send the output back to variant[[1]]
  start.coding.temp <- as.numeric(variants[[2]][,16])-as.numeric(variants[[2]][,7])+1
  variants[[1]] <- rbind(variants[[1]],variants[[2]][substr(variants[[2]][,14],start.coding.temp,start.coding.temp+peplength*3-1) == variants[[2]][,15],])
  variants[[2]] <- variants[[2]][substr(variants[[2]][,14],start.coding.temp,start.coding.temp+peplength*3-1) != variants[[2]][,15],]
  variants[[2]] <- as.matrix(variants[[2]])
  if(dim(variants[[2]])[2] == 1){
    variants[[2]] <- t(variants[[2]])
  }
  
  #test if the SNV changes the predicted peptide sequence & send back the ones that don't to variant[[1]]
  aa.table <- list()
  aa.table[["TTC"]] <- 'F'
  aa.table[["TTT"]] <- 'F'
  aa.table[['CTA']] <- 'L'
  aa.table[['CTC']] <- 'L'
  aa.table[['CTG']] <- 'L'
  aa.table[['CTT']] <- 'L'
  aa.table[['TTA']] <- 'L'
  aa.table[['TTG']] <- 'L'
  aa.table[['ATA']] <- 'I'
  aa.table[['ATC']] <- 'I'
  aa.table[['ATT']] <- 'I'
  aa.table[['ATG']] <- 'M'
  aa.table[['GTA']] <- 'V'
  aa.table[['GTC']] <- 'V'
  aa.table[['GTG']] <- 'V'
  aa.table[['GTT']] <- 'V'
  aa.table[['TCA']] <- 'S'
  aa.table[['TCC']] <- 'S'
  aa.table[['TCG']] <- 'S'
  aa.table[['TCT']] <- 'S'
  aa.table[['AGT']] <- 'S'
  aa.table[['AGC']] <- 'S'
  aa.table[['CCA']] <- 'P'
  aa.table[['CCC']] <- 'P'
  aa.table[['CCG']] <- 'P'
  aa.table[['CCT']] <- 'P'
  aa.table[['ACA']] <- 'T'
  aa.table[['ACC']] <- 'T'
  aa.table[['ACG']] <- 'T'
  aa.table[['ACT']] <- 'T'
  aa.table[['GCA']] <- 'A'
  aa.table[['GCC']] <- 'A'
  aa.table[['GCG']] <- 'A'
  aa.table[['GCT']] <- 'A'
  aa.table[['TAC']] <- 'Y'
  aa.table[['TAT']] <- 'Y'
  aa.table[['CAT']] <- 'H'
  aa.table[['CAC']] <- 'H'
  aa.table[['CAA']] <- 'Q'
  aa.table[['CAG']] <- 'Q'
  aa.table[['AAT']] <- 'N'
  aa.table[['AAC']] <- 'N'
  aa.table[['AAA']] <- 'K'
  aa.table[['AAG']] <- 'K'
  aa.table[['GAT']] <- 'D'
  aa.table[['GAC']] <- 'D'
  aa.table[['GAA']] <- 'E'
  aa.table[['GAG']] <- 'E'
  aa.table[['TGT']] <- 'C'
  aa.table[['TGC']] <- 'C'
  aa.table[['TGG']] <- 'W'
  aa.table[['CGA']] <- 'R'
  aa.table[['CGC']] <- 'R'
  aa.table[['CGG']] <- 'R'
  aa.table[['CGT']] <- 'R'
  aa.table[['AGA']] <- 'R'
  aa.table[['AGG']] <- 'R'
  aa.table[['GGA']] <- 'G'
  aa.table[['GGC']] <- 'G'
  aa.table[['GGG']] <- 'G'
  aa.table[['GGT']] <- 'G'
  if(dim(variants[[2]])[1] > 0){
    start.coding.temp <- as.numeric(variants[[2]][,16])-as.numeric(variants[[2]][,7])+1
    for(i in 1:dim(variants[[2]])[1]){
      temp.bases <- substr(variants[[2]][i,14],start.coding.temp[i],start.coding.temp[i]+peplength*3-1)
      temp.peptide <- vector(mode='character',0)
      #translate to find reference peptide sequence
      for(j in 1:peplength){
        temp.peptide <- c(temp.peptide,aa.table[[substr(temp.bases,(j-1)*3+1,j*3)]])
      }
      variants[[2]][i,18] <- paste(temp.peptide,collapse='')
      #extract mutation position(s)
      temp.bases.split <- strsplit(temp.bases,'')[[1]]
      mut.locations.temp <- setdiff(c(1:(peplength*3)),c(1:(peplength*3))[c(temp.bases.split == strsplit(variants[[2]][i,15],'')[[1]])])
      mut.locations.temp1 <- mut.locations.temp+start.coding.temp[i]-1+as.numeric(variants[[2]][i,9])
      for(j in 1:length(mut.locations.temp)){
        variants[[2]][i,(j-1)*5+19] <- mut.locations.temp1[j]
        variants[[2]][i,(j-1)*5+20] <- mut.locations.temp1[j]
        variants[[2]][i,(j-1)*5+21] <- 'SNV'
        variants[[2]][i,(j-1)*5+22] <- temp.bases.split[mut.locations.temp[j]]
        variants[[2]][i,(j-1)*5+23] <- strsplit(variants[[2]][i,15],'')[[1]][mut.locations.temp[j]]
      }
    }
    ##send silent mutations back to the variants[[1]]
    variants[[1]] <- rbind(variants[[1]],variants[[2]][variants[[2]][,18] == variants[[2]][,1],])
    ##remove silent variants
    variants[[2]] <- variants[[2]][variants[[2]][,18] != variants[[2]][,1],]
  }
  
  #now analyze sequences with gap openings ONLY into list[[3]]
  variants[[3]] <- blastout.temp[as.numeric(blastout.temp[,5]) == 0 & as.numeric(blastout.temp[,6]) > 0,]
  variants[[3]] <- as.matrix(variants[[3]])
  if(dim(variants[[3]])[2] == 1){
    variants[[3]] <- t(variants[[3]])
  }
  
  #not this only works if your max.gaps = 1!!!!
  if(dim(variants[[3]])[1] > 0){
    ##annotate with ins vs del, location of start & end, and reference vs altered sequence
    ##separate into insertions vs deletions based on length of contigs vs reference
    variants[[3]][(nchar(variants[[3]][,14])-(as.numeric(variants[[3]][,8])-as.numeric(variants[[3]][,7])+1)) < 0,21] <- 'INS'
    variants[[3]][(nchar(variants[[3]][,14])-(as.numeric(variants[[3]][,8])-as.numeric(variants[[3]][,7])+1)) > 0,21] <- 'DEL'
    
    ##find gap insertion and return the start & end locations, and reference and altered sequences
    for(i in 1:dim(variants[[3]])[1]){
      if(variants[[3]][i,21] == 'INS'){
        #try adding each base and see if it fixes the alignment
        gap.length <- (as.numeric(variants[[3]][i,8])-as.numeric(variants[[3]][i,7])+1)-nchar(variants[[3]][i,14])
        temp.ref <- variants[[3]][i,14]
        temp.rna <- substr(variants[[3]][i,17],as.numeric(variants[[3]][i,7]),as.numeric(variants[[3]][i,8]))
        for(j in 1:(nchar(temp.ref)-1)){
          temp.compare <- paste(substr(temp.ref,1,j),substr(temp.rna,j,j+gap.length-1),substr(temp.ref,j+gap.length,nchar(temp.ref)),sep='')
          if(temp.compare == temp.rna){
            variants[[3]][i,19] <- j+as.numeric(variants[[3]][i,9])
            variants[[3]][i,20] <- j+as.numeric(variants[[3]][i,9])+gap.length
            variants[[3]][i,22] <- ''
            variants[[3]][i,23] <- substr(variants[[3]][i,17],as.numeric(variants[[3]][i,7])+j-1,as.numeric(variants[[3]][i,7])+j-2+gap.length)
            break
          }
        }
        
        start.coding.temp <- as.numeric(variants[[3]][i,16])-as.numeric(variants[[3]][i,7])+1
        ref.peptide <- vector(mode='character',0)
        #translate to find reference peptide sequence
        for(j in 1:peplength){
          ref.peptide <- c(ref.peptide,aa.table[[substr(temp.ref,(j-1)*3+start.coding.temp,start.coding.temp+j*3-1)]])
        }
        variants[[3]][i,18] <- paste(ref.peptide,collapse='')
      
      } else if(variants[[3]][i,21] == 'DEL'){
        #try deleting each base and see if it fixes the alignment
        gap.length <- nchar(variants[[3]][i,14])-(as.numeric(variants[[3]][i,8])-as.numeric(variants[[3]][i,7])+1)
        temp.ref <- variants[[3]][i,14]
        temp.rna <- substr(variants[[3]][i,17],as.numeric(variants[[3]][i,7]),as.numeric(variants[[3]][i,8]))
        for(j in 1:(nchar(temp.ref)-1)){
          temp.compare <- paste(substr(temp.ref,1,j),substr(temp.ref,j+gap.length+1,nchar(temp.ref)),sep='')
          if(temp.compare == temp.rna){
            variants[[3]][i,19] <- j+as.numeric(variants[[3]][i,9])
            variants[[3]][i,20] <- j+as.numeric(variants[[3]][i,9])+gap.length
            variants[[3]][i,22] <- substr(variants[[3]][i,14],j+1,j+gap.length)
            variants[[3]][i,23] <- ''
            break
          }
        }
        
        start.coding.temp <- as.numeric(variants[[3]][i,16])-as.numeric(variants[[3]][i,7])+1
        ref.peptide <- vector(mode='character',0)
        #translate to find reference peptide sequence
        for(j in 1:peplength){
          ref.peptide <- c(ref.peptide,aa.table[[substr(temp.ref,(j-1)*3+start.coding.temp,start.coding.temp+j*3-1)]])
        }
        variants[[3]][i,18] <- paste(ref.peptide,collapse='')
      }
    }
    
    ##confirm that the reference sequence and peptide coding sequence do not produce the same peptides, if so, send back to variant[[1]]
    variants[[1]] <- rbind(variants[[1]],variants[[3]][variants[[3]][,1] == variants[[3]][,18],])
    variants[[3]] <- variants[[3]][variants[[3]][,1] != variants[[3]][,18],]
  }
  
  #now analyze sequences with gap opening AND mutations into list[[4]]
  #variants[[4]] <- blastout.temp[as.numeric(blastout.temp[,5]) > 0 & as.numeric(blastout.temp[,6]) > 0,]
  
  ##separate into insertions vs deletions based on length of contigs vs reference
  #variants[[4]][(nchar(variants[[4]][,14])-(as.numeric(variants[[4]][,8])-as.numeric(variants[[4]][,7])+1)) < 0,21] <- 'INS'
  #variants[[4]][(nchar(variants[[4]][,14])-(as.numeric(variants[[4]][,8])-as.numeric(variants[[4]][,7])+1)) > 0,21] <- 'DEL'
  
  ##find gap insertion and return the start & end locations, and reference and altered sequences
  ##Currently not annotating variants[[4]]
  
  
  
  #now combine all of your types of variants into a single blastout object and return!
  blastout <- rbind(variants[[1]],variants[[2]],variants[[3]])
  blastout <- blastout[order(blastout[,1],decreasing=FALSE),]
  
  return(blastout)
  
}

framed <- function(){
  #Determines if peptide is in or out of frame. Additionally handles the frameshift mutations from variant caller
  ccds <- fread('references/ccds/CCDS.current.txt')
  ccds <- ccds[ccds[[6]] == 'Public',]
  ccds.exon <- ccds[[10]]
  ccds.exon <- substr(ccds.exon,2,nchar(ccds.exon)-1)
  
}
  





