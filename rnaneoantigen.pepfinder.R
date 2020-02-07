################################################ rnaneoantigen.pepfinder.R ########################################################################

#Initialize
library(data.table)

#Load in neoantigen peptides built from exomeSeq from TCIA.at
tcia.loader = function(file_location,patients,pep_lengths){
  setwd(paste(file_location,'corrollary_files/',sep=''))
  tcia = fread('TCIA-NeoantigensData_luad.txt') #load in TCIA data
  tcia = tcia[tcia[[1]] %in% patients,]  #only look for patients you care about
  tcia = tcia[,c(1,3:4)] #clean up columns
  tcia = tcia[nchar(tcia[[3]]) %in% pep_lengths,] #subset to peptide lengths that you measured
  return(tcia)
}

#Load in TCGA maf file
tcga.loader = function(file_location,patients){
  setwd(paste(file_location,'corrollary_files/LUAD.maf',sep='')) #change to directory containing tcga maf files
  maf = lapply(dir(),fread) #load in all files in this directory
  names(maf) = substr(dir(),1,12) #label the maf by patient ID
  maf = maf[names(maf) %in% patients] #subset to only patients that we're considering
  return(maf)
}

translator = function(codon){
  #translated codons to amino acids
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
  return(aa.table[[codon]])
}

antisenser = function(sequences){
  antisense = list()
  antisense[['A']] = 'T'
  antisense[['C']] = 'G'
  antisense[['G']] = 'C'
  antisense[['T']] = 'A'
  for(i in 1:length(sequences)){
    substr(sequences,i,i) = antisense[[substr(sequences,i,i)]]
  }
  sequences = paste(rev(strsplit(sequences,'')[[1]]),collapse='')
  return(sequences)
}

seq.to.pep = function(sequences,pep_length,frame){
  #creates all possible peptides of length pep_length from a sequence that contain the middle base, note input sequence must be of odd length
  peptides = vector(mode='character',pep_length) #output variable 
  for(i in 1:pep_length){
    peptides.temp = ''
    for(j in 1:pep_length){
      peptides.temp = paste(peptides.temp,translator(substr(sequences,(i-1)*3+frame+(j-1)*3+1,(i-1)*3+frame+(j-1)*3+3)),sep='') #walk along sequence and translate each codon into peptides of length pep_length
    }
    peptides[i] = peptides.temp
  }
  return(peptides)
}

#convert mutations into neoantigenic peptides of a given length
tcga.to.pep = function(file_location,patients,maf,pep_lengths){
  peps = list()
  for(i in 1:length(maf)){
    print(paste('analyzing patient',patients[i],sep=' '))
    setwd(file_location)
    bed = maf[[i]][,c(5,6,7,11,9,8,13,37)] #first make the bed file to get sequence from reference for each patient
    maf.out.temp1 = maf[[i]][bed[[5]] == 'Missense_Mutation' & nchar(bed[[4]]) == 1,]
    bed = bed[bed[[5]] == 'Missense_Mutation' & nchar(bed[[4]]) == 1,] #currently only works for singlet missense mutations due to the difficulty of creating peptides from other alterations
    bed[[1]] = paste('chr',bed[[1]],sep='')
    peps.temp = list()
    maf.out = list()
    for(j in 1:length(pep_lengths)){
      bed.write = bed #create a dummy variable for this loop
      bed.write[[2]] = as.numeric(bed.write[[2]])-1-(pep_lengths[j]*3-1) #open your bed window based on the peptide lengths
      bed.write[[3]] = as.numeric(bed.write[[3]])+(pep_lengths[j]*3-1) #note that the base we care about is on the end location
      write.table(bed.write,file = paste('bedfiles_maf/',names(maf)[i],pep_lengths[j],sep=''),col.names = FALSE,row.names = FALSE,sep = '\t',quote = FALSE) #write your bed file
      setwd('references/') #change wd to the location of your bedtools executable, here I've put it in the same folder as the reference genome
      fasta = system(paste('./bedtools getfasta -fi ',file_location,'references/hg19.fa -bed ',file_location,'bedfiles_maf/',names(maf)[i],pep_lengths[j],sep=''),intern=TRUE) #get your fasta from bed
      setwd('..')
      bed.write = bed[bed[[4]] == substr(fasta[seq(2,length(fasta),2)],pep_lengths[j]*3,pep_lengths[j]*3),] #test that the reference base in the maf file matches the reference genome and remove those that don't
      maf.out.temp2 = maf.out.temp1[bed[[4]] == substr(fasta[seq(2,length(fasta),2)],pep_lengths[j]*3,pep_lengths[j]*3),]
      fasta = fasta[seq(2,length(fasta),2)][bed[[4]] == substr(fasta[seq(2,length(fasta),2)],pep_lengths[j]*3,pep_lengths[j]*3)]  #update fasta file
      substr(fasta,pep_lengths[j]*3,pep_lengths[j]*3) = bed.write[[7]] #mutate the reference file
      for(l in 1:length(fasta)){ #invert any sequences with antisense coding, also puts it in reverse orientation
        if(bed.write[[8]][l] == '-'){fasta[l] = antisenser(toupper(fasta[l]))}
      }
      frame.temp = maf.out.temp2[[41]]  #extract the coding frame for each mutation
      frame.temp = lapply(strsplit(frame.temp,')'),strsplit,">")
      frame1 = vector(mode='numeric',length(frame.temp))
      for(count in 1:length(frame.temp)){
        temp = strsplit(frame.temp[[count]][[2]],'')
        frame1[count] = which(temp[[1]] != temp[[2]])
      }
      frame = frame1
      frame[frame1==1] = 2
      frame[frame1==2] = 1
      frame[frame1==3] = 0
      peps.temp[[j]] = mapply(seq.to.pep,fasta,pep_lengths[j],frame)  #translate fasta sequences into peptides
      maf.out[[j]] = maf.out.temp2
    }
    peps[[i]] = list(peps.temp,maf.out[[1]])
  }
  names(peps) = names(maf)
  return(peps)
}

#the function that loads peptides of choice
peptide.loader = function(peps_patient,patient,your_file_location,pep_length){
  filesize = as.integer(system(paste('wc -l < ',your_file_location,sep = ''),intern = TRUE)[[1]]) #have to load the peptide file in pieces of size, lines_per
  lines_per = 10000000
  imax = ceiling(filesize/lines_per)
  peps.found = matrix(NA,nrow = 0,ncol = 2) #initialize matrix storing the peptides that you found in your file
  peps.expr = 0
  for(i in 1:imax){
    #print(i)
    if(i != imax){
      peps.temp = fread(your_file_location,nrows = lines_per,skip = (i-1)*lines_per,header=FALSE)
      peps.found = rbind(peps.found,peps.temp[peps.temp[[2]] %in% peps_patient,])
      peps.expr = peps.expr+sum(as.integer(peps.temp[[1]]))
    } else {
      peps.temp = fread(your_file_location,nrows=filesize-(lines_per*(i-1)),skip=(i-1)*lines_per,header=FALSE)
      peps.found = rbind(peps.found,peps.temp[peps.temp[[2]] %in% peps_patient,])
      peps.expr = peps.expr+sum(as.integer(peps.temp[[1]]))
    }
  }
  return(list(peps.found,peps.expr))
}


#Intersect a list of peptides in your generated peptide file
peptide.finder = function(peps,patients,your_file_dir,pep_lengths){  #load peptide data from your files
  peps.found.tumor.kmer = list()
  peps.found.normal.kmer = list()
  for(i in 1:length(pep_lengths)){
    print(paste('analysing peptide length ',pep_lengths[i]))
    peps.found.tumor <- list()
    peps.found.normal <- list()
    for(j in 1:2){
      
      if(j == 1){type = 'normal'} else {type = 'tumor'}
      print(paste('analysing ',type,' patients'))
      setwd(paste(your_file_dir,'TCGA_neoantigens_',pep_lengths[i],'mer/',type,sep='')) #change directory to the location of the peptide files
      peps.found.temp = list()
      for(l in 1:length(dir())){
        patients_pep = dir()
        print(paste('analysing patient ',l,' of ',length(patients_pep)))
        peps_patient = as.vector(peps[which(patients == substr(patients_pep[l],1,12))][[1]][[1]][[i]]) #peptides for the patient you're looking for
        peps.found.temp[[l]] = list(peps_patient,peptide.loader(peps_patient,patient,paste(your_file_dir,'TCGA_neoantigens_',pep_lengths[i],'mer/',type,'/',dir()[l],sep=''),pep_lengths[i]))
      }
      if(j == 1){
        peps.found.normal = peps.found.temp
      } else {
        peps.found.tumor = peps.found.temp
      }
    }
    peps.found.tumor.kmer[[i]] = peps.found.tumor
    peps.found.normal.kmer[[i]] = peps.found.normal
  }
  return(list(peps.found.tumor.kmer,peps.found.normal.kmer))
}

#protein blast your peptides to see if 'mutated' peptides are present in non-mutated proteins
pep.blast <- function(peptides,file_location,file){
  setwd(paste(file_location,'references/nr.142',sep=''))
  filesize = as.integer(system(paste('wc -l < ',file_location,'peps.tumornormal/',file,sep = ''),intern = TRUE)[[1]]) #have to load the peptide file in pieces of size, lines_per
  if(filesize > 0){
    pblastout <- system(paste('./blastp -db nr.142 -task blastp-short -gilist /Users/Michael/Desktop/Hsapiens_gilist.gi -comp_based_stats 0 -gapopen 9 -gapextend 1 -window_size 40 -threshold 11 -matrix PAM30 -evalue 200000 -word_size 2 -query ',file_location,'peps.tumornormal/',file,' -outfmt 6',sep=''),intern=TRUE)
    pblastout <- system(paste('./blastp -db nr -remote -entrez_query "Homo sapiens [Organism]"  -query ',file_location,'peps.tumornormal/',file,' -outfmt 6',sep=''),intern=TRUE)
  }
  pblastout.temp <- strsplit(pblastout,'\t')
  pblast.mat <- matrix(NA,nrow=length(pblastout),ncol=12)
  for(i in 1:length(pblastout.temp)){
    pblast.mat[i,] <- pblastout.temp[[i]]
  }
  pblast.mat.temp <- pblast.mat[as.numeric(pblast.mat[,3]) == 100 & as.numeric(pblast.mat[,4]) == 9 & as.numeric(pblast.mat[,5]) == 0 & as.numeric(pblast.mat[,6]) == 0,]
  
}

#initialize variables
file_location = '/Volumes/Michael-backup3/rnaneoantigen/'
setwd(file_location) #set location to directory containing your corrolary files
patients = unique(fread('corrollary_files/LUAD_uuids.txt',header=FALSE)[[1]]) #load in the list of TCGA patients you're looking at
pep_lengths = 8:9 #what peptide lengths are you looking for?
your_file_location = '/Volumes/Michael-backup3/rnaneoantigen/'  #location of your peptide files


#run your code 
maf = tcga.loader(file_location,patients)
peps = tcga.to.pep(file_location,patients,maf,pep_lengths)
tcia = tcia.loader(file_location,patients,pep_lengths)
peps.found = peptide.finder(peps,patients,"/Volumes/Michael-backup3/rnaneoantigen/",pep_lengths)











