### MitoSonar: 12S metabarcoding-based Fish Taxonomy Identifier

### Obtaining input
arguments <- commandArgs(trailingOnly = TRUE)

## Filtering parameters
maxN = arguments[1]
truncQ = arguments[2]
truncLen = arguments[3]
trimLeft = arguments[4]
maxEE = arguments[5]

## Defaults for testing
#maxN=0
#truncQ=2
#truncLen=100
#trimLeft=18
#maxEE=2


#Every few months you should uninstall and reinstall the packages you will be using to ensure they are most recent versions.

# remove.packages("BiocManager")
# 
# remove.packages("dada2")
# 
# remove.packages("ShortRead")
# 
# remove.packages("ggplot2")
# 
# remove.packages("phylosec")
# 
# remove.packages("Biostrings")
# 
# remove.packages("dplyr")
# 
# remove.packages("easycsv")
# 
# remove.packages("tidyverse")
# 
# remove.packages("devtools")
# 
# remove.packages("rBLAST")



### Installing Packages

# install.packages("BiocManager", force = TRUE)

# BiocManager::install("ShortRead", force = TRUE)

# BiocManager::install("dada2", force = TRUE)

# BiocManager::install("ggplot2", force = TRUE)

# BiocManager::install("phyloseq", force = TRUE)

# BiocManager::install("Biostrings", force = TRUE)

# BiocManager::install("dplyr", type = "binary", force = TRUE) !!!!!>>>> type 'binary' is not supported on this platform
# install.packages("dplyr")

# BiocManager::install("easycsv", force = TRUE)

# BiocManager::install("tidyverse", force = TRUE)

# install.packages("devtools")

# devtools::install_github("mhahsler/rBLAST", force = TRUE)



### Loading Packages

system(paste("echo Loading R Packages..."))

library(BiocManager); #packageVersion("BiocManager") 

library(dada2); #packageVersion("dada2")

library(ShortRead); #packageVersion("ShortRead") 

library(ggplot2); #packageVersion("ggplot2") 

library(phyloseq); #packageVersion("phyloseq") 

library(Biostrings); #packageVersion("Biostrings") 

library(dplyr); #packageVersion("dplyr") 

library(easycsv)

library(devtools)

library(rBLAST); #packageVersion("rBLAST")

library(tidyverse); #packageVersion("tidyverse")

### Loading Data and Setting Paths and Filenames

path0 = system("pwd", intern = TRUE)
#FOR TESTING:
#path0 = "/home/carloslima/projects/MitoSonar"
path = paste0(path0,"/work-dir/")
setwd(path)

path1 = "data-raw/fastqs/"

blastout <- "data-raw/blast.out" 

blastout2 <- "data-raw/blast.out.unflitered" 

tax_sequences <- "data/taxsequences.fna"

tax_sequences2 <- "data/taxsequences.fna.unfiltered"

vert_fasta <- "data-raw/MiFish_all_mitogenomes.fasta"

fns <- list.files(paste(path, path1, sep = ""))

fastqs <- fns[grepl(".fastq$", fns)]

fastqs <- sort(fastqs) # Sorting fastqs so that Forward (R1) and Reverse (R2) reads correspond to each other.

fnFs <- fastqs[grepl("_R1", fastqs)]

fnRs <- fastqs[grepl("_R2", fastqs)]



### Visualize the quality profile of the forward reads: 

system(paste("echo Examining quality profiles of forward reads..."))

for(fnF in fnFs) {
  fastq_name <- sub("_R1_001.fastq$", "", fnF)
  sample_dir <- paste0(path,"data/images/plots/",fastq_name,"_quality_report")
  if(!dir.exists(sample_dir)){dir.create(sample_dir)}
  png(filename = paste0(sample_dir,"/",fastq_name,"_forward_quality.png"), width = 800, height = 600)
  print(fnF)
  qqF <- qa(paste0(path, path1, fnF))[["perCycle"]]$quality
  print(ShortRead:::.plotCycleQuality(qqF, main="Forward reads quality profile"))
  dev.off()
} 

### Visualize the quality profile of the reverse reads: 

system(paste("echo Examining quality profiles of reverse reads..."))

for(fnR in fnRs) {
  fastq_name <- sub("_R2_001.fastq$", "", fnR)
  sample_dir <- paste0(path,"data/images/plots/",fastq_name,"_quality_report")
  if(!dir.exists(sample_dir)){dir.create(sample_dir)}
  png(filename = paste0(sample_dir,"/",fastq_name,"_reverse_quality.png"), width = 800, height = 600)
  print(fnR)
  qqR <- qa(paste0(path, path1, fnR))[["perCycle"]]$quality
  print(ShortRead:::.plotCycleQuality(qqR, main="Reverse reads quality profile")) 
  dev.off()
} 

system(paste("echo Quality reports saved into 'work-dir/data/images/plots'"))

### Filtering and trimming 

### Filter the forward and reverse reads: 

system(paste("echo Filtering and Trimming..."))

filtFs <- paste0(path, sapply(strsplit(fnFs, "\\."), `[`, 1), "_filt.fastq.gz") 

filtRs <- paste0(path, sapply(strsplit(fnRs, "\\."), `[`, 1), "_filt.fastq.gz") 

for(i in seq_along(fnFs)) { #Adjust parameters according to quality profiles
  fastqPairedFilter(paste0(path, path1, c(fnFs[i], fnRs[i])), c(filtFs[i], filtRs[i]), maxN=maxN, maxEE=maxEE, truncQ=as.numeric(truncQ), trimLeft=as.numeric(trimLeft), truncLen=c(truncLen,truncLen), compress=TRUE, verbose=TRUE)
}

### Visualize the quality profile of the forward reads after filtering: 

system(paste("echo Examining quality profiles after filtering..."))

for(fnF in filtFs) {
  fastq_name <- sub(paste0("^", "/home/carloslima/projects/MitoSonar/work-dir/"), "", fnF)
  fastq_name <- sub("_R1_001_filt.fastq.gz$", "", fastq_name)
  sample_dir <- paste0(path,"data/images/plots/",fastq_name,"_quality_report")
  if(!dir.exists(sample_dir)){dir.create(sample_dir)}
  png(filename = paste0(sample_dir,"/",fastq_name,"_forward_filtered_quality.png"), width = 800, height = 600)
  print(fnF)
  qqF <- qa(paste0(fnF))[["perCycle"]]$quality
  print(ShortRead:::.plotCycleQuality(qqF, main="Forward reads quality profile after filtering"))
  dev.off()
} 

### Visualize the quality profile of the reverse reads after filtering: 

for(fnR in filtRs) {
  fastq_name <- sub(paste0("^", "/home/carloslima/projects/MitoSonar/work-dir/"), "", fnR)
  fastq_name <- sub("_R2_001_filt.fastq.gz$", "", fastq_name)
  sample_dir <- paste0(path,"data/images/plots/",fastq_name,"_quality_report")
  if(!dir.exists(sample_dir)){dir.create(sample_dir)}
  png(filename = paste0(sample_dir,"/",fastq_name,"_reverse_filtered_quality.png"), width = 800, height = 600)
  print(fnR)
  qqR <- qa(paste0(fnR))[["perCycle"]]$quality
  print(ShortRead:::.plotCycleQuality(qqR, main="Reverse reads quality profile after filtering"))
  dev.off()
} 

system(paste("echo Quality reports have been sent to 'work-dir/data/images/plots'"))



### Dereplication 

### Dereplicate the filtered fastq files:

system(paste("echo Derreplicating FASTQs..."))

derepFs <- lapply(filtFs, derepFastq, verbose=TRUE) 
derepRs <- lapply(filtRs, derepFastq, verbose=TRUE)

### Name the derep-class objects by the sample names

sam_names <- sapply(strsplit(fnFs, "/"), tail, n=1)
sam_names <- sapply(strsplit(sam_names, "_"), `[`, 1)

names(derepFs) <- sam_names
names(derepRs) <- sam_names

### Sample Inference by DADA2

system(paste("echo Initiating inference phase..."))

dadainfer <- function(derepFastqs){
  dada(derepFastqs, err=inflateErr(tperr1,3), errorEstimationFunction=loessErrfun, selfConsist = TRUE)
}

dadaFs <- dadainfer(derepFs)
dadaRs <- dadainfer(derepRs)



### Visualize estimated error rates:

system(paste("echo Drawing estimated error rates plots..."))

estimErrPng <- function(dadaObj,samId,sample_dir,acgtBase) {
  png(filename = paste0(sample_dir,"/",samId,"_estimErr_",acgtBase,".png"), width = 800, height = 600)
  plot <- plotErrors(dadaObj, acgtBase, nominalQ=TRUE)
  print(plt)
  dev.off()
}

if (length(dadaFs) > 1) {
  i=1
  while (i <= length(sam_names)) {
    sample_dir <- paste0(path,"data/images/plots/",sam_names[i],"_estimated_errors")
    if(!dir.exists(sample_dir)){dir.create(sample_dir)}
    estimErrPng(dadaFs[[i]],sam_names[i],sample_dir,"A")
    estimErrPng(dadaFs[[i]],sam_names[i],sample_dir,"C")
    estimErrPng(dadaFs[[i]],sam_names[i],sample_dir,"G")
    estimErrPng(dadaFs[[i]],sam_names[i],sample_dir,"T")
    i=i+1
  }
}else {
  sample_dir <- paste0(path,"data/images/plots/",sam_names,"_estimated_errors")
  if(!dir.exists(sample_dir)){dir.create(sample_dir)}
  estimErrPng(dadaFs,sam_names,sample_dir,"A")
  estimErrPng(dadaFs,sam_names,sample_dir,"C")
  estimErrPng(dadaFs,sam_names,sample_dir,"G")
  estimErrPng(dadaFs,sam_names,sample_dir,"T")
}

system(paste("echo Plots saved into 'work-dir/data/images/plots'"))



### Identify chimeric sequences:

system(paste("echo Removing chimeric sequences..."))

bimFs <- sapply(dadaFs, isBimeraDenovo, verbose=TRUE)
bimRs <- sapply(dadaRs, isBimeraDenovo, verbose=TRUE)

print(unname(sapply(bimFs, mean)), digits=2)
print(unname(sapply(bimRs, mean)), digits=2)



### Merge paired reads

mergers <- mapply(mergePairs, dadaFs, derepFs, dadaRs, derepRs, SIMPLIFY=FALSE)



### Remove chimeras 

mergers.nochim <- mapply(function(mm, bF, bR) mm[!bF[mm$forward] & !bR[mm$reverse],], mergers, bimFs, bimRs, SIMPLIFY=FALSE) 



### Constructing the sequence table

system(paste("echo Building tables..."))

seqtab <- makeSequenceTable(mergers.nochim)
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(80,120)]



### Make otu_table

seqs <- colnames(seqtab2)
otab <- otu_table(seqtab2, taxa_are_rows=FALSE)
colnames(otab) <- paste0("Seq_", seq(ncol(otab)))



### BLAST+
### Get Taxonomy Based on Blasting against a known reference 

system(paste("echo Preparing database for taxonomy annotation..."))

writeFasta <- function(seqs, output) { 
  seqsout <- mapply( function(idx, sequence) paste0(">Seq_",idx,"\n",sequence,"\n"), seq(length(seqs)), seqs)
  write(paste0(seqsout), file = output, sep = "")
} 

seqs_for_blast <- DNAStringSet(seqs)

names(seqs_for_blast) <- sapply(seq(length(seqs)),function(x) {paste0("Seq_",x)})

writeFasta(seqs, paste(path, tax_sequences, sep = ""))

transotab = t(otab)

write.table(transotab, paste0(path, "data/otutable.txt"))

write.csv(transotab, paste0(path, "data/otutable.csv"))


### Making BLAST database

system(paste("/home/carloslima/tools/ncbi-blast-2.15.0+/bin//makeblastdb -dbtype nucl -in", paste0(path,vert_fasta)))

system(paste("echo Blasting..."))

system(paste("/home/carloslima/tools/ncbi-blast-2.15.0+/bin//blastn -query", paste0(path, tax_sequences),  "-db", paste0(path, vert_fasta), "-outfmt '6 qseq qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore' -max_target_seqs 1 -out", paste0(path, blastout)))
# Additional column headers can be found on https://www.metagenomics.wiki/tools/blast/blastn-output-format-6.
# Can rename columns but must be in same order as above.

blasttable <- read.table(paste0(path, blastout))

colnames(blasttable) <- c("qseq", "qid", "sid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

paste("/home/carloslima/tools/ncbi-blast-2.15.0+/bin//blastn -query", paste0(path, tax_sequences),  "-db", paste0(path, vert_fasta), "-outfmt 6 -max_target_seqs 1 -out", paste0(path, blastout))


blastresults <- blasttable %>%
  filter(pident == 100) %>%
  select(qid, sid) %>%
  data.frame()


taxtab <- data.frame(seqs = colnames(otab), stringsAsFactors = FALSE) %>%
  left_join(blastresults, by=c('seqs'='qid'))


taxtab$taxa <- taxtab$sid
rownames(taxtab) <- taxtab$seqs

taxtab <- taxtab %>%
  select(-seqs, -sid) %>%
  as.matrix() %>%
  tax_table


ps <- phyloseq(otab, taxtab)



### Plotting top 10 taxa by abundance

# system(paste("echo Done! Now plotting top 10 most abundant taxa..."))
# 
# top10 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:10]
# 
# ps.top10 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
# 
# ps.top10 <- prune_taxa(top10, ps.top10)
# 
# png(filename = paste0(path,"data/images/plots/top10taxa.png"), width = 800, height = 600)
# plot_bar(ps.top10, fill="taxa")
# dev.off()



### Plotting most abundant taxonomies by sample

system(paste("echo Done! Now plotting most abundant taxonomies by sample..."))
toptaxa <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:5]
ps.toptaxa <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.toptaxa <- prune_taxa(toptaxa, ps.toptaxa)

for (sample_name in row.names(otu_table(ps.toptaxa))){
  ps.toptaxa_sample <- prune_samples(sample_name, ps.toptaxa)
  png(filename = paste0(path, "data/images/plots/", sample_name, "_toptaxa.png"), width = 800, height = 600)
  plot <- plot_bar(ps.toptaxa_sample, x = "taxa", y = "Abundance", fill = "taxa", title = paste0("Most abundant taxa for sample ",sample_name))
  print(plot)
  dev.off()
}



### Saving taxonomy table into a file

system(paste("echo Saving taxonomy table and blast output..."))

tax_table1<-"data/taxtable.txt"
write.table(taxtab, file=paste0(path,tax_table1))

tax_tablecsv<-"data/taxtable.csv"
write.csv(taxtab, file=paste0(path,tax_tablecsv)) # just so it can be easily read by excel

# same thing for blasttable
blasttable <- read.table(paste0(path, blastout))

colnames(blasttable) <- c("qseq", "qid", "sid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

blast_out<-"data-raw/blastout.csv"

write.csv(blasttable,file=paste0(path,blast_out))

system(paste("echo ##### Analysis Complete #####"))



### Calling Rmd for generating reports

create_report <- function(sample,top,sec){
  rmarkdown::render(input = paste0(path0,"/MitoSonar_report.Rmd"),
                    output_format = rmarkdown::pdf_document(),
                    output_file = paste0(sample,"_MSReport"),
                    output_dir = paste0(path,"reports"),
                    intermediates_dir = paste0(path,"reports"),
                    clean = TRUE,
                    params = list(sample = sample,
                                  maxN = maxN,
                                  truncQ = truncQ,
                                  truncLen = truncLen,
                                  trimLeft = trimLeft,
                                  maxEE = maxEE,
                                  species = top,
                                  species2 = sec))
}

rep_dir <- paste0(path,"reports")
if(!dir.exists(rep_dir)){dir.create(rep_dir)}

for (sample_name in sam_names) {
  ps.toptaxa_sample <- prune_samples(sample_name, ps.toptaxa)
  otu_data <- as.data.frame(otu_table(ps.toptaxa_sample))
  topseq <- names(otu_data)[which.max(otu_data)]
  toptaxa <- tax_table(ps.toptaxa_sample)[topseq]
  sectaxas <- c()
  for (seq in (setdiff(names(otu_data), topseq))) {
    if (otu_data[seq] >= 0.1) {
      secseq <- seq
      sectaxas <- append(sectaxas,tax_table(ps.toptaxa_sample)[secseq])
    }
  }
  #print(sample_name)
  #print(paste0("Top taxa: ",toptaxa))
  #print(paste0("Other relevant taxas: ",paste(sectaxas, collapse = ", ")))
  #print("-----")
  create_report(sample_name,toptaxa,sectaxas)
}
system(paste0("rm -r ",path,"reports/work-dir"))

