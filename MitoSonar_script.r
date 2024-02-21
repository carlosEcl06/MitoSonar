### MitoSonar: 12S metabarcoding-based Fish Taxonomy Identifier



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

library(BiocManager); packageVersion("BiocManager") 

library(dada2); packageVersion("dada2")

library(ShortRead); packageVersion("ShortRead") 

library(ggplot2); packageVersion("ggplot2") 

library(phyloseq); packageVersion("phyloseq") 

library(Biostrings); packageVersion("Biostrings") 

library(dplyr); packageVersion("dplyr") 

library(easycsv)

library(devtools)

library(rBLAST); packageVersion("rBLAST")

library(tidyverse); packageVersion("tidyverse")

### Loading Data and Setting Paths and Filenames

path = "/home/carloslima/projects/MitoSonar/work-dir/"
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

for(fnF in fnFs[1:1]) {
  print(fnF)
  qqF <- qa(paste0(path, path1, fnF))[["perCycle"]]$quality
  print(ShortRead:::.plotCycleQuality(qqF, main="Forward")) 
} 



### Visualize the quality profile of the reverse reads: 

for(fnR in fnRs[1:1]) {
  qqR <- qa(paste0(path, path1, fnR))[["perCycle"]]$quality
  print(ShortRead:::.plotCycleQuality(qqR, main="Reverse")) 
} 



### Filtering and trimming 

### Filter the forward and reverse reads: 

filtFs <- paste0(path, sapply(strsplit(fnFs, "\\."), `[`, 1), "_filt.fastq.gz") 

filtRs <- paste0(path, sapply(strsplit(fnRs, "\\."), `[`, 1), "_filt.fastq.gz") 

for(i in seq_along(fnFs)) { #Adjust parameters according to quality profiles
  fastqPairedFilter(paste0(path, path1, c(fnFs[i], fnRs[i])), c(filtFs[i], filtRs[i]), maxN=0, maxEE=2, truncQ=2, trimLeft=c(18, 18), truncLen=c(100,100), compress=TRUE, verbose=TRUE)
}



### Dereplication 

### Dereplicate the filtered fastq files:

derepFs <- lapply(filtFs, derepFastq, verbose=TRUE) 
derepRs <- lapply(filtRs, derepFastq, verbose=TRUE)

### Name the derep-class objects by the sample names

sam_names <- sapply(strsplit(fnFs, "/"), tail, n=1)
sam_names <- sapply(strsplit(sam_names, "_"), `[`, 1)

names(derepFs) <- sam_names
names(derepRs) <- sam_names



### Sample Inference by DADA2

dadainfer <- function(derepFastqs){
  dada(derepFastqs, err=inflateErr(tperr1,3), errorEstimationFunction=loessErrfun, selfConsist = TRUE)
}

dadaFs <- dadainfer(derepFs)
dadaRs <- dadainfer(derepRs)



### Visualize estimated error rates:

plotErrors(dadaFs[[1]], "A", nominalQ=TRUE)
plotErrors(dadaFs[[1]], "C", nominalQ=TRUE)
plotErrors(dadaFs[[1]], "G", nominalQ=TRUE)
plotErrors(dadaFs[[1]], "T", nominalQ=TRUE)



### Identify chimeric sequences:

bimFs <- sapply(dadaFs, isBimeraDenovo, verbose=TRUE)
bimRs <- sapply(dadaRs, isBimeraDenovo, verbose=TRUE)

print(unname(sapply(bimFs, mean)), digits=2)
print(unname(sapply(bimRs, mean)), digits=2)



### Merge paired reads

mergers <- mapply(mergePairs, dadaFs, derepFs, dadaRs, derepRs, SIMPLIFY=FALSE)



### Remove chimeras 

mergers.nochim <- mapply(function(mm, bF, bR) mm[!bF[mm$forward] & !bR[mm$reverse],], mergers, bimFs, bimRs, SIMPLIFY=FALSE) 



### Constructing the sequence table

seqtab <- makeSequenceTable(mergers.nochim)
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(80,120)]



### Make otu_table

seqs <- colnames(seqtab2)
otab <- otu_table(seqtab2, taxa_are_rows=FALSE)
colnames(otab) <- paste0("Seq_", seq(ncol(otab)))



### BLAST+
### Get Taxonomy Based on Blasting against a known reference 

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

top10 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:10]

ps.top10 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))

ps.top10 <- prune_taxa(top10, ps.top10)

plot_bar(ps.top10, fill="taxa")



### Saving taxonomy table into a file

tax_table1<-"data/taxtable.txt"
write.table(taxtab, file=paste0(path,tax_table1))

tax_tablecsv<-"data/taxtable.csv"
write.csv(taxtab, file=paste0(path,tax_tablecsv)) # just so it can be easily read by excel

# same thing for blasttable
blasttable <- read.table(paste0(path, blastout))

colnames(blasttable) <- c("qseq", "qid", "sid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

blast_out<-"data-raw/blastout.csv"

write.csv(blasttable,file=paste0(path,blast_out))

