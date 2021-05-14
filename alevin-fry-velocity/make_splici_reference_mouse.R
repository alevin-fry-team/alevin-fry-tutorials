# load the packages we will use
library(eisaR)
library(Biostrings)
library(BSgenome)
library(stringr)

gtf <- file.path("genes/genes.gtf")
# fl is the flank length, here we set it to 
# the read length - 5 (151 - 5bp = 146bp) 
grl <- getFeatureRanges(
  gtf = gtf,
  featureType = c("spliced", "intron"), 
  intronType = "separate", 
  flankLength = 146L, 
  joinOverlappingIntrons = FALSE, 
  verbose = TRUE
)

# load the genome sequence
x <- Biostrings::readDNAStringSet(file.path("10X-mouse-2.1.0/mm10/fasta/genome.fa"))
# fix the names
names(x) <- sapply(strsplit(names(x), " "), .subset, 1)

seqlevels(grl) <- seqlevels(x)
seqlengths(grl) <- seqlengths(x)

seqs <- GenomicFeatures::extractTranscriptSeqs(
  x = x,
  transcripts = trim(grl)
)

df <- getTx2Gene(grl)
write.table(df, file.path("t2g.tsv"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
df[, "status"] = sapply(strsplit(df$transcript_id, "-"), function(x) if(length(x) == 2){"U"} else {"S"})
df[, "gene_id"] = sapply(strsplit(df$gene_id, "-"), function(x) x[1])

writeXStringSet(seqs, file.path("transcriptome_splici.fa"), format = "fasta")
write.table(df, file.path("t2g_3col.tsv"), sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
