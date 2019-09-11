# This file is generated from the corresponding .Rmd file



options(width=90)


# Install BiocManager from CRAN package repository
install.packages("BiocManager")

# BiocManager can then install packages from Bioconductor package repository
BiocManager::install("Biostrings")
#or
library(BiocManager)
install("Biostrings")


vignette()
vignette(package="Biostrings")
vignette("BiostringsQuickOverview", package="Biostrings")


# __________________________________________
# ==== DNA sequences and genomic ranges ====

library(Biostrings)     # Provides DNAString, DNAStringSet, etc
library(BSgenome)       # Provides getSeq()
library(GenomicRanges)  # Provides GRanges, etc
library(rtracklayer)    # Provides import() and export()


## ____________________
## ----> DNAString ----

myseq <- DNAString("ACCATTGATTAT")
myseq

class(myseq)

reverseComplement(myseq)
translate(myseq)

subseq(myseq, 3,5)
myseq[3:5]

as.character(myseq)


methods(class="DNAString")


?"DNAString-class"


## _______________________
## ----> DNAStringSet ----

myset <- DNAStringSet( list(chrI=myseq, chrII=DNAString("ACGTACGT")) )
myset

# A DNAStringSet is list-like
myset$chrII
# or myset[["chrII"]]
# or myset[[2]]


## __________________
## ----> GRanges ----

range1 <- GRanges("chrI", IRanges(start=3,end=5), strand="+")
range1
getSeq(myset, range1)

range2 <- GRanges("chrI", IRanges(start=3,end=5), strand="-")
getSeq(myset, range2)


seqnames(range1)
start(range1)
end(range1)
strand(range1)
as.data.frame(range1)


# GRanges are sometimes like vectors:
c(range1, range2)

# GRanges can have metadata columns, so they are also like data frames:
mcols(range1)$wobble <- 10
range1
mcols(range1)
range1$wobble

# A handy way to create a GRanges
as("chrI:3-5:+", "GRanges")


## ___________________
## ----> Question ----
# 
# Based on what we saw for `DNAString`, how can we learn more about
# using `GRanges` and `IRanges` objects?
# 
# 
#
## ____________________
## ----> Challenge ----
# 
# Reverse complement the following DNA sequence and then translate to an
# amino acid sequence:
# 

TTCCATTTCCAT

# 
# 
#
# _______________________
# ==== Loading files ====

## ____________________________
## ----> Loading sequences ----

### The start of the .fa file looks like this:
# >Chromosome dna:chromosome chromosome:GCA_000800765.1:Chromosome:1:4558660:1
# AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTC
# TGATAGCAGCTTCTGAACTGGTTACCTGCCGTGAGTAAATTAAAATTTTATTGACTTAGG
# TCACTAAATACTTTAACCAATATAGGCATAGCGCACAGACAGATAAAAATTACAGAGTAC
# ACAACATCCATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGT
# AACGGTGCGGGCTGACGCGTACAGGAAACACAGAAAAAAGCCCGCACCTGACAGTGCGGG
# CTTTTTTTTTCGACCAAAGGTAACGAGGTAACAACCATGCGAGTGTTGAAGTTCGGCGGT
# ...


seqs <- readDNAStringSet("r-bioc-files/Escherichia_coli_k_12.GCA_000800765.1.29.dna.genome.fa")
seqs

# Our chromosome name is too verbose.
# Remove everything from the name after the first space.
names(seqs)
names(seqs) <- sub(" .*","",names(seqs))
names(seqs)


## ___________________________
## ----> Loading features ----

### The start of the .gtf file looks like this:
# #!genome-build ASM80076v1
# #!genome-version GCA_000800765.1
# #!genome-date 2014-12
# #!genome-build-accession GCA_000800765.1
# #!genebuild-last-updated 2014-12
# Chromosome      ena     gene    190     255     .       +       .       gene_id "ER3413_4519"; gene_version "1"; gene_name "thrL"; gene_source "ena"; gene_biotype "protein_coding";
# Chromosome      ena     transcript      190     255     .       +       .       gene_id "ER3413_4519"; gene_version "1"; transcript_id "AIZ54182"; transcript_version "1"; gene_name "thrL"; gene_source "ena"; gene_biotype "protein_coding"; transcript_name "thrL-1"; transcript_source "ena"; transcript_biotype "protein_coding";
# Chromosome      ena     exon    190     255     .       +       .       gene_id "ER3413_4519"; gene_version "1"; transcript_id "AIZ54182"; transcript_version "1"; exon_number "1"; gene_name "thrL"; gene_source "ena"; gene_biotype "protein_coding"; transcript_name "thrL-1"; transcript_source "ena"; transcript_biotype "protein_coding"; exon_id "AIZ54182-1"; exon_version "1";
# Chromosome      ena     CDS     190     252     .       +       0       gene_id "ER3413_4519"; gene_version "1"; transcript_id "AIZ54182"; transcript_version "1"; exon_number "1"; gene_name "thrL"; gene_source "ena"; gene_biotype "protein_coding"; transcript_name "thrL-1"; transcript_source "ena"; transcript_biotype "protein_coding"; protein_id "AIZ54182"; protein_version "1";
# ...


features <- import("r-bioc-files/Escherichia_coli_k_12.GCA_000800765.1.29.gtf")

# Optional: just retain the columns of metadata we need
mcols(features) <- mcols(features)[,c("type","gene_name","gene_id")]

features


feat <- features[4,]
feat
feat_seq <- getSeq(seqs, feat)
feat_seq
translate(feat_seq)


subset(features, gene_name == "lacA")
# Equivalently:
#   features[features$gene_name == "lacA" & !is.na(features$gene_name),]


cds <- subset(features, type == "CDS")
cds
# Equivalently:
#   features[features$type == "CDS",]


# _______________________________________
# ==== Further operations on GRanges ====

## ______________________
## ----> Intra-range ----

?"intra-range-methods"


feat <- features[4,]
feat_stop <- resize(feat, width(feat)+3)
seq_stop <- getSeq(seqs, feat_stop)
translate(seq_stop)


# input
#                         (----)
#                         .    .
# resize                  (--------)
# resize(fix="end")   (--------)
#                         .    .
# flank              (---).    .
# flank(start=F)          .    .(---)
#                         .    .
# promoters          (------)  .
#                         .    .
# narrow                  .(--).
#                         .    .
# shift (ignores strand!) .  (----)


## ______________________
## ----> Inter-range ----

?"inter-range-methods"
?"setops-methods"


query <- as("Chromosome:9500-10000:+", "GRanges")
hits <- findOverlaps(query, features, ignore.strand=TRUE)
hits
subjectHits(hits)
features[subjectHits(hits),]

findOverlaps(query, features, ignore.strand=FALSE)


# input
#         (--------)
#              (----------------)
#                      (-----)
#                                     (---)
#
# range   (-------------------------------)
#
# reduce  (----------------------)    (---)
#
# disjoin (---)(---)(-)(-----)(--)    (---)
#
# GenomicRanges::setdiff(range(input),input)
#                                 (--)


## ____________________
## ----> Challenge ----
# 
# What are E. coli's most common start and stop codons?
# 
# The start codon is the first three bases of the CDS, and the stop
# codon is the three bases following the end of the CDS.
# 
# Hint: Recall that we could get all CDS ranges with:
# 

cds <- subset(features, type == "CDS")

# 
# Hint: Use `flank()` and `resize()` to manipulate these ranges.
# 
# 
#
# _______________________________________
# ==== Further data types to explore ====

# _______________________________
# ==== Finding a known motif ====

vmatchPattern("AGGAGGT", seqs)


vmatchPattern(reverseComplement(DNAString("AGGAGGT")), seqs)


query <- DNAString("AGGAGGT")
max.mismatch <- 1

fwd <- vmatchPattern(query, seqs, max.mismatch=max.mismatch)
fwd <- as(fwd, "GRanges")
strand(fwd) <- "+"
rev <- vmatchPattern(reverseComplement(query), seqs, max.mismatch=max.mismatch)
rev <- as(rev, "GRanges")
strand(rev) <- "-"

complete <- c(fwd, rev)
complete

# Write to GFF file
export(complete, "motif-matches.gff")


# _______________________________
# ==== De novo motif finding ====

# Note: bacteria do not have introns
# In a eukaryote, you would need to merge CDS by transcript

size <- 20

initiation_regions <- flank(cds, size, start=TRUE)
initiation_seqs <- getSeq(seqs, initiation_regions)
names(initiation_seqs) <- initiation_regions$gene_id

# Look for any composition bias
library(seqLogo)
letter_counts <- consensusMatrix(initiation_seqs)
probs <- prop.table(letter_counts[1:4,], 2)
seqLogo(probs, ic.scale=FALSE)
seqLogo(probs)

# Generate a background set of sequences by shuffling
shuffle <- function(dna) {
    # Convert to a vector of single bases
    charvec <- strsplit(as.character(dna),"")[[1]]
    # Shuffle the vector
    shuffled_charvec <- sample(charvec)
    # Convert back to a DNA string
    DNAString( paste(shuffled_charvec, collapse="") )
}

background_seqs <- DNAStringSet( lapply(initiation_seqs, shuffle) )
names(background_seqs) <- paste0(names(background_seqs), "-shuffled")


library(motifRG)

results <- findMotifFgBg(
    initiation_seqs, background_seqs,
    both.strand=FALSE, start.width=4)

summaryMotif(results$motif, results$category)

motifHtmlTable(results)

refined <- refinePWMMotif(results$motifs[[1]]@match$pattern, initiation_seqs)
seqLogo(refined$model$prob)


writeXStringSet(initiation_seqs, "fg.fa")
writeXStringSet(background_seqs, "bg.fa")


system("meme -dna -maxsize 1000000 fg.fa")


# ____________________
# ==== Next steps ====

### _______________________________________
### ---->> BSgenome.Hsapiens.UCSC.hg38 ----

### _____________________________________________
### ---->> TxDb.Hsapiens.UCSC.hg38.knownGene ----

### ________________________
### ---->> org.Hs.eg.db ----

### ___________________
### ---->> biomaRt ----

### _________________________
### ---->> AnnotationHub ----

library(AnnotationHub)
ah <- AnnotationHub()

# ah contains a large collection of records that can be retrieved
ah
length(ah)
colnames( mcols(ah) )
table( ah$rdataclass )

# query() searches for terms in an unstructured way
records <- query(ah, c("Ensembl", "85", "Saccharomyces cerevisiae"))
records

mcols(records)
mcols(records)[,c("title","rdataclass")]

# Having located records of interest,
# your R script can refer to the specific AH... record,
# so it always uses the same version of the data.
ah[["AH51399"]]
sc_genome <- import( ah[["AH51399"]] )
sc_granges <- ah[["AH51088"]]

# More recent versions of Bioconductor also allow you to
# retrieve TxDb (and similar EnsDb) objects.


query(ah, c("OrgDb", "Saccharomyces cerevisiae"))
sc_orgdb <- ah[["AH49589"]]

# An OrgDb contains information about genes in an organism
columns(sc_orgdb)
head( keys(sc_orgdb, "ORF") )
select(sc_orgdb, "YFL039C", c("GENENAME", "DESCRIPTION"), keytype="ORF")


sessionInfo()
