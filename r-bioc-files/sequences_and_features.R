# This file is generated from the corresponding .Rmd file



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

# or browse https://bioconductor.org/packages/release/BiocViews.html#___Software




#////////////////////////////////////////
# 1 DNA sequences and genomic ranges ----

library(Biostrings)      # Provides DNAString, DNAStringSet, etc
library(BSgenome)        # Provides getSeq()
library(GenomicRanges)   # Provides GRanges containing genomic ranges, etc
library(GenomicFeatures) # Provides TxDb objects containing genes/transcripts/exons
library(rtracklayer)     # Provides import() and export()


# 1.1 DNAString ----

myseq <- DNAString("CCGCGCACCAAC")
myseq

class(myseq)

reverseComplement(myseq)
translate(myseq)

subseq(myseq, 3,5)
myseq[3:5]

as.character(myseq)


methods(class="DNAString")


?"DNAString-class"


# 1.2 DNAStringSet ----

myset <- DNAStringSet( list(chrI=myseq, chrII=DNAString("ACGTACGT")) )
myset

# A DNAStringSet is list-like
myset$chrII
# or myset[["chrII"]]
# or myset[[2]]


# 1.3 Challenge ----
# 
# Reverse complement the following DNA sequence and then translate to an
# amino acid sequence:
# 

GCTTTCGTTTTCGCC

# 
# 
#
# 1.4 GRanges ----

range1 <- GRanges("chrI", IRanges(start=3,end=5), "+")
range1
getSeq(myset, range1)

range2 <- GRanges("chrI", IRanges(start=3,end=5), "-")
getSeq(myset, range2)


seqnames(range1)
start(range1)
end(range1)
width(range1)
strand(range1)
as.data.frame(range1)


# Look at completions for
# range1@


# GRanges are like vectors:
c(range1, range2)

# GRanges can have metadata columns, and are often used like data frames:
mcols(range1)$wobble <- 10
range1
mcols(range1)
range1$wobble

# A handy way to create a GRanges
as("chrI:3-5:+", "GRanges")


# 1.5 Question ----
# 
# Based on what we saw for `DNAString`, how can we learn more about
# using `GRanges` and `IRanges` objects?
# 
# 
# 
#


#/////////////////////
# 2 Loading files ----

# 2.1 Loading sequences ----

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


# 2.2 Loading features ----

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
mcols(features) <- mcols(features)[,c("type","gene_name","gene_id","transcript_id","exon_id")]

features


# 2.3 Seqinfo and genome versions ----

seqinfo(features)
seqinfo(seqs)

myseqinfo <- seqinfo(seqs)
isCircular(myseqinfo) <- c(TRUE)

seqinfo(features) <- myseqinfo


Seqinfo(genome="hg19")
Seqinfo(genome="hg38")


# 2.4 Getting sequences ----

feat <- features[4]
feat
feat_seq <- getSeq(seqs, feat)
feat_seq
translate(feat_seq)


subset(features, gene_name == "lacA")
# Equivalently:
#   features[features$gene_name == "lacA" & !is.na(features$gene_name)]


cds <- subset(features, type == "CDS")
cds
# Equivalently:
#   features[features$type == "CDS"]




#/////////////////////////////////////
# 3 Further operations on GRanges ----

# 3.1 Intra-range ----

?"intra-range-methods"


feat <- features[4]
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


# 3.2 Inter-range ----

?"inter-range-methods"
?"nearest-methods"
?"setops-methods"


query <- as("Chromosome:9500-10000:+", "GRanges")
hits <- findOverlaps(query, features, ignore.strand=TRUE)
hits
subjectHits(hits)
features[subjectHits(hits)]

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
# setdiff(range(input),input)     (--)


# 3.3 Challenge ----
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
# Hint: Use `table()` to count occurrences of different strings.
# 
# 
#
# 3.4 Transcript database objects ----

subset(features, gene_name == "lacY")


# --------------------------------------------------> gene
#
# -------------------------------------------->       transcript
# ---------->         --->    ---------------->       exon
#       ---->         --->    ---------->             CDS
#
#
#                -----------------------------------> transcript
#                -------->       ---->    ----------> exon
#                     --->       ---->    -->         CDS


txdb <- makeTxDbFromGRanges(features)
txdb


genes(txdb)
transcriptsBy(txdb, by="gene")
exonsBy(txdb, use.names=TRUE)
cdsBy(txdb, use.names=TRUE)

cds_ranges <- cdsBy(txdb, use.names=TRUE)
cds_ranges$AIZ54182
cds_ranges[[1]]
unlist(cds_ranges)


extractTranscriptSeqs(seqs, cds_ranges)




#///////////////////////////////////////
# 4 Genome and annotation resources ----

# (return to slideshow)


# 4.1 Example packages ----

# 4.2 biomaRt ----

# 4.3 AnnotationHub ----

library(AnnotationHub)
ah <- AnnotationHub()

# ah contains a large collection of records that can be retrieved
ah
length(ah)
colnames( mcols(ah) )
table( ah$rdataclass )

# query() searches for terms in an unstructured way
records <- query(ah, c("Ensembl", "96", "Saccharomyces cerevisiae"))
records

mcols(records)
mcols(records)[,c("title","rdataclass")]

# Having located records of interest,
# your R script can refer to the specific AH... record,
# so it always uses the same version of the data.
sc_genome <- ah[["AH70449"]]
sc_granges <- ah[["AH69700"]]
sc_txdb <- ah[["AH69265"]]

# sc_genome is a TwoBitFile
# Can use getSeq on it without loading everything into memory
seqinfo(sc_genome)
getSeq(sc_genome, as("IV:1-10","GRanges"))
import(sc_genome)

# An OrgDb contains information about genes in an organism, and lets you map between different identifiers
query(ah, c("OrgDb", "Saccharomyces cerevisiae"))
sc_orgdb <- ah[["AH84129"]]
columns(sc_orgdb)
head( keys(sc_orgdb, "ENSEMBL") )
select(sc_orgdb, "YFL039C", c("GENENAME", "DESCRIPTION"), keytype="ENSEMBL")




#///////////////////////////////////////////////
# 5 Example: finding and discovering motifs ----

# 5.1 Finding a known motif ----

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


# 5.2 De novo motif finding ----

# Note: bacteria do not have introns
# In a eukaryote we would need to work with transcript
# sequences (extractTranscriptSeqs) and work out where
# the CDSs start within them.

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


writeXStringSet(initiation_seqs, "seqs.fasta")


### in a terminal ###
# conda activate
# conda install -c bioconda meme


### in a terminal ###
# meme -dna seqs.fasta


# 5.3 Challenge ----
# 
# Which genes have close matches to AGGAGGT (as we found earlier)
# immediately upstrand of their CDS?
# 
# Check some of the genes you find using IGV. (Once the E. coli FASTA
# file is loaded as the genome and the E. coli GTF file is loaded, you
# can search for a gene by typing its name in the location box.)
# 
# 
#

sessionInfo()

