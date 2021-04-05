
# ---- shortreads and rsamtoos ----
#' The yeastRNASeq experiment data package contains FASTQ files from an RNA seq experiment
#' in yeast. When the package is installed, you can access one of the FASTQ files by the path 
#' Question: What fraction of reads in this file has an A nucleotide in the 5th base of the read?
#0.3638

library(Rsamtools)
library(Biostrings)

#BiocManager::install("yeastRNASeq")
library(ShortRead)
library(yeastRNASeq)
fastqFilePath <- system.file("reads", "wt_1_f.fastq.gz", package = "yeastRNASeq")
reads <- readFastq(fastqFilePath)
reads #1000000 reads; width: 36 cycles
aset <- sread(reads) # DNAstringSet
b5 <- subseq(aset, 5,5)
sum(b5 == DNAString('A'))/length(aset)

#' This is a continuation of Question 1.
#' Question: What is the average numeric quality value of the 5th base of these reads?
#28.93

mean(as(quality(reads),'matrix')[,5])


#' Question 3 
#' The leeBamViews experiment data package contains aligned BAM files from an RNA seq experiment in yeast 
#' (the same experiment as in Questions 1 and 2, but that is not pertinent to the question). 
#' You can access one of the BAM files by the path given by
BiocManager::install('leeBamViews')
library(leeBamViews)
bamFilePath <- system.file("bam", "isowt5_13e.bam", package="leeBamViews")

#' These reads are short reads (36bp) and have been aligned to the genome using a standard aligner, 
#' ie. potential junctions have been ignored (this makes some sense as yeast has very few junctions 
#' and the reads are very short). A read duplicated by position is a read where at least one more 
#' read shares the same position. We will focus on the interval from 800,000 to 801,000 on yeast chromosome 13.
#' Question: In this interval, how many reads are duplicated by position?
#129.00


library(Rsamtools)
bamFile <- BamFile(bamFilePath)
bamFile 
seqinfo(bamFile)

gr <- GRanges(seqnames = 'Scchr13', ranges = IRanges(start = 800000, end = 801000))
params <- ScanBamParam(which = gr, what = scanBamWhat() )
aln <- scanBam(bamFile, param = params)
aln
names(aln)
counts <- sort(table(aln[[1]]$pos),decreasing = T)
sum(counts[counts>1])


#' This is a continuation of Question 3.
#' The package contains 8 BAM files in total, representing 8 different samples from 4 groups. 
#' A full list of file paths can be had as

bpaths <- list.files(system.file("bam", package="leeBamViews"), pattern = "bam$", full=TRUE)

#' An objective of the original paper was the discovery of novel transcribed regions in yeast. 
#' One such region is Scchr13:807762-808068.
#' Question: What is the average number of reads across the 8 samples falling in this interval?

#90.25

summaryFunction <- function(bamFile, params){
    aln <- scanBam(bamFile, param=params)[[1]]
    n_reads <- length(aln$seq)
    return(n_reads)
}

gr <- GRanges(seqnames = 'Scchr13', ranges = IRanges(start = 807762, end = 808068))
params <- ScanBamParam(which = gr, what = scanBamWhat() )
N=c()
for( i in 1:8) {
    bamFile <- BamFile(bpaths[i])
    n = summaryFunction(bamFile, params)
    N=c(N,n)
}
mean(N)

# ---- oligo ----
#' In the lecture on the oligo package an ExpressionSet with 18 samples is constructed, 
#' representing normalized data from an Affymetrix gene expression microarray. 
#' The samples are divided into two groups given by the group variable.
#' What is the average expression across samples in the control group for the “8149273” 
#' probeset (this is a character identifier, not a row number).
# 7.0218

library(oligo)
library(GEOquery)

# get data 
getGEOSuppFiles('GSE38792')
list.files('GSE38792')
untar('GSE38792/GSE38792_RAW.tar', exdir = 'GSE38792/CEL' )
list.files('GSE38792/CEL')

# read data 
celfiles <- list.files('GSE38792/CEL', full.names = TRUE)
rawData <- read.celfiles(celfiles)
rawData # GeneFeatureSet. very similar to a GeneExpressionSet
getClass('GeneFeatureSet')

pData(rawData)
# build pData
filename <- sampleNames(rawData)
filename
pData(rawData)$filename <- filename
sampleNames <- sub(".*_", "", filename) # remove GSM***_
sampleNames
sampleNames <- sub(".CEL.gz$", "", sampleNames )
sampleNames(rawData) <- sampleNames
pData(rawData)$group <- ifelse(grepl('^OSA', sampleNames(rawData)), 'OSA', 'Control')

pData(rawData)

# norm 
normData <- rma(rawData) # rma always works well
normData # ExpressionSet

featureNames(normData)[1:10]
mean(exprs(normData)["8149273", 1:8])


#' Question 6 This is a continuation of Question 5.
#' Use the limma package to fit a two group comparison between the control group and the OSA group, 
#' and borrow strength across the genes using eBayes(). Include all 18 samples in the model fit.
#'  What is the absolute value of the log foldchange (logFC) of the gene with the lowest P.value
#0.7126

library(limma)
normData$group <- factor(normData$group, levels = c('Control','OSA'))
normData$group

design <- model.matrix( ~ normData$group) # refer to the user guide
head(design)

fit <- lmFit(normData, design)
fit <- eBayes(fit)
topTable(fit)
topTable(fit, n=1)

#' Question 7 This is a continuation of Question 6.
#' Question: How many genes are differentially expressed between the two groups at an adj.P.value cutoff of 0.05?
# 0


#' Question 8
#' An example 450k dataset is contained in the minfiData package. This dataset contains 6 samples; 
#' 3 cancer and 3 normals. Cancer has been shown to be globally hypo-methylated (less methylated) 
#' compared to normal tissue of the same kind. 
#' Take the RGsetEx dataset in this package and preprocess it with the preprocessFunnorm function. 
#' For each sample, compute the average Beta value (percent methylation) across so-called OpenSea loci.
#' Question: What is the mean difference in beta values between the 3 normal samples and the 3 cancer samples,
#' across OpenSea CpGs?
0.1811
0.0034
0.0846
0.1732

library(minfi)
library(minfiData)
data("RGsetEx")
RGsetEx # RGChannelSet 
pd <- pData(RGsetEx)

MsetEx <- preprocessFunnorm(RGsetEx)
'''
[preprocessFunnorm] Background and dye bias correction with noob
[preprocessFunnorm] Mapping to genome
[preprocessFunnorm] Quantile extraction
[preprocessFunnorm] Normalization
'''
MsetEx #GenomicRatioSet  ready for analysis
colnames(MsetEx)
beta <- getBeta(MsetEx)
colnames(beta)
pd

annotation(MsetEx)
head(Manifest)

sts <- getIslandStatus(MsetEx) #OpenSea
openseas <- rownames(MsetEx)[sts=='OpenSea']

#normal 125, cancer 346

mean(beta[openseas,c(1,2,5)]) - mean(beta[openseas,c(3,4,6)])
0.08863657

colMeans(beta[openseas,])
'''
5723646052_R02C02 5723646052_R04C01 5723646052_R05C02 5723646053_R04C02 5723646053_R05C02 5723646053_R06C02 
0.6969593         0.7094411         0.6039988         0.6433666         0.7062855         0.5994108 
'''

#' Question 9 This is a continuation of Question 8.
#' The Caco2 cell line is a colon cancer cell line profiled by ENCODE. Obtain the narrowPeak DNase hyper 
#' sensitive sites computed by the analysis working group (AWG).
#' Question: How many of these DNase hypersensitive sites contain one or more CpGs on the 450k array?

# 40151

library(AnnotationHub)
ah <- AnnotationHub()
qh <- query(ah, c('caco2','encode','awg'))
dhs <- qh[[1]]
dhs # GRanges object
cpg_gr <- granges(MsetEx)

findOverlaps(dhs,cpg_gr)
sum(countOverlaps(dhs,cpg_gr) >0)


#' Question 10 The zebrafishRNASeq package contains summarized data from an RNA-seq experiment in zebrafish 
#' in the form of a data.frame called zfGenes. The experiment compared 3 control samples to 3 treatment samples.
#' Each row is a transcript; the data.frame contains 92 rows with spikein transcripts; these have a rowname 
#' starting with “ERCC”. Exclude these rows from the analysis.
#' Use DESeq2 to perform a differential expression analysis between control and treatment. 
#' Do not discard (filter) genes and use the padj results output as the p-value.
#' Question: How many features are differentially expressed between control and treatment (ie. padj <= 0.05)?

87
984
0
531

BiocManager::install("zebrafishRNASeq")
library(zebrafishRNASeq)
data("zfGenes")
head(zfGenes)

# exclude spikein
spikein <- grep('^ERCC', rownames(zfGenes))
zfGenes <- zfGenes[-spikein, ]

# coldata 
colnames(zfGenes)

coldata <- data.frame('condition' = c(rep('control',3), rep('treated',3)),
                      'sample_name' = colnames(zfGenes))
rownames(coldata) = coldata$sample_name
coldata$condition = factor(coldata$condition)

# deg
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = zfGenes,
                              colData = coldata,
                              design = ~ condition)
dds

dds <- DESeq(dds)
res <- results(dds)
res
sum(res$padj < 0.05, na.rm = T) # 116
