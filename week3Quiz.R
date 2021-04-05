setwd('~/Desktop/Bioconductor/')

library(biomaRt)
library(GenomicRanges)
library(ALL)
data(ALL)


# ---- Q1 ----
'''
Question: What is the mean expression across all features for sample 5 in the
ALL dataset (from the ALL package)?

#5.629627
'''
library(ALL)
data(ALL)
ALL # ExpressionSet
sampleNames(ALL)[5]
exprs(ALL)[1:4,1:5]
mean(exprs(ALL)[,5])
#5.629627

# ---- Q2 ----
"""
Question 2
We will use the biomaRt package to annotate an Affymetrix microarray. 
We want our results in the hg19 build of the human genome and we therefore
need to connect to Ensembl 75 which is the latest release on this genome version. 
How to connect to older versions of Ensembl is described in the biomaRt package vignette; 
it can be achived with the command 
mart <- useMart(host='feb2014.archive.ensembl.org', biomart = 'ENSEMBL_MART_ENSEMBL')
Question: Using this version of Ensembl, annotate each feature of the ALL dataset 
with the Ensembl gene id. How many probesets (features) are annotated with more 
than one Ensembl gene id?

#1045
"""

mart <- useMart(host='feb2014.archive.ensembl.org', biomart = 'ENSEMBL_MART_ENSEMBL')
listDatasets(mart)
searchDatasets(mart, 'hsapiens')
#                  dataset                     description    version
# 26 hsapiens_gene_ensembl Homo sapiens genes (GRCh37.p13) GRCh37.p13
ensemble75 <- useDataset(dataset = 'hsapiens_gene_ensembl', mart = mart)

# alternative way: 
listEnsembl(version = 75)
ensemble75 <- useEnsembl(biomart = 'genes', 
                         dataset = 'hsapiens_gene_ensembl',
                         version = 75)
ALL # Annotation hgu95av2
searchFilters(ensemble75, 'hg_u95av2')
'''             name                                         description
78  with_affy_hg_u95av2 with Affymetrix Microarray hg u95av2 ID(s) probeset
164      affy_hg_u95av2       Affy hg u95av2 probeset ID(s) [e.g. 32647_at]
'''
listAttributes(ensemble75)[1:4,]

ids <- featureNames(ALL)
mapping <- getBM(attributes = c('ensembl_gene_id', 'affy_hg_u95av2'),
                 filters = 'affy_hg_u95av2',
                 values = ids,
                 mart = ensemble75)

counts <- table(mapping$affy_hg_u95av2)
sum(counts > 1)
#1045


# ---- Q3 ----
'''
Question 3
Question: How many probesets (Affymetrix IDs) are annotated with one or more genes 
on the autosomes (chromosomes 1 to 22).

#11016
'''
searchAttributes(ensemble75, 'chr')[1:5,]
"chromosome_name"
autosomes <- as.character(1:22)

mapping <- getBM(attributes = c('affy_hg_u95av2', 'ensembl_gene_id', 'chromosome_name'),
                 filters = c('affy_hg_u95av2', 'chromosome_name'), 
                 values = list(ids, autosomes),
                 mart = ensemble75)
length(unique(mapping$affy_hg_u95av2))
#11016

# ---- Q4 ----
'''
Question 4
Use the MsetEx dataset from the minfiData package. Part of this question is to use 
the help system to figure out how to address the question.
Question: What is the mean value of the Methylation channel across the features 
for sample “5723646052_R04C01”?

# 7228.277
'''
BiocManager::install("minfiData")
data('MsetEx', package = 'minfiData')
MsetEx
# annotation ilmn12.hg19
class(MsetEx)
class?MethylSet
mean(getMeth(MsetEx)[, '5723646052_R04C01'])


# ---- Q5 ----
'''
Question 5
Question: Access the processed data from NCBI GEO Accession number GSE788. 
What is the mean expression level of sample GSM9024?

# 756.432

'''
library(GEOquery)
edata <- getGEO('GSM9024')  #GSM dataset
class(edata)
class?GSM
head(Meta(edata))
Table(edata)[1:5,]
Columns(edata)

mean(Table(edata)$VALUE)


# ---- Q6 ----
'''
Question 6
We are using the airway dataset from the airway package.
Question: What is the average of the average length across the samples in the expriment?

# 113.75
'''
BiocManager::install("airway")
data('airway',package = 'airway')
show(airway)
?airway # null
class?RangedSummarizedExperiment

rowRanges(airway) # GRangesList object of length 64102:
metadata(airway)
rowData(airway)
colData(airway)
mean(colData(airway)$avgLength)


# ---- Q7 ----
'''
Question 7
We are using the airway dataset from the airway package. 
The features in this dataset are Ensembl genes.
Question: What is the number of Ensembl genes which have a count of 1 read or more 
in sample SRR1039512?

#  25699
'''
edata <- assay(airway)
edata[1:5, ]
sum(edata[,'SRR1039512'] >=1)


# ---- Q8 ----
'''
Question 8
Question: The airway dataset contains more than 64k features. 
How many of these features overlaps with transcripts on the autosomes (chromosomes 1-22) 
as represented by the TxDb.Hsapiens.UCSC.hg19.knownGene package?
Clarification: A feature has to overlap the actual transcript, 
not the intron of a transcript. 
So you will need to make sure that the transcript representation does not contain introns.

#  26276
'''

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
?TxDb.Hsapiens.UCSC.hg19.knownGene
# no introns
rowRanges(airway) # has exon_id, exon_name

tx <- exonsBy(txdb, by='tx')
tx
# transcript range 
tx_gr <- reduce(unlist(tx))

# get GRanges on autosomes 
autosomes <- paste0('chr', 1:22)
tx_gr <- tx_gr[seqnames(tx_gr) %in% autosomes,]

subsetByOverlaps(rowRanges(airway), tx_gr)
# 0
seqlevels(rowRanges(airway))  # re format  chr1 -> 1 
?renameSeqlevels
seqlevels(tx_gr) <- sub('chr', '',seqlevels(tx_gr))

ov <- subsetByOverlaps(rowRanges(airway), tx_gr)
length(ov)

# ---- Q9 ----
'''
Question 9
The expression measures of the airway dataset are the number of reads mapping to each feature. 
In the previous question we have established that many of these features do not overlap 
autosomal transcripts from the TxDb.Hsapiens.UCSC.hg19.knownGene. 
But how many reads map to features which overlaps these transcripts?

Question: For sample SRR1039508, how big a percentage (expressed as a number between 0 and 1) 
of the total reads in the airway dataset for that sample, 
are part of a feature which overlaps an autosomal TxDb.Hsapiens.UCSC.hg19.knownGene transcript?

#0.9004193
'''

sr <- sum(edata[names(ov), 'SRR1039508'])
sr2 <- sum(edata[,'SRR1039508'])
sr/sr2

# ---- Q10 ----
'''
Question 10
Consider sample SRR1039508 and only consider features which overlap autosomal transcripts 
from TxDb.Hsapiens.UCSC.hg19.knownGene. We should be able to very roughly divide these transcripts 
into expressed and non expressed transcript. 

Expressed transcripts should be marked by H3K4me3 at their promoter. 
The airway dataset have assayed “airway smooth muscle cells”. 
In the Roadmap Epigenomics data set, the E096 is supposed to be “lung”. 
Obtain the H3K4me3 narrowPeaks from the E096 sample using the AnnotationHub package.

Question: What is the median number of counts per feature (for sample SRR1039508) containing a 
H3K4me narrowPeak in their promoter (only features which overlap autosomal transcripts 
from TxDb.Hsapiens.UCSC.hg19.knownGene are considered)?

Clarification: We are using the standard 2.2kb default Bioconductor promotor setting.

Conclusion Compare this to the median number of counts for features without a H3K4me3 peak. 
Note that this short analysis has not taken transcript lengths into account and 
it compares different genomic regions to each other; this is highly suscepticle to bias 
such as sequence bias.

243
224
232
240
'''
# my answer is 163 

library(AnnotationHub)
ah <- AnnotationHub()
ah <- subset(ah, species == 'Homo sapiens')
qh <- query(ah, c('E096','H3K4me3'))
h3k4me3 <- qh[[2]] # second one is the narrow peak

tx <- transcripts(txdb)
tx <- subset(tx, seqnames(tx) %in% autosomes)
prom <- promoters(tx)
prom_peak <- subsetByOverlaps(prom, h3k4me3) 
expressed_tx <- subset(tx[tx$tx_id %in% prom_peak$tx_id],)
seqlevels(expressed_tx) <- sub('chr', '', seqlevels(expressed_tx))


# expressed features 
expressed_features <- names(subsetByOverlaps(rowRanges(airway), expressed_tx)) 

# total 
sr <- edata[names(ov), 'SRR1039508']

# expressed: 

median( sr[names(sr) %in% expressed_features]) 
#163
median(sr[! (names(sr) %in% expressed_features) ] )
# 0







