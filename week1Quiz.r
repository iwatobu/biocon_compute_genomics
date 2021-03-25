#' week1 quiz
library(GenomicRanges)
library(IRanges)
library(AnnotationHub)
ah <- AnnotationHub()

#' Use the AnnotationHub package to obtain data on "CpG Islands" in the human genome.
#' Question: How many islands exists on the autosomes?
#' Question: How many CpG Islands exists on chromosome 4.
ah$species
#query(ah, c('Homo sapiens','hg19', 'CpG Islands'))
ah <- subset(ah, species=='Homo sapiens')
unique(ah$genome)
query(ah, 'CpG Islands')
query(ah, 'CpG Islands')$genome

# get the hg19
cpg <- query(ah, 'CpG Islands')[[1]]
show(cpg)
seqlevels(cpg)
autosomes <- paste0('chr',1:22)
#keepStandardChromosomes(cpg)
length(subset(cpg, seqnames(cpg) %in% autosomes))  # 26641
length(subset(cpg, seqnames(cpg) == 'chr4'))  # 1031
cpg <- subset(cpg, seqnames(cpg) %in% autosomes)

#' Obtain the data for the H3K4me3 histone modification for the H1 cell line (E003) from Epigenomics Roadmap, 
#' using AnnotationHub. Subset these regions to only keep regions mapped to the autosomes (chromosomes 1 to 22).
#' use narrow peak
#' Question: How many bases does these regions cover?
test <- query(ah, c('H3K4me3','E003','Roadmap')) 
View(data.frame(mcols(test))) 

h3k4me3 <- test[[2]]
h3k4me3 <- subset(h3k4me3, seqnames(h3k4me3) %in% autosomes)
print(sum(width(h3k4me3)))

#' Obtain the data for the H3K27me3 histone modification for the H1 cell line from Epigenomics Roadmap, 
#' using the AnnotationHub package. Subset these regions to only keep regions mapped to the autosomes. 
#' In the return data, each region has an associated "signalValue".
#' Question: What is the mean signalValue across all regions on the standard chromosomes?
test <- query(ah, c('H3K27me3','E003','Roadmap')) 
View(data.frame(mcols(test))) 
h3k27me3 <- test[[2]]
h3k27me3 <- subset(h3k27me3, seqnames(h3k27me3) %in% autosomes)
print(mean(h3k27me3$signalValue))



#' Bivalent regions are bound by both H3K4me3 and H3K27me3.
#' we will say a region is bivalent if it is enriched in both H3K4me3 and H3K27me3.
#' Note that histone modification marks does not have a strand.
#' Question: Using the regions we have obtained above, how many bases on the standard chromosomes are bivalently marked?
strand(h3k4me3)
strand(h3k27me3)

#subsetByOverlaps(h3k4me3, h3k27me3)
#subsetByOverlaps(h3k27me3, h3k4me3)

biv <- intersect(h3k4me3, h3k27me3)
sum(width(biv)) #10289096

#' We will examine the extent to which bivalent regions overlap CpG Islands.
#' Question: how big a fraction (expressed as a number between 0 and 1) of the 
#' bivalent regions, overlap one or more CpG Islands?

length(subsetByOverlaps(biv,cpg)) / length(biv)

#' Question: How big a fraction (expressed as a number between 0 and 1) of the bases 
#' which are part of CpG Islands, are also bivalent marked.

sum(width(intersect(biv,cpg)))/sum(width(cpg))


#' Question: How many bases are bivalently marked within 10kb of CpG Islands?
#' Tip: consider using the "resize()"" function.
#' 
sum(width(intersect(biv, resize(cpg, width = 20000 + width(cpg) , fix = 'center'))))


#' Question: How big a fraction (expressed as a number between 0 and 1) of the human genome is 
#' contained in a CpG Island?
#' Tip 1: the object returned by AnnotationHub contains "seqlengths".
#' Tip 2: you may encounter an integer overflow. As described in the session on R Basic Types, 
#' you can address this by converting integers to numeric before summing them, "as.numeric()".

sum(as.numeric(seqlengths(cpg)[1:22]))
sum(width(cpg))

20304026/2881033286

#' Question: Compute an odds-ratio for the overlap of bivalent marks with CpG islands.
#2x2
inOut = matrix(0, ncol = 2, nrow = 2)
colnames(inOut) = c('In', 'Out')
rownames(inOut) = c('In', 'Out')
inOut

inOut[1,1] = sum(width(intersect(biv, cpg))) # both biv and cpg
inOut[1,2] = sum(width(setdiff(biv, cpg, ignore.strand = TRUE ))) # biv but no cpg # anyway all the strands are * 
inOut[2,1] = sum(width(setdiff(cpg, biv, ignore.strand = TRUE )))# cpg but not biv
inOut[2,2] = sum(seqlengths(cpg)[1:22]) - sum(inOut)  # not cpg, not biv

inOut
fisher.test(inOut)

OR= inOut[1,1] * inOut[2,2] / (inOut[1,2] * inOut[2,1])

#169.0962





