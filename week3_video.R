BiocManager::install("GEOquery")

library(biomaRt)

## ---- expression set ----

library(ALL)
data(ALL)
show(ALL)

experimentData(ALL)

?ALL 
exprs(ALL)[1:4, 1:4]  # retrieve expression data from eSet

head(sampleNames(ALL))
head(featureNames(ALL))


pData(ALL)[1:4,] # p for pheno

ALL$sex

#first five samples
ALL[,1:5]

# first 5 features 
ALL[1:5, ]

ALL[1:5, 1:5]
# Annotation: hgu95av2

featureData(ALL) # none for this dataset, common.

## ---- annotation ------
ids = featureNames(ALL)[1:5]
ids

library(hgu95av2.db)


phenoData(ALL)
"""
An object of class 'AnnotatedDataFrame'
sampleNames: 01005 01010 ... LAL4 (128 total)
varLabels: cod diagnosis ... date last seen (21 total)
varMetadata: labelDescription
"""


## ---- GEO ----
# https://www.bioconductor.org/packages/release/bioc/vignettes/GEOquery/inst/doc/GEOquery.html

library(GEOquery)
eList <- getGEO('GSE11675') #GEO_Accession_number
length(eList)
names(eList)
eData <- eList[[1]]


## ---- biomaRt ---------
#https://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/biomaRt.html
#https://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/accessing_ensembl.html


library(biomaRt)
head(listMarts()) # listEnsemble()
mart <- useMart('ENSEMBL_MART_ENSEMBL')  # useEnsemble()
searchDatasets(mart = mart, pattern = "hsapiens")

mart
head(listDatasets(mart))
ensembl <- useDataset('hsapiens_gene_ensembl', mart)




values <- c('202763_at', '209310_s_at', '207500_at')

getBM(attributes = c('ensembl_gene_id', 'affy_hg_u133_plus_2'), 
      filters = 'affy_hg_u133_plus_2', values = values, mart = ensembl)

#the function 'listFilters' to get valid filter names
listFilters(ensembl)

att <- listAttributes(ensembl) #sometime need to go to ensembl website to understand the meanning
nrow(att)  #3183
tail(att) 

## very long list. difficult to find info effectively.
# organized in pages 
attributePages(ensembl)
#"feature_page" "structure"    "homologs"     "snp"          "snp_somatic"  "sequences"
att <- listAttributes(ensembl, page = 'feature_page')
att 
nrow(att) #204 #more tractable 

# there are featurs belong to mulitple page . can link things together. 
# query page individually then merge 

# read the veggnet to learn more 

## ---- S4 class -----
library(ALL)
library(GenomicRanges)

# in R , 3 class, S3 , S4 and S5(reference class)
# object of S3 class can change the classes , without error (not means it makes sense ) 
# S4 more validity check. Define an object 

data(ALL)
ALL
isS4(ALL)
# how to get help on class 
class?ExpressionSet
?"ExpressionSet-class"

# some convention of bioconducter class
# class starts with Capital letter
# class has a constructer
# eg : 
ExpressionSet() # to construct a ExpressionSet 
?ExpressionSet
new("ExpressionSet") # not recommended anymore. ()

# definition of a class
getClass('ExpressionSet')
'''
Class "ExpressionSet" [package "Biobase"]

Slots:
                                                                                  
Name:      experimentData          assayData          phenoData        featureData
Class:              MIAME          AssayData AnnotatedDataFrame AnnotatedDataFrame
                                                               
Name:          annotation       protocolData  .__classVersion__
Class:          character AnnotatedDataFrame           Versions

Extends: 
Class "eSet", directly
Class "VersionedBiobase", by class "eSet", distance 2
Class "Versioned", by class "eSet", distance 3
'''

# slot
@
slot()
#eg
ALL@annotation
slot(ALL, 'annotation')

opar = options(width = 80) #display 
getClass('ExpressionSet')
'''
Class "ExpressionSet" [package "Biobase"]

Slots:
  
  Name:      experimentData          assayData          phenoData
Class:              MIAME          AssayData AnnotatedDataFrame

Name:         featureData         annotation       protocolData
Class: AnnotatedDataFrame          character AnnotatedDataFrame

Name:   .__classVersion__
Class:           Versions

Extends: 
  Class "eSet", directly
Class "VersionedBiobase", by class "eSet", distance 2
Class "Versioned", by class "eSet", distance 3
'''

# BUT, you are not supposte to access data using @ sign in S4 class. because some slots are not to users 
# supposed to use accessor functions. 
# names are same as the slot, or like "getXX()"
annotation(ALL)

# class may get updated in bioconductor later versions . the object saved in older version may throw an error. 
# way to solve:
NEW_OBJECT = updateObject(OLD_OBJECT)
# not granteed to work. as ppl who write the pkg may not write the updateObject function. Should compalin loudly in the support forum 

validObject(ALL) # sometime takes long time to run on a big object. 

## ---- S4 methods -----

## ---- mimic method -- generic function 
mimicMethod <- function(x){
    if (is(x, 'matrix'))
        method1(x)
    if (is(x, 'data.frame'))
        method2(x)
    if (is(x, 'IRanges'))
        method3(x)
}

# eg 
library(GenomicRanges)
as.data.frame
base::as.data.frame

showMethods("as.data.frame")
getMethod('as.data.frame', 'GenomicRanges')
getMethod('as.data.frame', signature (x='GenomicRanges'))

# get help on methods
method?"as.data.frame, DataFrame"  # not supported anymore.
?"as.data.frame, GenomicRanges-method" # no documentation

showMethods('findOverlaps')
getMethod('findOverlaps', signature (query="GRangesList", subject="GenomicRanges"))
?'findOverlaps,GRangesList,GenomicRanges-method'























































