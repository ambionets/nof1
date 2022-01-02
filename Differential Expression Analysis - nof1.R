# Differential Expression Analysis - N-of-1-Analysis
# 2021-12-21
# Andrew Miller

setwd("/Users/andrew/Documents/Work & Volunteer/BME PhD/Jeyapalina Lab/01 - Lower_Limb_Implant_RNA_seq/03 - Explore/code/nof1")

## Get Data

### Sample data
library(tibble)
sampleData <- readr::read_tsv('../../data/SampleData.txt')
names(sampleData ) <- c('sample', 'surgery01_date', 'surgery02_date', 'blood_draw_date',
                        'time_period_3months', 'patient', 'surgery',
                        'time_since', 'dob', 'rin', 'batch', 'lane01_id', 'comment')

sampleData$time_period_3months <- factor(sampleData$time_period_3months,
                                         labels = c('PrS1-W1', 'PoS1-W1', 'PrS2-W1', 'PoS2-W1',
                                                    'PoS2-W2', 'PoS2-M1', 'PoS2-M3', 'PoS2-M6',
                                                    'PoS2-M9', 'PoS2-M12'))
sampleData$patient <- factor(sampleData$patient, levels = c("P1", "P2", "P3", "P4", "P5",
                                                            "P6", "P7", "P8", "P9", "P10"))
sampleData$batch <- factor(sampleData$batch, levels = c(1, 2))

#### Calculate age of patients
library(lubridate)
sampleData$surgery01_date = stringr::str_replace(sampleData$surgery01_date, "[0-9]+$", "20\\0")
sampleData$surgery01_date = mdy(sampleData$surgery01_date)
sampleData$dob = stringr::str_replace(sampleData$dob, "[0-9]+$", "19\\0")
sampleData$dob =  mdy(sampleData$dob)
#sampleData$age = floor_date(sampleData$surgery01_date - sampleData$dob, "day")
sampleData$age = year(sampleData$surgery01_date) - year(sampleData$dob) # 30-62

### Count data
countData <- readr::read_tsv('../../data/CountMatrix.txt')
countData <- column_to_rownames(countData, 'Geneid')
originalCountData <- countData

### Summarize data
#skimr::skim(sampleData)
#skimr::skim(countData)


## Sample Preprocessing
### Remove sample duplicates and filter low quality RIN.
library(dplyr)
sampleData <- sampleData %>%
    filter(rin > 3) %>%
    filter(!(is.na(time_period_3months)))
countData <- countData[, sampleData$sample]
sum(sampleData$sample == names(countData))

### Find and quantify the number of hemoglobin-related counts in the count data
library(biomaRt)
#### NOTE: It turns out that hemoglobin transcripts compose ~30% (29.36%) of the total number of transcripts.
ensembl <- biomaRt::useMart('ensembl', dataset = 'hsapiens_gene_ensembl')
#ensembl <- useEnsembl('ensembl', dataset = 'hsapiens_gene_ensembl', mirror = 'uswest')  # Wrapper for the useMart function. Useful if the default server/website is going slow.
countAnnotations <- biomaRt::getBM(attributes = c('ensembl_gene_id', 'external_gene_name', 'entrezgene_description'),
                                   filters = 'ensembl_gene_id',
                                   values = rownames(countData),
                                   mart = ensembl)
hemoglobinGenes <- countAnnotations %>%
    filter(stringr::str_detect(entrezgene_description, 'hemoglobin'))

#### Ensemble hemoglobin gene IDs: ENSG00000196565, ENSG00000213931, ENSG00000213934, ENSG00000223609, ENSG00000244734, ENSG00000086506, ENSG00000130656, ENSG00000188536, ENSG00000206172, ENSG00000206177
hemoglobinGenes <- hemoglobinGenes [2:11, ]    # Only using hemoglobin subunit genes.
#hemoglobin_genes <- c("ENSG00000263961", "ENSG00000196565", "ENSG00000213934", "ENSG00000223609", "ENSG00000244734", "ENSG00000086506", "ENSG00000130656", "ENSG00000169877", "ENSG00000188536", "ENSG00000206172", "ENSG00000206177")
transcriptTotals <- rowSums(countData)
hemoglobinTotals <- transcriptTotals[hemoglobinGenes$ensembl_gene_id]
hemoglobinPortion <- sum(hemoglobinTotals) / sum(transcriptTotals)
hemoglobinPortion


####  Visualize hemoglobin makeup
##### NOTE: Hemoglobin-related transcripts constituted between 39.7% at PoS2-W2 to 20.5% at PoS2-M9 when looking at hemoglobin percentages by time point.
library(ggplot2)

hgSampleCounts <- countData[hemoglobinGenes$ensembl_gene_id, ]
hgSampleTotals <- colSums(hgSampleCounts)
hgSamplePercentages <- hgSampleTotals / colSums(countData)

hgData <- tibble(sample = names(hgSampleTotals),
                 hg_counts = hgSampleTotals,
                 hg_percent = hgSamplePercentages * 100)

##### Join
hgSampleTime <- sampleData[, c('sample', 'time_period_3months')]
hgData <- left_join(hgData, hgSampleTime, by = c('sample' = 'sample'))


##### NOTE: total counts is not a good measure because of the variability total reads for sequencing.
hgSummaryData <- hgData %>%
    group_by(time_period_3months) %>%
    summarise(total = sum(hg_counts), avg = mean(hg_percent)) %>%
    arrange(desc(avg)) %>%
    print()

ggplot(hgSummaryData, aes(x = time_period_3months, y = avg)) +
    geom_point()


### Coefficient of variation (COV) distribution
countData = countData[rowSums(countData) != 0, ]    # Remove transcripts with no counts. 50917
oringalCountData = countData

#### Using apply
covDist = apply(countData, 1, function(x) sd(x)/mean(x) * 100)

#### Using purrr
##### Get the same values as the apply method.
#covDist2 = countData %>%
#    purrr::pmap(~ c(...)) %>%
#    purrr::map(~ (sd(.x)/mean(.x) * 100)) %>%
#    do.call(rbind, .)

library(ggplot2)
ggplot(data.frame(cov = covDist), aes(x = cov)) +
    geom_histogram() +
    #geom_density(alpha = 0.5) +
    labs(title = "Coefficient of Variation Distribution")
hist(covDist)

summary(covDist)

### Affect of filtering 50% and 75% of genes out with the lowest COV.
cov50CountData = countData[covDist > 168.85, ]  # Median, 25451 genes
cov75CountData = countData[covDist > 471.61, ]  # Median, 13009 genes
countData = cov50CountData

### Hemoglobin makes up between 20.5% to 39.7% depending on the time point. We
### will analyse the data without the hemoglobin included and requiring at least
### 2 counts in 5 samples for the gene to be included. Note that DESeq2 performs
### it's own filtering based the the p-value distribution.
####  Remove hemoglobin subunit transcripts and apply minimal filtering.
#countData <- countData[setdiff(rownames(countData), hemoglobinGenes$ensembl_gene_id), ] # 58,735 -> 58,725
#countData <- countData[rowSums(countData >= 2) >= 5, ]   #29242 genes

## N-of-1 Analysis
### Installation
#library(nof1)
#
#normCountData = normalize_multi_gene_data

## Differential Expression Analysis using DESeq2
library(DESeq2)
modelColumns <- c('sample', 'time_period_3months', 'patient', 'rin', 'batch', 'age')
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = sampleData[, modelColumns],
                              design = ~ rin + patient + time_period_3months)

### Incorporating patient id as a covariate.
start <- Sys.time()
design(dds) <- ~ patient + time_period_3months
dds = DESeq(dds, BPPARAM = MulticoreParam())
deseqDuration <- Sys.time() - start
deseqDuration

ddsPatient = dds

#### Summary
comparisons = c('PoS1-W1', 'PrS2-W1', 'PoS2-W1',
                'PoS2-W2', 'PoS2-M1', 'PoS2-M3',
                'PoS2-M6', 'PoS2-M9', 'PoS2-M12')

sigGeneCountsPatient = tibble(time_points = comparisons,
                              sig_genes = integer(9))

for (index in 1:length(sigGeneCountsPatient$time_points)) {
    res <- results(dds, contrast = c('time_period_3months', sigGeneCountsPatient$time_points[index], 'PrS1-W1'), alpha = 0.05, BPPARAM = MulticoreParam())
    sigGeneCountsPatient$sig_genes[index] <- sum(res$padj < 0.05, na.rm = TRUE)
}

### Incorporate RIN and patient id as covariates.
start <- Sys.time()
design(dds) <- ~ rin + patient + time_period_3months
dds = DESeq(dds, BPPARAM = MulticoreParam())
deseqDuration <- Sys.time() - start
deseqDuration

ddsRinPatient = dds

#### Summary
comparisons = c('PoS1-W1', 'PrS2-W1', 'PoS2-W1',
                'PoS2-W2', 'PoS2-M1', 'PoS2-M3',
                'PoS2-M6', 'PoS2-M9', 'PoS2-M12')

sigGeneCountsRinPatient = tibble(time_points = comparisons,
                                 sig_genes = integer(9))

for (index in 1:length(sigGeneCountsRinPatient$time_points)) {
    res <- results(dds, contrast = c('time_period_3months', sigGeneCountsRinPatient$time_points[index], 'PrS1-W1'), alpha = 0.05, BPPARAM = MulticoreParam())
    sigGeneCountsRinPatient$sig_genes[index] <- sum(res$padj < 0.05, na.rm = TRUE)
}