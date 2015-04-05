#### 5.2 the available data
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

# source("http://bioconductor.org/biocLite.R")
# biocLite()
# biocLite("ALL")

library(biobase)
library(ALL)
data(ALL)

pD <- phenoData(ALL)
varMetadata(pD)

dta <- tbl_df(pD@data)

# B- or T-cell frequencies (the type of acute lymphoblastic leukemia)

dta %>%
    group_by(BT) %>%
    summarise(N = n()) %>%
    mutate(rfrq = N / sum(N)) %>%
    arrange(-rfrq)

# BT by mol.biol cross-tabulation
# mol.biol codes the type of cytogenetic abnormality found on each sample
# a value of 'NEG' codes 'no abnormality'
dta %>%
    group_by(BT, mol.biol) %>%
    tally() %>%
    spread(key = mol.biol, value = n, fill = 0, drop = FALSE)

# some information on the genes and samples
featureNames(ALL)[1:10]
sampleNames(ALL)[1:10]

# focus of analysis will be on B-cell ALL cases
# in particular, samples with a subset of mutations (which will be target class)

dtab <- dta %>%
    filter(BT %in% levels(BT)[1:5], mol.biol %in% levels(mol.biol)[1:4]) %>%
    mutate(BT = factor(BT), mol.biol = factor(mol.biol))

tgt.cases <- which(ALL$BT %in% levels(ALL$BT)[1:5] &
                   ALL$mol.biol %in% levels(ALL$mol.biol)[1:4])
ALLb <- ALL[, tgt.cases]
ALLb$BT <- factor(ALLb$BT)
ALLb$mol.biol <- factor(ALLb$mol.biol)

saveRDS(ALLb, file = 'ALLb.rds')

#### 5.2.1 exploring the dataset
es <- exprs(ALLb)
dim(es)

# 12,625 variables
# 94 observations
# way too many variables
# first goal is to reduce the cardinality of the set of variables

quantile(as.vector(es), seq(0, 1, 0.10))

# Histogram of Overall Gene Expression values
hst <- ggplot(data.frame(V1 = as.vector(es)), aes(x = V1)) +
       geom_histogram()

# TO DO: 
#    1. histogram density instead of count
#    2. (slightly) smaller bins
#    3. title, y-axis, and x-axis labeling
#    4. vertical, dashed lines indicating 
#       1Q, median, 3Q, mean of values between 1Q and 3Q (shorth)

# are the distributions of the gene expression levels of the samples with a 
# particular mutation different from each other?

df <- data.frame(sapply(X = levels(ALLb$mol.biol), 
             FUN = function(x) {summary(as.vector(es[, which(ALLb$mol.biol == x)]))}))
df$rng <- apply(X = df[, 1:4], MARGIN = 1, FUN = function(x) {max(x) - min(x)})
df$sd  <- apply(X = df[, 1:4], MARGIN = 1, FUN = sd)
df

#### 5.3 Gene (Feature) Selection

# Feature selection is the process of selecting a subset of features 
# (i.e. variables) in order to 'focus' analysis.
# A more general approach is to decide the weight attached to each feature.
# Two types of appraoches to feature selection:
#     1. filter - step-by-step use of statistical properties of features to 
#                 select final set
#     2. wrapper - iterative search process using data mining tools in selection
#                  process

#### 5.3.1 Simple Filters Based on Dustribution Properties
# let's get a view of the distribution of the expression levels of each gene
# across all individuals 

# recall that rows of es are features and columns are samples

rowIQRs <- function(em) {
    rowQ(imat = em, which = ceiling(0.75 * ncol(em))) -
    rowQ(imat = em, which = floor(0.25 * ncol(em)))
}
plot(x = rowMedians(es), y = rowIQRs(es),
     xlab = 'Median expression level',
     ylab = 'IQR expression level',
     main = 'Main characteristics of Gene Expression Levels')

# many genes have IQRs near zero (low variability), and these are safe(r) to 
# remove because they provide reason to believe that they will not be useful in 
# discriminating among the different mutations of B-cell ALL.
# CAVEAT: we're focusing on genes individually and not conjointly
#         there is risk that some of these soon-to-removed genes, when included,
#         could be useful for classification
#         see RELIEF method

# use heuristic threshold for removal - remove gene if IQR is less than 1 / 5 of
# global IQR

biocLite("genefilter")
biocLite("hgu95av2.db")
library(genefilter)

ALLb <- nsFilter(eset = ALLb,
                 var.func = IQR,
                 var.cutoff = IQR(as.vector(es)) / 5,
                 feature.exclude = "^AFFX")

ALLb <- ALLb$eset
es <- exprs(ALLb)
dim(es)

# alternative approach
# dtaes <- tbl_df(data.frame(es)) %>%
#          mutate(gene = rownames(es)) %>%
#          filter(str_sub(gene, 1, 2) != "AF") 
# 
# dtaes$iqr <- apply(dtaes[, 1:(ncol(dtaes) - 1)], MARGIN = 1, FUN = IQR)
# dtaes
# 
# global_iqr <- IQR(as.vector(es))
# #iqr_20p <- quantile(dtaes$iqr, seq(0, 1, 0.20))["20%"]
# 
# dtaes <- dtaes %>%
#          select(gene, iqr) %>%
#          filter(iqr > global_iqr * (1 / 5) )

#### 5.3.2 ANOVA filters