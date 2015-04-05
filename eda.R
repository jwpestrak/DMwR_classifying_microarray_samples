#### 5.2 the available data
library(dplyr)
library(tidyr)
library(ggplot2)

source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite("ALL")

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
