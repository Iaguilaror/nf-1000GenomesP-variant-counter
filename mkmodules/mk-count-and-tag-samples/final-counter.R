## Load libraries
library("dplyr")
library("tidyr")

## Read args from command line
args = commandArgs(trailingOnly=TRUE)

## For debugging only
# args[1] <- "test/results/sample.chr21_GRCh38.genotypes.20170504.count_block.tmp"
# args[2] <- "test/metadata/integrated_call_samples_v3.20130502.ALL.panel"
# args[3] <- "test/results/sample.chr21_GRCh38.genotypes.20170504.counts.tsv"
# args[4] <- "testchr"

## Load data
widedata.df <- read.table(file = args[1], header = T, sep = "\t", stringsAsFactors = F)

## Transform wide to long format
longdata.df <- gather(widedata.df, variant_type, numbers, nRefHom:nMissing, factor_key=TRUE)

## Load metadata for tagging
metadata.df <- read.table(file = args[2], header = T, sep = "\t", stringsAsFactors = F)
## Perform join
joint.df <- longdata.df %>% left_join(metadata.df, by = "sample")

## add chromosome info
joint.df$chromosome <- args[4]

#### Finally, save tagged dataframe
write.table(x = joint.df, file = args[3], append = F, quote = F, sep = "\t", row.names = F, col.names = T)