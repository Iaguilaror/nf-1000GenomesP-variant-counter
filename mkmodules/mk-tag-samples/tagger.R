## Load libraries
library("dplyr")

## Read args from command line
args = commandArgs(trailingOnly=TRUE)

## For debugging only
# args[1] <- "test/data/chibiALL.chrY.phase3_integrated_v1b.20130502.genotypes.allstats.tsv"
# args[2] <- "test/metadata/integrated_call_samples_v3.20130502.ALL.panel"
# args[3] <- "test/data/chibiALL.chrY.phase3_integrated_v1b.20130502.genotypes.allstats.tagged.tsv"

data.df <- read.table(file = args[1], header = F, sep = "\t", stringsAsFactors = F)
  
metadata.df <- read.table(file = args[2], header = T, sep = "\t", stringsAsFactors = F)

## Perform join
joint.df <- data.df %>% left_join(metadata.df, by = c("V1" = "sample") )

#### Finally, save tagged dataframe
write.table(x = joint.df, file = args[3], append = F, quote = F, sep = "\t", row.names = F, col.names = F)