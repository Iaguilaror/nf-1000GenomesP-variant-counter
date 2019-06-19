## Load libraries
library("dplyr")
library("ggplot2")
library("viridis")

## Read args from command line
args = commandArgs(trailingOnly=TRUE)

## For debugging only
# args[1] <- "test/data/sample.allchrom_counts.tsv"
# args[2] <- "test/data/sample.allchrom_counts.basicQC.pdf"

## Load data
variants.df <- read.table(file = args[1], header = T, sep = "\t")

## get type of number of variants
register_types.v <- levels(variants.df$variant_type)
## set color pallete
mycolors <- viridis(length(register_types.v))

## start PDF device to recive multi plots
pdf(args[2])
## start loop for multiplotting
for (i in 1:length(register_types.v)) {

  ## automated plotting of numbers per sample
  plotable.df <- variants.df %>% filter(variant_type == register_types.v[i])
  plotable.df <- plotable.df %>% group_by(sample) %>% summarize(total = sum(numbers))

  # propper plot
  ## pre-calculate values for axis
  max_y <- round(max(plotable.df$total*1.1))
  axis_breaks_y <- round(seq(0,max_y, by = max_y/10))    
  
  plot.tmp <- ggplot(data = plotable.df, aes(x=sample, y=total)) +
    geom_bar(stat = "identity", fill = mycolors[i], color = "black", lwd = 0.5) +
    ggtitle(label = args[1]) +
    scale_y_continuous(name = paste0("total positions (", register_types.v[i],")"),
                       limits = c(0, max_y*1.1),
                       breaks = axis_breaks_y, 
                       labels = axis_breaks_y, 
                       expand = c(0.01, 0)) +
    scale_x_discrete(name = "",
                     expand = c(0.01, 0) ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, size = 3),
          plot.title = element_text(size=5))
  
  print(plot.tmp)
  
}

## create plot for total variants per sample
plotable_SNPs.df <- variants.df %>%
  filter(variant_type == "nNonRefHom" | variant_type == "nHets") %>%
  group_by(sample) %>%
  summarize(total_SNP = sum(numbers))

plotable_indels.df <- variants.df %>%
  filter(variant_type == "nIndels") %>%
  group_by(sample) %>%
  summarize(total_indel = sum(numbers))

## merge indel and snp data
plotable.df <- left_join(x = plotable_SNPs.df, y = plotable_indels.df, "sample")
## create 2 dummy columns for downstream legen setup
plotable.df$tag1 <- "SNP"
plotable.df$tag2 <- "indel"

## pre-calculate axis values
max_y <- max(plotable.df$total_SNP)
axis_breaks_y <- round(seq(0,max_y, by = max_y/10))    

## plot total SNPs and indels
plot.tmp <- ggplot() +
  geom_point(data = plotable.df,
             aes(x=sample, y=total_SNP, color = tag1, shape = tag1),
             size = 5, alpha = 0.8) +
  geom_point(data = plotable.df,
             aes(x=sample, y=total_indel, color = tag2, shape = tag2),
             size = 5, alpha = 0.8) +
  ggtitle(label = args[1]) +
  scale_color_manual(name = "", values = viridis(6)) +
  scale_shape_manual(values = c(18,16) ) +
  scale_y_continuous(name = paste0("total variants"),
                     limits = c(0, max_y*1.1),
                     breaks = axis_breaks_y,
                     labels = axis_breaks_y,
                     expand = c(0.0, 0)) +
  scale_x_discrete(name = "",
                   expand = c(0.02, 0.02) ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 3),
             plot.title = element_text(size=5))

print(plot.tmp)

## plot data per chromosome
plotable_SNPs.df <- variants.df %>%
  filter(variant_type == "nNonRefHom" | variant_type == "nHets") %>%
  group_by(sample, chromosome) %>%
  summarize(total = sum(numbers))
plotable_SNPs.df$type <- "SNP"

plotable_indels.df <- variants.df %>%
  filter(variant_type == "nIndels") %>%
  group_by(sample, chromosome) %>%
  summarize(total = sum(numbers))
plotable_indels.df$type <- "indel"

## concat dataframes
plotable.df <- rbind(plotable_SNPs.df, plotable_indels.df)
## reorder chromosomes
plotable.df$chromosome <- factor(plotable.df$chromosome, levels = c(paste0("chr",seq(1,22)),"chrX","chrY"))
## reorder type of variant
plotable.df$type <- factor(plotable.df$type, levels = c("SNP","indel"))

# pre-calculate axis values
max_y <- max(plotable.df$total)
axis_breaks_y <- round(seq(0,max_y, by = max_y/10))   

plot.tmp <- ggplot(data = plotable.df, aes(x=chromosome, y = total)) +
  geom_boxplot(aes(fill=type)) +
  ggtitle(label = "Total variants per chromosome", subtitle = args[1]) +
  scale_y_continuous(name = paste0("total variants"),
                     limits = c(0, max_y*1.1),
                     breaks = axis_breaks_y,
                     labels = axis_breaks_y) +
  scale_fill_manual(values = viridis(16)[c(4,8)]) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 10),
        plot.title = element_text(size=10),
        plot.subtitle = element_text(size=5))

print(plot.tmp)

dev.off() ## close PDF device
