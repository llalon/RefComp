# Libs
library(ggplot2)
library(ggpubr)

# Global vars
ref.bb <- "GCA_900302385.1_ASM90030238v1_genomic"
ref.cod <- "GCF_902167405.1_gadMor3"
refs <- c(ref.bb, ref.cod)
method.bwa <- "bwa"
method.bt <- "bowtie2"
methods <- c(method.bwa, method.bt)

# Project directories
dir.root <- "C:/Users/liaml/git/toastedfrog"
dir.stats <- paste(dir.root, "stats", sep = "/")

# Set wd to root dir
setwd(dir.root)

# Source helper functions
source("utils.R")

# Read files
files.aln <- list.files(path = dir.stats, pattern = "*stats_raw.tsv", full.names = FALSE, recursive = FALSE)
# remove extensions
files.aln.cleaned <- unlist(lapply(files.aln, function(x) gsub(".sorted.bam.stats_raw.tsv", "", x)))

# Build list of stats
sam.stats <- lapply(files.aln.cleaned, function(x) read_stats(x, dir.stats))

# Change names to cleaned file names
files.aln.cleaned <- unlist(lapply(files.aln, function(x) gsub(".sorted.bam.reads_per_chunk.tsv", "", x)))
names(sam.stats) <- files.aln.cleaned

# Names without methods
names.no.methods <- gsub("*.bowtie2", "", files.aln.cleaned)
names.no.methods <- names.no.methods[-grep("bwa", names.no.methods)]

# Split the data into combinations of each method/reference
stats.all <- sam.stats
stats.bwa <- sam.stats[which(sam.stats %l% 3 == method.bwa)]
stats.bt <- sam.stats[which(sam.stats %l% 3 == method.bt)]
stats.cod <- sam.stats[which(sam.stats %l% 1 == ref.cod)]
stats.bb <- sam.stats[which(sam.stats %l% 1 == ref.bb)]
stats.bwa.bb <- stats.bb[which(stats.bb %l% 3 == method.bwa)]
stats.bwa.cod <- stats.cod[which(stats.cod %l% 3 == method.bwa)]
stats.bt.bb <- stats.bb[which(stats.bb %l% 3 == method.bt)]
stats.bt.cod <- stats.cod[which(stats.cod %l% 3 == method.bt)]

# Comparing %mapped between each method
map.bt <- data.frame(unname(stats.bt %l% 4))
map.bwa <- data.frame(unname(stats.bwa %l% 4))

names(map.bt) <- c("mapped")
names(map.bwa) <- c("mapped")

map.bt$method = method.bt
map.bwa$method = method.bwa

# Comparing %mapped between each reference
map.cod <- data.frame(unname(stats.cod %l% 4))
map.bb <- data.frame(unname(stats.bb %l% 4))

names(map.cod) <- c("mapped")
names(map.bb) <- c("mapped")

map.cod$reference = ref.cod
map.bb$reference = ref.bb

# Built the plots

# Stacked histograms of each comparison
hista <- ggplot(rbind(map.bt, map.bwa), aes(mapped, fill = method)) +
  geom_histogram(binwidth = 1) +
  ylab("Frequency") +
  xlab("% Mapped") +
  ggtitle("A.") +
  labs(fill = 'Alignment Method')

histb <- ggplot(rbind(map.cod, map.bb), aes(mapped, fill = reference)) +
  geom_histogram(binwidth = 1) +
  ylab("Frequency") +
  xlab("% Mapped") +
  ggtitle("B.") +
  labs(fill = 'Reference Genome')

# Violin plots showing differences between each sample. Split between reference, and method. Y axis logged to normalize range.

# bwa/burbot
pa <- ggplot(t_ggp(stats.bwa.bb), aes(x = seq, y = log(map), color = seq)) +
  geom_violin() +
  scale_x_discrete(name = "Sample Sequence", labels = c(1:10)) +
  ylab("log(% mapped)") +
  ggtitle("A. bwa/GCA_900302385.1_ASM90030238v1_genomic") +
  labs(color = 'Sample Sequence (1-10)')

# bwa/cod
pb <- ggplot(t_ggp(stats.bwa.cod), aes(x = seq, y = log(map), color = seq)) +
  geom_violin() +
  scale_x_discrete(name = "Sample Sequence", labels = c(1:10)) +
  ylab("log(% mapped)") +
  ggtitle("B. bwa/GCF_902167405.1_gadMor3") +
  labs(color = 'Sample Sequence (1-10)')

# bowtie2/burbot
pc <- ggplot(t_ggp(stats.bt.bb), aes(x = seq, y = log(map), color = seq)) +
  geom_violin() +
  scale_x_discrete(name = "Sample Sequence", labels = c(1:10)) +
  ylab("log(% mapped)") +
  ggtitle("C. bowtie2/GCA_900302385.1_ASM90030238v1_genomic") +
  labs(color = 'Sample Sequence (1-10)')

# bowtie2/cod
pd <- ggplot(t_ggp(stats.bt.cod), aes(x = seq, y = log(map), color = seq)) +
  geom_violin() +
  scale_x_discrete(name = "Sample Sequence", labels = c(1:10)) +
  ylab("log(% mapped)") +
  ggtitle("D. bowtie2/GCF_902167405.1_gadMor3") +
  labs(color = 'Sample Sequence (1-10)')

# Arrange plots in grid with common legend and axis
vgrid <- ggarrange(pa + theme(axis.title.x = element_blank()),
                   pb + theme(axis.title.x = element_blank(), axis.title.y = element_blank()),
                   pc,
                   pd + theme(axis.title.y = element_blank()),
                   ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom")


t.test(as.vector(stats.bwa %l% 4), as.vector(stats.bt %l% 4))
t.test(as.vector(stats.bb %l% 4), as.vector(stats.cod %l% 4))

df.summ <- do.call(rbind.data.frame, list(
  bwa = summary(as.vector(stats.bwa %l% 4)),
  bowtie2 = summary(as.vector(stats.bt %l% 4)),
  cod = summary(as.vector(stats.cod %l% 4)),
  burbot = summary(as.vector(stats.bb %l% 4))
))
row.names(df.summ) <- c("bwa", "bowtie2", "cod", "burbot")
names(df.summ) <- c("Min", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max")