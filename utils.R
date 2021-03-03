# Helper functions
library(sumsamstats)

# Takes in file name, returns a list with 7 items, the seq name, reference name, %mapped, and 3 data frames containing the outputs of samtools idxstats, and samtools stats (COV and SN).
read_stats <- function(fname, dir.stats) {
  
  # file names for each commands output (idxstats, and stats)
  fname.stats <- paste(fname, "sorted.bam.stats_raw.tsv", sep = ".")
  fname.idx <- paste(fname, "sorted.bam.reads_per_chunk.tsv", sep = ".")
  
  # Parse info
  s <- unlist(strsplit(fname, ".", fixed = TRUE))
  ref <- paste(s[4], s[5], sep = ".")
  seq <- s[1]
  method <- s[3]

  # idxstats
  df1 <- read.delim(paste(dir.stats, fname.idx, sep = "/"), header = FALSE)
  names(df1) <- c("name", "length", "mapped", "unmapped")
  
  # Calculate the % of reads from idxstats.
  total <- sum(df1[, 4]) + sum(df1[, 3])
  reads.idx <- 100*(sum(df1[, 3]) / total)

  # stats
  stats <- readSamtoolsStats(paste(dir.stats, fname.stats, sep = "/"), section = c("SN", "COV"))
  total <- as.numeric(stats$SN[1, 2])
  reads.stats <- 100*(as.numeric(stats$SN[7,2]) / total)
  
  # the %mapped found from idxstats and stats are identical
  
  # Return as list
  return(list(ref = ref, seq = seq, method = method, reads = reads.stats, data.idx = df1, data.SN = stats$SN, data.COV = stats$COV))
}

# Reformat for ggplot2
t_ggp <- function(datalist) {
  df <- data.frame()
  for (x in datalist) {
    met <- x$method
    ref <- x$ref
    seq <- x$seq
    reads <- x$reads
    map <- as.numeric(as.vector(x$data.idx$mapped))
    
    df <- rbind(df, data.frame(ref, seq, met, reads, map))
  }
  return(df)
}

# List indexing operator shortcut
'%l%' <- function(x, y) {
  sapply(x, "[[", y)
}
