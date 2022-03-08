suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg19))
genome <- BSgenome.Hsapiens.UCSC.hg19
chrs <- c(paste0('chr', 1:22), 'chrX', 'chrY')
res <- sapply(chrs, function(x){
  chr <- genome[[x]]
  res.extended <- matchPattern(DNAString('GCGC'), chr)
  length(res.extended)
})