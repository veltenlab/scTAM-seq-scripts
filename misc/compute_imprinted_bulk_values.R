library(RnBeads)
load('/users/mscherer/mount/abianchi/WGBS.data.bulk.Renee/WGBS_57_10reads_NAs_grch38.RData')
map <- read.table('/users/mscherer/mount/abianchi/WGBS.data.bulk.Renee/CpGs.wgbs.hg38.blood.B.cell.panel.mission.bio.bed')
ampli_inf <- read.table('/users/lvelten/project/Methylome/infos/BCells/Blood.Bone.Marrow.Amplicons.design.dropout.added.selected.tsv')
map <- map[grepl('imprinted', map$V4), ]
map.gr <- makeGRangesFromDataFrame(map, seqnames.field='V1', start.field='V2', end.field='V3')
dat.gr <- GRanges(meth.final$Contig, IRanges(start=meth.final$Pos, end=meth.final$Pos))
op <- findOverlaps(map.gr, dat.gr)
sel_dat <- meth.final[subjectHits(op), ]
ampli_inf$CpG_Names <- sapply(strsplit(ampli_inf$CpG_Names, ','), function(x)x[1])
ampli_inf <- ampli_inf[match(map$V4, ampli_inf$CpG_Names), ]

