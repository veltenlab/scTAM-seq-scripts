library(argparse)
ap <- ArgumentParser()
ap$add_argument("-f", "--file", action="store", help="The barcode distribution file as output of the MissionBio pipeline (tsv)")
ap$add_argument("-a", "--ampli", action="store", help="A file containing information about the amplicons (tsv)")
ap$add_argument("-o", "--output", action="store", help="The output file name")
ap$add_argument("-c", "--cutoff", action="store", type="double", default=0.6, help="Integer number specifiying the number of amplicons that need to fulfill the coverage threshold")
opt <- ap$parse_args()

# V1 tends to perform better for us
vers.2 <- FALSE
val <- as.numeric(opt$cutoff)
resemble.cellfinder <- function(d){
  cut <- ifelse(vers.2,ncol(d) * 8,ncol(d) * 10)
  reads.per.cell <- apply(d,1,sum)
  print(paste("Removing",sum(reads.per.cell<cut),"cells due to amplicon cutoff"))
  d <- d[reads.per.cell>=cut,]
  m.ean <- mean(apply(d,1,mean))
  thres <- ifelse(vers.2,0.2*m.ean,min(0.2*m.ean,10))
}
ampli.info <- read.table(opt$ampli,header=TRUE)
sel.amplis <- row.names(ampli.info)[grepl('NonHhaI',ampli.info$Type)&ampli.info$good.performance]
dat.s3 <- read.table(opt$file,
                     header = T)
thres <- resemble.cellfinder(dat.s3[,sel.amplis])
ful.cutoff <- dat.s3[,sel.amplis]>thres
ful.perc <- apply(ful.cutoff,1,function(x)sum(x)/length(x))
ful.perc.s3 <- ful.perc
dat.s3 <- dat.s3[ful.perc.s3>val,]
write.table(dat.s3,opt$output,sep='\t',quote = FALSE)
