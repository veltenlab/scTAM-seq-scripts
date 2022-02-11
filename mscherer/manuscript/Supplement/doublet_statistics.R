####################### doublet_statistics.R ####################### 
#' With this script, I simply extract the information from the doublet detection files to 
#' include them as supplementary information.

doublet_s6 <- read.csv('/users/mscherer/cluster/project/Methylome/analysis/missionbio/re_sequencing/Sample6_70_percent_good_performance/tsv/doublet_scores_DoubletDetection.csv')
doublet_s7 <- read.csv('/users/mscherer/cluster/project/Methylome/analysis/missionbio/re_sequencing/Sample7_70_percent_good_performance/tsv/all_doublets.csv')
doublet_s8 <- read.table('/users/mscherer/cluster/project/Methylome/analysis/missionbio/re_sequencing/Sample8_70_percent_good_performance/tsv/all_doublets.tsv')
dim(doublet_s6)
dim(doublet_s7)
dim(doublet_s8)
length(which(doublet_s6$DoubletDetectionLabel==1))
length(which(doublet_s7$DoubletDetectionLabel==1))
length(which(doublet_s8$DoubletDetectionLabel==1))
length(which(doublet_s7$DoubletDetectionLabel==0&doublet_s7$Doublet==1))
length(which(doublet_s8$DoubletDetectionLabel==0&doublet_s8$Doublet==1))
rowinfo_s6 <- read.csv('/users/mscherer/cluster/project/Methylome/analysis/missionbio/re_sequencing/Sample6_70_percent_good_performance/tsv/rowinfo.csv')
rowinfo_s7 <- read.csv('/users/mscherer/cluster/project/Methylome/analysis/missionbio/re_sequencing/Sample7_70_percent_good_performance/tsv/rowinfo.csv')
rowinfo_s8 <- read.csv('/users/mscherer/cluster/project/Methylome/analysis/missionbio/re_sequencing/Sample8_70_percent_good_performance/tsv/rowinfo.csv')
dim(rowinfo_s6)
dim(rowinfo_s7)
dim(rowinfo_s8)
