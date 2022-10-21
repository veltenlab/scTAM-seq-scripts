####################### Supp_Fig_8.R ####################### 
#' this file determines why the performance of the always methylated
#' amplicons is the worst

library(ggplot2)
plot_theme <- theme(panel.background = element_rect(color='black',fill='white'),
                    panel.grid=element_blank(),
                    text=element_text(color='black',size=6),
                    axis.text=element_text(color='black',size=5),
                    axis.ticks=element_line(color='black', size=.25),
                    axis.text.x=element_text(angle=90, hjust=1, vjust=1),
                    strip.background = element_blank(),
                    strip.text.x = element_blank(),
                    legend.key=element_rect(color=NA, fill=NA),
                    legend.position='none')
plot_path <- '~'
type_clear <- c("Other"="Other",
                "NonHhaI"="No HhaI cutsite",
                "CpG.B.cell.diff"="Differential CpG Bcells",
                "CpG.always.unmeth.B"="Always unmethylated CpG in Bcells",
                "CpG.always.meth.B"="Always methylated CpG in Bcells",
                "CpG.imprinted.multiple.cutsites"="Imprinted CpG",
                "CpG.imprinted"="Imprinted CpG")
colors_amplicons <- c("Other"="#e5c494",
                      "Differential CpG Bcells"="#fc8d62",
                      "Always unmethylated CpG in Bcells"="#a6d854",
                      "Always methylated CpG in Bcells"="#66c2a5",
                      "Imprinted CpG"="#e78ac3",
                      "Imprinted CpG multiple"="#e78ac3",
                      "No HhaI cutsite"="#b3b3b3")
ampli.info <- read.csv('../../misc/amplicons_CpG_GC_annotated.csv', row.names=1)
ampli.info$Type <- type_clear[as.character(ampli.info$Type)]
ampli.info <- ampli.info[!is.na(ampli.info$Type), ]
plot <- ggplot(ampli.info, aes(x=Type, y=GC_Content, fill=Type))+geom_boxplot()+plot_theme+
  scale_fill_manual(values=colors_amplicons)
ggsave(file.path(plot_path, 'amplicon_performance_GCContent.pdf'),
       width=100, height=75, unit='mm')

plot <- ggplot(ampli.info, aes(x=Type, y=CpG_count, fill=Type))+geom_boxplot()+plot_theme+
  scale_fill_manual(values=colors_amplicons)
ggsave(file.path(plot_path, 'amplicon_performance_CpG_count.pdf'),
       width=100, height=75, unit='mm')

sample <- 'GSM5935921_BM_HhaI'
dat <- read.table(paste0('../../data/', sample, '.tsv.gz'),
                  header=TRUE)
rowinfo <- read.csv(paste0('../../misc/', sample, '/tsv/rowinfo.csv'),
                    row.names=1)
selected_data <- dat[row.names(rowinfo), row.names(ampli.info[ampli.info$Type%in%'Always methylated CpG in Bcells', ])]
dropout <- apply(selected_data, 2, function(x){
  sum(x==0)/length(x)
})
to_plot <- data.frame(GC_Content=ampli.info[ampli.info$Type%in%'Always methylated CpG in Bcells', 'GC_content'],
                      Type=ampli.info[ampli.info$Type%in%'Always methylated CpG in Bcells', 'Type'],
                      Dropout=dropout)
plot <- ggplot(to_plot, aes(x=GC_Content, y=Dropout, color=Type))+geom_point()+geom_abline(slope=1, intercept = 0)+plot_theme+
  xlim(0, 1)+ylim(0,1)+scale_color_manual(values=colors_amplicons)
ggsave(file.path(plot_path, 'dropout_GCContent_methylated.pdf'),
       width=100, height=75, unit='mm')

to_plot <- data.frame(CpG_count=ampli.info[ampli.info$Type%in%'Always methylated CpG in Bcells', 'CpG_count'],
                      Type=ampli.info[ampli.info$Type%in%'Always methylated CpG in Bcells', 'Type'],
                      Dropout=dropout)
plot <- ggplot(to_plot, aes(x=CpG_count, y=Dropout, color=Type))+geom_point()+geom_abline(slope=1, intercept = 0)+plot_theme+
  ylim(0,1)+geom_smooth(method='lm', se=FALSE)+scale_color_manual(values=colors_amplicons)
ggsave(file.path(plot_path, 'dropout_CpG_count_methylated.pdf'),
       width=100, height=75, unit='mm')

selected_data <- dat[row.names(rowinfo), row.names(ampli.info[ampli.info$Type%in%'No HhaI cutsite', ])]
dropout <- apply(selected_data, 2, function(x){
  sum(x==0)/length(x)
})
to_plot <- data.frame(GC_Content=ampli.info[ampli.info$Type%in%'No HhaI cutsite', 'GC_content'],
                      Type=ampli.info[ampli.info$Type%in%'No HhaI cutsite', 'Type'],
                      Dropout=dropout)
plot <- ggplot(to_plot, aes(x=GC_Content, y=Dropout, color=Type))+geom_point()+geom_abline(slope=1, intercept = 0)+plot_theme+
  xlim(0, 1)+ylim(0,1)+scale_color_manual(values=colors_amplicons)
ggsave(file.path(plot_path, 'dropout_GCContent_nonHhaI.pdf'),
       width=100, height=75, unit='mm')
