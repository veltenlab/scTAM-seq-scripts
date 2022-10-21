#### F2_E_chromatin_states.R ###################################
#' With this file, you can generate the chromatin state pie charts for the different CpG
#' cluster and their corresponding significance tests

#packages
require(tidyverse)
require(ggplot2)
require(ggfortify)

# read data
load("../../misc/colinfo_diff_cells.chrom.states.RData")

# read labels
colors <- read.table("../../misc/Colors_ChromHMM.txt",
                     sep = "\t",
                     header=T)

Background <- CpGs.states.2 %>% select(4:18) %>% 
  sapply(table) %>% t() %>% 
  as.data.frame() %>% 
  transmute(E1E3 = E1 + E3,
            E2E6 = E6 + E6,
            E7E8E9 = E7 + E8 + E9,
            E4 = E4,
            E5 = E5,
            E10 = E10,
            E11 = E11,
            E12 = E12) %>% 
  colSums()/15

Cluster1 <- CpGs.states.2 %>%
  filter(cluster == 1) %>% 
  select(4:18) %>%
  sapply(table) %>% t() %>% 
  as.data.frame() %>% 
  transmute(E1E3 = E1 + E3,
            E2E6 = E6 + E6,
            E7E8E9 = E7 + E8 + E9,
            E4 = E4,
            E5 = E5,
            E10 = E10,
            E11 = E11,
            E12 = E12) %>% 
  colSums()/15

Cluster2 <- CpGs.states.2 %>%
  filter(cluster == 2) %>% 
  select(4:18) %>%
  sapply(table) %>% t() %>% 
  as.data.frame() %>% 
  transmute(E1E3 = E1 + E3,
            E2E6 = E6 + E6,
            E7E8E9 = E7 + E8 + E9,
            E4 = E4,
            E5 = E5,
            E10 = E10,
            E11 = E11,
            E12 = E12) %>% 
  colSums()/15

Cluster3 <- CpGs.states.2 %>%
  filter(cluster == 3) %>% 
  select(4:18) %>%
  sapply(table) %>% t() %>% 
  as.data.frame() %>% 
  transmute(E1E3 = E1 + E3,
            E2E6 = E6 + E6,
            E7E8E9 = E7 + E8 + E9,
            E4 = E4,
            E5 = E5,
            E10 = E10,
            E11 = E11,
            E12 = E12) %>% 
  colSums()/15

Cluster4 <- CpGs.states.2 %>%
  filter(cluster == 4) %>% 
  select(4:18) %>%
  sapply(table) %>% t() %>% 
  as.data.frame() %>% 
  transmute(E1E3 = E1 + E3,
            E2E6 = E6 + E6,
            E7E8E9 = E7 + E8 + E9,
            E4 = E4,
            E5 = E5,
            E10 = E10,
            E11 = E11,
            E12 = E12) %>% 
  colSums()/15

Cluster5 <- CpGs.states.2 %>%
  filter(cluster == 5) %>% 
  select(4:18) %>%
  sapply(table) %>% t() %>% 
  as.data.frame() %>% 
  transmute(E1E3 = E1 + E3,
            E2E6 = E6 + E6,
            E7E8E9 = E7 + E8 + E9,
            E4 = E4,
            E5 = E5,
            E10 = E10,
            E11 = E11,
            E12 = E12) %>% 
  colSums()/15

Cluster6 <- CpGs.states.2 %>%
  filter(cluster == 6) %>% 
  select(4:18) %>%
  sapply(table) %>% t() %>% 
  as.data.frame() %>% 
  transmute(E1E3 = E1 + E3,
            E2E6 = E6 + E6,
            E7E8E9 = E7 + E8 + E9,
            E4 = E4,
            E5 = E5,
            E10 = E10,
            E11 = E11,
            E12 = E12) %>% 
  colSums()/15


chr.states.plot <- rbind(Background,
                         Cluster1, 
                         Cluster2,
                         Cluster3,
                         Cluster4,
                         Cluster5,
                         Cluster6) %>%
  t()%>% 
  reshape2::melt() 

plot_theme <- theme(panel.background = element_blank(),
                    panel.grid=element_blank(),
                    text=element_text(color='black',size=6),
                    axis.text=element_text(color='black',size=5),
                    axis.ticks=element_blank(),
                    strip.background = element_blank(),
                    legend.key=element_rect(color=NA, fill=NA),
                    axis.text.x=element_blank(),
                    axis.title.x=element_blank(),
                    axis.ticks.x=element_blank(),
                    legend.title=element_blank(),
                    panel.spacing = unit(.1, "lines"),
                    legend.position = 'none',
                    title = element_blank())

chr.states.plot$Var1 <- c('E1E3'='active/weak promoter',
                          'E2E6'='strong enhancer',
                          'E7E8E9'='transcription',
                          'E4'='poised promoter',
                          'E5'='weak enhancer',
                          'E10'='H3K9me3-repressed',
                          'E11'='Heterochromatin-low',
                          'E12'='H3K9me3-repressed')[as.character(chr.states.plot$Var1)]
cols <- c('active/weak promoter'='#ff678c',
          'strong enhancer'='#ffdc64',
          'transcription'='#008c64',
          'poised promoter'='#6e1e8c',
          'weak enhancer'='#ffff00',
          'H3K9me3-repressed'='#787878',
          'Heterochromatin-low'='#aaaaaa',
          'H3K9me3-repressed'='#f0f0f0')
to_plot <- subset(chr.states.plot, subset = Var2=='Background')
plot <- ggplot(to_plot , 
               aes(x="", y=value, fill=Var1))+
  geom_bar(stat = "identity",)+xlab('')+ylab('')+ggtitle('')+
  coord_polar("y", start=0) +plot_theme+scale_fill_manual(values=cols)+scale_color_manual(values=cols)
ggsave('Background.pdf', plot, width = 25, height = 25, unit='mm')

to_plot <- subset(chr.states.plot, subset = Var2=='Cluster1')
plot <- ggplot(to_plot , 
               aes(x="", y=value, fill=Var1))+
  geom_bar(stat = "identity",)+xlab('')+ylab('')+ggtitle('')+
  coord_polar("y", start=0) +plot_theme+scale_fill_manual(values=cols)+scale_color_manual(values=cols)
ggsave('Cluster1.pdf', plot, width = 25, height = 25, unit='mm')

to_plot <- subset(chr.states.plot, subset = Var2=='Cluster2')
plot <- ggplot(to_plot , 
               aes(x="", y=value, fill=Var1))+
  geom_bar(stat = "identity",)+xlab('')+ylab('')+ggtitle('')+
  coord_polar("y", start=0) +plot_theme+scale_fill_manual(values=cols)+scale_color_manual(values=cols)
ggsave('Cluster2.pdf', plot, width = 25, height = 25, unit='mm')

to_plot <- subset(chr.states.plot, subset = Var2=='Cluster3')
plot <- ggplot(to_plot , 
               aes(x="", y=value, fill=Var1))+
  geom_bar(stat = "identity",)+xlab('')+ylab('')+ggtitle('')+
  coord_polar("y", start=0) +plot_theme+scale_fill_manual(values=cols)+scale_color_manual(values=cols)
ggsave('Cluster3.pdf', plot, width = 25, height = 25, unit='mm')

to_plot <- subset(chr.states.plot, subset = Var2=='Cluster4')
plot <- ggplot(to_plot , 
               aes(x="", y=value, fill=Var1))+
  geom_bar(stat = "identity",)+xlab('')+ylab('')+ggtitle('')+
  coord_polar("y", start=0) +plot_theme+scale_fill_manual(values=cols)+scale_color_manual(values=cols)
ggsave('Cluster4.pdf', plot, width = 25, height = 25, unit='mm')

to_plot <- subset(chr.states.plot, subset = Var2=='Cluster5')
plot <- ggplot(to_plot , 
               aes(x="", y=value, fill=Var1))+
  geom_bar(stat = "identity",)+xlab('')+ylab('')+ggtitle('')+
  coord_polar("y", start=0) +plot_theme+scale_fill_manual(values=cols)+scale_color_manual(values=cols)
ggsave('Cluster5.pdf', plot, width = 25, height = 25, unit='mm')

to_plot <- subset(chr.states.plot, subset = Var2=='Cluster6')
plot <- ggplot(to_plot , 
               aes(x="", y=value, fill=Var1))+
  geom_bar(stat = "identity",)+xlab('')+ylab('')+ggtitle('')+
  coord_polar("y", start=0) +plot_theme+scale_fill_manual(values=cols)+scale_color_manual(values=cols)
ggsave('Cluster6.pdf', plot, width = 25, height = 25, unit='mm')

consensus <- as.data.frame(t(apply(CpGs.states.2, 1, function(x){
  vec <- x[-(1:3)]
  vec <- c('E1'='E1E3',
           'E3'='E1E3',
           'E2'='E2E6',
           'E6'='E2E6',
           'E7' = 'E7E8E9',
           'E8' = 'E7E8E9',
           'E9' = 'E7E8E9',
           'E4' = 'E4',
           'E5' = 'E5',
           'E10' = 'E10',
           'E11' = 'E11',
           'E12' = 'E12')[vec]
  cu <- plyr::count()
  cu <- cu$x[order(cu$freq, decreasing = TRUE)][1]
  c(x[1:3], cu)
})))
p.vals.cluster <- sapply(unique(consensus$cluster), function(clust){
  fr <- consensus[consensus$cluster==clust, ]
  bg <- consensus[consensus$cluster!=clust, ]
  p.vals.state <- sapply(unique(fr$V4), function(chr){
    tp <- sum(fr$V4==chr)
    fp <- sum(bg$V4==chr)
    tn <- sum(bg$V4!=chr)
    fn <- sum(fr$V4!=chr)
    p.vals <- fisher.test(matrix(c(tp, fp, fn, tn), ncol=2))$p.value
  })
})
