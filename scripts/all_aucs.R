library(ggplot2)  

# Parse arguments
args = commandArgs(TRUE)
filepath = normalizePath(args[1])
filepath = '/Users/oliviermartin/Dropbox/VU/Fundamentals of Bioinformatics/Project/blast_project/results/all_aucs.out'
directory = dirname(filepath)
setwd(directory)

# Plot
df = read.table(file = filepath, header = FALSE, sep = ',', dec='.')
colnames(df) = c('t1', 't2', 'AUC')

ggplot(df, aes(t1, t2, fill = AUC)) + 
  geom_raster() + 
  scale_fill_gradient2(low = "blue", 
                       mid = "white", 
                       high = "red", 
                       midpoint = 0.75, 
                       limits = c(0.5, 1)) +
  xlab('Lower GO Threshold') +
  ylab('Higher GO Threshold') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  theme(text = element_text(size = 25)) +
  geom_vline(xintercept = 0.14) +
  geom_hline(yintercept = 0.27)
ggsave('./all_aucs.png')

