library(ggplot2)

args = commandArgs(TRUE)
datafile = normalizePath(args[1])
directory = dirname(datafile)
setwd(directory)

go_data = read.table(datafile, sep='\t')
colnames(go_data) = c('p1', 'p2', 'homology', 'score')
quant = quantile(go_data$score, c(0.95, 0.99))
write.table(quant, file = './go_score_distribution_quantiles.txt')

ggplot(go_data, aes(score)) +
  geom_density() +
  xlab("GO Score") +
  ylab("Density") +
  geom_vline(xintercept = quant)
ggsave('./go_density.png')
