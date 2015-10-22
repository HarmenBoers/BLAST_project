library(ggplot2)
df = data.frame(Method = c(rep('BLAST', 12), rep('PSI-BLAST', 12)),
  Database = c(rep('GO', 4), rep('SCOP', 4), rep('Pfam', 4)),
  AUC = c(0.676, 0.740, 0.778, 0.689, 0.860, 0.859, 0.760, 0.860, 0.929, 0.932, 0.949, 0.928, 0.679, 0.750, 0.791, 0.685, 0.740, 0.889, 0.780, 0.850, 0.931, 0.932, 0.95, 0.93))
ggplot(df, aes(factor(Database, levels = c('GO', 'SCOP', 'Pfam')), AUC)) +
  geom_boxplot(aes(fill = Method)) +
  ylim(0.6, 1) +
  xlab('Database') +
  ylab('AUC') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  theme(text = element_text(size = 25))
ggsave('./compare_methods.png', width=10, height=6)