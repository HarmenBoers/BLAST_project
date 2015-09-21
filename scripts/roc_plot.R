library(ggplot2)

pdf("/Users/harmen/PycharmProjects/blast_project/results/roc_plots.pdf")
########### PSI-BLAST VS GO ########## AUC = 0.91
coor <- read.table("/Users/harmen/PycharmProjects/blast_project/roc_plot.txt", header=T)

roc <-ggplot(data=coor, aes(x=x, y=y, group=1)) +
      geom_line() +
      geom_abline() +
      expand_limits(y=0) +
      xlab("False Positive Rate") + ylab("True Positive Rate") +
      ggtitle("PSI-BLAST benchmarked using GO AUC = 0.91")
plot(roc)


########## BLAST VS GO ############ AUC = 0.86
coor <- read.table("/Users/harmen/PycharmProjects/blast_project/blast_rocplot.txt", header=T)

roc <-ggplot(data=coor, aes(x=x, y=y, group=1)) +
  geom_line() +
  geom_abline() +
  expand_limits(y=0) +
  xlab("False Positive Rate") + ylab("True Positive Rate") +
  ggtitle("BLAST benchmarked using GO AUC = 0.86")
plot(roc)
dev.off()

