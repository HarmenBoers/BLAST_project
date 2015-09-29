library(ggplot2)

pdf("/Users/harmen/PycharmProjects/blast_project/results/roc_plot_psiblast.pdf")
########### PSI-BLAST VS GO ########## AUC = 0.91
psi.coor <- read.table("/Users/harmen/PycharmProjects/blast_project/psiblast_rocplot.txt", header=T)
blast.coor <- read.table("/Users/harmen/PycharmProjects/blast_project/blast_rocplot.txt", header=T)

plot(psi.coor$x, psi.coor$y, type="l", col="red",
     main="ROC plot of PSI-BLAST and BLAST")
lines(blast.coor$x, blast.coor$y, type="l", col="blue")
abline(0,1)
legend("topleft",legend=c("PSI-BLAST", "BLAST"), lty=c(1,1), lwd=c(2.5,2.5),col=c("red"," blue"))
#PSI BLAST  AUC: 0.679
#BLAST  AUC: 0.676

roc <-ggplot(data=coor, aes(x=x, y=y, group=1)) +
      geom_line() +
      geom_abline() +
      expand_limits(y=0) +
      xlab("False Positive Rate") + ylab("True Positive Rate") +
      ggtitle("PSI-BLAST benchmarked using GO AUC = 0.91")
plot(roc)
dev.off()
pdf("/Users/harmen/PycharmProjects/blast_project/results/roc_plot_blast.pdf")
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

