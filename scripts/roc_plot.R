# Parse arguments
args = commandArgs(TRUE)
blastfile = normalizePath(args[1])
psiblastfile = normalizePath(args[2])
directory = dirname(blastfile)
setwd(directory)

# Create and save plot
psi.coor <- read.table(psiblastfile, header=T)
blast.coor <- read.table(blastfile, header=T)
pdf("./roc_plot_blast.pdf")
plot(psi.coor$x, psi.coor$y, type="l", col="red",
     main="ROC plot of PSI-BLAST and BLAST",
     ylab="Coverage",
     xlab="Error")
lines(blast.coor$x, blast.coor$y, type="l", col="blue")
abline(0,1)
legend("topleft",legend=c("PSI-BLAST", "BLAST"), lty=c(1,1), lwd=c(2.5,2.5),col=c("red"," blue"))
dev.off()
