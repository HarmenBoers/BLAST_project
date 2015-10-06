# Parse arguments
args = commandArgs(TRUE)
blastfile = normalizePath(args[1])
psiblastfile = normalizePath(args[2])
directory = dirname(blastfile)
setwd(directory)

###BLAST Scores
blast_scores <- read.table(blastfile)
psi_scores <- read.table(psiblastfile)

#PLOT AUC on the evalue cutoffs that we compared - to see the difference
pdf("./blast_evalue_cutoffs.pdf", height=3)
par(mfrow=c(1,2))
hist_blast <- hist(log10(blast_scores$V3[!is.na(blast_scores$V3) & blast_scores$V3 < 10]), breaks=30, 
                   main="a) BLAST log10 e-values < 10", cex.main=.9,
                   xlab="log10 evalue",
                   ylab="Frequency",
                   ylim=c(0,200))
hist_psiblast <-hist(log10(psi_scores$V3[!is.na(psi_scores$V3) & psi_scores$V3 <10]), breaks=30, 
                     main="b) PSI-BLAST log10 e-values < 10", cex.main=.9,
                     xlab="log10 evalue",
                     ylab="Frequency")
dev.off()