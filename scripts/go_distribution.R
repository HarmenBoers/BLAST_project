go_scores <- read.table("/Users/harmen/PycharmProjects/blast_project/scripts/test_go.out", sep="\t")

plot(density(go_scores$V4[go_scores$V4 > 0]))
#First threshold
quantile(go_scores$V4,probs=0.95)
#Second threshold
quantile(go_scores$V4,probs=0.99)

human_uniprot <- read.table("/Users/harmen/PycharmProjects/blast_project/reviewed_human_uniprot.txt")

random_human <- sample(human_uniprot$V1, 1000)

write(as.vector(random_human), file="/Users/harmen/PycharmProjects/blast_project/scripts/random_1000_human_ids.txt")

human_go <- read.table("/Users/harmen/PycharmProjects/blast_project/results/random_go/random_human_go.out")
plot(density(human_go$V4[human_go$V4 > 0]),
     title="Density of GO terms in 1000 random proteins")

#First threshold
quantile(human_go$V4,probs=0.95)
#Second threshold
quantile(human_go$V4,probs=0.99)



###BLAST Scores
blast_scores <- read.table("/Users/harmen/PycharmProjects/blast_project/results/blast_results.out")
psi_scores <- read.table("/Users/harmen/PycharmProjects/blast_project/results/psi_blast.out")
plot(density(blast_scores$V3[blast_scores$V3 != "NA"]))
blast_scores$V3
plot(density(blast_scores$V3[!is.na(blast_scores$V3)]))


#PLOT AUC on the evalue cutoffs that we compared - to see the difference
pdf("/Users/harmen/PycharmProjects/blast_project/results/blast_evalue_cutoffs.pdf", height=6)
hist(log10(blast_scores$V3[!is.na(blast_scores$V3)]), breaks=30, ylim=c(0,600), 
     main="Histogram of log10 BLAST E-values",
     xlab="log10 evalue",
     ylab="Frequency")
arrows(0,0,0,550)
text(0, 575, "0.87", col = "red", cex=0.8)
arrows(-15  ,0,-15,550)
text(-15, 575, "0.80", col = "red", cex=0.8)
arrows(-30  ,0,-30,550)
text(-30, 575, "0.89", col = "red", cex=0.8)
arrows(-45  ,0,-45,550)
text(-45, 575, "0.54", col = "red", cex=0.8)
text(-60,580,"AUC:", col="red", cex=0.9)
dev.off()
?arrows
#hist(log10(psi_scores$V3[!is.na(psi_scores$V3)]), breaks=30)
plot(density(log10(blast_scores$V3[!is.na(blast_scores$V3)])))
psi_scores<-!is.na(psi_scores)
min(psi_scores$V3)
abline(v=-4, )

######### FALSE POSITIVES ###########
fp <- read.table("/Users/harmen/PycharmProjects/blast_project/false_positives.txt")
hist(log10(fp$V2))
