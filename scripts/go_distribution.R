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
plot(density(human_go$V4[human_go$V4 > 0]))

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

hist(log10(blast_scores$V3[!is.na(blast_scores$V3)]), breaks=100)
hist(log10(psi_scores$V3[!is.na(psi_scores$V3)]), breaks=100)

psi_scores<-!is.na(psi_scores)
min(psi_scores$V3)


#!/usr/bin/python