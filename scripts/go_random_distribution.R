go_scores <- read.table("/Users/harmen/PycharmProjects/blast_project/scripts/test_go.out", sep="\t")

plot(density(go_scores$V4[go_scores$V4 > 0]))
