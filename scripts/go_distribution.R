# Parse arguments
args = commandArgs(TRUE)
randomfile = normalizePath(args[1])
directory = dirname(randomfile)
setwd(directory)

#Assign the random selection GO similarity score output to a matrix to make a plot from.
random_go <- read.table(randomfile)
#First threshold
one <- unname(quantile(random_go$V4,probs=0.95))
#Second threshold
two <- unname(quantile(random_go$V4,probs=0.99))

#random_go <- read.table("/Users/harmen/PycharmProjects/blast_project/results/random_go/random_human_go.out")
pdf("./go_random_density.pdf", height=5)
plot(density(random_go$V4[random_go$V4 > 0]), 
     main="Density of GO terms in 1000 random human proteins")
#Add lines and text to density plot
abline(v=one, col="red")
abline(v=two, col="red")
text(one+0.04, 14, "95%\n0.14", cex=0.8)
text(two+0.04, 14, "99%\n0.27", cex=0.8)
dev.off()
