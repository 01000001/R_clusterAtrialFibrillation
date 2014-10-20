#CODE to cluster and maintain 10 clusters of the dataset with divisive analysis clustering (DIANA):
  
dmed <- read.csv("Dmed.csv", header=T, sep=",")
dmedAF <- subset(dmed, AF==1) #take only the AF patients
dmed.diana <- diana(dmedAF[,6:20]) #create variable for diana clustering of small dataset, DC: 0.70
dmed.diana.dist <- dist(dmedAF[,6:20])
pltree(dmed.diana)
hcl.dmed.diana <- as.hclust(dmed.diana)
plot(dist(dmedAF[,6:20]),cophenetic(hcl.dmed.diana)) # validate
cor(dist(dmedAF[,6:20]),cophenetic(hcl.dmed.diana)) #calculate the correlation: 0.64

dmed.diana.cut10 <- cutree(hcl.dmed.diana,10)
cluster.d=cutree(hcl.dmed.diana, 10)
#for( i in 1 : 926) { dmedAF$cluster[i] = cluster.d[i]}