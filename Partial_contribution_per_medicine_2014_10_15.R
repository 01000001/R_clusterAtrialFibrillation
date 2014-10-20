med <- read.csv("mergemed.csv", header=T, sep=",")


afmed2 <- subset (med, med$cluster== 2, select = c("pid", "M01AB05","A10BA02","C07AB02","M01AE01","N02BE01","N02BE51","C03AA03","B01AC08","N05BA04","C09AA03"))
afmed3 <- subset (med, med$cluster== 3, select = c("pid", "M01AB05","A10BA02","C07AB02","M01AE01","N02BE01","N02BE51","C03AA03","B01AC08","N05BA04","C09AA03"))
#matrix for medication related with AF for cluster 2 and 3

afmed2.d <- dist(afmed2)
afmed3.d <- dist(afmed3)
#euclidean distances

afmed2.d.mt <- as.matrix (afmed2.d)
afmed3.d.mt <- as.matrix (afmed3.d)
#coerce to a normal matrix

NON AF MEDICATION
#same for non AF medications
nafmed2 <- subset (med, med$cluster== 2, select = c("pid", "C10AA01","A02BC01","R05DA04","J01CA04","S01XA20","D02AX","D07AB09","J01AA02","D07AA02","J01CR02"))
nafmed3 <- subset (med, med$cluster== 3, select = c("pid", "C10AA01","A02BC01","R05DA04","J01CA04","S01XA20","D02AX","D07AB09","J01AA02","D07AA02","J01CR02"))
#matrix for medication related with AF for cluster 2 and 3

nafmed2.d <- dist(nafmed2)
nafmed3.d <- dist(nafmed3)
#euclidean distances

nafmed2.d.mt <- as.matrix (nafmed2.d)
nafmed3.d.mt <- as.matrix (nafmed3.d)
#coerce to a normal matrix


#Compare matrixes
#cor(dist(nafmed2.d.mt), dist(afmed2.d.mt))
#cor(dist(nafmed3.d.mt), dist(afmed3.d.mt))
#help(cor)

#diff.mt <- nafmed2.d.mt - afmed2.d.mt
#for (i in 1:164){
#  diff.mt$mean = rowMeans(diff.mt,dim(diff.mt))
#}

##append the mean of rows:
#diff.mt <- cbind(diff.mt, rowMeans(diff.mt))

#make a summed matrix
summed2.d.mt <- afmed2.d.mt + nafmed2.d.mt

#afmed2.d.mt[1,2] + nafmed2.d.mt[1,2]
#summed2.d.mt[1,2]

#make binary variable of 0 (CHAD < 2) and 1 (Chad > 1) for "high" chad score
med$high_chads <- med$chads > 1

#mean(med$A10BA02==TRUE & !med$high_chads == TRUE)

#mean(med$cluster == 1)

##get the total number of patients in particular clusters and append it to dataframe
for (i in 1:nrow(med)){med$total_n_same_cluster[i] <- nrow(subset(med, cluster == med$cluster[i]))}

#all people in cluster 1 using medication A
nrow(subset(med, cluster == 1 & A10BA02==TRUE))

# high_chads / total number in claster x using medication y
#ecample
nrow(subset(med, cluster == 1 & A10BA02==TRUE & high_chads == FALSE))/nrow(subset(med, cluster == 1 & A10BA02==TRUE))

# nrow(subset(med, cluster == 1 & med[,8:27][3]==TRUE & high_chads == FALSE))/nrow(subset(med, cluster == 1 & A10BA02==TRUE))

# Magic:
for (i in 1:10){
  for (m in 1:20){
  cat(( sprintf( "n with chads > 1 / n in cluster %i using med %s: %f \n", i,colnames(med[,8:27])[m], nrow(subset(med, cluster == i & med[,8:27][m]==TRUE & high_chads == TRUE))/nrow(subset(med, cluster == i & med[,8:27][m]==TRUE)))))
  cat(( sprintf( "n with chads < 2 / n in cluster %i using med %s: %f \n", i, colnames(med[,8:27])[m], nrow(subset(med, cluster == i & med[,8:27][m]==TRUE & high_chads == FALSE))/nrow(subset(med, cluster == i & med[,8:27][m]==TRUE)))))}}