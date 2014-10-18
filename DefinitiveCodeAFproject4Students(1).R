# Using the Cluster Package of R for Various Cluster Techniques
library(cluster)
library(stats)
library(rms)
library(MASS)
library(FSelector)


dmed <- read.csv("~/Documents/R_clusterAtrialFibrillation/Dmed.csv", header=T, sep=",")
dmedAF <- subset(dmed, AF==1) #take only the AF patients
dmed.s <- head(dmedAF,100) #create smaller dataset
dmed.s.diana <- diana(dmed.s) #create variable for diana clustering of small dataset

dmed.s100 <- head(dmedAF,100) #take 100 samle
dmed.s100ok <- head(dmed,100) #take 100 samle without the AF
dmed.s100.dist <- dist(dmed.s100[,6:20]) #calculate the distance
dmed.s100.hclust <- hclust(dmed.s100.dist) #do clustering
plot(dmed.s100.hclust, labels=dmed.s100$pid) #plot data with labels

dmed.s1000 <- head(dmedAF,1000) #take 100 samle
dmed.s1000.dist <- dist(dmed.s1000[,6:20]) #calculate the distance
dmed.s1000.hclust <- hclust(dmed.s1000.dist) #do clustering
plot(dmed.s1000.hclust, labels=dmed.s1000$pid) #plot data with labels

dmed.s1000.diana = diana(dmed.s1000.dist)
pltree(dmed.s1000.diana)
dmed.s1000.hclust <- as.hclust(dmed.s1000.diana)
plot(dmed.s1000.hclust)

counts = sapply(1:10,function(ncl)table(cutree(dmed.s100.hclust,ncl))) #cut the tree on different levels
names(counts) = 1:10
counts

dmed.3 <- cutree(dmed.s100.hclust, 3) #variable dmed.3 for only 3 cluster
sapply(unique(dmed.3),function(g)dmed.s100$AF[dmed.3 == g]) #show how many AV patients are in different groups

table(dmed.3,dmed.s100$AF) #how on a table how many AF patients are in different clusters

#for (i in 1:20) { dmed.s100[i] <- cutree(dmed.s100.hclust, i) } #create automaticly groups with clusters

### for all
dmed <- read.csv("~/Documents/R_clusterAtrialFibrillation/Dmed.csv", header=T, sep=",")
dmedAF <- subset(dmed, AF==1) #take only the AF patients
dmed.diana <- diana(dmedAF[,6:20]) #create variable for diana clustering of small dataset
dmed.diana.dist <- dist(dmedAF[,6:20])

##chads
cheds <- read.csv("~/Documents/R_clusterAtrialFibrillation/chads.csv", header=T, sep=",")


dlab <- read.csv("~/Documents/R_clusterAtrialFibrillation/Dlab.csv", header=T, sep=",")
dicp <- read.csv("~/Documents/R_clusterAtrialFibrillation/Dicp.csv", header=T, sep=",")

#Introduction
#In this practical and project you will learn to apply and analyze data using different clustering techniques available in R. Cluster algorithms are used to analyze data not containing a pre-specified class attribute. This is why clustering is referred to as unsupervised learning, in contrast to classification in which a relationship is sought between a class attribute and the other attributes. So although most of the datasets we use have a class attribute these will be ignored except for evaluation purposes. We will use two main types of cluster methods in this lab: hierarchical methods (mainly Diana) and a partitioning method (Pam). 


# Preparation. Before we start we need to include the relevant 
# libraries. For example for the practical we will use the packages
# `cluster??? and `stats???. The `stats??? package is part of 
# R and does not need installation. If the `cluster??? package is 
# not yet installed you need to install it. To install, go to 
# ???Install package(s)??? in the "Packages" menu. The first time 
# you are asked to select a CRAN mirror. Choose a mirror that is 
# close to you, for example Netherlands(Amsterdam 2). Select from the list the package you would like to install (in this case ???cluster???) and press ???OK???. The packages include data and functions that we need in the practical and the project. Installing a package creates a library with the same name but the library still needs to be loaded. After installation you need to load the libraries (the installation is needed only once, loading the library is needed for each session):

library(cluster)
library(stats)
library(rms)
library(MASS)
library(FSelector)


# We provide you with a function that yo might use later. 
# The function calculates silhouette widths for pam clustering 
# with k clusters
pam.silhouetteWidths <- function(data, k){
  for (i in 2:k){
    print(paste("k=",i," s=",pam(data,i)$silinfo$avg.width,sep=""))
  }
}

# Now we are ready to start.. We start with the pratcical and later we will move on to the
# project.

# ------Clusterting the animals example data set-----------------------

# We will use the animals database (which is part of the cluster package).

#To get help on the data type
?animals 

# The rows indicate the animals (this was not shown in the help file but provided for you here):
# ant, bee, cat, caterpillar, chimpanzee, cow, duck, eagle, elephant, fly, frog, herring, lion, lizard, lobster, man, rabbit, salmon, spider, and whale.

# The variables and their meanings are:
#  WAR : warm-blooded  FLY :  can fly 		VER :  vertebrate
# END : endangered 		GRO :  live in groups 	HAI :  have hair


# Let's take a look at the dataset
animals

# In principle our variables (WAR, FLY etc) are binary. But they are defined as integers:

class(animals$war)

# for simplicity in this practical we will consider them as integers.

# You can see that there are missing values (denoted with NA). For example we do not know if the lion is an endangered animal or not. Also note that each variable can take one of two values: 1 and 2. As it is coded now "2" means "Yes", but it is easier to think in terms of values of 0 ("No") and 1 ("Yes"). So let us translate the 1 and 2 to 0 and 1. This is simple: just substract "1" from "animals" which will work on each cell. We can do this because the variables are of type integer.

myanimals <- animals -1

# Now look at "myanimals"

myanimals?dist

# You see that when you substract a number from a missing value it stays a missing value.
# The codes mean:
# WAR : 1 = warm-blooded
# FLY : 1 = can fly
# VER : 1 = vertebrate
# END : 1 = endangered
# GRO : 1 = live in groups
# HAI : 1 = have hair

# You see for example that the cat is warm blooded, unlike the ant and the bee.

#---------------------------------------
# Calculating distances


# We use the function dist (part of the cluster library) to calculate all the distances between the animals.
# "dist" calculates, by default, the eucledian distance. Look at "?dist" for other options.
animals.d <- dist(myanimals)

# => Look at the distance "matrix".
animals.d

# animals.d forms a compact representation of the distance matrix. We can coerce (force) it into a normal matrix
animals.d.mt <- as.matrix(animals.d) 

# Inspect the distances between the first 4 animals
animals.d.mt[1:4, 1:4]

# What is the ditance between an ant and a bee?
# 1.414214
# What is the distance between an ant and a cat?
# 2
# Note that the matrix is symmetric (the distance between ant and bee is the same as between bee and ant, so you can get the answers by looking at rows or columns first)

#----------------------------------------
Clustering - using the hclust algorithm

# Look at the available functions in the ???cluster??? package.
help(package="cluster")

# What is hclust?
?hclust

# As you can read this is an agglomorative algorithm: each instance (animal) is first considered as a separate cluster and then clusters are merged. This is different from divisive algorithms such as Diana where one begins with all animals as one cluster and then partition the cluster recursively.

# Let us now cluster the dataset with the hclust algorithm.
# The default method (see ?hclust) is complete linking. If you can't rememeber what this means then look at the slides of the clustering lecture. 
# complete linking: Largest distance between observations

cl.hclust <- hclust(animals.d)

# Let us plot the resulting dendrogram
plot(cl.hclust)

# In which subcluster is man? Is this surprising?
# man is clustered together with the chimpanzee and the lion, indeed, it is somewhat surprising

# Let us now use single linking. Make sure you understand what this means.
# single linking - link entities using the smallest distance

cl.hclust.single <- hclust(animals.d, method="single")

# Let's create a new graphical device and plot the result.
plot(cl.hclust.single)

# Are there differences? Where is man now?
# Man is now clustered with rabbit, cow, chimansee and lion


#----------------------
# Clustering - using Diana too
# Let us now move to cluster with the Diana algoritm (which is a divisive algorithm) and analyse the results.

### Agglomerative: Start with the points as individual clusters and, at each step, merge the closest pair of clusters. This requires defining a notion of cluster proximity.

### Divisive: Start with one, all-inclusive cluster and, at each step, split a cluster until only singleton clusters of individual points remain. In this case, we need to decide which cluster to split at each step and how to do the splitting.

cl.diana <- diana(animals)

# We only specified the dataset so Diana will use the default parameters. If you type ????diana??? for help you will see that the call to the function looks like this:
  
#  diana(x, diss = inherits(x, "dist"), metric = "euclidean", stand = FALSE,
#        keep.diss = n < 100, keep.data = !diss)
# You see that the default metric (distance measure) is ???euclidean??? and that ???stand??? is False meaning that the variable values are not standardized.

pltree(cl.diana)


# Show one agreements between the results of clustering using Diana and using hclust.

# You may be wondering what the height on the dendrogram means when two clusters meet. For example consider the height just above "1.5" in the Diana dendrogram where the cluster "ant, lob" and the cluster "cpl, spi" meet. This is the mean of distances between all the elements in the two clusters.

# Note the agreement and differences between the Diana clustering with the former cl.hclust. Where is man now?
# man is clustered now with the chipansee only

# We can get a summary of the results
summary(cl.diana)

# Look at the divisive coefficient (read Chapter 6 in the "Finding Groups" book). Because the summary itself (like most things in R) is an object we can also ask for the divisive coefficient directly:
## DC aslo called agglomerative coefficient, DC = 0 No cluster structure, DC = 1 Clear cluster structure

summary(cl.diana)$dc

# To know which names belong to an object try: 
names(summary(cl.diana))
# now you see that "dc" is one of the variables that is used in the summary. You can of course also first save the summary and then ask about the names

mysum <- summary(cl.diana)
names(mysum)
mysum$dc

# We can in fact coerce the diana representation to an hclust representation to allow us performing operations defined for hclust. 

# The original type (class) is:
class(cl.diana)

# We see that "twins" is also reported. This means that operation defined for the class "twins" are also defined for "diana". In other words "diana" inherts operations from "twins".
# Let's make our object cl.diana of type hclust (of course the original cl.diana does not change but we create a new variable with the hclust type)
hcl.diana <- as.hclust(cl.diana)

# After coercion:
class(hcl.diana)

# now we can cut the tree at a given height
clusters.d <- cutree(hcl.diana, h=1.8)

# Look at the result
clusters.d

# Now you see that there are 5 clusters (the numbers from 1 to 5). The animals ant, lob, cpl and spi belong to one cluster (cluster number 1), and bee and fly form another cluster (clutsre 2), etc. Everything below the cut will be put in one cluster.

# Let us find out the elements of cluster 2
clusters.d[clusters.d==2]

## Instead of saying where to cut you may specify how many clusters you want by using the argument "k" in cutree

da2 <- cutree(hcl.diana, k = 2)

# We cal also plot rectangles around these clusters
plot(hcl.diana)
rect.hclust(hcl.diana, 2)

# we can make a table to see the number of objects in each cluster 
table(da2)
# We see that there are 16 objects in cluster 1 and 4 in cluster 2.

# Get the animals in cluster 1
rownames(myanimals)[da2 == 1]

# And in cluster 2
rownames(myanimals)[da2 == 2]

##-----------
# Different dissimilarity metrics

# For two groups, does the metric matter ?
da.manhattan <- diana(animals, metric="manhattan")
da.manhattan2 <- cutree(as.hclust(da.manhattan), k = 2)
da2
da.manhattan2
da2 == da.manhattan2
table(da2 == da.manhattan2)## identical group assignments
myanimals[da.manhattan2==1,] # give the rows in the data frame for the first cluster

#cross tabulation of two cluster methods
table(cutree(cl.hclust, k = 3), cutree(hcl.diana, k = 3))

# We see for example that there are 4 common animals in cluster 1 but no common elements in cluster 3. 

#--------------------------------------------------
# Quality of a cluster

# Any clustering algorithm provides clusters, the important question is whether the clusters represent the "natural structure" in the data. One simple way to check the quality of the cluster is the following: we know the distances in the data set between any two objects (animals), but we also can calculate a "dendrogrammatic distance" between any two objects as well. This distance is defined as the height of the node at which these two points are first joined together (see earlier remark on height). Intuitively we expect from a good clustering algorithm that if the distance between two objects in the original data set is large then the distance between them in the dendrogram is also large. We can quantify the association between the two types of distances by using the correlation between them. This is called the cophenetic correlation, which is just the correlation between the distances within the objects in the data and within the dendrogam. The function "cophenetic" in the stats library (which you have uploaded at the beginning of this session) calculates the dendrogrammatic distance.

# Calculate the cophenetic distance
cophenetic(hcl.diana)

# plot the cophenetic distance as function of the (regular) distance.
# If the two are perfectly correlated we will get the points on the
# line y=x
plot(dist(myanimals),cophenetic(hcl.diana))
abline(0,1) # plot y = x (intercept = 0 and slope = 1)

# Not completely perfect but a good result

# Calculate the correlation
cor(dist(myanimals), cophenetic(hcl.diana))

# High correlation indeed.

# Calculate characteristics of clusters
aggregate(myanimals, by=list(da2), mean, na.rm=T)
aggregate(myanimals, by=list(da2), sd, na.rm=T)

#################################################################
## NOW we move to the PROJECT 
# -------------Cluster on AF Data------------------------------
# # Description: we have 3 ANONYMIZED datasets. The source of data is the 6 GAZO centres (Gezondheidscentra Amsterdam ZuidOost) where about 45 General Practitioners work. The AF patients has ICPC code "K78". The patients that do not have ICPC codes were selected at random. The data corresponds to the years 2000 to 2012. There was no restriction on age.
# 
# 1. The dataset Dmed.csv, which includes patient identifier, sex AF status (no, yes), number of medications that the patient is using, and whether the patient is using any of the most frequently used 20 medication types in the sample (for each type the patient will have a "TRUE" or "FALSE" value indicating whether the patient uses this medication or not):
# pid,age,sex,AF,PatMedNum,C10AA01,M01AB05,A10BA02,A02BC01,C07AB02,M01AE01,R05DA04,N02BE01,J01CA04,S01XA20,N02BE51,D02AX,D07AB09,C03AA03,J01AA02,B01AC08,D07AA02,N05BA04,C09AA03,J01CR02
# 
# Definition of meds:
# SIMVASTATINE SDZ 20MG T FO;C10AA01
# DICLOFEN NA ACT 50MG T MSR;M01AB05
# METFORMINE HCL PCH 500MG T;A10BA02
# OMEPRAZOL ACT 20MG CAPS MSR;A02BC01
# METOPROLOLTAR ACT 100MG TAB;C07AB02
# IBUPROFEN KRING 200MG TABL;M01AE01
# CODEINEFOSFAAT A 10MG TABL;R05DA04
# PARACETAMOL ACT 500MG TABL;N02BE01
# AMOXICILLINE RP 500MG CAPS;J01CA04
# DURATEARS OOGDR FL 15ML;S01XA20
# PARACE/CODEIN PCH 500/10 T;N02BE51
# CETOMACROGOLZALF FNA;D02AX
# CR TRIAMCINOLONI ACET 1MG/G;D07AB09
# HYDROCHL THIAZ PCH 12,5MG T;C03AA03
# DOXYCYCL ACTA 100MG DISP TA;J01AA02
# CARBASAL CALC CARDIO 100MG;B01AC08
# HYDROCORTISON 10MG/G ZALF;D07AA02
# OXAZEPAM ACTAVIS 10MG TABL;N05BA04
# LISINOPRIL ACTAVIS 10MG TAB;C09AA03
# AMOXI/CL ACTAV 500/125 TABL;J01CR02
# 
# 2. The datset Dlab.csv which includes patient identifier, sex AF status (no, yes), number of Lab results that the patient has, and whether he had any of the most frequently used 20 lab result types (for each type the patient will have a "TRUE" or "FALSE" value indicating whether the patient had a lab result of that type or not):
#   
# pid,age,sex,AF,PatLabNum,KREA,GLUC,RRSY,RRDI,CHOL,HDL,K,TRIG,KREM,CHHD,GEW,HBAC,LNGP,GLHB,QUET,ROOK,LDL,LIBW,ALBK,ALCO,DMHB,COHZ,HB,OMVA,MOFV,INSP,ERY,HT,LEUK,DBLO
# 
# Definition of Lab orders inclusive example of their possible value:
# KREA,kreatinine,5/31/12,140
# GLUC,glucose nuchter. draagbare meter,9/6/12,5
# RRSY,systolische bloeddruk,9/6/12,130
# RRDI,diastolische bloeddruk,9/6/12,81
# CHOL,cholesterol totaal,5/31/12,3
# HDL,HDL-cholesterol,5/31/12,0
# K,kalium,5/31/12,4
# TRIG,triglyceriden,5/31/12,1
# KREM,eGFR volgens MDRD formule,6/15/12,46
# CHHD,cholesterol/HDL-cholesterol ratio,5/31/12,4
# GEW,gewicht pati_??nt,9/6/12,68
# HBAC,HbA1c (glycohemoglobine) IFCC,5/31/12,50
# LNGP,lengte pati_??nt,9/6/12,175
# GLHB,glycohemoglobine (HbA1c) DCCT,3/14/11,7
# QUET,Quetelet-index (BMI) pati_??nt,9/6/12,22.2
# ROOK,roken,9/6/12,1
# LDLD,LDL-cholesterol direct,5/31/12,1
# LIBW,lichaamsbeweging,9/6/12,2
# ALBK,albumine (micro-) /kreatinine urine,5/31/12,1
# ALCO,alcoholgebruik,9/6/12,0
# DMHB,hoofdbehandelaar diabetes,9/6/12,1
# COHZ,coron HZ in naaste familie <60jr,6/15/12,2
# HB,hemoglobine (Hb),5/31/12,7
# OMVA,buikomvang (middelomtrek),9/6/12,96
# MOFV,monofilament linkervoet,6/15/12,1
# INSP,inspectie linkervoet (diabetes),6/15/12,1
# ERY,erytrocyten,5/31/12,4
# HT,hematocriet (Ht),5/31/12,0
# LEUK,leukocyten,5/31/12,9
# DBLO,doorbloeding rechtervoet,6/15/12,1
# 
# 3. The dataset Dicp.csv, which includes patient identifier, sex AF status (no, yes), number of ICPC codes that the patient has, and whether he had any of the most frequent 20 ICPC codes (for each code the patient will have a "TRUE" or "FALSE" value indicating whether the patient had this ICPC code or not):
# 
# pid,age,sex,AF,PatIcpNum,T90.02,K86.00,K78.00,U71.00,L03.00,L90.00,F92.00,T93.01,T82.00,U99.01,R05.00,D12.00,K74.00,R97.00,K87.00,K77.00,L08.00,L92.00,R95.00,A00.00,L15.00,L04.00,R74.00,H84.00,P76.00,T92.00,D02.00,K90.00,K75.00,S88.00
# 
# The ICPC codes can be obtained online from https://www.nhg.org


# There are 2 ideas that you may want to pursue: 1) if you ignore the "AF" attribute then you can cluster and then consider the clusters with a large number of AF patients. You may then inspect the non-AF patients in such clusters and try to see if they may in fact be underdiagnosed (they seem to be AF patients but they are not diagnosed); 2) cluster only the AF patients and try to see whether there are interesting subgroups within this group.
# Tip: For the experiments it might be useful to use a small sample of the data first before you embark on an analysis of the complete data sets. For example if you want to work with 200 patients you could use: smallDataset <- pat[sample(1:nrow(pat), size=200),].
# The assignmet in short is:
# A. Decide on the goal of your analysis: do you want to inspect the differences between AF and non AF patients in clusters? or do you want to try to identify subclusters of only AF patients.
# B. Decide on a clustering method, try to motivate for yourself why.
# C. There are probably too many variables in the datasets to be practically useful and they may hurt the analysis as well: In addition to sex and age, there are 20 meds, 20 lab types and 20 icpc codes, along with the numbers of meds, the number of lab orders and the number of icpc codes a patient has. Decide on how you want to decide which attributes are the most relevant for your problem.
# D. perform the clustering.
# E. Provide information about the quality of the clusters obtained.
# F. Answer your research question and discuss strength and limitations of your approach.
# G. Write a scientific report of max 10 pages.

#################
# Here is just some code to give you some ideas of using clustering
# reading a dataset
dmed <- read.csv("Dmed.csv", header=T, sep=",")

1. Inspecting the dataset
# How does my dataset look like. Use structure
str(dmed)
# What is the mean age of the patients
mean(dmed$age)
# What is the mean age of AF patients and non AF patients. The command tapply(Var1, Var2, Fun) is an extremely powerful command. It applies the function "Fun" (here mean) for Var1 (here age) for each value of Var2 (here AF status). You might also be interested in the function aggregate.
tapply(dmed$age, dmed$AF, mean)

# What is the percentage of patients using medication "A10BA02"?
mean(dmed$A10BA02==TRUE)

# What is the percentage of males?
mean(dmed$sex=="M")

# make a subset of AF patients
AFpats <- subset(dmed, AF==1)

#how many AF patients are there?
nrow(AFpats)

# 2. Clustering
# For a given set of variables clustering can be performed (you will see below that it might pay to select specific variables insetad of using all variables). Different algorithms are possible. For the binary variables one may use MONA. But in general one may use DAISY to calculate distances among patients using any variable type. For example you may want to use the ???daisy??? function to calculate distance for all of your variables (binary and continuous) and then proceed with the various common hierarchical algorithms. You may also consider the Mona (Monothetic Analysis) algorithm, which is designed to work with binary data and considers only one variable at a time, while other hierarchical approaches (e.g. classical HC, Agnes, Diana) use all variables at each step. To analyze your results look at the summaries, silhouette widths and divisive/agglomerative coefficients. 

# Let's use MONA to cluster. MONA works with binary data
mona(dmed[,6:20]) # assuming the 6th til the 20th variables are binary

# Run Diana on same dataset (can take some time)
Med.diana <- diana(dmed)

pltree(Med.diana)
Med.hclust <- as.hclust(Med.diana)
rect.hclust(Med.hclust, 2)



summary(Med.diana)

# Let's make 5 clusters
clusterD <- cutree(as.hclust(Med.diana), k=5)

# Selecting attributes associated with AF. We can use information gain for that.
# We will use "apply" with "2" to work on each row and activate "information.gain"
# For each row of Med.sub we call the function "information gain" from FSelector. Now weigths will be a vector of infomration gain per variable.

weights <- information.gain(AF ~ . -pid, dmed)
print(weights)
vars.4 <- cutoff.k(weights, 4) # top 4 variables
Med4 <- subset(dmed,select=vars.4)


#Run various clustering algorithms
Med4.diana <- diana(Med4, stand=TRUE)
pltree(Med4.diana)
bannerplot(Med4.diana)
summary(Med4.diana)
clusterD4.3 <- cutree(as.hclust(Med4.diana), k=3)
plot(Med4, col = clusterD4.3)

# For those who want to work with PAM (Partioning around Medoids)
# here is some code to get you started
pam.silhouetteWidths(Med4,10)
Med4.pam <- pam(Med4, 2, stand=TRUE)
Med4.pam$silinfo$avg.width
clusplot(Med4.pam)
plot(silhouette(Med4.pam))
summary(Med4.pam)
plot(Med4, col = Med4.pam$clustering)


# Look at AF in each cluster
# Describe each cluster in terms of its descriptive statistics. When you inspect clustering results try to discern differences between them. Because AF is important, inspect differences in AF occurrence among them. See if there are clusters of patients with AF that is much higher or lower than the mean AF of the whole sample. Which attribute distributions are very different within the cluster?

# Provide various measures to inspect the quality of the clusters. For example, the cophenetic correlation. If you have an outcome such as AF then you may want to provide performance measures such as the Brier score (Google it up).

Med.cl.inf <- cbind(dmed,"cluster.pam"=Med4.pam$clustering, "cluster.diana"=clusterD4.3)
mean(Med.cl.inf$AF)
tapply(Med.cl.inf$AF, Med.cl.inf$cluster.pam, mean)
tapply(Med.cl.inf$AF, Med.cl.inf$cluster.diana, mean)

# Looking at age and med#
mean(as.numeric(Med.cl.inf$age), na.rm=T)
tapply(as.numeric(Med.cl.inf$age), Med.cl.inf$cluster.pam, mean, na.rm=T)
mean(as.numeric(Med.cl.inf$PatMedNum), na.rm=T)
tapply(as.numeric(Med.cl.inf$PatMedNum), Med.cl.inf$cluster.pam, mean, na.rm=T)

# Attribute/cluster plots
boxplot(as.numeric(Med.cl.inf$age)~as.numeric(Med.cl.inf$cluster.pam))
boxplot(as.numeric(Med.cl.inf$PatMedNum)~as.numeric(Med.cl.inf$cluster.pam))

### REPORT
# Your report should have a maximum of 10 pages. Your report should clearly describe the following elements: What is the problem that you are addressing? What is your research goal? What are the research questions? What was your approach? Which results did you obtain? Discuss the results (what could they mean? What are the strengths/limitations etc?). Provide a conclusion that includes the answer to the research questions. Which further research should be done? It is important that you reflect on the experiment: you have performed clustering, but what did you really learn from it? Did you attain your goal? Do the clusters provide adequate answers to the questions you had?
