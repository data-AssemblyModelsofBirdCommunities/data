install.packages('phytools')
library(phytools)

library(readxl)
trait1 <- read_excel("E:/xxx/Functional characteristic matrix of birds.xlsx")
trait1 <- as.data.frame(trait1)
col_names <- trait1[, 2]
trait <- trait1[, 6]
names(trait) <- col_names

library(ape)
setwd("")
tree <- read.nexus("Bird pedigree tree.tre")


x1<-log(setNames(trait,
                 col_names))

lambda.x1<-phylosig(tree,x1,
                    method="lambda",test=TRUE)
print(lambda.x1)



install.packages('rtrees')
library(rtrees)
library(ape)
library(FD)
library(readxl)
library(phytools)
plot(tree, type = "fan", cex=0.35,
         use.edge.length = TRUE,font = 4,label.offset=1.5,edge.width= 2.4)



library(FD)
library(picante)

trait1 <- as.data.frame(trait1)
rownames(trait1) <- trait1[,2]
trait <- trait1[,c(3:6)]

GOWer_distance <- gowdis(trait)
Hclust_trait <- hclust(GOWer_distance,"average")
FD_TREE <- as.phylo(Hclust_trait)

abon1 <- read_excel("E:/111/222/Bird Diversity Matrix.xlsx")
abon1 <- as.data.frame(abon1)
abon <- abon1[-1,]
rownames(abon) <- abon[,1]
abon <- abon[,-1]

FD <- pd(abon,FD_TREE)



SES.FD <- ses.pd(abon,FD_TREE,null.model = 'richness',
                 runs=999,iterations=1000)




GOWer_distance <- as.matrix(GOWer_distance)
abon[] <- lapply(abon, as.numeric)
MFD <- mpd(abon,GOWer_distance,abundance.weighted = TRUE)



library(dplyr)
library(magrittr)

cbind(MFD,row.names(abon))  %>%as_tibble() %>%mutate(MFD=as.numeric(MFD))

SES.MFD <- ses.mpd(abon,GOWer_distance,null.model = 'taxa.labels', abundance.weighted = TRUE,
                   runs=999,iterations=1000)




PD <- pd(abon,tree,include.root = TRUE)



SES.PD <- ses.pd(abon,tree,null.model = 'richness',
                 runs=999,iterations=1000)




pd_dist <- cophenetic(tree)
abon <- as.data.frame(abon)

MPD <- mpd(abon,pd_dist,abundance.weighted = TRUE)





cbind(MPD,row.names(abon))  %>%as_tibble() %>%mutate(MPD=as.numeric(MPD))

SES.MPD <- ses.mpd(abon,pd_dist,null.model = 'taxa.labels', abundance.weighted = TRUE,
                   runs=999,iterations=1000)



install.packages("performance")
library(performance)



library(readxl)
data1 <- read_excel("E:/111/222/Polynomial regression matrix.xlsx")
data1 <- as.data.frame(data1)
rownames(data1) <- data1[,1]
x <- c(2050, 2350, 2650, 2950, 3250, 3550, 3850)
y <- data1[,2]

mode1 <- lm(y ~ poly(x, degree = 1, raw = TRUE))
mode2 <- lm(y ~ poly(x, degree = 2, raw = TRUE))
mode3 <- lm(y ~ poly(x, degree = 3, raw = TRUE))
#
plot(x, y, xlab = "Elevation (m)", ylab = "Species richness", 
     lwd = 2, xlim = c(1900, 4000), xaxt = "n",
     pch = 16, col = "blue")
# 
axis(1, at = seq(1900, 4000, by = 300))

# 
smooth_line <- predict(mode2, newdata = data.frame(x = seq(min(x), max(x), length.out = 1000)))

lines(seq(min(x), max(x), length.out = 1000), smooth_line, col = "green", type = "l", lwd = 2,
      lty = 1)







#
x <- c(2050, 2350, 2650, 2950, 3250, 3550, 3850)
y <- data1[,3]
mode1 <- lm(y ~ poly(x, degree = 1,raw = TRUE))
mode2 <- lm(y ~ poly(x, degree = 2,raw = TRUE))
mode3 <- lm(y ~ poly(x, degree = 3,raw = TRUE))
plot(x, y, xlab = "Elevation (m)", ylab = "Functional diversity", xlim = c(1900, 4000)
     ,xaxt = "n",
     pch = 16, col = "blue")
# 
axis(1, at = seq(1900, 4000, by = 300))

# 
smooth_line <- predict(mode3, newdata = data.frame(x = seq(min(x), max(x), length.out = 1000)))

lines(seq(min(x), max(x), length.out = 1000), smooth_line, col = "green", type = "l", lwd = 2,
      lty = 1)













# 
x <- c(2050, 2350, 2650, 2950, 3250, 3550, 3850)
y <- data1[,4]
mode1 <- lm(y ~ poly(x, degree = 1,raw = TRUE))
mode2 <- lm(y ~ poly(x, degree = 2,raw = TRUE))
mode3 <- lm(y ~ poly(x, degree = 3,raw = TRUE))
plot(x,y, xlab = "Elevation (m)", ylab = "Phylogenetic diversity", lwd = 2, xlim = c(1900, 4000)
     ,xaxt = "n",
     pch = 16, col = "blue")
# 
axis(1, at = seq(1900, 4000, by = 300))

# 
smooth_line <- predict(mode2, newdata = data.frame(x = seq(min(x), max(x), length.out = 1000)))

lines(seq(min(x), max(x), length.out = 1000), smooth_line, col = "green", type = "l", lwd = 2,
      lty = 1)




##
# SES_FD
x <- c(2050, 2350, 2650, 2950, 3250, 3550, 3850)
y <- data1[,5]
mode1 <- lm(y ~ poly(x, degree = 1,raw = TRUE))
mode2 <- lm(y ~ poly(x, degree = 2,raw = TRUE))
mode3 <- lm(y ~ poly(x, degree = 3,raw = TRUE))
plot(x, y, xlab = "Elevation (m)", ylab = "Standardized effect size of FD"
     , xlim = c(1900, 4000),ylim = c(-1.8,1),
     xaxt = "n",
     pch = 16, col = "blue")
# 
axis(1, at = seq(1900, 4000, by = 300))

# 
smooth_line <- predict(mode3, newdata = data.frame(x = seq(min(x), max(x), length.out = 1000)))

lines(seq(min(x), max(x), length.out = 1000), smooth_line, col = "green", type = "l", lwd = 2,
      lty = 1)







# SES_PD
x <- c(2050, 2350, 2650, 2950, 3250, 3550, 3850)
y <- data1[,6]
mode1 <- lm(y ~ poly(x, degree = 1,raw = TRUE))
mode2 <- lm(y ~ poly(x, degree = 2,raw = TRUE))
mode3 <- lm(y ~ poly(x, degree = 3,raw = TRUE))

plot(x, y, xlab = "Elevation (m)", ylab = "Standardized effect size of PD"
     , xlim = c(1900, 4000),xaxt = "n",
     pch = 16, col = "blue")
# 
axis(1, at = seq(1900, 4000, by = 300))

abline(h = 0, col = "red", lty = 2)








# MFD
x <- c(2050, 2350, 2650, 2950, 3250, 3550, 3850)
y <- data1[,7]
mode1 <- lm(y ~ poly(x, degree = 1,raw = TRUE))
mode2 <- lm(y ~ poly(x, degree = 2,raw = TRUE))
mode3 <- lm(y ~ poly(x, degree = 3,raw = TRUE))
plot(x, y, xlab = "Elevation (m)", ylab = "Mean pairwise functional distance"
     , xlim = c(1900, 4000),ylim = c(0.3,0.6)
     ,xaxt = "n",
     pch = 16, col = "blue")

# 
axis(1, at = seq(1900, 4000, by = 300))










# MPD
x <- c(2050, 2350, 2650, 2950, 3250, 3550, 3850)
y <- data1[,8]
mode1 <- lm(y ~ poly(x, degree = 1,raw = TRUE))
mode2 <- lm(y ~ poly(x, degree = 2,raw = TRUE))
mode3 <- lm(y ~ poly(x, degree = 3,raw = TRUE))
plot(x, y, xlab = "Elevation (m)", ylab = "Mean pairwise phylogenetic distance"
     , xlim = c(1900, 4000),ylim = c(74,88)
     ,xaxt = "n",
     pch = 16, col = "blue" )


# 
axis(1, at = seq(1900, 4000, by = 300))

# 
smooth_line <- predict(mode3, newdata = data.frame(x = seq(min(x), max(x), length.out = 1000)))

lines(seq(min(x), max(x), length.out = 1000), smooth_line, col = "green", type = "l", lwd = 2,
      lty = 1)







#SES_MFD
x <- c(2050, 2350, 2650, 2950, 3250, 3550, 3850)
y <- data1[,9]
mode1 <- lm(y ~ poly(x, degree = 1,raw = TRUE))
mode2 <- lm(y ~ poly(x, degree = 2,raw = TRUE))
mode3 <- lm(y ~ poly(x, degree = 3,raw = TRUE))
plot(x, y, xlab = "Elevation (m)", ylab = "Standardized effect size of MFD"
     , xlim = c(1900, 4000),ylim = c(-5,2)
     ,xaxt = "n",
     pch = 16, col = "blue" )
# 
axis(1, at = seq(1900, 4000, by = 300))
abline(h = 0, col = "red", lty = 2)



# SES_MPD
x <- c(2050, 2350, 2650, 2950, 3250, 3550, 3850)
y <- data1[,10]

mode1 <- lm(y ~ poly(x, degree = 1,raw = TRUE))
mode2 <- lm(y ~ poly(x, degree = 2,raw = TRUE))
mode3 <- lm(y ~ poly(x, degree = 3,raw = TRUE))
plot(x, y, xlab = "Elevation (m)", ylab = "Standardized effect size of MPD"
     , xlim = c(1900, 4000),,ylim = c(-2.8,0)
     ,xaxt = "n",
     pch = 16, col = "blue" )
# 
axis(1, at = seq(1900, 4000, by = 300))
abline(h = 0, col = "red", lty = 2)

# 
smooth_line <- predict(mode3, newdata = data.frame(x = seq(min(x), max(x), length.out = 1000)))

lines(seq(min(x), max(x), length.out = 1000), smooth_line, col = "green", type = "l", lwd = 2,
      lty = 1)

####
compare_performance(mode1, mode2,mode3)





    



####
install.packages('devtools')
library(devtools)
devtools::install_github("cjbwalsh/hier.part")
library(hier.part)


?hier.part

x1 <- read_excel('E:/111/222//Standardized matrix of various diversity indices and environmental factors.xlsx')


s1 <- x1$sr
NDVI <- x1$NDVI
pre <- x1$pre
tem<- x1$tem
int <- x1$int
area <- x1$area
hab <- x1$hab
data <- data.frame(NDVI, pre, tem,int , area, bab) 


s2 <- x1$FD

FD <- hier.part(s2, data, fam = "gaussian", gof = "Rsqu")