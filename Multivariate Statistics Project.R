###############scatterplots with gplot2################
#ggplot 2 wigh multivariate statistics and geometry of principle components
#https://web.stanford.edu/class/bios221/labs/multivariate/lab_5_multivariate.html
install.packages("ade4")
install.packages("ggplot2")
install.packages("grid")
#install.packages("phyloseq")
source("https://raw.githubusercontent.com/joey711/phyloseq/master/inst/scripts/installer.R",
       local = TRUE)
require(ade4) # multivariate analysis
require(ggplot2) # fancy plotting
require(grid) # has the viewport function
require(phyloseq) # GlobalPatterns data

#Fun with SVD and Eigendecomposition

#x <- matrix(rnorm(20), ncol=4)

#Generating a rank one matrix
u=seq(2,30,by=2)
v=seq(3,12,by=3)
X1=u%*%t(v) # transpose so the dimensions match in the multiplication 
X1
#Note that this is the exact same as:
X1=outer(u,v)
X1

# Make a noise matrix to add to X1
Matex <- matrix(rnorm(60,1),nrow=15,ncol=4)
Matex
X <- X1+Matex
X
#Let's see what \(X\) actually looks like. Remember, 
#\(X\) is four dimensional - so to try to visualize 
#this it is easiest to do one of two things. We can 
#look at lots of plots in two dimensions and even make
#a movie where we rotate which two dimensions we're 
#looking from: this is the approach taken in ggobi 
#which you can learn about on your own if you want. 
#Another method is to plot the data in two dimensions 
#and use plotting aesthetics such as point color and point
#size to try to visualize the other dimensions. 
#When using plot aesthetics like this, I think about 
#big points as being closer to me (so I can imagine 
#3 dimensions relatively easily), and for me color
#is the next easiest way to represent a dimension 
#(I struggle with this for more than 2 colors 
#though - the default in ggplot2 ranges from black to blue).
ggplot(data=data.frame(X), aes(x=X1, y=X2, col=X3, size=X4)) + geom_point()

#Here we see that the data looks linear in all four dimensions. 
#This is what it means to be rank one.

#Now let's consider a rank 2 matrix.
set.seed(0)
n <- 100
p <- 4
Y2 <- outer(rnorm(n), rnorm(p)) + outer(rnorm(n), rnorm(p))
head(Y2)

ggplot(data=data.frame(Y2), aes(x=X1, y=X2, col=X3, size=X4)) + geom_point()

#Now there are obviously at least two dimensions because if we project 
#the data onto the first two coordinates (by default called X1 and X2 
#when you convert a matrix into a data frame in R), 
#then the data varies in both dimensions. 
#So the next step is to try to decide if there are more than two dimensions. 
#The top right points are the closest to you (they're biggest) 
#and as you go down and left in the plot those points are farther away. 
#In the left are the bluest points and they seem to get darker linearly as you move right.

#As you can probably tell, it is very hard to visually 
#discover a low dimensional space in higher dimensions, 
#even when "high dimensions" only means 4! This is one 
#reason why we rely on the singular value decomposition.

svd(Y2)$d # two non-zero eigenvalues
## [1] 1.986e+01 1.025e+01 2.330e-15 1.120e-15
Y <- Y2 + matrix(rnorm(n*p, sd=0.01),n,p) # add some noise to Y2
svd(Y)$d # four non-zero eigenvalues (but only 2 big ones)
## [1] 19.85669 10.24460  0.10767  0.08972

#Here we have two dimensions which are non-zero 
#and two dimensions which are approximately 0 
#(for Y2, they are within square root of computer tolerance of 0).


# Principal Component Analysis
# Relating PCA to SVD
# Let \(X\) be a centered but unscaled matrix. 
# We will show that there is a matrix \(X_r\) 
# whose principal component output (without rescaling the columns)
#is the same as the eigendecomposition of \(X'X\).
# 
# The first \(k\) principal components of \(X\) are the first \(k\) 
#directions explaining maximum variance. 
#This is equivalent to the first \(k\) eigenvectors of the covariance matrix. 
#We estimate the sample covariance matrix as \(S = X'X/N\). 
#Hence, the principal component analysis of \(X\) gives the first \(k\) eigenvectors
#of \(X'X/N\).
# 
# So if we let \(X_r = X*\sqrt{N}\), 
#then the pca output will be the first \(k\) eigenvectors 
#of \((X*\sqrt{N})'(X*\sqrt{N}) / N = X'X\). 
#That is the eigendecomposition of (the centered) \(X\).
# 
# Let's apply this to the \(X\) above
#to verify that our calculations are correct. 
#First, we need to center it, and we will call that centered version Xc. 
#In our case, \(N=15\).

Xmeans <- apply(X,2,mean)
Xc <- sweep(X,2,Xmeans)
Sc <- crossprod(Xc)
Sce <- eigen(Sc)
# here is Xr
Xr <- Xc * sqrt(15)
Xr.pca <- dudi.pca(Xr, scale=F, scannf=F, nf=4)
Xr.pca$eig
## [1] 3.036e+05 1.472e+01 1.019e+01 1.724e+00
Sce$values
## [1] 3.036e+05 1.472e+01 1.019e+01 1.724e+00
# eigenvectors are the same (up to sign)
Xr.pca$c1
##        CS1     CS2      CS3       CS4
## V1 -0.1821  0.1380  0.87966 -0.417089
## V2 -0.3648 -0.9282  0.07263  0.005302
## V3 -0.5463  0.1748 -0.46144 -0.676784
## V4 -0.7316  0.2980  0.08938  0.606607
Sce$vectors
##         [,1]    [,2]     [,3]      [,4]
## [1,] -0.1821  0.1380  0.87966 -0.417089
## [2,] -0.3648 -0.9282  0.07263  0.005302
## [3,] -0.5463  0.1748 -0.46144 -0.676784
## [4,] -0.7316  0.2980  0.08938  0.606607

# SVD can be used to determine the direction 
# of the most variance (and next most variance, 
#                       and next most variance, .) 
# and how much of the variation is explained by each of those directions.
# This is exactly the goal of PCA.
# 
# When we use PCA to plot data, we only plot
# the directions in which there is the most change in the data. 
# For example in two dimensional data Y, we can easily plot that
# in two dimensions now and there is very little (actually 0) 
# variation in all other dimensions. By default (using dudi.pca),
# we center the data and then rescale it so each column has a Euclidean norm of 1. 
# Here we show an example and use the default plotting function
# of the package ade4 and then a fancy plot from ggplot2.
Y.pca <- dudi.pca(Y, scannf=F, nf=4) 
scatter(Y.pca) # default quick plot

#To save time later, we'll save a default plot 
#and a screeplot making function.


# set up a plot we'll use later
ppp <- ggplot() + coord_fixed() + 
  labs(x="Comp1, Axis1", y="Comp2, Axis2") +
  geom_hline(yintercept=0, col="darkgrey") + 
  geom_vline(xintercept=0, col="darkgrey")
# make the scree plot in a viewport
myscree <- function(eigs, x=0.8, y=0.1, just=c("right","bottom")){
  vp <- viewport(x=x, y=y, width=0.2, height=0.2, just=just)
  sp <- qplot(factor(1:length(eigs)), eigs, 
              geom="bar", stat="identity") +  
    labs(x = NULL, y = NULL)
  print(sp, vp=vp)
}
#Now actually make the plot

ppp + geom_point(data=Y.pca$li, aes(x=Axis1, y=Axis2, size=Axis3, col=Axis4)) +
  geom_text(data=Y.pca$co, aes(x=Comp1, y=Comp2, label=paste("V",1:4)), col="red")
myscree(Y.pca$eig / sum(Y.pca$eig))

# Here we can see that there is lots of variation 
# in the first two axes (horizontal axis has the
# most variation, vertical axis has second most variation). 
# There is very little variation in the other axes
# (notice the scale for axis 3 only goes from -0.02 to 0.04
#   compared to the scale for axis 1 which goes from -4 to 2).


#Centering and Scaling:
#Why is the default to center and to scale?
  
#Suppose that we did not center. 
# Recall (above) that we can relate PCA to directions
# with highest covariance. When we calculate sample covariance,
# we subtract the mean from each observation. 
# If we skip this step (not centering), 
# then the first axis of the PCA would always
# be pointing towards the center of mass.
# 
# Some functions in R that calculate the PCA do not center by default.
# There might be a good reason to not center 
# (e.g., you centered a large dataset already and 
#   you are only looking at a subsample), 
# but in general, you should always center your data when doing a PCA.
# 
# Why is the default to rescale the data?
#   
# Recall the difference between correlation and covariance. 
# In correlation you rescale by dividing by the norm of each dimension. 
# This is more in line with what we're interested in.
# If one of our variable is measured in inches and then
# we decide to change that measurement to feet, 
# the variance decreases by a factor of \(12^{-2}\). 
# We don't want the result of our PCA to change based
# on the units a dimension is measured in. 
# To avoid problems like this, we rescale our data 
# so that each dimension has variance 1.
# 
# Real example for PCA:
# The dataset deug contains data on 104 French students' 
# scores in 9 subjects: Algebra, Analysis, Proba, Informatic, 
# Economy, Option1, Option2, English, Sport. We will look at the PCA of this data.
data(deug)
pca1 <- dudi.pca(deug$tab, scan = FALSE)
biplot(pca1)

#The default plotting functions in ade4 are very limitted.
#We use ggplot2 here to show what's going on.
grade <- factor(deug$result, levels=c("A+", "A", "B", "B-", "C-", "D"))
pca1.dfs <- data.frame(pca1$li, grade)
# multiply the loadings by 5 so they are more spread out
subject <- names(deug$tab)
pca1.dfl <- data.frame(5*pca1$co[,1:2], subject)
ppp + geom_point(data=pca1.dfs, aes(x=Axis1, y=Axis2, col=grade)) + 
geom_text(data=pca1.dfl, aes(x=Comp1, y=Comp2, label=subject)) 
#myscree(pca1$eig / sum(pca1$eig))

#Here we see what is called a "size effect".
#The first principal component has the largest variation, 
#and it corresponds to an overall effect of how 
#well the students did. Since we also know the students' 
#overall grades, it makes sense that those students who 
#are very far in the direction of all subjects should be
#the ones who got an "A+" and the ones on the opposite extreme got a "D".

#The second component seems to break up "Analysis" on the
#one end versus "English" on the other. 
#This component gives us an idea of what the students were good at. 
#Looking at the screeplot though, it is evident that this dimension 
#is not very well defined since there is a small jump in variance
#explained from this direction to the direction with next most variance.

#The screeplot suggests that there is one very
#important dimension (which corresponds to a size effect). 
#There is not a low rank structure left after accounting for this effect, 
#and plotting this in two dimenions tells us little more than plotting 
#only in one dimension.

#Correspondence Analysis
#Correspondence analysis takes a different sort of approach 
#to figuring out where data changes the most. 
#Instead of looking for a direction with a high variance, 
#correspondence analysis looks for the directions where the
#data is "most surprising" from a chi-squared test perspective. 
#This means that correspondence analysis is best suited for count data.

#Let's look at an example with the parathyroid data from before. 
#First we load this data and clean it some.
install.packages("parathyroid")
library(parathyroid)
data("parathyroidGenes")
# This time we won't look at the estrogen genes because there are very
# few counts of them (no observation has more than 20 counts)
gene.names <- matrix(ncol=3, byrow=T, data=c(
  #  "ESR1",  "ENSG00000091831", "estrogen",
  #  "ESR2",  "ENSG00000140009", "estrogen",
  "CASR",  "ENSG00000036828", "parathyroid",
  "VDR",   "ENSG00000111424", "parathyroid",
  "JUN",   "ENSG00000177606", "parathyroid",
  "CALR",  "ENSG00000179218", "parathyroid",
  "ORAI2", "ENSG00000160991", "parathyroid"))
gene.counts <- t(counts(parathyroidGenes)[gene.names[,2],])
colnames(gene.counts) <- gene.names[,1]
gene.dat <- parathyroidGenes@phenoData@data
# change dat$time to a numeric 24 and 48 instead of factor 24h and 48h
gene.dat$time <- as.numeric(gsub("h","",(as.character(gene.dat$time))))
dim(gene.counts)
## [1] 27  5
gene.coa <- dudi.coa(gene.counts, scannf=F, nf=2)
scatter(gene.coa)


#Or you can make it very fancy in ggplot2.
# main plot
ppp + geom_point(data=data.frame(gene.coa$li, gene.dat), 
                 aes(x=Axis1, y=Axis2, col=patient, shape=treatment, size=time)) +
  geom_text(data=gene.coa$co, aes(x=Comp1, y=Comp2, label=gene.names[,1])) +
  scale_size_area(breaks=c(24,48)) +
  scale_shape(solid=F) +
  labs(title="Correspondence analysis: select genes in Parathyroid data")
myscree(gene.coa$eig / sum(gene.coa$eig))

#Here we see what is called a "size effect".
#The first principal component has the largest variation, 
#and it corresponds to an overall effect of how well the students did. 
#Since we also know the students' overall grades, 
#it makes sense that those students who are very far in the direction
#of all subjects should be the ones who got an "A+" 
#and the ones on the opposite extreme got a "D".

#The second component seems to break up "Analysis" 
#on the one end versus "English" on the other. 
#This component gives us an idea of what the students were good at.
#Looking at the screeplot though, it is evident that this dimension 
#is not very well defined since there is a small jump in variance 
#explained from this direction to the direction with next most variance.

#The screeplot suggests that there is one very important
#dimension (which corresponds to a size effect). 
#There is not a low rank structure left after accounting
#for this effect, and plotting this in two dimenions tells
#us little more than plotting only in one dimension.



###############################################################
###############################################################
############scatterplot matrix, profile plots, screeplot#######
#https://little-book-of-r-for-multivariate-analysis.readthedocs.io/en/latest/src/multivariateanalysis.html

wine <- read.table("http://archive.ics.uci.edu/ml/machine-learning-databases/wine/wine.data",
                   sep=",")
wine

#tie plot
#A Matrix Scatterplot
#One common way of plotting multivariate data is to make a "matrix scatterplot",
#showing each pair of variables plotted against each other. 
#We can use the "scatterplotMatrix()" function from the "car" R package to do this. 
#To use this function, we first need to install the "car" R package 
#(for instructions on how to install an R package, see How to install an R package).

#Once you have installed the "car" R package, 
#you can load the "car" R package by typing:
install.packages("car")
library("car")
#You can then use the "scatterplotMatrix()" function to plot the multivariate data.
#To use the scatterplotMatrix() function, you need to give it as its 
#input the variables that you want included in the plot. 
#Say for example, that we just want to include the variables 
#corresponding to the concentrations of the first five chemicals. 
#These are stored in columns 2-6 of the variable "wine". 
#We can extract just these columns from the variable "wine" by typing:
wine[2:6]
#To make a matrix scatterplot of just these 13 variables
#using the scatterplotMatrix() function we type:
scatterplotMatrix(wine[2:6])
#In this matrix scatterplot, the diagonal cells show histograms 
#of each of the variables, in this case the concentrations of 
#the first five chemicals (variables V2, V3, V4, V5, V6).
#Each of the off-diagonal cells is a scatterplot of two of 
#the five chemicals, for example, the second cell in the 
#first row is a scatterplot of V2 (y-axis) against V3 (x-axis).
#Also see multivariate textbook page 191. 



#Profile plots
#Another type of plot that is useful is a "profile plot", 
#which shows the variation in each of the variables, 
#by plotting the value of each of the variables for each of the samples.
#The function "makeProfilePlot()" below can be used to make a profile plot. 
#This function requires the "RColorBrewer" library. 
#To use this function, we first need to install the "RColorBrewer" 
#R package (for instructions on how to install an R package, 
#see How to install an R package).
install.packages("RColorBrewer")
makeProfilePlot <- function(mylist,names)
{
  require(RColorBrewer)
  # find out how many variables we want to include
  numvariables <- length(mylist)
  # choose 'numvariables' random colours
  colours <- brewer.pal(numvariables,"Set1")
  # find out the minimum and maximum values of the variables:
  mymin <- 1e+20
  mymax <- 1e-20
  for (i in 1:numvariables)
  {
    vectori <- mylist[[i]]
    mini <- min(vectori)
    maxi <- max(vectori)
    if (mini < mymin) { mymin <- mini }
    if (maxi > mymax) { mymax <- maxi }
  }
  # plot the variables
  for (i in 1:numvariables)
  {
    vectori <- mylist[[i]]
    namei <- names[i]
    colouri <- colours[i]
    if (i == 1) { plot(vectori,col=colouri,type="l",ylim=c(mymin,mymax)) }
    else         { points(vectori, col=colouri,type="l")}
    lastxval <- length(vectori)
    lastyval <- vectori[length(vectori)]
    text((lastxval-10),(lastyval),namei,col="black",cex=0.6)
  }
}
#To use this function, you first need to copy and paste it into R. 
#The arguments to the function are a vector containing the names 
#of the varibles that you want to plot, and a list variable 
#containing the variables themselves.

#For example, to make a profile plot of the concentrations of the 
#first five chemicals in the wine samples (stored in columns 
#V2, V3, V4, V5, V6 of variable "wine"), we type:

library(RColorBrewer)
names <- c("V2","V3","V4","V5","V6")
mylist <- list(wine$V2,wine$V3,wine$V4,wine$V5,wine$V6)
makeProfilePlot(mylist,names)




#screeplots
#Principal Component Analysis
#The purpose of principal component analysis is to find the best low-dimensional 
#representation of the variation in a multivariate data set. 
#For example, in the case of the wine data set, we have 13 chemical 
#concentrations describing wine samples from three different cultivars. 
#We can carry out a principal component analysis to investigate whether
#we can capture most of the variation between samples using a smaller 
#number of new variables (principal components), 
#where each of these new variables is a linear combination of all 
#or some of the 13 chemical concentrations.
#To carry out a principal component analysis (PCA) on a multivariate data set,
#the first step is often to standardise the variables under study 
#using the "scale()" function (see above). 
#This is necessary if the input variables have very different variances, 
#which is true in this case as the concentrations of the 13 chemicals have
#very different variances (see above).
#Once you have standardised your variables, 
#you can carry out a principal component analysis using the "prcomp()" function in R.
#For example, to standardise the concentrations of the 
#13 chemicals in the wine samples, and carry out a principal 
#components analysis on the standardised concentrations, we type:
# standardise the variables
standardisedconcentrations <- as.data.frame(scale(wine[2:14])) 
# do a PCA
wine.pca <- prcomp(standardisedconcentrations)                
#You can get a summary of the principal component analysis results 
#using the "summary()" function on the output of "prcomp()":
summary(wine.pca)

# Importance of components:
#   PC1   PC2   PC3    PC4    PC5    PC6    PC7    PC8    PC9   PC10
# Standard deviation     2.169 1.580 1.203 0.9586 0.9237 0.8010 0.7423 0.5903 0.5375 0.5009
# Proportion of Variance 0.362 0.192 0.111 0.0707 0.0656 0.0494 0.0424 0.0268 0.0222 0.0193
# Cumulative Proportion  0.362 0.554 0.665 0.7360 0.8016 0.8510 0.8934 0.9202 0.9424 0.9617
# PC11   PC12    PC13
# Standard deviation     0.4752 0.4108 0.32152
# Proportion of Variance 0.0174 0.0130 0.00795
# Cumulative Proportion  0.9791 0.9920 1.00000
#This gives us the standard deviation of each component, 
#and the proportion of variance explained by each component. 
#The standard deviation of the components is stored in a named 
#element called "sdev" of the output variable made by "prcomp":
wine.pca$sdev
#[1] 2.1692972 1.5801816 1.2025273 0.9586313 0.9237035 0.8010350 0.7423128 0.5903367
#[9] 0.5374755 0.5009017 0.4751722 0.4108165 0.3215244
#The total variance explained by the components is the sum of the 
#variances of the components:
sum((wine.pca$sdev)^2)
#[1] 13
#In this case, we see that the total variance is 13,
#which is equal to the number of standardised variables (13 variables). 
#This is because for standardised data, the variance of each 
#standardised variable is 1. The total variance is equal to the sum
#of the variances of the individual variables, 
#and since the variance of each standardised variable is 1, 
#the total variance should be equal to the number of variables (13 here).
#Deciding How Many Principal Components to Retain
#In order to decide how many principal components should be retained, 
#it is common to summarise the results of a principal components 
#analysis by making a scree plot, which we can do in R using the "scree
screeplot(wine.pca, type="lines")
#The most obvious change in slope in the scree plot occurs at component 4, 
#which is the "elbow" of the scree plot. 
#Therefore, it cound be argued based on the basis of the scree plot that 
#the first three components should be retained.

#Another way of deciding how many components to retain is to use Kaiser's 
#criterion: that we should only retain principal components for which the 
#variance is above 1 (when principal component analysis was applied to standardised data).
#We can check this by finding the variance of each of the principal components:
(wine.pca$sdev)^2
#We see that the variance is above 1 for principal components 1, 2, and 3 
#(which have variances 4.71, 2.50, and 1.45, respectively). 
#Therefore, using Kaiser's criterion, we would retain the first three principal components.






###############################################################
###############################################################
####################parallel line plot#########################
#parallel line plot


#Another source:
#https://stat.ethz.ch/R-manual/R-devel/library/MASS/html/parcoord.html
install.packages("MASS")
library(MASS)
#Parallel Coordinates Plot
parcoord(state.x77[, c(7, 4, 6, 2, 5, 3)])

ir <- rbind(iris3[,,1], iris3[,,2], iris3[,,3])
parcoord(log(ir)[, c(3, 4, 2, 1)], col = 1 + (0:149)%/%50)


#install.packages(("plotly"))
library(plotly)

#Adding Dimensions
# pa1 <- plot_ly(type = 'parcoords', line = list(color = 'blue'),
#              dimensions = list(
#                list(range = c(1,5),
#                     constraintrange = c(1,2),
#                     label = 'A', values = c(1,4)),
#                list(range = c(1,5),
#                     tickvals = c(1.5,3,4.5),
#                     label = 'B', values = c(3,1.5)),
#                list(range = c(1,5),
#                     tickvals = c(1,2,4,5),
#                     label = 'C', values = c(2,4),
#                     ticktext = c('text 1', 'text 2', 'text 3', 'text 4')),
#                list(range = c(1,5),
#                     label = 'D', values = c(4,2))
#              )
# )
# pa1

#Basic Parallel Cordinates Plot
# df <- read.csv("https://raw.githubusercontent.com/bcdunbar/datasets/master/iris.csv")
# 
# pa2 <- df %>%
#   plot_ly(type = 'parcoords',
#           line = list(color = ~species_id,
#                       colorscale = list(c(0,'red'),c(0.5,'green'),c(1,'blue'))),
#           dimensions = list(
#             list(range = c(2,4.5),
#                  label = 'Sepal Width', values = ~sepal_width),
#             list(range = c(4,8),
#                  constraintrange = c(5,6),
#                  label = 'Sepal Length', values = ~sepal_length),
#             list(range = c(0,2.5),
#                  label = 'Petal Width', values = ~petal_width),
#             list(range = c(1,7),
#                  label = 'Petal Length', values = ~petal_length)
#           )
#   )
# pa2





###############################################################
###############################################################
#######################bubbleplots#############################
#bubble plots
#https://plot.ly/r/bubble-charts/
#install.packages(("plotly"))
library(plotly)
#packageVersion('plotly')

data <- read.csv("https://raw.githubusercontent.com/plotly/datasets/master/school_earnings.csv")

p <- plot_ly(data, x = ~Women, y = ~Men, text = ~School, type = 'scatter', mode = 'markers',
             marker = list(size = ~Gap, opacity = 0.5)) %>%
  layout(title = 'Gender Gap in Earnings per University',
         xaxis = list(showgrid = FALSE),
         yaxis = list(showgrid = FALSE))

p

#Setting Markers Color
p2 <- plot_ly(data, x = ~Women, y = ~Men, text = ~School, type = 'scatter', mode = 'markers',
             marker = list(size = ~Gap, opacity = 0.5, color = 'rgb(255, 65, 54)')) %>%
  layout(title = 'Gender Gap in Earnings per University',
         xaxis = list(showgrid = FALSE),
         yaxis = list(showgrid = FALSE))
p2

#Setting Multiple Colors
colors <- c('rgba(204,204,204,1)', 'rgba(222,45,38,0.8)', 'rgba(204,204,204,1)', 'rgba(204,204,204,1)', 'rgba(204,204,204,1)',
            'rgba(204,204,204,1)', 'rgba(204,204,204,1)', 'rgba(204,204,204,1)', 'rgba(204,204,204,1)', 'rgba(204,204,204,1)',
            'rgba(204,204,204,1)', 'rgba(204,204,204,1)', 'rgba(204,204,204,1)', 'rgba(204,204,204,1)', 'rgba(204,204,204,1)',
            'rgba(204,204,204,1)', 'rgba(204,204,204,1)', 'rgba(204,204,204,1)', 'rgba(204,204,204,1)', 'rgba(204,204,204,1)',
            'rgba(204,204,204,1)')
# Note: The colors will be assigned to each observations based on the order of the observations in the dataframe.

p3 <- plot_ly(data, x = ~Women, y = ~Men, text = ~School, type = 'scatter', mode = 'markers',
             marker = list(size = ~Gap, opacity = 0.5, color = colors)) %>%
  layout(title = 'Gender Gap in Earnings per University',
         xaxis = list(showgrid = FALSE),
         yaxis = list(showgrid = FALSE))
p3

#Mapping a Color Variable (Continuous)
p4 <- plot_ly(data, x = ~Women, y = ~Men, text = ~School, type = 'scatter', mode = 'markers', color = ~Gap, colors = 'Reds',
             marker = list(size = ~Gap, opacity = 0.5)) %>%
  layout(title = 'Gender Gap in Earnings per University',
         xaxis = list(showgrid = FALSE),
         yaxis = list(showgrid = FALSE))
p4

#Mapping a Color Variable (Categorical)
data$State <- as.factor(c('Massachusetts', 'California', 'Massachusetts', 'Pennsylvania', 'New Jersey', 'Illinois', 'Washington DC',
                          'Massachusetts', 'Connecticut', 'New York', 'North Carolina', 'New Hampshire', 'New York', 'Indiana',
                          'New York', 'Michigan', 'Rhode Island', 'California', 'Georgia', 'California', 'California'))

p5 <- plot_ly(data, x = ~Women, y = ~Men, text = ~School, type = 'scatter', mode = 'markers', size = ~Gap, color = ~State, colors = 'Paired',
             marker = list(opacity = 0.5, sizemode = 'diameter')) %>%
  layout(title = 'Gender Gap in Earnings per University',
         xaxis = list(showgrid = FALSE),
         yaxis = list(showgrid = FALSE),
         showlegend = FALSE)
p5


#Scaling the Size of Bubble Charts
data$State <- as.factor(c('Massachusetts', 'California', 
                          'Massachusetts', 'Pennsylvania', 'New Jersey', 
                          'Illinois', 'Washington DC',
                          'Massachusetts', 'Connecticut',
                          'New York', 'North Carolina', 'New Hampshire', 
                          'New York', 'Indiana',
                          'New York', 'Michigan', 
                          'Rhode Island', 'California', 'Georgia', 
                          'California', 'California'))

p6 <- plot_ly(data, x = ~Women, y = ~Men, text = ~School, type = 'scatter', 
              mode = 'markers', size = ~Gap, color = ~State, colors = 'Paired',
             #Choosing the range of the bubbles' sizes:
             sizes = c(10, 50),
             marker = list(opacity = 0.5, sizemode = 'diameter')) %>%
  layout(title = 'Gender Gap in Earnings per University',
         xaxis = list(showgrid = FALSE),
         yaxis = list(showgrid = FALSE),
         showlegend = FALSE)
p6

#Styled Bubble Chart
data2 <- read.csv("https://raw.githubusercontent.com/plotly/datasets/master/gapminderDataFiveYear.csv")
data_2007 <- data2[which(data2$year == 2007),]
data_2007 <- data_2007[order(data_2007$continent, data_2007$country),]
slope <- 2.666051223553066e-05
data_2007$size <- sqrt(data_2007$pop * slope)
colors <- c('#4AC6B7', '#1972A4', '#965F8A', '#FF7070', '#C61951')

p7 <- plot_ly(data_2007, x = ~gdpPercap, y = ~lifeExp, color = ~continent, size = ~size, colors = colors,
             type = 'scatter', mode = 'markers', sizes = c(min(data_2007$size), max(data_2007$size)),
             marker = list(symbol = 'circle', sizemode = 'diameter',
                           line = list(width = 2, color = '#FFFFFF')),
             text = ~paste('Country:', country, '<br>Life Expectancy:', lifeExp, '<br>GDP:', gdpPercap,
                           '<br>Pop.:', pop)) %>%
  layout(title = 'Life Expectancy v. Per Capita GDP, 2007',
         xaxis = list(title = 'GDP per capita (2000 dollars)',
                      gridcolor = 'rgb(255, 255, 255)',
                      range = c(2.003297660701705, 5.191505530708712),
                      type = 'log',
                      zerolinewidth = 1,
                      ticklen = 5,
                      gridwidth = 2),
         yaxis = list(title = 'Life Expectancy (years)',
                      gridcolor = 'rgb(255, 255, 255)',
                      range = c(36.12621671352166, 91.72921793264332),
                      zerolinewidth = 1,
                      ticklen = 5,
                      gridwith = 2),
         paper_bgcolor = 'rgb(243, 243, 243)',
         plot_bgcolor = 'rgb(243, 243, 243)')
p7





###############################################################
###############################################################
####################Chernoff faces#############################
#Chernoff faces
install.packages("aplpack")
library(aplpack)

crime <- read.csv("http://datasets.flowingdata.com/crimeRatesByState-formatted.csv")
crime[1:6,]
par(4, 4)
#faces(crime[,2:8])
crime_filled <- cbind(crime[,1:6], rep(0, length(crime$state)), crime[,7:8])
crime_filled[1:6,]
#faces(crime_filled[,2:8])
faces(crime_filled[,2:8], labels=crime_filled$state)





###############################################################
###############################################################
################(spider) web plots#############################
#(spider) web plots
#https://www.r-graph-gallery.com/142-basic-radar-chart.html
#install.packages("fmsb")
library(fmsb)

#https://personality-project.org/r/html/spider.html
#install.packages("psych")
library(psych)
#op <- par(mfrow=c(3,2))
spider(y=1,x=2:9,data=Thurstone,connect=FALSE) #a radar plot
spider(y=1,x=2:9,data=Thurstone) #same plot as a spider plot
spider(y=1:3,x=4:9,data=Thurstone,overlay=TRUE)
#make a somewhat oversized plot
par(op)
#much more detail with the variales in bfi dataset
spider(y=26:28,x=1:25,data=cor(bfi,use="pairwise"),fill=TRUE,scale=2) 


names(bfi)

