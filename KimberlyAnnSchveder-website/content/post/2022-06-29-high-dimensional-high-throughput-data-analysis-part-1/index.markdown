---
title: "High Dimensional High Throughput Data Analysis part 1"
author: "Kimberly Schveder"
date: '2022-06-29'
output: pdf_document
categories:
- R
- statistics
- big data
tags: []
slug: high-dimensional-high-throughput-data-analysis-part-1
---
In Spring of 2020 I was finishing my last semester of my masters in statistics program. My colleagues and I did not get to continue to learn this subject material in the classroom due to the lock down for the COVID-19 pandemic in March. 

I have decided to re-post the notes I took from Dr. Datta's class on "High-Dimensional, High-Throughput Data Analysis" - an important topic for data scientists to learn of. It was a nice summary course of the masters in statistics program. I have written out my notes in a Word Doc below and will break them into a series of blog posts.

The class went over: 
- [x] Reviewing of univariate and multivariate linear regression. 
- [x] Multiple Hypotheses testing and multiplicity adjustment (family-wise error rate control for independent and dependent p-values, false discovery rate control, etc). 
- [x] High Dimensional data visualization. 
- [x] Dimension reduction methods (sufficiency and equivalence, principal components and principal components regression, penalized regression, sliced inverse regression, etc.)
- [x] Model selection and estimator selection (feature selection, cross-validation, complexity selection)
- [x] Machine learning methods (introduction to come classification techniques and clustering techniques)
- [x] Aggregation of estimators and classifiers (boosting, bagging, other ensemble methods)
- [x] Graphical and network models (random graphs, graphs as statistical models, Gaussian graphical models)

Our first lecture, we went over what was covered in the article "Challenges of Big Data analysis", by Fan, Jianqing et. al.  in the National Science Review 1:293-314, 2014, doi: 10.1093/nsr/nwt032. 

The highlights I have about this article are the following: 
- [x] Goals and challenges of analyzing big data and high dimensional data analysis are to develop effective methods that can predict accurately future observations and to gain insight to relationship(s) between the features and responses for science. 
- [x] Due to large sample size, big data gives rise to additional goals: to understand heterogeneity and commonality across different sub-populations. That is, big data gives promises for: 
- (i) exploring hidden structures of each sub-population of the data that is not feasible and might even be treated as 'outliers' and 
- (ii) extracting important common features across many sub-populations even when those are large individual variations.
- [x] However, with big data there can be big challenges to deal with. There are 3 main challenges: 
- (i) high dimensionality brings noise accumulation and spurious correlations and incidental homogeneity; 
- (ii) high dimensionality combined with large samples size creates issues such as heavy computational cost and algorithmic instability;
- (iii) the massive samples in Big Data are typically aggregated from multiple sources at different time points using different technologies. This creates issues of heterogeneity, experimental variations and statistical biases, and requires us to develop more adaptive and robust procedures. 



There are a variety of linear regression modeling that were covered in this class. Here, I summarize the basics.


**Part 1: Introduction and Regression Basics**
A simple linear regression model is defined as 

`$y_i = \beta_0 + \beta_1 x_i + \epsilon_i$`

where i = 1, ..., n, 
`$x_i$` is fixed, `$\epsilon_i	\sim  N(0, \sigma^2)$`, `$y_i	\sim  N(\beta_0 + \beta_1 X_i, \sigma^2)$` (independently and identically normally distributed about the regression line equation). 

We can demonstrate a simple linear model by looking at the dataset cars programmed into Rstudio. What is the distance traveled in miles based on the speed of the car? 


```r
head(mtcars)
```

```
##                    mpg cyl disp  hp drat    wt  qsec vs am gear carb
## Mazda RX4         21.0   6  160 110 3.90 2.620 16.46  0  1    4    4
## Mazda RX4 Wag     21.0   6  160 110 3.90 2.875 17.02  0  1    4    4
## Datsun 710        22.8   4  108  93 3.85 2.320 18.61  1  1    4    1
## Hornet 4 Drive    21.4   6  258 110 3.08 3.215 19.44  1  0    3    1
## Hornet Sportabout 18.7   8  360 175 3.15 3.440 17.02  0  0    3    2
## Valiant           18.1   6  225 105 2.76 3.460 20.22  1  0    3    1
```

```r
summary(mtcars)
```

```
##       mpg             cyl             disp             hp       
##  Min.   :10.40   Min.   :4.000   Min.   : 71.1   Min.   : 52.0  
##  1st Qu.:15.43   1st Qu.:4.000   1st Qu.:120.8   1st Qu.: 96.5  
##  Median :19.20   Median :6.000   Median :196.3   Median :123.0  
##  Mean   :20.09   Mean   :6.188   Mean   :230.7   Mean   :146.7  
##  3rd Qu.:22.80   3rd Qu.:8.000   3rd Qu.:326.0   3rd Qu.:180.0  
##  Max.   :33.90   Max.   :8.000   Max.   :472.0   Max.   :335.0  
##       drat             wt             qsec             vs        
##  Min.   :2.760   Min.   :1.513   Min.   :14.50   Min.   :0.0000  
##  1st Qu.:3.080   1st Qu.:2.581   1st Qu.:16.89   1st Qu.:0.0000  
##  Median :3.695   Median :3.325   Median :17.71   Median :0.0000  
##  Mean   :3.597   Mean   :3.217   Mean   :17.85   Mean   :0.4375  
##  3rd Qu.:3.920   3rd Qu.:3.610   3rd Qu.:18.90   3rd Qu.:1.0000  
##  Max.   :4.930   Max.   :5.424   Max.   :22.90   Max.   :1.0000  
##        am              gear            carb      
##  Min.   :0.0000   Min.   :3.000   Min.   :1.000  
##  1st Qu.:0.0000   1st Qu.:3.000   1st Qu.:2.000  
##  Median :0.0000   Median :4.000   Median :2.000  
##  Mean   :0.4062   Mean   :3.688   Mean   :2.812  
##  3rd Qu.:1.0000   3rd Qu.:4.000   3rd Qu.:4.000  
##  Max.   :1.0000   Max.   :5.000   Max.   :8.000
```

```r
fit <- lm(mpg ~ wt, data = mtcars)
fit
```

```
## 
## Call:
## lm(formula = mpg ~ wt, data = mtcars)
## 
## Coefficients:
## (Intercept)           wt  
##      37.285       -5.344
```

```r
plot(mpg ~ wt, data = mtcars)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/cars-1.png" width="672" />
Diagnostics and multiple linear regression next.



The next blog post will contain information on the gist of multivariate statistical methods. 

**Coming up.....**

For a multivariate regression model is defined with matrices organizing the information, 
`$\underline{Y} = \begin{pmatrix}
    y_{1}\\
    ...\\
    y_{n}
  \end{pmatrix} = X_n \underline{\beta} + \underline{\epsilon_i} $`



