# NCETS

NCETS(Negative Control Exposure based on time-series studies),is an R package for causal effect estimation of environmental exposure on health outcome. Just taking advantage of a post-outcome exposure as an auxiliary variable, the NCETS  can obtain an unbiased and robust causal effect estimation of exposure on outcome. 


# Installation
It is easy to install the development version of NCETS package using the 'devtools' package. The typical install time on a "normal" desktop computer is less than one minute.

```
# install.packages("devtools")
library(devtools)
install_github("yuyy-shandong/NCETS")
```


# Usage
There are three main functions in NCETS package. One is diff_test for testing causal effect among time-varying Exposure. The second is ncets_cat which could eliminate the unmeasured confounders and estimate causal effect for categorical outcome. And the last one is ncets_con which eliminate the unmeasured confounders and estimate causal effect for continuous outcome.
You can find the instructions by '?diff_test', '?ncets_cat'  and  '?ncets_con'.

library(NCETS)

?diff_test

?ncets_cat

?ncets_con


# Example


```
u <- rnorm(1000,0,1)
x1 <- 0.5*u +rnorm(1000,0,1)
x3 <- 0.5*u +rnorm(1000,0,1)
y <- 2*x1 + 1*u +rnorm(1000,0,1)
data <- data.frame (u,x1,x3,y)
model <- ncets_con (y~x1+x3, data = data, sdmethod ="normal",x1_name = "x1",x3_.name = "x3",boots_no = NULL)
model
```


# Development
This R package is developed by Yuanyuan Yu, HongKai Li and Fuzhong Xue.
