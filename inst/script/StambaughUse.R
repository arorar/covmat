rm(list = ls())

library(factorAnalytics)
library(robust)
library(ggplot2)
library(RCurl)
library(VIM)

source("Stambaugh.R")

#test
data <- read.csv("../../data/hfunds5.ue.ts.csv",header=TRUE)
plot.missing(data,3)
plot.missing(data,4)


data <- xts(data[,-1],order.by=as.Date(data[,1]))
data <- data[,-which(colnames(data) == "EVENT.1")]
cov.control <- covRob.control(estim="mcd",alpha=0.9)

models <- stambaugh.fit(data, cov.control=cov.control)
plot(models,1)
plot(models,2,0.975)

# #test
data <- read.csv("../../data/normal.hectic.csv",header=TRUE)
data <- data[,-which(colnames(data) == "EMU")]
data <- xts(data[,-1],order.by=as.Date(data[,1],format="%m/%d/%Y"))
data <- data[-(1:60),]

models <- stambaugh.fit(data)
plot(models,1)
plot(models,2)


