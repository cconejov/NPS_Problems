#libraries

library(tidyverse)

# load data

problemphones <- read.table("temps-7.txt", header = TRUE)
pastphones <- read.delim("temps-other.txt", header = TRUE)


#rename variables

names(problemphones)[names(problemphones) == "x"] <- "problemtemps"
names(pastphones)[names(pastphones) == "x"] <- "pasttemps"

## maintain this order of temperatures of problematic phones and the temperatires of past phones 


#attaching for ease of use

attach(problemphones)
attach(pastphones)

#sumamries of data

summary(problemtemps)
summary(pastphones)


#histograms

p1 <- hist(problemtemps)
p2 <- hist(pasttemps)

plot( p1, col=rgb(0,0,1,1/4), xlim=c(0,65))  # first histogram
plot( p2, col=rgb(1,0,0,1/4), xlim=c(0,65), add=T)  

#pastemps in red, so lower, makes since they didn't burn


