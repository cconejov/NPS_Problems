## Exercise 1 - Description --------------------------------------------------

## Exercise 1: Assess S claims of battery failure from temperatures from a sample of previous batteries that experienced failure and from a sample
## in general 

## a. Perform a kernel density estimation for temps-7 and temps-other using what you consider is the most adequate bandwidth. Since the temperatures are positive,
## is it required to perform any transformation?.

## b. Is there any important difference on the results from considering the LSCV selector over the DPI selector?

## c. It seems that in temps-7 there is a secondary mode. Compute a kernel derivative estimation for temps-7 and temps-other using what you consider are the 
## most adequate bandwidths.

## d. Precisely determine the location of the extreme points

## e. Check with a kernel second derivative that the extreme points are actually modes.

## Preliminary Work--------------------------------------------------

# Load data

problemphones <- read.table("temps-7.txt", header = TRUE)
pastphones <- read.delim("temps-other.txt", header = TRUE)

# Attaching data for ease of use

attach(problemphones)
attach(pastphones)

# Rename variables

names(problemphones)[names(problemphones) == "x"] <- "problemtemps"
names(pastphones)[names(pastphones) == "x"] <- "pasttemps"


# Summaries of data

summary(problemtemps)

# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 6.85   12.51   15.68   17.34   19.52   62.84 

summary(pastphones)

# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 5.39   12.69   17.22   17.76   22.01   46.86


# Histograms - default bandwidths

p1 <- hist(problemtemps)
p2 <- hist(pasttemps)

plot( p1, col=rgb(0,0,1,1/4), xlim=c(0,65))  # first histogram
plot( p2, col=rgb(1,0,0,1/4), xlim=c(0,65), add=T)  

# Density plots - default bandwidths

par(mfrow = c(1,2))
plot(density(x = temps_7), xlim = c(-5,65), main = 'Problematic battery')
plot(density(x = temps_other), xlim = c(-5,65), main = 'Past battery')
par(mfrow = c(1,1))

## what is it using sas the defualt bandwidth



