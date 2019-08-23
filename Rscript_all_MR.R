##########################################################
#	EVOLUTION OF MUTATION RATE MR
#
##########################################################

# to avoid conflict with decimal values separated by points
# read as factors, use dec = "." inside read.csv function
# or convert to numeric 
#using myData$fitness <- as.numeric(levels(myData$fitness))[myData$fitness]
# another alternative is to read csv as follows:
#myData <- read.csv("data_all.csv.csv", stringsAsFactors = FALSE, header = T)

#convert MR to numeric:
myData$MR <- as.numeric(levels(myData$MR))[myData$MR]

t.limit <- 150
a <- 0.7

plot(subset(myData$MR, myData$X.step. == t.limit &
myData$level.autocorr == a)~subset(interaction
(myData$rate.change.of.optimum), 
myData$X.step. == t.limit & myData$level.autocorr == a), las = 1)

##########################################################
#	BOXPLOT USING GGPLOT2
##########################################################

library(ggplot2)

# Evolution of mutation rate and mutation effect size under scenarios
# of rate of environmental change for each treatment of noise color

#type of environment:
t.limit <- 300

# condition for data selection from table
condition <- myData$X.step. == t.limit & 
#myData$beneficial.mutations == 0.2 &
myData$mut.effect.size == 0.2 &
myData$probability.recombination == 0.5 &
myData$level.autocorr == 0

# d.type would be the levels of mut-effect-size or prob-recombination
#d.type = myData$mut.effect.size[condition]
d.type = myData$beneficial.mutations[condition]
#d.type = myData$probability.recombination

# set d.type as factor
d.type <- as.factor(d.type)

# sub.type would be the levels of rate of environmental change
sub.type = myData$rate.change.of.optimum[condition]

# val would be my response variable MR
#val = subset(myData$MR, condition)

# base 10 logarithm
# val would be my response variable MR
val = subset(log10(myData$MR), condition)

#plotting
d <- data.frame(d.type, sub.type, val)
p <- ggplot(d, aes(factor(sub.type), val)) 
p + geom_boxplot() + facet_grid(. ~ d.type) + ylim(-5,0) +
labs(x = "rate of environmental change", y = "mutation rate")

# count runs that did not end in extinction
x <- myData$N[myData$beneficial.mutations == 0.5 &
myData$mut.effect.size == 0.2 &
myData$probability.recombination == 1 &
myData$rate.change.of.optimum == 0.04 &
myData$level.autocorr == 0]

length(x[x > 0])
