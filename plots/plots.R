#' ---
#' title: "Plots in R"
#' author: "Brandon Gonzalez"
#' date: "January 24, 2018"
#' ---

# Class 05 R graphics intro 

#My first boxplot
x <- rnorm(1000,0)
boxplot(x, horizontal = TRUE) 

#get summary of x
summary(x)

hist(x)

#Hands on session 1
df_weight <- read.table("./bimm143_05_rstats/weight_chart.txt", header=TRUE)
plot(df_weight, typ="b", pch=15, cex=1.5, lwd=2, xlab="Age (months)", ylab="Weight (kg)", main="Age and Wt", 
     ylim=c(2,10), col="blue" )

#Barplot with feature counts
df_featureCounts <- read.table("./bimm143_05_rstats/feature_counts.txt",sep='\t', header=TRUE)
barplot(height = df_featureCounts$Count,names.arg =df_featureCounts$Feature, horiz=TRUE,
        main ="Feature Counts", las=1, xlab="Expression", mar=c(7,5,5,3))

#histogram?
hist(rnorm(10000), rnorm(10000)+4, breaks=20, col=heat.colors(20))

#
df_maleFemaleCounts <- read.table("./bimm143_05_rstats/male_female_counts.txt",
                                  header = TRUE, sep='\t')
barplot(df_maleFemaleCounts$Count, col = cm.colors(nrow(df_maleFemaleCounts)),
        ylab="Counts", las=2, names.arg = df_maleFemaleCounts$Sample)

#changing margin so we can see labels
par(mar=c(5.1,10,4.1,2.1))

#Coloring by value
df_upDown <- read.table("./bimm143_05_rstats/up_down_expression.txt", header = TRUE)
palette(c("red","gray","green"))
plot(df_upDown$Condition1, df_upDown$Condition2, col=df_upDown$State)



#' coloring by **density**
meth <- read.table("./bimm143_05_rstats/expression_methylation.txt",header=TRUE,sep='\t')
#get indexes of expr greater than 0
indexes <- meth$expression > 0


#' create custom colring based on **density** of the *points*, total number of points: `r length(indexes)`

dcols.custom = densCols(meth$gene.meth[indexes], meth$expression[indexes], colramp=colorRampPalette(c("green2", "yellow","red2","purple2")))

plot(meth$gene.meth[indexes], meth$expression[indexes], pch=20, col=dcols.custom)


