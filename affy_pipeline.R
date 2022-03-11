# Analysis of affymetrix expression data 
# Platform: 
# 

# Install packages to load and analyze data
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("affy")
#BiocManager::install("limma")
install.packages("scatterplot3d", repo="http://cran.ma.imperial.ac.uk")
library(affy)
library(limma)
library(scatterplot3d)

# Load the target file into adf, AnnotatedDataFrame object
adf <- read.AnnotatedDataFrame("targets.txt", header=TRUE, row.names=1, as.is=TRUE)

# Load the expression values of all CEL files listed in the target file
mydata <- ReadAffy(filenames=pData(adf)$Filename, phenoData=adf)

# Make quality control plots of the raw data
# Density plot
png("QC_histogram")
hist(mydata, main='Density plot of signal intensity', col=pData(adf)$Color) # Insert title and label samples by color
dev.off()
# Boxplot for quantile distribution
png("QC_boxplot")
boxplot(mydata, col=pData(adf)$Color, las=2, names=rownames(pData(adf)), main='Boxplot of samples', xlab='Samples', ylab='Log Intensity')
dev.off()
# Double check present and absent calls, difference between groups vs sample quality
calls <- exprs(mas5calls(mydata))
called <- data.frame(colSums(calls=='A'), colSums(calls=='P'))

# Normalize data using RMA, returns log2-transformed normalized values
# eset holds normalized data
eset <- rma(mydata)
# 'values' hold matrix of normalized expression values
values <- exprs(eset)
# Boxplot again but for normalized values
png("RMA_boxplot")
boxplot(values, col=pData(adf)$Color, las=2, names=rownames(pData(adf)), main='Boxplot of normalized expression values', xlab='Samples', ylab='Log Intensity')
dev.off()

# MA plot to check dependencies between the log-ratio of two variables and their mean intensity level
# x-axis = mean intensity, y-axis = log-ratio
# Except most points (genes) in MA ~0 as majority of genes don't change
png("MA_normalized_all", width = 4.5, height = 4, units = 'in', res = 300)
mva.pairs(values, labels=rownames(pData(adf)), main="MVA plot of CPC vs Control samples", cex.main=1, cex.axis=0.5, cex=0.5)
dev.off()

# Determine relationship between features (samples) by hierarchical clustering
# Change column header from filename to sample name
colnames(values) <- rownames(pData(adf))
# Use Pearson's correlation coefficient to perform HC with average linkage
HC <- hclust(as.dist(1-cor(values, method="pearson")), method="average")
png("HC_pearsonavg_normalized")
plot(HC)
dev.off()

# Determine the artificial feats in dataset that can account for observed variability in the observed features
# Transpose values matrix and scales variables to have unit variance before analysis
pca <- prcomp(t(values), scale=T)
png("PCA_normalized")
s3d <- scatterplot3d(pca$x[,1:3], pch=19, color=rainbow(nrow(adf)))
s3d.coords <- s3d$xyz.convert(pca$x[,1:3])
text(s3d.coords$x, s3d.coords$y, labels=colnames(values), pos=3, offset=0.5)
dev.off()

# Convert normalized expression values from log2
values10 <- 2**values
# Check conversion
values[1:20,]
values10[1:20,]

# Build a fold change table
# Check order of Sample Name and Probe IDs
mysamples <- sampleNames(eset)
probesets <- probeNames(mydata)
mysamples
probesets[1:10]
CPC.mean <- apply(values10[, c("CPC.1","CPC.2","CPC.3")], 1, mean)
CTL.mean <- apply(values10[, c("Ctl.1","Ctl.2","Ctl.3")], 1, mean)
CPC_CTL_FC <- CPC.mean/CTL.mean
all_data <- cbind(CPC.mean, CTL.mean, CPC_CTL_FC)
colnames(all_data)
write.table(all_data, file="group_means.txt", quote=F, sep="/t", col.names=NA)

