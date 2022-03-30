# Analysis of affymetrix expression data, Platform: MOE430v2
# Script modified from original tutorial script from the Functional Genomic Technologies course 2022 by Simon Tomlinson

# Install packages to load and analyze data
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("affy")
#BiocManager::install("limma")
#install.packages("scatterplot3d", repo="http://cran.ma.imperial.ac.uk")
library(affy)
library(limma)
library(scatterplot3d)
library(mouse4302.db) #load library for the chip
library(annotate)
BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)
library(pheatmap)

# Load targets file as ADF object
adf <- read.AnnotatedDataFrame("targets.txt", header=TRUE, row.names=1, as.is=TRUE)

# Load CEL files listed in ADF
mydata <- ReadAffy(filenames=pData(adf)$Filename, phenoData=adf)

# QC plots of raw data
# Density plot
png("QC_histogram")
hist(mydata, main='Density plot of signal intensity', col=pData(adf)$Color)
dev.off()
# Boxplot for quantile distribution
png("QC_boxplot")
boxplot(mydata, col=pData(adf)$Color, 
	las=2, names=rownames(pData(adf)), 
	main='Boxplot of samples', 
	xlab='Samples', ylab='Log Intensity')
dev.off()
# Confirm abs/pres calls, difference between groups vs sample quality
calls <- exprs(mas5calls(mydata))
called <- data.frame(colSums(calls=='A'), colSums(calls=='P'))
write.table(called, file="calls.txt", 
	quote=F, sep="\t", col.names=NA)

# RMA-normalization outputs log2-transformed values
# eset holds normalized data
eset <- rma(mydata)
# 'values' hold matrix of normalized expression values
values <- exprs(eset)
# Boxplot of normalized values
png("RMA_boxplot")
boxplot(values, col=pData(adf)$Color, 
	las=2, names=rownames(pData(adf)), 
	main='Boxplot of normalized expression values',
	xlab='Samples', ylab='Log Intensity')
dev.off()

# Dependency of log-ratio & mean intensity level of variables
# x-axis = mean intensity, y-axis = log-ratio
# Except most points (genes) in MA ~0 as majority of genes don't change
png("MA_normalized_all", width = 4.5, height = 4, units = 'in', res = 300)
mva.pairs(values, labels=rownames(pData(adf)), 
	main="MVA plot of CPC vs Control samples", 
	cex.main=1, cex.axis=0.5, cex=0.5)
dev.off()

# Features relationship by hierarchical clustering
# Change column header from filename to sample name
colnames(values) <- rownames(pData(adf))
# Use Pearson's correlation coefficient to perform HC with average linkage
HC <- hclust(as.dist(1-cor(values, method="pearson")), method="average")
png("HC_pearsonavg_normalized")
plot(HC)
dev.off()

# Determine dataset artificial feats explaining variability
# Transpose values matrix and scales variables to have unit variance
pca <- prcomp(t(values), scale=T)
png("PCA_normalized")
s3d <- scatterplot3d(pca$x[,1:3], pch=19, color=rainbow(nrow(adf)))
s3d.coords <- s3d$xyz.convert(pca$x[,1:3])
text(s3d.coords$x, s3d.coords$y, 
	labels=colnames(values), 
	pos=3, offset=0.5)
dev.off()

# Convert normalized expression values from log2 and check
values10 <- 2**values
values[1:20,]
values10[1:20,]

# Rename samples in eset
sampleNames(eset) <- rownames(pData(adf))

# Annotating eset with gene names using MOE430v2 annotations
eset@annotation
# ls("package:mouse4302.db")
# Get transcript cluster IDs from eset, gene symbol, & gene name
ID <- featureNames(eset)
Symbol <- getSYMBOL(ID, "mouse4302.db")
Name <- as.character(lookUp(ID, "mouse4302.db", "GENENAME"))

# Build annotation table, fix padding with 'NA' characters
tmp <- data.frame(ID=ID, Symbol=Symbol, Name=Name, stringsAsFactors=F)
tmp[tmp=="NA"] <- NA
# Assign annotation table as feature data of eset
fData(eset) <- tmp

# Check order of Sample Name and Probe IDs
mysamples <- sampleNames(eset)
probesets <- probeNames(mydata)
mysamples
probesets[1:10]
# Build FC table
CPC.mean <- apply(values10[, c("CPC.1","CPC.2","CPC.3")], 1, mean)
CTL.mean <- apply(values10[, c("Ctl.1","Ctl.2","Ctl.3")], 1, mean)
CPC_CTL_FC <- CPC.mean/CTL.mean
FC_GSE49448 <- cbind(Symbol, Name, CPC.mean, CTL.mean, CPC_CTL_FC)
colnames(FC_GSE49448)
write.table(FC_GSE49448, file="FC_GSE49448.txt", 
	quote=F, sep="\t", col.names=NA)

# Construct design matrix, response~model
# response: gene profiles
# model: the formula used to predict the response
design <- model.matrix(~-1+factor(c(1,1,1,2,2,2)))
colnames(design) <- unique(pData(adf)$Group)
# Construct contrast matrix
# Specifies which comparisons are to be made between groups
contrastmatrix <- makeContrasts(CPC-CTL, levels=design)
# Fit the model lmFit(response, model)
fit <- lmFit(eset,design)
# Do contrasts w/borrowed variance to moderate the t-stats
fit2 <- eBayes(contrasts.fit(fit, contrastmatrix))

# Top DEGs as specified in contrastmatrix
# 'coef' is contrastmatrix column
# 'fdr' is FDR statistical adjustment
# 'number' is maximum number of genes to list
limmaresults <- topTable(fit2, coef=1, adjust="fdr", number=nrow(eset))
write.table(limmaresults, "limmaresults.txt", sep="\t")

# Volcano plot of DEGs with cutoffs
limma_df <- as.data.frame(limmaresults)
png("volcano_DEgenes_padj.png", width=1000, height=800)
EnhancedVolcano(
	limma_df, x="logFC", y="adj.P.Val", lab=limma_df$Symbol, 
	title="GSE49448: CPC versus Control", subtitle="Limma results", 
	labSize=6.0, pointSize=1.0, pCutoff=0.001, FCcutoff=4, 
	legendPosition='right', legendLabSize=12, legendIconSize=4.0, 
	xlab=bquote(~Log[2]~ 'fold change'), ylab=bquote(~Log[10]~ 'adj.P.val'), 
	ylim=c(0,-log10(0.00001))
	)
dev.off()

# Overlapping DEGs between groups
classify <- classifyTestsF(fit2)
png("vennDiagram_classified.png", width=1000, height=1000)
vennDiagram(classify)
dev.off()

# Load molecular signature db
Mm.H <- readRDS("Mm.h.all.v7.1.entrez.rds")

# Show annotation keys in the database
#keytypes(mouse4302.db)
# Annotate eset with ENTREZID to map to gene signatures
# Select some keys from the package, with PROBEID as primary key
res <- select(mouse4302.db, keys=rownames(eset), 
columns=c("ENTREZID","ENSEMBL","SYMBOL"), 
	keytype="PROBEID")

# Do a table join of eset and res
# Match is done for exact mapping, as gene list may be scrambled
idx <- match(rownames(eset), res$PROBEID)

# Use the index to set phenotypic data in eset
# fData(eset) before: (probe)ID, Symbol, and Name
# fData(eset) after: PROBEID, ENTREZID, ENSEMBL, SYMBOL
fData(eset) <- res[idx,]

# Remove rows without ENTREZID annotation ie. mask false values
eset_t <- eset[is.na(fData(eset)$ENTREZID)==0,]

# Reuse the design matrix and contrast matrix
# for functional enrichment analysis

# Convert ID list in mol sig db into an index in eset_t
# ids2indices(gene.sets, identifiers, remove.empty=TRUE)
H.indices <- ids2indices(Mm.H, fData(eset_t)$ENTREZID)

# CAMERA competitive test for Functional Enrichment
camera_results <- camera(eset_t, index=H.indices, 
	design=design, 
	contrast=contrastmatrix[,1])
write.table(camera_results, 
	"func.enrichment_camera.txt", 
	sep="\t")

# Plot heatmap of CAMERA results
cm_matrix <- apply(camera_results[,c("PValue","FDR")], c(1,2), "*",-1)
cm_annotate <- data.frame("Direction"=camera_results$Direction)
rownames(cm_annotate) <- rownames(cm_matrix)
png("heatmap_camera.png", height=1000, width=800)
pheatmap(as.matrix(cm_matrix), 
	main="GSE49448: CAMERA Enriched Genes in CPC group over control",
	cluster_cols=F, annotation_row=cm_annotate, 
	cutree_rows=4, fontsize=10)
dev.off()

# Other methods
mroast_results <- mroast(eset_t, index=H.indices, design=design,
contrast=contrastmatrix[,1], adjust.method="BH")
write.table(mroast_results, "func.enrichment_mroast.txt", sep="\t")
#
romer_results <- romer(eset_t, index=H.indices, design=design,
contrast=contrastmatrix[,1])
write.table(romer_results, "func.enrichment_romer.txt", sep="\t")

# Optional: Extract the model from the fit if using the same one in both the differential gene analysis and the enrichment analysis
#sv <- squeezeVar(fit$sigma^2, df=fit$df.residual)

#save.image("~/FGT/ICA1/workspace.RData")
