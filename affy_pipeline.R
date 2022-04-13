# Analysis of affymetrix expression data, Platform: MOE430v2
# By student exam number B195515
# Script modified from FGT course 2022 by Simon Tomlinson
############
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("affy")
#BiocManager::install("limma")
#install.packages("scatterplot3d", repo="http://cran.ma.imperial.ac.uk")
#BiocManager::install('EnhancedVolcano')
library(affy)
library(limma)
library(scatterplot3d)
library(mouse4302.db) # chip
library(annotate)
library(EnhancedVolcano)
library(pheatmap)

GSE <- "GSE49448"

# Load targets file as ADF object
adf <- read.AnnotatedDataFrame("targets.txt", 
	header=TRUE, row.names=1, as.is=TRUE)
# Load CEL files listed in ADF
allsamples <- ReadAffy(filenames=pData(adf)$Filename, phenoData=adf)

# QC plots of raw data
# Density plot
png("images/QC_histogram.png", width=600, height=600)
hist(allsamples, 
	main=paste0(GSE,": Density plot of signal intensity"), 
	col=pData(adf)$Color, cex=2, lwd=3)
dev.off()

# Boxplot for quantile distribution
png("images/QC_boxplot.png", width=600, height=600)
boxplot(allsamples, 
	col=pData(adf)$Color, las=2, names=rownames(pData(adf)), 
	main=paste0(GSE,": Boxplot of signal intensity"), 
	xlab='Samples', ylab='Log Intensity', cex=2)
dev.off()

# Confirm abs/pres calls
# Distinguish biological vs quality differences
calls <- exprs(mas5calls(allsamples))
called <- data.frame(colSums(calls=='A'), colSums(calls=='P'))
write.table(called, file="MAS5calls.txt", 
	quote=F, sep="\t", col.names=NA)

# RMA-norm outputs log2 values
# eset holds normalized data
eset <- rma(allsamples)
# Matrix of normalized values
RMA_values <- exprs(eset)
colnames(RMA_values) <- rownames(pData(adf))
# Boxplot normalized values
png("images/RMA_boxplot.png", width=600, height=600)
boxplot(RMA_values, 
	col=pData(adf)$Color, las=2, names=rownames(pData(adf)), 
	main=paste(GSE,":Boxplot of normalized expression values"),
	xlab='Samples', ylab='Log Intensity', cex=2)
dev.off()

# Dependency of log-ratio & mean intensity level of variables
# x-axis = mean intensity, y-axis = log-ratio
# Most points at 0 as majority of genes don't change
png("images/MA_normalized_all.png", width=900, height=800)
mva.pairs(RMA_values, 
	main=paste(GSE,":MVA plot of CPC vs Control samples"), 
	ylim=c(-1.5,1.5))
dev.off()

# Features relationship by hierarchical clustering
# Pearson's corr coeff for HC with average linkage
HC <- hclust(as.dist(1-cor(RMA_values, method="pearson")), method="average")
png("images/HC_pearsonavg_normalized.png", width=600, height=600)
plot(HC)
dev.off()

# Determine artificial feats explaining variability
# Transpose values matrix and scale variables to have unit variance
pca <- prcomp(t(RMA_values), scale=T)
png("images/PCA_normalized.png", width=600, height=600)
s3d <- scatterplot3d(pca$x[,1:3], pch=19, color=rainbow(nrow(adf)))
s3d.coords <- s3d$xyz.convert(pca$x[,1:3])
text(s3d.coords$x, s3d.coords$y, 
	labels=colnames(RMA_values), 
	pos=3, offset=0.5)
dev.off()

# Rename samples in eset
sampleNames(eset) <- rownames(pData(adf))

#-------------SAMPLE ANNOTATION
eset@annotation
ls("package:mouse4302.db")
ID <- featureNames(eset)
Symbol <- getSYMBOL(ID, "mouse4302.db")
Name <- as.character(lookUp(ID, "mouse4302.db", "GENENAME"))
tmp <- data.frame(ID=ID, Symbol=Symbol, Name=Name, stringsAsFactors=F)
tmp[tmp=="NA"] <- NA
# Annotation table as feature data of eset
fData(eset) <- tmp
keytypes(mouse4302.db) #annotation keys
# Annotate eset with ENTREZID - map to gene signatures
# Key selection, PROBEID as primary key
res <- select(mouse4302.db, keys=rownames(eset), 
	columns=c("ENTREZID","PFAM", "GENENAME", "ENSEMBL","SYMBOL"), 
	keytype="PROBEID")
# Table join of eset and res
# Exact mapping as gene list may be scrambled
idx <- match(rownames(eset), res$PROBEID)
# Index sets phenotypic data in eset
fData(eset) <- res[idx,]
# Mask false ENTREZID rows
eset_t <- eset[is.na(fData(eset)$ENTREZID)==0,]

#-------------FOLD CHANGE

# log2 values to actual values
RMA_values10 <- 2**RMA_values
# Simple FC table
treated_mean <- apply(RMA_values10[, c("CPC.1","CPC.2","CPC.3")], 1, mean)
control_mean <- apply(RMA_values10[, c("Ctl.1","Ctl.2","Ctl.3")], 1, mean)
fc_mean <- treated_mean/control_mean
foldchange <- cbind(Symbol, Name, treated_mean, control_mean, fc_mean)
write.table(foldchange, file="foldchange.txt", 
	quote=F, sep="\t", col.names=NA)

#-------------LIMMA ANALYSIS

# Construct design matrix, response~model
# (gene profiles~formula used to predict the response)
design <- model.matrix(~-1+factor(c(1,1,1,2,2,2)))
colnames(design) <- unique(pData(adf)$Group)
# Contrast matrix: which group comparison(s)
contrastmatrix <- makeContrasts(CPC-CTL, levels=design)

# Fit the model lmFit(response, model)
fit <- lmFit(eset_t, design)
# Contrasts w/borrowed variance to moderate t-statistics
fit2 <- eBayes(contrasts.fit(fit, contrastmatrix))

# Top DEGs as specified in contrastmatrix
# 'coef' is contrastmatrix column
# 'fdr' is FDR statistical adjustment
# 'number' is maximum number of genes to list
limmaresults <- topTable(fit2, coef=1, adjust="fdr", number=nrow(eset_t))
write.table(limmaresults, "limmaresults.csv", 
	sep=",", col.names=NA, row.names=TRUE)

# DEGs Volcano plot with cutoffs
png("images/volcano_DEgenes_padj.png", width=1000, height=800)
EnhancedVolcano(
	limmaresults, x="logFC", y="adj.P.Val", 
	lab=limmaresults$SYMBOL, 
	title=paste0(GSE,": CPC versus Control"), 
	subtitle="Limma results", 
	labSize=6.0, pointSize=1.0, 
	pCutoff=0.001, FCcutoff=4, 
	legendPosition='right', legendLabSize=12, legendIconSize=4.0, 
	xlab=bquote(~Log[2]~ 'fold change'), 
	ylab=bquote(~Log[10]~ 'adj.P.val'), 
	ylim=c(0,-log10(0.00001)))
dev.off()

# Overlapping DEGs between groups
#classify <- classifyTestsF(fit2)
classify <- decideTests(fit2)
ded <- unname(colSums(classify!=0))
png("images/vennDiagram_classified.png", width=600, height=600)
vennDiagram(classify)
dev.off()

#-----------FUNCTIONAL ENRICHMENT ANALYSIS

# Load molecular signature db
Mm.H <- readRDS("Mm.h.all.v7.1.entrez.rds")
# Filter log2 values significant DEGs from decideTests
# Convert Mm.H ID list into eset_t index
# ids2indices(gene.sets, identifiers, remove.empty=TRUE)
# CAMERA competitive test for Functional Enrichment
# Reuse the design & contrast matrix
RMA_values_filt <- RMA_values[rownames(limmaresults)[1:ded],]
H.indices <- ids2indices(Mm.H, limmaresults[1:ded,]$ENTREZID)
camera_results <- camera(RMA_values_filt, index=H.indices,
	design=design, contrast=contrastmatrix[,1],
	weights=-log10(limmaresults[1:ded,"adj.P.Val"]))
write.table(camera_results, "func.enrichment_camera.txt", sep="\t")

# Heatmap of CAMERA results
cm_annotate <- data.frame("Direction"=camera_results$Direction)
rownames(cm_annotate) <- rownames(camera_results)
png("images/heatmap_camera.png", height=1000, width=800)
pheatmap(as.matrix(camera_results[,c("PValue","FDR")]),
	main=paste0(GSE,": CAMERA Enriched Genes in CPC group over control"),
	cluster_cols=F, annotation_row=cm_annotate, 
	cutree_rows=4, fontsize=10)
dev.off()

# Barplot modified from https://github.com/YuLab-SMU/DOSE/issues/20
library(ggplot2)
library(forcats)
cameraa <- camera_results
cameraa$type <- "upregulated"
cameraa$type[cameraa$Direction=="Down"] <- "downregulated"
cameraa$term <- rownames(cameraa)
png("images/barplot_camera.png", height=1000, width=800)
p <- ggplot(cameraa, aes(x=NGenes/ded, y=fct_reorder(term, NGenes))) +
	geom_point(aes(size=NGenes/ded, color=FDR)) +
	theme_bw(base_size=14) +
	scale_colour_gradient(limits=c(0,1), low="red") +
	ylab(NULL) + xlab("GeneRatio") + labs(size="GeneRatio") +
	ggtitle(paste0(GSE,": CAMERA functional enrichment"))
p + facet_grid(.~type)
dev.off()

# Other methods
# mroast_results <- mroast(eset_t, index=H.indices, design=design,
# contrast=contrastmatrix[,1], adjust.method="BH")
# write.table(mroast_results, "func.enrichment_mroast.txt", sep="\t")
# romer_results <- romer(eset_t, index=H.indices, design=design,
# contrast=contrastmatrix[,1])
# write.table(romer_results, "func.enrichment_romer.txt", sep="\t")

# Optional: Extract the model from the fit if using the same one 
# in both the DEG analysis and the enrichment analysis
#sv <- squeezeVar(fit$sigma^2, df=fit$df.residual)

#---------------SHINY DATA

# Top 50 DEGs by Limma
library(dplyr)
PROBEID <- data.frame(limmaresults$PROBEID)
colnames(PROBEID) <- "PROBEID"
expression <- left_join(
	limmaresults[1:50,], 
	cbind(PROBEID, RMA_values10[rownames(limmaresults),]), 
	by="PROBEID", keep=TRUE, copy=TRUE)
rownames(expression) <- with(expression, paste(SYMBOL, PROBEID.x, sep="/"))
# Select only cols with sample names, convert to num
expression <- expression[,(ncol(expression)-6+1):ncol(expression)]
expression[] <- lapply(expression, as.numeric)
save(expression, file="expression.Rdata")

# Read targets file as a DF
experiment <- read.table("targets.txt", 
	header=T, as.is=T, row.names=1)
save(experiment, file="experiment.Rdata")

# Volcano plot (-log10(FDR) vs log2FC) (csv, RData)
differential <- limmaresults[, c("ENSEMBL","SYMBOL","logFC","P.Value","adj.P.Val","B")]
differential$minus_log10_Pval <- -log10(limmaresults$adj.P.Val)
#differential$sig <- as.factor(
#	abs(differential$logFC) > 2 & differential$adj.P.Val < 0.01)
write.table(differential, "differential.csv", 
	sep=",", col.names=NA, row.names=TRUE)
save(differential, file="differential.Rdata")

#--------------SAVE

#save.image("~/FGT/ICA1/workspace.RData")
#writeLines(capture.output(sessionInfo()), "sessionInfo.txt")