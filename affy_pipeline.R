# Analysis of affymetrix expression data, Platform: MOE430v2
# Student B195515
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

adf <- read.AnnotatedDataFrame("targets.txt", 
	header=TRUE, row.names=1, as.is=TRUE) # load targets file
allsamples <- ReadAffy(filenames=pData(adf)$Filename, phenoData=adf) # load CEL files

png("images/QC_histogram.png", width=600, height=600)
hist(allsamples, 
	main=paste0(GSE,": Density plot of signal intensity"), 
	col=pData(adf)$Color, cex=2, lwd=3) # density plot
dev.off()

png("images/QC_boxplot.png", width=600, height=600)
boxplot(allsamples, 
	col=pData(adf)$Color, las=2, names=rownames(pData(adf)), 
	main=paste0(GSE,": Boxplot of signal intensity"), 
	xlab='Samples', ylab='Log Intensity', cex=2) # quantile plot
dev.off()

calls <- exprs(mas5calls(allsamples)) # abs/pres calls for QC
called <- data.frame(colSums(calls=='A'), colSums(calls=='P'))
write.table(called, file="MAS5calls.txt", 
	quote=F, sep="\t", col.names=NA)

eset <- rma(allsamples) # normalize intensity
RMA_values <- exprs(eset) # expression data (log2)
colnames(RMA_values) <- rownames(pData(adf))
sampleNames(eset) <- rownames(pData(adf))

png("images/RMA_boxplot.png", width=600, height=600)
boxplot(RMA_values, 
	col=pData(adf)$Color, las=2, names=rownames(pData(adf)), 
	main=paste(GSE,":Boxplot of normalized expression values"),
	xlab='Samples', ylab='Log Intensity', cex=2) # quantile plot normalized
dev.off()

png("images/MA_normalized_all.png", width=900, height=800)
mva.pairs(RMA_values, 
	main=paste(GSE,":MVA plot of CPC vs Control samples"), 
	ylim=c(-1.5,1.5)) # MA plot: mean intensity level vs log-ratio (mostly 0)
dev.off()

# Hierarchical clustering using Pearson's correlation coefficient
HC <- hclust(as.dist(1-cor(RMA_values, method="pearson")), method="average")
png("images/HC_pearsonavg_normalized.png", width=600, height=600)
plot(HC)
dev.off()

pca <- prcomp(t(RMA_values), scale=T) # Transpose matrix, scale variables for unit variance
png("images/PCA_normalized.png", width=600, height=600)
s3d <- scatterplot3d(pca$x[,1:3], pch=19, color=rainbow(nrow(adf)))
s3d.coords <- s3d$xyz.convert(pca$x[,1:3])
text(s3d.coords$x, s3d.coords$y, 
	labels=colnames(RMA_values), 
	pos=3, offset=0.5) # plot artificial features
dev.off()

#-------------SAMPLE ANNOTATION

#eset@annotation
#ls("package:mouse4302.db")
ID <- featureNames(eset)
Symbol <- getSYMBOL(ID, "mouse4302.db")
Name <- as.character(lookUp(ID, "mouse4302.db", "GENENAME"))
tmp <- data.frame(ID=ID, Symbol=Symbol, Name=Name, stringsAsFactors=F)
tmp[tmp=="NA"] <- NA 
fData(eset) <- tmp # annotation as feature data of eset
keytypes(mouse4302.db) # annotation keys
res <- select(mouse4302.db, keys=rownames(eset), 
	columns=c("ENTREZID","PFAM", "GENENAME", "ENSEMBL","SYMBOL"), 
	keytype="PROBEID") # entrezID maps to gene signatures
idx <- match(rownames(eset), res$PROBEID) # table join with exact mapping
fData(eset) <- res[idx,] # index as feature data of eset
eset_t <- eset[is.na(fData(eset)$ENTREZID)==0,] # Mask false ENTREZID rows

#-------------FOLD CHANGE

RMA_values10 <- 2**RMA_values # get actual expr values
treated_mean <- apply(RMA_values10[, c("CPC.1","CPC.2","CPC.3")], 1, mean)
control_mean <- apply(RMA_values10[, c("Ctl.1","Ctl.2","Ctl.3")], 1, mean)
fc_mean <- treated_mean/control_mean # simple fold change
foldchange <- cbind(Symbol, Name, treated_mean, control_mean, fc_mean)
write.table(foldchange, file="foldchange.txt", 
	quote=F, sep="\t", col.names=NA)

#-------------LIMMA ANALYSIS

design <- model.matrix(~-1+factor(c(1,1,1,2,2,2))) # (response ~ model)
colnames(design) <- unique(pData(adf)$Group)
contrastmatrix <- makeContrasts(CPC-CTL, levels=design) # which group comparison(s)

fit <- lmFit(eset_t, design) # lmFit(response, model)
fit2 <- eBayes(contrasts.fit(fit, contrastmatrix)) # moderate t-stats w/borrowed variance

# Top DEGs as specified in contrastmatrix
# 'coef'=contrastmatrix column, 'fdr'=statistical adjustment, 'number'=max genes
limmaresults <- topTable(fit2, coef=1, adjust="fdr", number=nrow(eset_t))
write.table(limmaresults, "limmaresults.csv", 
	sep=",", col.names=NA, row.names=TRUE)

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
	ylim=c(0,-log10(0.00001))) # DEGs Volcano plot with cutoffs
dev.off()

#classify <- classifyTestsF(fit2)
classify <- decideTests(fit2)
ded <- unname(colSums(classify!=0)) # count DEGs
png("images/vennDiagram_classified.png", width=600, height=600)
vennDiagram(classify)
dev.off()

#-----------FUNCTIONAL ENRICHMENT ANALYSIS

Mm.H <- readRDS("Mm.h.all.v7.1.entrez.rds") # load signature
RMA_values_filt <- RMA_values[rownames(limmaresults)[1:ded],] # expr values of DEGs
H.indices <- ids2indices(Mm.H, limmaresults[1:ded,]$ENTREZID) # index the signature
camera_results <- camera(RMA_values_filt, index=H.indices,
	design=design, contrast=contrastmatrix[,1],
	weights=-log10(limmaresults[1:ded,"adj.P.Val"])) # competitive test
write.table(camera_results, "func.enrichment_camera.txt", sep="\t")

cm_annotate <- data.frame("Direction"=camera_results$Direction) # add annotation
rownames(cm_annotate) <- rownames(camera_results)
png("images/heatmap_camera.png", height=1000, width=800)
pheatmap(as.matrix(camera_results[,c("PValue","FDR")]),
	main=paste0(GSE,": CAMERA Enriched Genes in CPC group over control"),
	cluster_cols=F, annotation_row=cm_annotate, 
	cutree_rows=4, fontsize=10) # heatmap of enriched terms
dev.off()

#---Barplot modified from https://github.com/YuLab-SMU/DOSE/issues/20
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

# mroast_results <- mroast(eset_t, index=H.indices, design=design,
# contrast=contrastmatrix[,1], adjust.method="BH")
# write.table(mroast_results, "func.enrichment_mroast.txt", sep="\t")
# romer_results <- romer(eset_t, index=H.indices, design=design,
# contrast=contrastmatrix[,1])
# write.table(romer_results, "func.enrichment_romer.txt", sep="\t")

# Optional: Extract the model from the fit if using the same one 
# in both the DEG analysis and the enrichment analysis
# sv <- squeezeVar(fit$sigma^2, df=fit$df.residual)

#---------------SHINY DATA

library(dplyr)
PROBEID <- data.frame(limmaresults$PROBEID)
colnames(PROBEID) <- "PROBEID"
expression <- left_join(
	limmaresults[1:50,], 
	cbind(PROBEID, RMA_values10[rownames(limmaresults),]), 
	by="PROBEID", keep=TRUE, copy=TRUE) # expr values of top 50 DEGs
rownames(expression) <- with(expression, paste(SYMBOL, PROBEID.x, sep="/"))
expression <- expression[, (ncol(expression)-6+1):ncol(expression)] # only sample names
expression[] <- lapply(expression, as.numeric) # convert to numeric
save(expression, file="expression.Rdata")

experiment <- read.table("targets.txt", 
	header=T, as.is=T, row.names=1)
save(experiment, file="experiment.Rdata")

differential <- limmaresults[, c("ENSEMBL","SYMBOL","logFC","P.Value","adj.P.Val","B")]
differential$minus_log10_Pval <- -log10(limmaresults$adj.P.Val)
#differential$sig <- as.factor(
#	abs(differential$logFC) > 2 & differential$adj.P.Val < 0.01)
write.table(differential, "differential.csv", 
	sep=",", col.names=NA, row.names=TRUE)
save(differential, file="differential.Rdata")

#--------------SAVE

#save.image("workspace.RData")
#writeLines(capture.output(sessionInfo()), "sessionInfo.txt")