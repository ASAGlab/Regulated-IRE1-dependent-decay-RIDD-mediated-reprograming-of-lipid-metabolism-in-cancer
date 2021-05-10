# RNAseq - edgeR v.final
library("edgeR")
library("limma")
library("sva")
library("biomaRt")
library("RColorBrewer")
library("ggplot2")
library("stringr")
library("Glimma")

## load raw counts from edgeR
eset <- read.delim("C:/Users/alman/Desktop/NUIG/RNAseq/mRNA/counts.txt", comment.char="#")
counts_short <- eset[,7:30] # first 6 columns is annotation
row.names(counts_short) <- eset$Geneid

## metadata
samples <- data.frame( groups = c("DMSO 8", "DMSO 24", "DMSO 24", "DMSO 8", "DMSO 8", "DMSO 24", "MKC8866 24", "MKC8866 8", "MKC8866 8", "MKC8866 8", "MKC8866 24", "MKC8866 24", "PACLITAXEL 8", "PACLITAXEL 24", "PACLITAXEL 24", "PACLITAXEL 24", "PACLITAXEL 8", "PACLITAXEL 8", "COTREATMENT 8", "COTREATMENT 8", "COTREATMENT 24", "COTREATMENT 24", "COTREATMENT 8", "COTREATMENT 24"),
                       repeats = c("R1", "R3", "R2", "R3", "R2", "R1", "R3", "R2", "R1", "R3", "R2", "R1", "R1", "R2", "R3", "R1", "R2", "R3", "R3", "R1", "R2", "R1", "R2", "R3"),
                       times = c("8", "24", "24", "8", "8", "24", "24", "8", "8", "8", "24", "24", "8", "24", "24", "24", "8", "8", "8", "8", "24", "24", "8", "24"),
                       treatment = c("DMSO", "DMSO", "DMSO", "DMSO", "DMSO", "DMSO", "MKC8866", "MKC8866", "MKC8866", "MKC8866", "MKC8866", "MKC8866", "PACLITAXEL", "PACLITAXEL", "PACLITAXEL", "PACLITAXEL", "PACLITAXEL", "PACLITAXEL", "COTREATMENT", "COTREATMENT", "COTREATMENT", "COTREATMENT", "COTREATMENT", "COTREATMENT"))

row.names(samples) <- substring(colnames(counts_short), 51, 56)
colnames(counts_short) <- rownames(samples)

# remove paclitaxel samples
ind <- samples$treatment == "DMSO" | samples$treatment == "MKC8866"
counts_short <- counts_short[,ind]
samples <- samples[ind,]

# Annotate genes 
listEnsembl(GRCh=38)
listEnsembl(version=91)
ensembl = useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
todo_genes <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol','chromosome_name'), mart=ensembl)

genetable <- data.frame(gene.id=rownames(counts_short))
genes <- merge(genetable, todo_genes, by.x= "gene.id", by.y="ensembl_gene_id", all.x=TRUE)

dup <- genes$gene.id[duplicated(genes$gene.id)]
mat <- match(genetable$gene.id, genes$gene.id)
genes <- genes[mat,]

# create DGEobject
y <- DGEList(counts=counts_short, samples=samples, genes=genes)
names(y)

# REMOVE LOW EXPRESSED GENES
keep.exprs <- filterByExpr(y)
y <- y[keep.exprs, keep.lib.sizes=FALSE]

# COUNT PER MILLION, LOG2 COUNTS. for plotting
cpm <- cpm(y) 
lcpm_pre <- cpm(y, log=TRUE)

## Normalisation by the method of trimmed mean of M-values (TMM)
nsamples <- ncol(counts_short)
col <- brewer.pal(nsamples, "Paired")

wee <- log2(y$counts)
boxplot(wee, las=2, col=col, main="")
title(main="Log2 Raw data",ylab="Log-cpm")

boxplot(lcpm_pre, las=2, col=col, main="")
title(main="Log2 CPM data",ylab="Log-cpm")

y <- calcNormFactors(y)
lcpm <- cpm(y, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="Normalised data",ylab="Log-cpm")

############################### BATCH CORRECTION. SVA estimation of surrogate variables
mod1 <- model.matrix(~0 + groups, y$samples)
mod0 <- model.matrix(~1, y$samples)

svobj <- svaseq(cpm(y), mod1, mod0) 
design <- cbind(mod1, svobj$sv)

# "Clean" gene expression data
cleanY = function(y, mod, svs) {
  X = cbind(mod, svs)
  Hat = solve(t(X) %*% X) %*% t(X)
  beta = (Hat %*% t(y))
  rm(Hat)
  gc()
  P = ncol(mod)
  return(y - t(as.matrix(X[,-c(1:P)]) %*% beta[-c(1:P),]))
}

cleaned_count <- cleanY(cpm(y), mod1, svobj$sv) #you can also specify to not use all sva, just 1,2, etc.
log_cleaned_count <- log2(cleaned_count)

### PCA to inspect batch correction
# no sva
pca <- prcomp(t(na.omit(lcpm_pre)), center = T, scale. = T)
summary(pca)
eigs <- pca$sdev^2
PCAi <- as.data.frame(pca$x)
ggplot(PCAi , aes(PC1, PC2, col= samples$groups)) + geom_point(aes(size=3)) + geom_text(aes(label = row.names(PCAi)), vjust = -1, nudge_y = 1) + 
  labs (y = paste("PC2", round(eigs[2] / sum(eigs), digits = 2)), x = paste("PC1", round(eigs[1] / sum(eigs), digits = 2))) + # substitute with % variances
  ggtitle("PCA. by treatment") +
  theme_bw() 

# after sva
pca <- prcomp(t(na.omit(log_cleaned_count)), center = T, scale. = T)
summary(pca)
eigs <- pca$sdev^2
variance <- eigs*100/sum(eigs)
cumvar <- cumsum(variance)
eig.veamos <- data.frame(eig = eigs, variance = variance, cumvariance = cumvar)

variances <- barplot(eig.veamos[, 2], names.arg=1:nrow(eig.veamos),main = "Variances", xlab = "Principal Components",ylab = "Percentage of variances", col ="steelblue")
lines(x = variances, eig.veamos[, 2], type="b", pch=19, col = "red")

PCAi <- as.data.frame(pca$x)
ggplot(PCAi , aes(PC1, PC2, col= samples$groups)) + geom_point(aes(size=3)) + geom_text(aes(label = row.names(PCAi)), vjust = -1, nudge_y = 1) + 
  labs (y = paste("PC2", round(eigs[2] / sum(eigs), digits = 2)), x = paste("PC1", round(eigs[1] / sum(eigs), digits = 2))) + # substitute with % variances
  ggtitle("PCA. by treatment") +
  theme_bw() 

############################ DIFFERENTIAL EXPRESSION ANALYSIS EdgeR
#each gene will get its own unique dispersion estimate, but the common dispersion is still used in the calculation.
design <- model.matrix(~0 + groups, y$samples)
design <- cbind(design, svobj$sv)

y <- estimateDisp(y, design)

# mean-variance plot. raw variances of the counts (grey dots), the variances using the tagwise

meanVarPlot <- plotMeanVar( y , show.raw.vars=TRUE ,
                            show.tagwise.vars=TRUE ,
                            show.binned.common.disp.vars=FALSE ,
                            show.ave.raw.vars=FALSE ,
                            dispersion.method = "qcml" , NBline = TRUE ,
                            nbins = 100 ,
                            pch = 16 ,
                            xlab ="Mean Expression (Log10 Scale)" ,
                            ylab = "Variance (Log10 Scale)" ,
                            main = "Mean-Variance Plot" )

#
fit <- glmFit(y, design)

## 24 hours
D24vsM24 <- glmLRT(fit, contrast = c(-1,0,1,0,0,0)) # -1 is the denominator (check design matrix)
res_D24vsM24 <- topTags(D24vsM24, n=nrow(y))
resultados_D24vsM24 <- as.data.frame(res_D24vsM24)
topTags(D24vsM24) # just the top 10 by default

table(resultados_D24vsM24$PValue < 0.05)

## 8 hours
D8vsM8 <- glmLRT(fit, contrast = c(0,-1,0,1,0,0)) # -1 is the denominator (check design matrix)
res_D8vsM8 <- topTags(D8vsM8, n=nrow(y))
resultados_D8vsM8 <- as.data.frame(res_D8vsM8)
topTags(D8vsM8) # just the top 10 by default

table(resultados_D8vsM8$PValue < 0.05)

write.csv(resultados_D24vsM24, "resultados_D24vsM24_ALL.csv")
write.csv(resultados_D8vsM8, "resultados_D8vsM8_ALL.csv")

########################################
### Volcano plots
#  with ggplot2
# 8 hours
resultados_D8vsM8$color <- "grey"
resultados_D8vsM8$color[resultados_D8vsM8$logFC > 0.25 & resultados_D8vsM8$PValue < 0.05] <- "red"
resultados_D8vsM8$color[resultados_D8vsM8$logFC < -0.25 & resultados_D8vsM8$PValue < 0.05] <- "green"
table(resultados_D8vsM8$color)

p <- ggplot(resultados_D8vsM8, aes(resultados_D8vsM8$logFC, -log10(resultados_D8vsM8$PValue)))
p + geom_point(aes(colour = resultados_D8vsM8$color, size = 1, alpha = 0.5)) +
  labs(x = "Log2 Fold Change MKC8866/DMSO", y = "-log10(p.value)", element_text(face = "bold", angle = 0)) +
  scale_colour_manual(values = c( "green", "grey50",  "red")) +  
  geom_vline(xintercept = c(-0.25, 0.25) ,linetype = 2, size = 1, color = "grey30") +
  ylim(0, 20) +
  geom_hline(yintercept = -log10(0.05) ,linetype = 2, size = 1, color = "grey30") +
  annotate( "text", x = 2, y = 17,label = "N == 783",parse = TRUE, size=5) +
  annotate( "text", x = -2, y = 17,label = "N == 100",parse = TRUE, size=5) +
  ggtitle("8 hours")
  
# 24 hours
resultados_D24vsM24$color <- "grey"
resultados_D24vsM24$color[resultados_D24vsM24$logFC > 0.25 & resultados_D24vsM24$PValue < 0.05] <- "red"
resultados_D24vsM24$color[resultados_D24vsM24$logFC < -0.25 & resultados_D24vsM24$PValue < 0.05] <- "green"
table(resultados_D24vsM24$color)

p <- ggplot(resultados_D24vsM24, aes(resultados_D24vsM24$logFC, -log10(resultados_D24vsM24$PValue)))
p + geom_point(aes(colour = resultados_D24vsM24$color, size = 1, alpha = 0.5)) +
  labs(x = "Log2 Fold Change MKC8866/DMSO", y = "-log10(p.value)", element_text(face = "bold", angle = 0)) +
  scale_colour_manual(values = c( "green", "grey50",  "red")) +  
  geom_vline(xintercept = c(-0.25, 0.25) ,linetype = 2, size = 1, color = "grey30") +
  geom_hline(yintercept = -log10(0.05) ,linetype = 2, size = 1, color = "grey30") +
  annotate( "text", x = 2, y = 12,label = "N == 273",parse = TRUE, size=5) +
  annotate( "text", x = -2, y = 12,label = "N == 122",parse = TRUE, size=5) +
  ggtitle("24 hours") +
  ylim(0, 15)


### glimmaVolcano
### it saves interactive volcano in your working directory
counts_8h <- cleaned_count[,samples$times == 8]
counts_24h <- cleaned_count[,samples$times == 24]
samples_8h <- samples[samples$times == 8,]
samples_24h <- samples[samples$times == 24,]

uuu <- D8vsM8$table

status_8h <-rep(0, 17268)
status_8h[uuu$logFC > 0.25 & uuu$PValue < 0.05] <- 1
status_8h[uuu$logFC < -0.25 & uuu$PValue < 0.05] <- -1
table(status_8h)

glimmaVolcano(D8vsM8,
              counts=counts_8h,
              anno=y$genes,
              status = status_8h,
              groups=samples_8h$groups,
              samples=colnames(counts_8h),
              transform.counts= "none",
              main = "MKC8866/DMSO 8 hours")

### 
status_24h <-rep(0, 17268)
status_24h[ooo$logFC > 0.25 & ooo$PValue < 0.05] <- 1
status_24h[ooo$logFC < -0.25 & ooo$PValue < 0.05] <- -1
table(status_24h)

glimmaVolcano(D24vsM24,
              counts=counts_24h,
              anno=y$genes,
              status = status_24h,
              groups=samples_24h$groups,
              samples=colnames(counts_24h),
              transform.counts= "none",
              main = "MKC8866/DMSO 24 hours")
