# RNAseq - edgeR v.final
library("edgeR")
library("limma")
library("sva")
library("biomaRt")
library("RColorBrewer")
library("ggplot2")
library("stringr")
library("reshape2")
library("pheatmap")

## load raw counts and sample info
raw_counts <- read.delim(url("https://raw.githubusercontent.com/ASAGlab/Regulated-IRE1-dependent-decay-RIDD-mediated-reprograming-of-lipid-metabolism-in-cancer/main/Data/raw_counts.txt"), sep="") 
sample_info_RNAseq <- read.delim(url("https://raw.githubusercontent.com/ASAGlab/Regulated-IRE1-dependent-decay-RIDD-mediated-reprograming-of-lipid-metabolism-in-cancer/main/Data/sample_info_RNAseq.txt"), sep="") 

eset <- raw_counts[,7:18] # first 6 columns is annotation

# Annotate genes 
listEnsembl(GRCh=38)
listEnsembl(version=105)
ensembl = useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
todo_genes <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol','chromosome_name'), mart=ensembl)

genetable <- data.frame(gene.id=rownames(eset))
genes <- merge(genetable, todo_genes, by.x= "gene.id", by.y="ensembl_gene_id", all.x=TRUE)

dup <- genes$gene.id[duplicated(genes$gene.id)]
mat <- match(genetable$gene.id, genes$gene.id)
genes <- genes[mat,]

# create DGEobject
y <- DGEList(counts=eset, samples=sample_info_RNAseq, genes=genes)
names(y)

# REMOVE LOW EXPRESSED GENES
keep.exprs <- filterByExpr(y)
y <- y[keep.exprs, keep.lib.sizes=FALSE]

# COUNT PER MILLION, LOG2 COUNTS. for plotting
cpm <- cpm(y) 
lcpm_pre <- cpm(y, log=TRUE)

## Normalisation by the method of trimmed mean of M-values (TMM)
nsamples <- ncol(eset)
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

cleaned_count <- as.data.frame(cleanY(cpm(y), mod1, svobj$sv)) #you can also specify to not use all sva, just 1,2, etc.
log_cleaned_count <- log2(cleaned_count)

#
boxplot(log_cleaned_count, las=2, col=col, main="")
title(main="SVA normalised data",ylab="Log-cpm")

### PCA to inspect batch correction
# no sva
pca <- prcomp(t(na.omit(lcpm_pre)), center = T, scale. = T)
summary(pca)
eigs <- pca$sdev^2
PCAi <- as.data.frame(pca$x)
ggplot(PCAi , aes(PC1, PC2, col= sample_info_RNAseq$groups)) + geom_point(aes(size=3)) + geom_text(aes(label = row.names(PCAi)), vjust = -1, nudge_y = 1) + 
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
ggplot(PCAi , aes(PC1, PC2, col= sample_info_RNAseq$groups)) + geom_point(aes(size=3)) + geom_text(aes(label = row.names(PCAi)), vjust = -1, nudge_y = 1) + 
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

## 8 hours
D8vsM8 <- glmLRT(fit, contrast = c(0,-1,0,1,0,0)) # -1 is the denominator (check design matrix)
res_D8vsM8 <- topTags(D8vsM8, n=nrow(y))
resultados_D8vsM8 <- as.data.frame(res_D8vsM8)
topTags(D8vsM8) # just the top 10 by default

write.table(resultados_D24vsM24, "results_D24vsM24_ALL.txt", sep = " ", dec = ".", row.names = TRUE, col.names = TRUE)
write.table(resultados_D8vsM8, "results_D8vsM8_ALL.txt", sep = " ", dec = ".", row.names = TRUE, col.names = TRUE)

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
  ggtitle("8 hours") +
  theme(legend.position = "none")
  
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
  ylim(0, 15) +
  theme(legend.position = "none")

### Functional analysis of hits at 8&24h done with Bioinfominer
### cutoff < 0.05 p.value & abs(log2 FC) > 0.25

hits_8h <- resultados_D8vsM8$hgnc_symbol[resultados_D8vsM8$color != "grey"]
hits_24 <- resultados_D24vsM24$hgnc_symbol[resultados_D24vsM24$color != "grey"]
hits_combined <- as.data.frame(unique(c(hits_8h, hits_24)))

### Identify all DEG annotated to lipid metabolism terms. At either time point
### Download parent term GO_0044255 "cellular lipid metabolic process" and child terms
### read GO signatures downloaded from MsigDB
file_list <- list.files(path= "~/GitHub/Regulated-IRE1-dependent-decay-RIDD-mediated-reprograming-of-lipid-metabolism-in-cancer/Data/",pattern='GO_', full.names = TRUE)
nombres <- list.files(path= "~/GitHub/Regulated-IRE1-dependent-decay-RIDD-mediated-reprograming-of-lipid-metabolism-in-cancer/Data/",pattern='GO_')

signatures <- list()
for (i in 1:length(file_list)){
  temp_data <- read.csv(file_list[i]) 
  signatures[i] <- temp_data #for each iteration, bind the new data to the building list
}

# create new columns and fill with NA
for(i in 1:9) {                                   
  new <- rep(NA, nrow(hits_combined))                       
  hits_combined[ , ncol(hits_combined) + 1] <- new                  
  colnames(hits_combined)[ncol(hits_combined)] <- paste0(nombres[i])
}

# fill columns with LOGICAL. Is gene present in each GO lipid list?
for (i in 1:9){
  hits_combined[,i +1] <- hits_combined$`unique(c(hits_8h, hits_24))` %in% signatures[[i]]
}

### for plotting, change TRUE/FALSE to 1/0
hits_combined <- data.frame(lapply(hits_combined, function(x) {gsub("TRUE", 1, x)}))
hits_combined <- data.frame(lapply(hits_combined, function(x) {gsub("FALSE", 0, x)}))
hits_combined[,2:10] <- lapply(hits_combined[,2:10], function(x) { as.numeric(x)})

hits_combined <- hits_combined[rowSums(hits_combined[,2:10])>1,] # remove non lipid genes
colnames(hits_combined) <- c("hgnc_symbol","fatty_acid_metabolic_process", "neutral_lipid_metabolic_process", "membrane_lipid_metabolic_process", "phospholipid_metabolic_process", "isoprenoid_metabolic_process", "lipid_modification", "cellular_lipid_catabolic_process","cellular_lipid_metabolic_process", "glycerolipid_metabolic_process" )
hits_combined <- hits_combined[,c(1,9,2,3,4,5,6,7,8,10)] # move parent GO term to column 2

# absent/present heatmap
my_palette <- colorRampPalette(c("grey90","grey0"))(n = 3) # or whatever colors

melt2 <- melt(hits_combined, id.vars = 1, measure.vars = 2:10)
melt2$value <- as.numeric(melt2$value)
melt2 <-melt2[order(melt2$hgnc_symbol, decreasing = TRUE),]
melt2$hgnc_symbol <- factor(melt2$hgnc_symbol, levels = c(unique(melt2$hgnc_symbol)))

margin(5,5,5,5)
p <- ggplot(melt2, aes(melt2$variable, melt2$hgnc_symbol, fill = factor(melt2$value)))
p + scale_fill_manual(values=my_palette) + geom_raster() + geom_tile(aes(fill = factor(melt2$value)), colour = "grey50") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), axis.text.y = element_text(hjust = 1, size = 10))

#### HEATMAP TG METABOLISM GENES
# selected from literature
cleaned_count_TG <- cleaned_count[ c("ENSG00000062282", "ENSG00000169692", "ENSG00000131408", "ENSG00000072310", "ENSG00000025434", "ENSG00000166035", "ENSG00000079435" ),]
row.names(cleaned_count_TG) <- c("DGAT2", "AGPAT2", "NR1H2", "SREBF1", "NR1H3", "LIPC", "LIPE")

# average eset for plotting heatmap
TG_eset_avg <- data.frame("DMSO 8" = rowMeans(cleaned_count_TG[,sample_info_RNAseq$groups == "DMSO 8"]), 
                          "MKC8866 8" = rowMeans(cleaned_count_TG[,sample_info_RNAseq$groups == "MKC8866 8"]), 
                          "DMSO 24" = rowMeans(cleaned_count_TG[,sample_info_RNAseq$groups == "DMSO 24"]), 
                          "MKC8866 24" =rowMeans(cleaned_count_TG[,sample_info_RNAseq$groups == "MKC8866 24"]))

my_palette <- colorRampPalette(c("navy","white","red"))(n = 299) # scale blue to red, or whatever colors
pheatmap((TG_eset_avg), 
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         scale = "row", 
         color = my_palette,
         fontsize_col = 11, 
         border_color = NA,
         show_rownames = TRUE,
         show_colnames = TRUE,
         cellwidth = 18,
         cellheight = 18,
         angle_col = 45,
         gaps_col = 2,
         gaps_row = 5)


### CONTROL GENES EXPRESSION
# DNAJB9 - ENSG00000128590 
plotting_genes <-  as.data.frame(t(y$counts))

# using ensembl identifier
sample_info_RNAseq$groups <- factor(sample_info_RNAseq$groups, levels = c("DMSO 8" , "MKC8866 8", "DMSO 24", "MKC8866 24"))

p <- ggplot(plotting_genes, aes(sample_info_RNAseq$groups, plotting_genes$ENSG00000128590, fill = sample_info_RNAseq$groups)) # change for any other gene here
p + geom_boxplot(outlier.shape = NA, size =2, alpha = 0.5) + 
  geom_point(aes(), size = 2, position = position_jitterdodge()) + 
  labs(x="", y = "DNAJB9 norm counts", element_text(face = "bold", angle = 0)) + 
  scale_fill_manual(values =  c("white", "red", "white", "red")) + # change colors here 
  theme(panel.background = element_rect( colour = "grey50"), axis.text.x = element_text(angle = 45, hjust = 1, size = 12)) + 
  theme_bw() +
  theme(text = element_text(size=20))



