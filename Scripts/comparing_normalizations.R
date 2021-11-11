# comparing metabolomics normalization methods
# MSTUS vs QUANTILE vs PARETO vs LOESS
library("limma")
library("RColorBrewer")
library("ggplot2")
library("sva")
library("GGally")

todo_treatments <- read.csv(url("https://raw.githubusercontent.com/ASAGlab/Regulated-IRE1-dependent-decay-RIDD-mediated-reprograming-of-lipid-metabolism-in-cancer/main/Data/Metasyx_MD.csv"))
Sample_Info <- read.csv(url("https://raw.githubusercontent.com/ASAGlab/Regulated-IRE1-dependent-decay-RIDD-mediated-reprograming-of-lipid-metabolism-in-cancer/main/Data/Metasyx_sample_info.csv"))

eset <- todo_treatments[,12:ncol(todo_treatments)]

# MSTUS NORMALIZATION (normalized to sum of peak areas)
MSTUS_norm <- sweep(eset ,2,colSums(eset)/100000000,`/`)
row.names(MSTUS_norm) <- todo_treatments$Peak_ID

# PARETO SCALING (normalized to  square-root of column-wise standard deviations)
PARETO_norm <- as.data.frame(lapply(eset, function(x) x/sqrt(sd(x))))
row.names(PARETO_norm) <- todo_treatments$Peak_ID

# QUANTILE NORMALIZATION (normalization based upon quantiles)
quantile_normalisation <- function(df){
  df_rank <- apply(df,2,rank,ties.method="min")
  df_sorted <- data.frame(apply(df, 2, sort))
  df_mean <- apply(df_sorted, 1, mean)
  index_to_mean <- function(my_index, my_mean){
    return(my_mean[my_index])
  }
  df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
  rownames(df_final) <- rownames(df)
  return(df_final)
}

QUANTILE_norm <- as.data.frame(quantile_normalisation(eset))
row.names(QUANTILE_norm) <- todo_treatments$Peak_ID

# SVA BATCH EFFECT CORRECTION
cleanY = function(y, mod, svs) {
  X = cbind(mod, svs)
  Hat = solve(t(X) %*% X) %*% t(X)
  beta = (Hat %*% t(y))
  rm(Hat)
  gc()
  P = ncol(mod)
  return(y - t(as.matrix(X[,-c(1:P)]) %*% beta[-c(1:P),]))
}

sva_corr = function(x){
mod1 <- model.matrix(~0 + Group, Sample_Info)
mod0 <- model.matrix(~1, Sample_Info)
svobj <- svaseq(as.matrix(x), mod1, mod0) 
cleaned_count <- cleanY(x, mod1, svobj$sv) 
}

# sva batch effect adjust each norm data
MSTUS_norm_c <- sva_corr(log2(MSTUS_norm))
PARETO_norm_c <- log2(sva_corr(PARETO_norm))
QUANTILE_norm_c <- sva_corr(log2(QUANTILE_norm))

# HISTOGRAMS
nsamples <- length(unique(Sample_Info$Group))
col <- brewer.pal(nsamples, "Paired")
col <- rep(col, each=5)

boxplot(log2(eset), las=2, col=col, main="")
title(main="Raw data",ylab="")

boxplot(MSTUS_norm_c, las=2, col=col, main="")
title(main="MSTUS + sva norm data",ylab="")

boxplot(PARETO_norm_c, las=2, col=col, main="")
title(main="PARETO + sva scaled data",ylab="")

boxplot(QUANTILE_norm_c, las=2, col=col, main="")
title(main="QUANTILE + sva normalised data",ylab="")

# PCA for each
pca_explore = function(x){ 
pca <- prcomp(t(na.omit(x)), center = T, scale. = F)
summary(pca)
eigs <- pca$sdev^2
PCAi <- as.data.frame(pca$x)
ggplot(PCAi , aes(PC1, PC2, col= Sample_Info$Group)) + geom_point(aes(size=3)) + geom_text(aes(label = row.names(PCAi)), nudge_y = 1) + 
  labs (y = paste("PC2", round(eigs[2] / sum(eigs), digits = 2)), x = paste("PC1", round(eigs[1] / sum(eigs), digits = 2))) + # substitute with % variances
  ggtitle(paste(deparse(substitute(x)))) +
  theme_bw() 
}

pca_explore(log2(eset))
pca_explore(MSTUS_norm_c)
pca_explore(PARETO_norm_c)
pca_explore(QUANTILE_norm_c)

# Correlation of logFC for metabolites
FC_todo = function(x){ 
final_avg <- as.data.frame(cbind(rowMeans(x[,1:5]), rowMeans(x[,6:10]), rowMeans(x[,11:15]), rowMeans(x[,16:20]), rowMeans(x[,21:25]), rowMeans(x[,26:30])))
DvsM24h <- final_avg[,2] - final_avg[,1]
DvsM48h <- final_avg[,4] - final_avg[,3]
DvsM72h <- final_avg[,6] - final_avg[,5]
todo <- c(DvsM24h, DvsM48h, DvsM72h)
}

todo_junto <- data.frame( MSTUS = FC_todo(MSTUS_norm_c),
                          PARETO = FC_todo(PARETO_norm_c),
                          QUANTILE = FC_todo(QUANTILE_norm_c),
                          time_point = c(rep("24h", 798), rep("48h", 798), rep("72h", 798)),
                          ID = rep(row.names(MSTUS_norm_c), 3))

ggpairs(todo_junto,
        columns = 1:3,
        aes(color=todo_junto$time_point, alpha = 0.5))


# TAKE A LOOK AT TAGs
annotation <- todo_treatments[,c(1,10)]
ind <- subset(annotation, annotation$Lipid_Class == "Triacylglyceride")

TAGs <- merge(ind, todo_junto, by.x =1, by.y = 5)


ggplot() +  geom_point(data=TAGs ,pch=21, size=3, colour="black", aes(fill= TAGs$Lipid_Class, group= TAGs$time_point, y= TAGs$MSTUS, x= TAGs$time_point), position=position_jitterdodge(jitter.width=0.4, dodge.width=0.9)) +
  theme_classic() +
  theme(text = element_text(size=20), axis.title.y = element_text(size=16), legend.position = "none") +
  ggtitle("MSTUS") +
  geom_hline(yintercept = 0 ,linetype = 2, size = 1, color = "black") +
  labs(x="", y = "Log2 FC (MKC8866/DMSO)", element_text(face = "bold", angle = 0,)) +
  ylim(-1, 1)


ggplot() +  geom_point(data=TAGs ,pch=21, size=3, colour="black", aes(fill= TAGs$Lipid_Class, group= TAGs$time_point, y= TAGs$PARETO, x= TAGs$time_point), position=position_jitterdodge(jitter.width=0.4, dodge.width=0.9)) +
  theme_classic() +
  theme(text = element_text(size=20), axis.title.y = element_text(size=16), legend.position = "none") +
  ggtitle("PARETO") +
  geom_hline(yintercept = 0 ,linetype = 2, size = 1, color = "black") +
  labs(x="", y = "Log2 FC (MKC8866/DMSO)", element_text(face = "bold", angle = 0,)) +
  ylim(-1, 1)


ggplot() +  geom_point(data=TAGs ,pch=21, size=3, colour="black", aes(fill= TAGs$Lipid_Class, group= TAGs$time_point, y= TAGs$QUANTILE, x= TAGs$time_point), position=position_jitterdodge(jitter.width=0.4, dodge.width=0.9)) +
  theme_classic() +
  theme(text = element_text(size=20), axis.title.y = element_text(size=16), legend.position = "none") +
  ggtitle("QUANTILE") +
  geom_hline(yintercept = 0 ,linetype = 2, size = 1, color = "black") +
  labs(x="", y = "Log2 FC (MKC8866/DMSO)", element_text(face = "bold", angle = 0,)) +
  ylim(-1, 1)
