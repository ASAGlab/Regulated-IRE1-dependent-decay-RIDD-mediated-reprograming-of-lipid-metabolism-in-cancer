# TAG lenght&saturation analysis
library("stringr")
library("ggplot2")
library("plyr")
library("reshape2")
library("svglite")
library("pheatmap")
library("ggpubr")

# indicate path to .csv file with average values for each replicate
todo_treatments <- read.delim(url("https://raw.githubusercontent.com/ASAGlab/Regulated-IRE1-dependent-decay-RIDD-mediated-reprograming-of-lipid-metabolism-in-cancer/main/Data/Beatson_MS.txt"), sep = " ")
sample_info <- read.delim(url("https://raw.githubusercontent.com/ASAGlab/Regulated-IRE1-dependent-decay-RIDD-mediated-reprograming-of-lipid-metabolism-in-cancer/main/Data/Beatson_sample_info.txt"), sep = " ")

# MSTUS NORMALIZATION (normalized to sum of peak areas)
MSTUS_norm <- sweep(todo_treatments[,3:ncol(todo_treatments)] ,2,colSums(todo_treatments[,3:ncol(todo_treatments)])/100000000,`/`)

# LOG2 scale (plotting, Limma)
log_MSTUS_norm <- log2(MSTUS_norm)

### PCA
pca <- prcomp(t(na.omit(log_MSTUS_norm)), center = T, scale. = T)
summary(pca)
eigs <- pca$sdev^2
PCAi <- as.data.frame(pca$x)
ggplot(PCAi , aes(PC1, PC2, col= sample_info$group)) + geom_point(aes(size=3)) + geom_text(aes(label = row.names(PCAi)), nudge_y = 1) + 
  labs (y = paste("PC2", round(eigs[2] / sum(eigs), digits = 2)), x = paste("PC1", round(eigs[1] / sum(eigs), digits = 2))) + # substitute with % variances
  ggtitle("PCA. by treatment") +
  theme_bw() 

#PCA by repeat
ggplot(PCAi , aes(PC1, PC2, col= as.character(sample_info$replicate))) + geom_point(aes(size=3)) + geom_text(aes(label = row.names(PCAi)), nudge_y = 1) + 
  labs (y = paste("PC2", round(eigs[2] / sum(eigs), digits = 2)), x = paste("PC1", round(eigs[1] / sum(eigs), digits = 2))) + # substitute with % variances
  ggtitle("PCA. by treatment") +
  theme_bw() 


### averaged eset for plotting
final_avg <- as.data.frame(cbind(rowMeans(log_MSTUS_norm[,1:3]), rowMeans(log_MSTUS_norm[,4:6]), rowMeans(log_MSTUS_norm[,7:9]), rowMeans(log_MSTUS_norm[,10:12])))
row.names(final_avg) <- row.names(log_MSTUS_norm)
colnames(final_avg) <- c("cotreatment", "DGAT2inh", "DMSO", "MKC8866")
final_avg <- cbind(todo_treatments[,1:2], final_avg)

# calculate FoldChanges, standard error of mean, error of the ratios
final_avg$DMSO.sd <- apply(log_MSTUS_norm[,1:3], 1, function(x) sd(as.numeric(x)))
final_avg$MKC.sd <- apply(log_MSTUS_norm[,4:6], 1, function(x) sd(as.numeric(x)))
final_avg$DGAT2inh.sd <- apply(log_MSTUS_norm[,7:9], 1, function(x) sd(as.numeric(x)))
final_avg$cotreat.sd <- apply(log_MSTUS_norm[,10:12], 1, function(x) sd(as.numeric(x)))

final_avg$DvsD2i <- final_avg$DGAT2inh - final_avg$DMSO
final_avg$Dvscotreat <- final_avg$cotreatment - final_avg$DMSO
final_avg$DvsM <- final_avg$MKC8866 - final_avg$DMSO

final_avg$DvsD2i.se <- abs(sqrt((final_avg$DGAT2inh.sd/final_avg$DGAT2inh)^2 + (final_avg$DMSO.sd/final_avg$DMSO)^2)*final_avg$DvsD2i)
final_avg$Dvscotreat.se <- abs(sqrt((final_avg$cotreat.sd/final_avg$cotreatment)^2 + (final_avg$DMSO.sd/final_avg$DMSO)^2)*final_avg$Dvscotreat)
final_avg$DvsM.se <- abs(sqrt((final_avg$MKC.sd/final_avg$MKC8866)^2 + (final_avg$DMSO.sd/final_avg$DMSO)^2)*final_avg$DvsM)

## formatting, lenght & saturation
final_avg$lenght <- sapply(str_extract_all(final_avg$C..Identification, "[0-9]{2}(?:)"), function(x) sum(as.numeric(x)))
final_avg$saturation <- sapply(str_extract_all(final_avg$C..Identification, "(?<=:)[0-9]+"), function(x) sum(as.numeric(x)))
final_avg$nombre <- paste(final_avg$C..Lipid.Class, final_avg$lenght, ":", final_avg$saturation)
annotation <- unique(final_avg[,c(1,19)])

# consolidate duplicated entries for same metabolites
final_avg <- ddply(final_avg,19,numcolwise(mean))
final_avg <- merge(annotation, final_avg, by.x = 2, by.y = 1)
row.names(final_avg) <- final_avg$nombre

### TAKE CLOSER LOOK AT TAGS
### HEATMAP
TAGs <- subset(final_avg, final_avg$C..Lipid.Class == "TG")
TAGs <- TAGs[order(TAGs$lenght),]
TAGs <- TAGs[order(TAGs$saturation),]

my_palette <- colorRampPalette(c("navy","white","red"))(n = 20) # or whatever colors and number
pheatmap(na.omit(TAGs[,c(5,6,4,3)]), 
         cluster_cols=FALSE, 
         cluster_rows = FALSE,
         color = my_palette, 
         border_color = NA, 
         fontsize_col = 14, 
         scale = "row",
         show_rownames = TRUE)

### Plot changes for all TAG species
numeros_magicos <- cumsum(rev(table(TAGs$saturation))) + 0.5
melt_uno <- melt(TAGs, id.vars = 1, measure.vars = 11:13)
melt_dos <- melt(TAGs, id.vars = 1, measure.vars = 14:16)
melt_todo <- cbind(melt_uno, melt_dos[,3])
melt_todo$nombre <- factor(melt_todo$nombre, levels = rev(unique(melt_todo$nombre)))
melt_todo$variable <- gsub("DvsM", "MKC8866", melt_todo$variable)
melt_todo$variable <- gsub("DvsD2i", "PF-06427878", melt_todo$variable)
melt_todo$variable <- gsub("Dvscotreat", "Cotreatment", melt_todo$variable)
melt_todo$variable <- factor(melt_todo$variable, levels = c("Cotreatment","PF-06427878","MKC8866"))

ggplot() +
  geom_bar(data=melt_todo, aes(fill=melt_todo$variable, y=melt_todo$value,x=melt_todo$nombre, ymin=melt_todo$value-melt_todo$`melt_dos[, 3]`,ymax=melt_todo$value+melt_todo$`melt_dos[, 3]`), position="dodge", stat="identity", color="black") +
  geom_errorbar(data=melt_todo, aes(fill=melt_todo$variable, y=melt_todo$value,x=melt_todo$nombre,ymin=melt_todo$value-melt_todo$`melt_dos[, 3]`,ymax=melt_todo$value+melt_todo$`melt_dos[, 3]`), position="dodge",stat="identity") +
  labs(x = "", y = "log2 FC (MKC8866 vs control)", element_text(face = "bold", angle = 0)) +
  coord_flip() +
  geom_vline(xintercept = numeros_magicos[-length(numeros_magicos)] ,linetype = 2, size = 1, color = "black") 
  

### saturation & elongation indexes TGs




