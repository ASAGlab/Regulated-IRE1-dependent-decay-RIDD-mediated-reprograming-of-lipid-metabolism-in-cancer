### METASYX FINAL
library("ggplot2")
library("sva")
library("colorRamps")
library("pheatmap")
library("stringr")
library("reshape2")
library("plyr")

# indicate path to .csv file with average values for each replicate
# Load data and metadata 
# data: samples in columns, features in rows
todo_treatments <- read.csv(url("https://raw.githubusercontent.com/ASAGlab/Regulated-IRE1-dependent-decay-RIDD-mediated-reprograming-of-lipid-metabolism-in-cancer/main/Data/Metasyx_MD.csv"))
Sample_Info <- read.csv(url("https://raw.githubusercontent.com/ASAGlab/Regulated-IRE1-dependent-decay-RIDD-mediated-reprograming-of-lipid-metabolism-in-cancer/main/Data/Metasyx_sample_info.csv"))

# MSTUS NORMALIZATION (normalized to sum of peak areas)
MSTUS_norm <- sweep(todo_treatments[,12:ncol(todo_treatments)] ,2,colSums(todo_treatments[,12:ncol(todo_treatments)])/100000000,`/`)
row.names(MSTUS_norm) <- todo_treatments$Peak_ID

# LOG2 scale (plotting)
log_MSTUS_norm <- log2(MSTUS_norm)

### PCA
pca <- prcomp(t(na.omit(log_MSTUS_norm)), center = T, scale. = T)
summary(pca)
eigs <- pca$sdev^2
PCAi <- as.data.frame(pca$x)
ggplot(PCAi , aes(PC1, PC2, col= Sample_Info$Group)) + geom_point(aes(size=3)) + geom_text(aes(label = row.names(PCAi)), nudge_y = 1) + 
  labs (y = paste("PC2", round(eigs[2] / sum(eigs), digits = 2)), x = paste("PC1", round(eigs[1] / sum(eigs), digits = 2))) + # substitute with % variances
  ggtitle("PCA. by treatment") +
  theme_bw() 

#PCA by repeat
ggplot(PCAi , aes(PC1, PC2, col= as.character(Sample_Info$replicate))) + geom_point(aes(size=3)) + geom_text(aes(label = row.names(PCAi)), nudge_y = 1) + 
  labs (y = paste("PC2", round(eigs[2] / sum(eigs), digits = 2)), x = paste("PC1", round(eigs[1] / sum(eigs), digits = 2))) + # substitute with % variances
  ggtitle("PCA. by treatment") +
  theme_bw() 

# repeat 1 looks particularly ugly...
# BATCH EFFECT CORRECTION WITH SVA
mod1 <- model.matrix(~0 + Group, Sample_Info)
mod0 <- model.matrix(~1, Sample_Info)
svobj <- svaseq(as.matrix(MSTUS_norm), mod1, mod0) 

# "Clean" dataset
cleanY = function(y, mod, svs) {
  X = cbind(mod, svs)
  Hat = solve(t(X) %*% X) %*% t(X)
  beta = (Hat %*% t(y))
  rm(Hat)
  gc()
  P = ncol(mod)
  return(y - t(as.matrix(X[,-c(1:P)]) %*% beta[-c(1:P),]))
}

cleaned_count <- cleanY(log_MSTUS_norm, mod1, svobj$sv[,1:4]) # how many sv to use (the more, more transformation)

### PCA after sva
pca <- prcomp(t(na.omit(cleaned_count)), center = T, scale. = T)
summary(pca)
eigs <- pca$sdev^2

variance <- eigs*100/sum(eigs)
cumvar <- cumsum(variance)
eig.veamos <- data.frame(eig = eigs, variance = variance, cumvariance = cumvar)
variances <- barplot(eig.veamos[, 2], names.arg=1:nrow(eig.veamos),main = "Variances", xlab = "Principal Components",ylab = "Percentage of variances", col ="steelblue")
lines(x = variances, eig.veamos[, 2], type="b", pch=19, col = "red")

PCAi <- as.data.frame(pca$x)
ggplot(PCAi , aes(PC1, PC2, col= Sample_Info$Group)) + geom_point(aes(size=3)) + geom_text(aes(label = row.names(PCAi)), nudge_y = 1) + 
  labs (y = paste("PC2", round(eigs[2] / sum(eigs), digits = 2)), x = paste("PC1", round(eigs[1] / sum(eigs), digits = 2))) + # substitute with % variances
  ggtitle("PCA. by treatment") +
  theme_bw() 

### averaged eset for plotting
final_avg <- as.data.frame(cbind(rowMeans(cleaned_count[,1:5]), rowMeans(cleaned_count[,6:10]), rowMeans(cleaned_count[,11:15]), rowMeans(cleaned_count[,16:20]), rowMeans(cleaned_count[,21:25]), rowMeans(cleaned_count[,26:30])))
colnames(final_avg) <- c( "DMSO24", "MKC24", "DMSO48", "MKC48", "DMSO72", "MKC72")
final_avg <- cbind(todo_treatments[,c(1,2,10)], final_avg)

# standard deviation (SD)
final_avg$DMSO24.sd <- apply(cleaned_count[,1:5], 1, function(x) sd(as.numeric(x)))
final_avg$MKC24.sd <- apply(cleaned_count[,6:10], 1, function(x) sd(as.numeric(x)))
final_avg$DMSO48.sd <- apply(cleaned_count[,11:15], 1, function(x) sd(as.numeric(x)))
final_avg$MKC48.sd <- apply(cleaned_count[,16:20], 1, function(x) sd(as.numeric(x)))
final_avg$DMSO72.sd <- apply(cleaned_count[,21:25], 1, function(x) sd(as.numeric(x)))
final_avg$MKC72.sd <- apply(cleaned_count[,26:30], 1, function(x) sd(as.numeric(x)))

# MKC/DMSO log2 fold changes
final_avg$DvsM24h <- final_avg$MKC24 - final_avg$DMSO24
final_avg$DvsM48h <- final_avg$MKC48 - final_avg$DMSO48
final_avg$DvsM72h <- final_avg$MKC72 - final_avg$DMSO72

# add lenght and saturation columns
final_avg$lenght <- sapply(str_extract_all(final_avg$Compound_Name, "[0-9]+(?:)"), function(x) sum(as.numeric(x)))
final_avg$saturation <- sapply(str_extract_all(final_avg$Compound_Name, "(?<=:)[0-9]+"), function(x) sum(as.numeric(x)))
final_avg$lenght <- final_avg$lenght - final_avg$saturation #bad regex 

# solo lipids
solo_lipids <- subset(final_avg, final_avg$Lipid_Class != "Amino Acids") #also removes NA

# HEATMAP OF FOLD CHANGES FOR ALL LIPIDS THROUGH 3 TIME POINTS
annotation_row <- as.data.frame(solo_lipids$Lipid_Class)
row.names(annotation_row) <- solo_lipids$Peak_ID
row.names(solo_lipids) <- solo_lipids$Peak_ID

newCols <- rev(primary.colors(length(unique(annotation_row$`solo_lipids$Lipid_Class`)))) #reverse just cause I didnt like colors
names(newCols) <- unique(annotation_row$`solo_lipids$Lipid_Class`)
newCols <- list(`solo_lipids$Lipid_Class` = newCols)

my_palette <- colorRampPalette(c("navy","white","red"))(n = 20) # or whatever colors and number
myBreaks <- c(seq(-1, 0, length.out=ceiling(length(my_palette)/2) + 1), 
              seq(1/length(my_palette), 1, length.out=floor(length(my_palette)/2))) # from -1 to 1 log2FC 

pheatmap(na.omit(solo_lipids[,16:18]), 
         cluster_cols=FALSE, 
         color = my_palette, 
         border_color = NA, 
         fontsize_col = 14, 
         show_rownames = FALSE,
         annotation_row = annotation_row,
         annotation_colors = newCols,
         breaks = myBreaks)


# lipid group through time
timeplot = function(lipido, color) {
  ooo <- subset(solo_lipids, solo_lipids$Lipid_Class == lipido)
  melt_ooo <- melt(ooo, id.vars = 3, measure.vars = 16:18)
  melt_ooo$variable <- gsub("DvsM", "", melt_ooo$variable)
  numero <- round(max(melt_ooo$value + 0.1), digits = 1)
  numero2 <- round(min(melt_ooo$value - 0.1), digits = 1)
  ggplot() +  geom_point(data=melt_ooo,pch=21, size=3, colour="black", aes(fill= `Lipid_Class`, group= `variable`, y= `value`,x= `variable`), position=position_jitterdodge(jitter.width=0.4, dodge.width=0.9)) +
    theme_classic() +
    theme(text = element_text(size=20), axis.title.y = element_text(size=16), legend.position = "none") +
    ggtitle(lipido) +
    scale_fill_manual(values = color) +
    geom_hline(yintercept = 0 ,linetype = 2, size = 1, color = "black") +
    ylim(numero2, numero) +
    labs(x="", y = "Log2 FC (MKC8866/DMSO)", element_text(face = "bold", angle = 0,)) + 
    annotate('text', x = c("24h", "48h", "72h"), y = numero ,label = c(paste("P-value==", t.test(ooo$DvsM24, mu = 0, alternative = "two.sided")$p.value),
                                                                       paste("P-value==", t.test(ooo$DvsM48, mu = 0, alternative = "two.sided")$p.value),
                                                                       paste("P-value==", t.test(ooo$DvsM72, mu = 0, alternative = "two.sided")$p.value)) ,parse = TRUE, size=3)
}


timeplot("Diacylglycerol", "green")

# plot all lipid groups. using same colors as in heatmap
for (i in unique(solo_lipids$Lipid_Class)) {
  g <- timeplot(i, newCols$`solo_lipids$Lipid_Class`[i])
  print(g) 
  #ggsave( paste( "time", i, ".pdf"), plot=last_plot(), width=8, height=5)
}

#or plot them all together to evaluate in 1 look
library(gridExtra)
p <- list()
for (i in unique(solo_lipids$Lipid_Class)) {
  p[[i]] <- timeplot(i, newCols$`solo_lipids$Lipid_Class`[i])
}
do.call(grid.arrange,p)

### saturation & elongation indexes TGs
# plot by saturation
# divided in 0, 1,2,3+ saturated bonds
satuplot = function(lipido, color) {
  ooo <- subset(final_avg, final_avg$Lipid_Class == lipido)
  uuu <- unique(ooo$saturation)
  ooo$saturange <- NA
  ooo$saturange[ooo$saturation == 0] <-  "0"
  ooo$saturange[ooo$saturation == 1] <-  "1"
  ooo$saturange[ooo$saturation == 2] <-  "2"
  ooo$saturange[ooo$saturation >= 3] <-  "3+"
  average_tri <- ddply(ooo,21,numcolwise(mean))
  average_tri$se24 <- c(sd(ooo$DvsM24[ooo$saturange == "0"]), sd(ooo$DvsM24[ooo$saturange == "1"]), sd(ooo$DvsM24[ooo$saturange == "2"]), sd(ooo$DvsM24[ooo$saturange == "3+"]))
  average_tri$se48 <- c(sd(ooo$DvsM24[ooo$saturange == "0"]), sd(ooo$DvsM48[ooo$saturange == "1"]), sd(ooo$DvsM48[ooo$saturange == "2"]), sd(ooo$DvsM48[ooo$saturange == "3+"]))
  average_tri$se72 <- c(sd(ooo$DvsM24[ooo$saturange == "0"]), sd(ooo$DvsM72[ooo$saturange == "1"]), sd(ooo$DvsM72[ooo$saturange == "2"]), sd(ooo$DvsM72[ooo$saturange == "3+"]))
  melt2 <- melt(average_tri, id.vars = 1, measure.vars = 14:16)
  melt3 <- melt(average_tri, id.vars = 1, measure.vars = 19:21)
  melt_todo <- cbind(melt2, melt3[,3])
  melt_todo$variable <- gsub("DvsM24", "24h", melt_todo$variable)
  melt_todo$variable <- gsub("DvsM48", "48h", melt_todo$variable)
  melt_todo$variable <- gsub("DvsM72", "72h", melt_todo$variable)
  melt_dots <- melt(ooo, id.vars = 21, measure.vars = 16:18)
  melt_dots$variable <- gsub("DvsM24", "24h", melt_dots$variable)
  melt_dots$variable <- gsub("DvsM48", "48h", melt_dots$variable)
  melt_dots$variable <- gsub("DvsM72", "72h", melt_dots$variable)
  numero <- round(max(melt_dots$value + 0.1), digits = 1)
  numero2 <- round(min(melt_dots$value - 0.1), digits = 1)
  ggplot() +
    geom_bar(data=melt_todo, aes(fill=melt_todo$variable, y=melt_todo$value,x=melt_todo$saturange,ymin=melt_todo$value-melt_todo$`melt3[, 3]`,ymax=melt_todo$value+melt_todo$`melt3[, 3]`), position="dodge", stat="identity", color="black") +
    geom_errorbar(data=melt_todo, aes(fill=melt_todo$variable, y=melt_todo$value,x=melt_todo$saturange,ymin=melt_todo$value-melt_todo$`melt3[, 3]`,ymax=melt_todo$value+melt_todo$`melt3[, 3]`), position="dodge",stat="identity") +
    geom_point(data=melt_dots,pch=21, colour="black", aes(fill=melt_dots$variable ,group=melt_dots$variable, y=melt_dots$value,x=melt_dots$saturange), position=position_jitterdodge(jitter.width=0.1, dodge.width=0.9)) +
    scale_fill_manual(values= c("white", "grey50", color)) +
    labs(x = "number of double bonds", y = "log2 FC (MKC8866 vs control)", element_text(face = "bold", angle = 0)) +
    ylim(numero2, numero)
}

# elongation
# plot by total chain lenght
# divided into ranges matching the 4 quartiles

lenghtplot = function(lipido, color) {
  ooo <- subset(final_avg, final_avg$Lipid_Class == lipido)
  uuu <- summary(ooo$lenght)
  ooo$saturange <- NA
  ooo$saturange[ooo$lenght <= round(uuu[[2]])] <- paste( "<" , round(uuu[[2]]))
  ooo$saturange[ooo$lenght > round(uuu[[2]])] <- paste(round(uuu[[2]]),"-", round(uuu[[3]]))
  ooo$saturange[ooo$lenght > round(uuu[[3]])] <- paste(round(uuu[[3]]),"-", round(uuu[[5]]))
  ooo$saturange[ooo$lenght >= round(uuu[[5]])] <- paste( ">" , round(uuu[[5]]))
  average_tri <- ddply(ooo,21,numcolwise(mean))
  average_tri$se24 <- c(sd(ooo$DvsM24[ooo$saturange == paste( "<" , round(uuu[[2]]))]), sd(ooo$DvsM24[ooo$saturange == paste(round(uuu[[2]]),"-", round(uuu[[3]]))]), sd(ooo$DvsM24[ooo$saturange == paste(round(uuu[[3]]),"-", round(uuu[[5]]))]), sd(ooo$DvsM24[ooo$saturange == paste( ">" , round(uuu[[5]]))]))
  average_tri$se48 <- c(sd(ooo$DvsM24[ooo$saturange == paste( "<" , round(uuu[[2]]))]), sd(ooo$DvsM48[ooo$saturange == paste(round(uuu[[2]]),"-", round(uuu[[3]]))]), sd(ooo$DvsM48[ooo$saturange == paste(round(uuu[[3]]),"-", round(uuu[[5]]))]), sd(ooo$DvsM48[ooo$saturange == paste( ">" , round(uuu[[5]]))]))
  average_tri$se72 <- c(sd(ooo$DvsM24[ooo$saturange == paste( "<" , round(uuu[[2]]))]), sd(ooo$DvsM72[ooo$saturange == paste(round(uuu[[2]]),"-", round(uuu[[3]]))]), sd(ooo$DvsM72[ooo$saturange == paste(round(uuu[[3]]),"-", round(uuu[[5]]))]), sd(ooo$DvsM72[ooo$saturange == paste( ">" , round(uuu[[5]]))]))
  melt2 <- melt(average_tri, id.vars = 1, measure.vars = 14:16)
  melt3 <- melt(average_tri, id.vars = 1, measure.vars = 19:21)
  melt_todo <- cbind(melt2, melt3[,3])
  melt_todo$variable <- gsub("DvsM24", "24h", melt_todo$variable)
  melt_todo$variable <- gsub("DvsM48", "48h", melt_todo$variable)
  melt_todo$variable <- gsub("DvsM72", "72h", melt_todo$variable)
  melt_dots <- melt(ooo, id.vars = 21, measure.vars = 16:18)
  melt_dots$variable <- gsub("DvsM24", "24h", melt_dots$variable)
  melt_dots$variable <- gsub("DvsM48", "48h", melt_dots$variable)
  melt_dots$variable <- gsub("DvsM72", "72h", melt_dots$variable)
  numero <- round(max(melt_dots$value + 0.1), digits = 1)
  numero2 <- round(min(melt_dots$value - 0.1), digits = 1)
  melt_todo$saturange <- factor(melt_todo$saturange, levels = c(paste( "<" , round(uuu[[2]])), paste(round(uuu[[2]]),"-", round(uuu[[3]])), paste(round(uuu[[3]]),"-", round(uuu[[5]])), paste( ">" , round(uuu[[5]]))))
  ggplot() +
    geom_bar(data=melt_todo, aes(fill=melt_todo$variable, y=melt_todo$value,x=melt_todo$saturange,ymin=melt_todo$value-melt_todo$`melt3[, 3]`,ymax=melt_todo$value+melt_todo$`melt3[, 3]`), position="dodge", stat="identity", color="black") +
    geom_errorbar(data=melt_todo, aes(fill=melt_todo$variable, y=melt_todo$value,x=melt_todo$saturange,ymin=melt_todo$value-melt_todo$`melt3[, 3]`,ymax=melt_todo$value+melt_todo$`melt3[, 3]`), position="dodge",stat="identity") +
    geom_point(data=melt_dots,pch=21, colour="black", aes(fill=melt_dots$variable ,group=melt_dots$variable, y=melt_dots$value,x=melt_dots$saturange), position=position_jitterdodge(jitter.width=0.1, dodge.width=0.9)) +
    scale_fill_manual(values= c("white", "grey50", color)) +
    labs(x = "number of carbons", y = "log2 FC (MKC8866 vs control)", element_text(face = "bold", angle = 0)) +
    ylim(numero2, numero)
}


# For TAGs
unique(final_avg$Lipid_Class)

satuplot("Phosphatidylcholine"  , "#00FF00"  )
lenghtplot("Phosphatidylcholine"  , "#00FF00"  )


