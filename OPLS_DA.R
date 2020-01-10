### OPLS - DA
library(ropls)
library(ggplot2)
library(RCurl)
library(plyr)
library(fgsea)

# load omics data
mfa_data <- read.csv("~/GitHub/OPLS-DA_Integration/data/mfa_data_NEW.csv", row.names = 1)

### PCA
data.pca <- opls(mfa_data, predI=10)

eig <- data.pca@pcaVarVn/sum(data.pca@pcaVarVn)
variances <- barplot(eig, names.arg=names(eig),main = "Variance explained by dimension", xlab = "Dimensions",ylab = "Percentage of variance", col ="steelblue")
lines(x = variances, eig, type="b", pch=19, col = "red")

eig
PCAi <- as.data.frame(data.pca@scoreMN)
ggplot(PCAi , aes(p1, p2, col= c(rep("DMSO", 5), rep("MKC", 5)))) + 
  geom_point(aes(size=2)) + 
  geom_text(aes(label = row.names(PCAi)), nudge_y = 1) + 
  labs (y = "PC2 (23%)", x = "PC1 (37%)") + # substitute with % variances
  theme_light()

plot(data.pca,
     typeVc = "x-score",
     parAsColFcVn = c(rep("DMSO", 5), rep("MKC", 5)))


### OPLS-DA
classes <- c(rep("DMSO", 5), rep("MKC", 5))
data.oplsda.pareto <- opls(mfa_data, classes,
                        predI = 1, orthoI = 1, scaleC="pareto")

scores_OPLS <- as.data.frame(cbind(data.oplsda.pareto@orthoScoreMN,data.oplsda.pareto@scoreMN)) 
colnames(scores_OPLS) <- c( "to[1]","t[1]")

ggplot(scores_OPLS , aes(scores_OPLS$`t[1]`, scores_OPLS$`to[1]`, col= c(rep("DMSO", 5), rep("MKC", 5)))) + 
  geom_point(aes(size=2)) + 
  geom_text(aes(label = row.names(scores_OPLS)), nudge_y = 1) + 
  stat_ellipse(geom = "path", alpha = 0.5, level = 0.9, linetype = 1, size=3) +
  labs (y = "to[1]", x = "t[1]") + # substitute with % variances
  theme_classic()

### extract loadings
loadings_OPLS_pareto <- as.data.frame(cbind(loadings.x = data.oplsda.pareto@loadingMN,
                                     weight.x = data.oplsda.pareto@weightMN,
                                     importance.x = data.oplsda.pareto@vipVn,
                                     loadings.y =data.oplsda.pareto@orthoLoadingMN,
                                     weight.y = data.oplsda.pareto@weightMN,
                                     importance.y = data.oplsda.pareto@orthoVipVn))

colnames(loadings_OPLS_pareto) <- c("loadings.x", "weight.x", "importance.x", "loadings.y", "weight.y", "importance.y")

# annotate loadings
annotacion_todo <- read.csv(text = getURL("https://raw.githubusercontent.com/xaitorx/OPLS-DA_Integration/data/var_dim2_annot.csv"))
annotacion_todo <- annotacion_todo[,c(1,2,6,7)]
annotacion_todo$X <- gsub("-", ".", annotacion_todo$X)

loadings_OPLS_annot <- cbind(annotacion_todo[match(row.names(loadings_OPLS_pareto), annotacion_todo$X),], loadings_OPLS_pareto)

#write file
write.csv(loadings_OPLS_annot, "loadings_OPLS_annot.csv", row.names = FALSE)

# compare with Foldhanges MKC/DMSO
FoldChanges <-  colMeans(mfa_data[6:10,], na.rm = T) - colMeans(mfa_data[1:5,], na.rm = T)
loadings_OPLS_annot$FoldChanges <- FoldChanges

# compare loadings.x with FC DMSO/MKC
ppp <- ggplot(loadings_OPLS_annot, aes(loadings_OPLS_annot$loadings.x, loadings_OPLS_annot$FoldChanges)) 
ppp + geom_point(aes(colour = loadings_OPLS_annot$type, size = 0.1, alpha = 0.5)) +
  labs (y = "FoldChanges", x = "loadings.x")

# Calculate contribution of each layer
total <- sum(loadings_OPLS_annot$importance.x)
contributions <- c(sum(loadings_OPLS_annot$importance.x[1:15445])/total,
                   sum(loadings_OPLS_annot$importance.x[15446:15657])/total,
                   sum(loadings_OPLS_annot$importance.x[15658:16455])/total)

names(contributions) <- c("mRNA", "miRNA", "metabolome")
contributions <- as.data.frame(contributions)


ggplot(data.frame(contributions),aes(seq_along(contributions),contributions, fill= c("mRNA", "miRNA", "metabolome")))+geom_bar(stat="identity")

### Inspect top variables for each group
loadings_OPLS_annot <- loadings_OPLS_annot[order(loadings_OPLS_annot$importance.x, decreasing = TRUE),]
loadings_OPLS_annot$Y <- factor(loadings_OPLS_annot$Y ,levels = unique(loadings_OPLS_annot$Y))

mRNA_comp1 <- subset(loadings_OPLS_annot, loadings_OPLS_annot$type== "mRNA")
miRNA_comp1 <- subset(loadings_OPLS_annot, loadings_OPLS_annot$type== "miRNA")
metabolome_comp1 <- subset(loadings_OPLS_annot, loadings_OPLS_annot$type== "metabolome")

# plot mRNA group
grafico <- ggplot(mRNA_comp1[1:20,], aes(y=loadings.x, x=Y))
grafico + geom_segment(aes(colour= mRNA_comp1$type[1:20], x=Y, xend=Y, y=0,size=1, yend=loadings.x)) +   
  geom_point( aes(colour=mRNA_comp1$type[1:20]), size=6, alpha=0.8, shape=19, stroke=2) +   
  theme_classic() +   
  theme(plot.margin = margin(1.5, 1.5, 1.5, 1.5, "cm") ,axis.text.x = element_text(angle = 45, hjust = 1)) +   
  labs (title= "Top mRNA contributors to component1",y = "loadings.x", x = "") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_colour_manual(values = c("steelblue","steelblue"))

# plot miRNA group
grafico <- ggplot(miRNA_comp1[1:20,], aes(y=loadings.x, x=Y))
grafico + geom_segment(aes(colour= miRNA_comp1$type[1:20], x=Y, xend=Y, y=0,size=1, yend=loadings.x)) +   
  geom_point( aes(colour=miRNA_comp1$type[1:20]), size=6, alpha=0.8, shape=19, stroke=2) +   
  theme_classic() +   
  theme(plot.margin = margin(1.5, 1.5, 1.5, 1.5, "cm") ,axis.text.x = element_text(angle = 45, hjust = 1)) +   
  labs (title= "Top miRNA contributors to component1",y = "loadings.x", x = "") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_colour_manual(values = c("green4","green4"))

# plot metabolome group. remove ambiguous variables
metabolome_comp1 <- metabolome_comp1[!is.na(metabolome_comp1$REACTOME),]

grafico <- ggplot(metabolome_comp1[1:20,], aes(y=loadings.x, x=Y))
grafico + geom_segment(aes(colour= metabolome_comp1$type[1:20], x=Y, xend=Y, y=0,size=1, yend=loadings.x)) +   
  geom_point( aes(colour=metabolome_comp1$type[1:20]), size=6, alpha=0.8, shape=19, stroke=2) +   
  theme_classic() +   
  theme(plot.margin = margin(1.5, 1.5, 1.5, 1.5, "cm") ,axis.text.x = element_text(angle = 45, hjust = 1)) +   
  labs (title= "Top metabolome contributors to component1",y = "loadings.x", x = "") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_colour_manual(values = c("indianred1","indianred1"))


### PATHWAY ENRICHMENT WITH FGSEA
## Download Reactome annotation from Reactome website
## All levels of the pathway hierarchy
# ChEBI to All pathways
# ENSEMBL to All pathways
# miRBase to All pathways
# Paste them together rbind()

### data frame with ENSEMBL/CHEBI/miRNA identifiers annotated to reactome pathways
### split into 1 list of IDs/pathway
Reactome_all <- read.csv(text = getURL("https://raw.githubusercontent.com/xaitorx/OPLS-DA_Integration/data/Reactome_all.csv"), stringsAsFactors=FALSE)
Reactome_all <- Reactome_all[order(Reactome_all$pathway),]

Reactome_all <- unique(Reactome_all)

#count numer of elements/pathway
counts <- as.data.frame(table(Reactome_all$pathway))
magic_number <- length(unique(counts$Var1))

breaks <- rep(1:magic_number, counts$Freq)

# break IDs into pathways, rename each list accordingly
lista1 <- split(Reactome_all$ID, breaks)
names(lista1) <- unique(counts$Var1)

# save for later
anotacion <- as.data.frame(unique(Reactome_all[,2:3]))

# save list of lists used in fgsea in .gmt format (for Cytoscape)
ooo <- plyr::ldply(lista1, rbind)
reactome_gmt <- cbind(anotacion, ooo)
reactome_gmt <- reactome_gmt[,-3]

reactome_gmt_unique = reactome_gmt[!duplicated(reactome_gmt$name),]
write.csv(reactome_gmt_unique, file = "reactome_gmt_unique.csv", row.names=FALSE)


write.table(reactome_gmt, file = "reactome.gmt", row.names=FALSE, col.names = FALSE, sep="\t")

### Perform GSEA with these lists on our preranked list of IDs
ranking <- loadings_OPLS_annot[,c(5,4)]
ranking <- ranking[!is.na(ranking$REACTOME),]

# consolidate duplicated IDs, order ranking
ranking <- ddply(ranking,2,numcolwise(mean))
ranking <- ranking[order(ranking$loadings.x, decreasing = TRUE),]

ranking_ooo <- ranking$loadings.x
names(ranking_ooo) <- ranking$REACTOME

###Perform GSEA with list of reactome pathways
fgseaRes <- fgsea(pathways = lista1, 
                  stats = ranking_ooo,
                  minSize=5,
                  maxSize=500,
                  nperm=100000)


fgseaRes_annot <- merge(anotacion, fgseaRes, by.x = 1, by.y = 1)

#save results
write.csv(fgseaRes_annot[,-9], "fgseaRes_annot.csv")

### plot desired term gsea
plotEnrichment(lista1[[(fgseaRes[447,])$pathway]],ranking_ooo) +
  ggtitle("XBP1(S) activates chaperone genes") 

# change for number of row
### plot
plotEnrichment(lista1[[(fgseaRes[59,])$pathway]],ranking_ooo) +
  ggtitle("Metabolism of proteins") 




https://www.bioconductor.org/packages/devel/bioc/vignettes/ropls/inst/doc/ropls-vignette.html#the-ropls-package
