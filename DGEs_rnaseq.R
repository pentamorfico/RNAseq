#Loading packages
install.packages("ggplot2")
devtools::install_github("kassambara/ggpubr")
devtools::install_github('alyssafrazee/RSkittleBrewer')

library(ggpubr)
library(ggplot2)
library(ballgown)
library(genefilter)
library(dplyr)
library(devtools)
library(plotly)
library(RSkittleBrewer)
library(AnnotationDbi)
library(org.Dm.eg.db)
library(DOSE)
library(GO.db)
library(org.Hs.eg.db)
library(topGO)
library(GSEABase)
library(clusterProfiler)
library(clusterProfiler)
library(pathview)
library(limma)
library(reshape2)

#Set Working Dir
setwd("~/Documents/samples")

#Read and process files
pheno.data <- read.csv("pheno_data.csv")
bg.data <- ballgown(dataDir = ".", samplePattern = "sample", pData=pheno.data)
gene.expression <- gexpr(bg.data)
colnames(gene.expression) <- c("WT1","WT2","Hp1a_1","Hp1a_2")

bg.data.filt <- subset(bg.data,"rowSums(texpr(bg.data)) >1 & rowVars(texpr(bg.data)) >1",genomesubset=TRUE)
results.genes <- stattest(bg.data,feature="gene",covariate="genotype",getFC=TRUE, meas="FPKM")

indices = match(results.genes$id, texpr(bg.data, 'all')$gene_id)
gene_names_for_result = texpr(bg.data, 'all')$gene_name[indices]
results.genes = data.frame(geneNames=gene_names_for_result, results.genes)

## Scatterplot replicas
d = data.frame( wt1 = log2(gene.expression[,1]+1) )
d$wt2 <- log2(gene.expression[,2]+1)
d$Hp1a_1 = log2(gene.expression[,3]+1)
d$Hp1a_2 <- log2(gene.expression[,4]+1)

p1 = ggplot(d,aes(wt1, wt2, color=wt1)) +
  geom_point(shape = 16, size = 3, show.legend = FALSE, alpha = .1) +
  theme_minimal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ scale_color_gradient(low = "#0091ff", high = "#f0650e")

p2 = ggplot(d,aes(Hp1a_1, Hp1a_2, color=Hp1a_1)) +
  geom_point(shape = 16, size = 3, show.legend = FALSE, alpha = .1) +
  theme_minimal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ scale_color_gradient(low = "#0091ff", high = "#f0650e")

ggarrange(p1, p2,
          labels = c("WT", "Hp1a"),
          ncol = 2, nrow = 1)

## Boxplot
ggplot(stack(d), aes(x=ind, y=values, fill=ind)) + geom_boxplot() + theme_minimal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank())+ scale_color_gradient(low = "#0091ff", high = "#f0650e")

## Mean expression matrix
head(gene.expression)
wt <- (gene.expression[,"WT1"] + gene.expression[,"WT2"])/2
hp1a <- (gene.expression[,"Hp1a_1"] + gene.expression[,"Hp1a_2"])/2
mean.expression <- matrix(c(wt,hp1a),ncol=2)
colnames(mean.expression) <- c("wt","hp1a")
rownames(mean.expression) <- rownames(gene.expression)

## Scatterplot condition1 vs condition2
d2 = data.frame( wt = log2(wt+1) )
d2$hp1a <- log2(hp1a+1)

p1 = ggplot(d2,aes(wt, hp1a, color=wt)) +
  geom_point(shape = 16, size = 3, show.legend = FALSE, alpha = .1) +
  theme_minimal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ scale_color_gradient(low = "#0091ff", high = "#f0650e")
p1

## Get DGEs by fold change
log2fc <- log2((hp1a+1)/(wt+1))
activated.genes <- names(which(log2fc > 2))
repressed.genes <- names(which(log2fc < -2))
length(activated.genes)
length(repressed.genes)

## Save DGEs
activated.genes.df <- subset(results.genes, id %in% activated.genes)
write.table(x = activated.genes.df[,1],file = "activated_genes.txt",quote = FALSE,row.names = FALSE,sep = "\t")
length(repressed.genes.df[,2])
repressed.genes.df <- subset(results.genes, id %in% repressed.genes)
write.table(x = repressed.genes.df[,1],file = "repressed_genes.txt",quote = FALSE,row.names = FALSE,sep = "\t")

## Scatterplot DGEs
unregulated=setdiff(names(wt), names(wt[activated.genes]))
unregulated=setdiff(unregulated, names(wt[repressed.genes]))

d1 <- data.frame(wt=log2(wt[unregulated]+1), hp1a=log2(hp1a[unregulated]+1), DEGs="unregulated")
d2 <- data.frame(wt=log2(wt[activated.genes]+1) , hp1a=log2(hp1a[activated.genes]+1), DEGs="activated")
d3 <- data.frame(wt=log2(wt[repressed.genes]+1) , hp1a=log2(hp1a[repressed.genes]+1), DEGs="repressed")

dataframe_plot=rbind(d1,d2,d3)

ggplot(dataframe_plot, aes(x=wt, y=hp1a, color=DEGs)) +
  geom_point(shape=16,size=3,alpha = .3) +
  scale_color_manual(values=c('#999999','#5ec335','#a90f40'))+
  theme(legend.position="top") + theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

## Gene Ontology term enrichment

### Activados ###

eg.activated = bitr(activated.genes.df[,1], fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Dm.eg.db")
ids.activated = bitr(activated.genes.df[,1], fromType="SYMBOL", toType=c("UNIPROT", "ENSEMBL"), OrgDb="org.Dm.eg.db")
ego.activated <- enrichGO(gene=eg.activated[,2],OrgDb=org.Dm.eg.db,ont="BP",pAdjustMethod="BH",pvalueCutoff=0.05,readable=TRUE)
summary.ego.activated <- as.data.frame(ego.activated)
write.table(summary.ego.activated,file="GO_enrichment_activated_genes.txt",quote=FALSE,row.names = FALSE,sep="\t")

barplot(dropGO(ego.activated,level = 1), drop=TRUE,showCategory=10)

kegg.activated <- enrichKEGG(gene= eg.activated[,2],keyType ='ncbi-geneid',organism = "dme", pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.05)
summary.kegg.activated <- as.data.frame(kegg.activated)
write.table(summary.kegg.activated,file="pathways_enrichment_activated.txt",quote=FALSE,sep = "\t",row.names = FALSE)
pathways.activated <- summary.kegg.activated[["Description"]]
pathways.activated.id <- summary.kegg.activated[["ID"]]
head(pathways.activated)
head(activated.genes.df[,4])

for(i in 1:length(pathways.activated.id))
{
  pathview(gene.data  = activated.genes.df[,4],
           pathway.id = pathways.activated.id[i],
           species    = "dme",
           limit = list(gene = max(activated.genes.df[,4]),cpd = 1),gene.idtype="KEGG")
}

### Reprimidos ###

eg.repressed = bitr(repressed.genes.df[,1], fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Dm.eg.db")
ids.repressed = bitr(repressed.genes.df[,1], fromType="SYMBOL", toType=c("UNIPROT", "ENSEMBL"), OrgDb="org.Dm.eg.db")
ego.repressed <- enrichGO(gene=eg.repressed[,2],OrgDb=org.Dm.eg.db,ont="BP",pAdjustMethod="BH",pvalueCutoff=0.05,readable=TRUE)
summary.ego.repressed <- as.data.frame(ego.repressed)
write.table(summary.ego.repressed,file="GO_enrichment_repressed_genes.txt",quote=FALSE,row.names = FALSE,sep="\t")

barplot(dropGO(ego.repressed,level = 1), drop=TRUE,showCategory=10)

kegg.repressed <- enrichKEGG(gene= eg.repressed[,2],keyType ='ncbi-geneid',organism = "dme", pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.05)
summary.kegg.repressed <- as.data.frame(kegg.repressed)
write.table(summary.kegg.repressed,file="pathways_enrichment_repressed.txt",quote=FALSE,sep = "\t",row.names = FALSE)
pathways.repressed <- summary.kegg.repressed[["Description"]]
pathways.repressed.id <- summary.kegg.repressed[["ID"]]
head(pathways.repressed)
head(repressed.genes.df[,4])

for(i in 1:length(pathways.repressed.id))
{
  pathview(gene.data  = repressed.genes.df[,4],
           pathway.id = pathways.repressed.id[i],
           species    = "dme",
           limit = list(gene = max(repressed.genes.df[,4]),cpd = 1),gene.idtype="KEGG")
}

## Getting barplots and transcript expression

tropical <- c('darkorange', 'dodgerblue', 'hotpink', 'limegreen', 'yellow')
palette(tropical)

#Getting the identifier of DGEs
activated.genes.df %>% arrange(fc) #Shows the activated genes by fc
identifier <- names(which(ballgown::geneIDs(bg.data) == "MSTRG.2592"))
identifier=identifier[1]

## Plot individual transcripts
fpkm <- texpr(bg.data, meas='FPKM')
fpkm <- log2(fpkm +1)
colnames(fpkm) <- c("WT1","WT2","Hp1a_1","Hp1a_2")

ballgown::transcriptNames(bg.data.filt)[identifier]
plot(fpkm[identifier,] ~ pheno.data$genotype, border=c(1,2),
     main=paste(ballgown::geneNames(bg.data)[identifier], ' : ',ballgown::transcriptNames(bg.data)[identifier]),
     pch=19, xlab="Genotype", ylab='log2(FPKM+1)')
points(fpkm[identifier,] ~ jitter(as.numeric(pheno.data$genotype)), col=as.numeric(pheno.data$genotype))

## Plot average expression
plotMeans(ballgown::geneIDs(bg.data)[identifier], bg.data, groupvar="genotype", legend=TRUE)
