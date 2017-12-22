###############################################################
## Tecnologías Ómicas y Bioinformática 2017/2018             ##
## Máster en Genética Molecular y Biotecnología              ##
## Universidad de Sevilla                                    ##
##                                                           ##
## Análisis transcriptómicos masivos basados en RNA-seq.     ##
##                                                           ##
## Prof. Francisco J. Romero-Campero fran@us.es              ##
###############################################################

## El paquete de bioconductor ballgown proporciona las funciones necesarias para
## realizar un análisis de expresión génica diferencial y visualizar los resultados
## a partir de procesamiento de los datos brutos de secuenciación realizados con
## hisat2 y stringtie.

## Para ejecutar con éxito este script es necesario descarga la carpeta samples
## completa a tu ordenador, mover este script a la carpeta samples y fijar el
## Working Directory To Source File Location.

## Instalación y carga de los paquetes necesarios. Sólo es necesario instalar los
## paquetes la primera vez que se ejecuta este script en un ordenador el resto de las
## veces bastará cargar los paquetes simplemente.
source("https://bioconductor.org/biocLite.R")

biocLite("ballgown")
biocLite("genefilter")
install.packages("dplyr")
install.packages("devtools")

library(ballgown)
library(genefilter)
library(dplyr)
library(devtools)

devtools::install_github('alyssafrazee/RSkittleBrewer')

library(RSkittleBrewer)

## Para cargar los datos es necesario crear previamente un fichero tabular
## que contenga como primera columna los nombres de las carpetas donde se guarda
## cada muesra típicamente sample1, sample2, etc ... El resto de columnas
## haran referencia al genotipo, tratamiento y demás caracteriśticas de cada muestra.

sewtd("/home/bms2017_08/dev/rnaseq_pipeline")
pheno.data <- read.csv("pheno_data.csv")

## La función ballgown se usa para cargar o leer los datos. Es necesario especificar
## el directorio donde se encuentra las muestras. En nuestro caso especificamos .
## para indicar que se encuentran en el actual directorio.
bg.data <- ballgown(dataDir = ".", samplePattern = "sample", pData=pheno.data)
bg.data
sampleNames(bg.data)

## La función gexpr extrae los niveles de expresión génicos en cada muestra
## medidos como FPKM (Fragments Per Kilobase of exon and Million of mapped reads)
gene.expression <- gexpr(bg.data)
head(gene.expression)

## Nombramos las columnas con los nombres de nuestras muestras.
colnames(gene.expression) <- c("WT1","WT2","Hp1a_1","Hp1a_2")

## Previsualizamos la similitud entre las réplicas
plot(log2(gene.expression[,1]+1),log2(gene.expression[,2]+1),pch=19,cex=0.7,xlab="Col-0",ylab=substitute(italic("atbmi1abc")),cex.lab=1.5)
plot(log2(gene.expression[,3]+1),log2(gene.expression[,4]+1),pch=19,cex=0.7,xlab="Col-0",ylab=substitute(italic("atbmi1abc")),cex.lab=1.5)

## Construimos un boxplot para comprobar que las distribuciones globales de las
## muestras son similares y comparables.
boxplot(log2(gene.expression + 1),col=rainbow(ncol(gene.expression)),ylab="log2(FPKM + 1)",cex.lab=1.5)

## Calculamos la matrix de expresión media.
col0 <- (gene.expression[,"col0_1"] + gene.expression[,"col0_2"])/2
abc <- (gene.expression[,"abc_1"] + gene.expression[,"abc_2"])/2

mean.expression <- matrix(c(col0,abc),ncol=2)
colnames(mean.expression) <- c("col0","abc")
rownames(mean.expression) <- rownames(gene.expression)

## Previsualizamos el efecto de la mutación en un scatterplot.
plot(log2(col0+1),log2(abc+1),pch=19,cex=0.7,xlab="Col-0",ylab=substitute(italic("atbmi1abc")),cex.lab=1.5)

## Pasamos a realizar un análisis de expresión génica diferencial.
## Primero debido a la alta sensibilidad del RNA-seq eliminamos genes con muy
## bajos niveles de expresión o que varían muy poco.
bg.data.filt <- subset(bg.data,"rowSums(texpr(bg.data)) >1 & rowVars(texpr(bg.data)) >1",genomesubset=TRUE)
bg.data.filt

## La función stattest calcula fold-change, p-valor y q-valor para realizar un análisis
## de expresión génica diferencial. El uso del paquete limma es recomendado cuando se tienen
## pocas réplicas.

results.genes <- stattest(bg.data.filt,feature="gene",covariate="genotype",getFC=TRUE, meas="FPKM")
class(results.genes)
head(results.genes)

## Extraemos los identidicaores de genes, fc y qvalor.
gene.ids <- as.vector(results.genes$id)
gene.fc <- as.vector(results.genes$fc)
gene.qval <- as.vector(results.genes$qval)

names(gene.fc) <- gene.ids
names(gene.qval) <- gene.ids

## Debido al gran efecto de la triple mutación elegimos un umbral para el fc de 4.
## En este caso los datos no están transformados por log2.
activated.genes <- gene.ids[gene.fc > 4]
repressed.genes <- gene.ids[gene.fc < 1/4]

length(activated.genes)
length(repressed.genes)

## Resaltar activados y reprimidos
plot(log2(col0+1),log2(abc+1),pch=19,cex=0.7,xlab="Col-0",ylab=substitute(italic("atbmi1abc")),cex.lab=1.5)
points(log2(col0[activated.genes]+1),log2(abc[activated.genes]+1),pch=19,cex=0.7,col="red",cex.lab=1.5)
points(log2(col0[repressed.genes]+1),log2(abc[repressed.genes]+1),pch=19,cex=0.7,col="blue",cex.lab=1.5)

## Guardamos en ficheros en formato txt los genes expresados de forma diferencial.
activated.genes.df <- subset(results.genes, id %in% activated.genes)
head(activated.genes.df)
activated.genes.df <- activated.genes.df[,2:5]
head(activated.genes.df)
write.table(x = activated.genes.df,file = "activated_genes.txt",quote = FALSE,row.names = FALSE,sep = "\t")

repressed.genes.df <- subset(results.genes, id %in% repressed.genes)
head(repressed.genes.df)
repressed.genes.df <- repressed.genes.df[,2:5]
head(repressed.genes.df)
write.table(x = repressed.genes.df,file = "repressed_genes.txt",quote = FALSE,row.names = FALSE,sep = "\t")

## Cuando se tiene pocas réplicas es más apropiado calcular FC directamente.
log2fc <- log2(abc / col0)

activated.genes <- names(which(log2fc > 2))
repressed.genes <- names(which(log2fc < -2))

length(activated.genes)
length(repressed.genes)

## Resaltar activados y reprimidos
plot(log2(col0+1),log2(abc+1),pch=19,cex=0.7,xlab="Col-0",ylab=substitute(italic("atbmi1abc")),cex.lab=1.5)
points(log2(col0[activated.genes]+1),log2(abc[activated.genes]+1),pch=19,cex=0.7,col="red",cex.lab=1.5)
points(log2(col0[repressed.genes]+1),log2(abc[repressed.genes]+1),pch=19,cex=0.7,col="blue",cex.lab=1.5)

plot(log2(col0),log2(abc),pch=19,cex=0.7,xlab="Col-0",ylab=substitute(italic("atbmi1abc")),cex.lab=1.5,xlim=c(0,15),ylim=c(0,15))
points(log2(col0[activated.genes]),log2(abc[activated.genes]),pch=19,cex=0.7,col="red",cex.lab=1.5)
points(log2(col0[repressed.genes]),log2(abc[repressed.genes]),pch=19,cex=0.7,col="blue",cex.lab=1.5)

## Las siguientes instrucciones realizan un análisis de enriquecimiento de términos GO
## y rutas metabólicas usando los paquetes clusterProfiler y pathview.
library(clusterProfiler)
library("pathview")

ego.activated <- enrichGO(gene = activated.genes, OrgDb="org.At.tair.db", ont = "BP", pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.05,universe=keys(org.At.tair.db),keytype="TAIR")

summary.ego.activated <- summary(ego.activated)
write.table(summary.ego.activated,file="GO_enrichment_activated_genes.txt",quote=FALSE,row.names = FALSE,sep="\t")


barplot(dropGO(ego.activated,level = 1), drop=TRUE,showCategory=10)

kegg.activated <- enrichKEGG(gene= activated.genes, organism = "arabidopsis", pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.05)
summary.kegg.activated <- as.data.frame(kegg.activated)
write.table(summary.kegg.activated,file="pathways_enrichment_activated.txt",quote=FALSE,sep = "\t",row.names = FALSE)
pathways.activated <- summary.kegg.activated[["Description"]]
pathways.activated.id <- summary.kegg.activated[["ID"]]

for(i in 1:length(pathways.activated.id))
{
  pathview(gene.data  = gene.fc[activated.genes],
           pathway.id = pathways.activated.id[i],
           species    = "ath",
           limit = list(gene = max(gene.fc[activated.genes]),cpd = 1),gene.idtype="KEGG")
}


ego.repressed <- enrichGO(gene = repressed.genes, OrgDb="org.At.tair.db", ont = "BP", pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.05,universe=keys(org.At.tair.db),keytype="TAIR")

summary.ego.repressed <- summary(ego.repressed)
write.table(summary.ego.repressed,file="GO_enrichment_repressed_genes.txt",quote=FALSE,row.names = FALSE,sep="\t")


barplot(dropGO(ego.repressed,level = 1), drop=TRUE,showCategory=10)

kegg.repressed <- enrichKEGG(gene= repressed.genes, organism = "arabidopsis", pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.05)
summary.kegg.repressed <- as.data.frame(kegg.repressed)
write.table(summary.kegg.repressed,file="pathways_enrichment_repressed.txt",quote=FALSE,sep = "\t",row.names = FALSE)
pathways.repressed <- summary.kegg.repressed[["Description"]]
pathways.repressed.id <- summary.kegg.repressed[["ID"]]

for(i in 1:length(pathways.repressed.id))
{
  pathview(gene.data  = gene.fc[repressed.genes],
           pathway.id = pathways.repressed.id[i],
           species    = "ath",
           limit = list(gene = max(gene.fc[repressed.genes]),cpd = 1),gene.idtype="KEGG")
}


## De forma similar podemos realizar un análisis de expresión diferencial de transcritos
## en lugar de genes.
transcript.expression <- texpr(bg.data, 'FPKM')
class(transcript.expression)
head(transcript.expression)

i <- 1
plotMeans(activated.genes[i], bg.data,groupvar="genotype",legend=FALSE)
i <- i + 1

i <- 16
plotMeans(activated.genes[i], bg.data,groupvar="genotype",legend=FALSE)

i <- 51
plotMeans(activated.genes[i], bg.data,groupvar="genotype",legend=FALSE)

i <- i + 1
gene.name <- activated.genes[i]
control.samples <- c("col0_1","col0_2")
treatment.samples <- c("abc_1","abc_2")
condition.names <- c("Col-0","AtBMI1abc")

head(gene.expression)
control.expression <- gene.expression[gene.name,control.samples]
treatment.expression <- gene.expression[gene.name,treatment.samples]

control.mean <- mean(control.expression)
treatment.mean <- mean(treatment.expression)

control.sd <- sd(control.expression)
treatment.sd <- sd(treatment.expression)

means <- c(control.mean,treatment.mean)
sds <- c(control.sd,treatment.sd)

ymax <- means + sds

pos <- barplot(means,col=rainbow(length(condition.names)),names.arg = condition.names,ylim=c(0,max(ymax)*1.1),main=gene.name,cex.main=2)
for(i in 1:length(condition.names))
{
  arrows(pos[i],means[i]-sds[i],pos[i],means[i]+sds[i],code=3,angle=90,length=0.1)
}
