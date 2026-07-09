#Required Packages
library(WGCNA)
library(flashClust)
library(gplots)
library(cluster)
library(igraph)
library(RColorBrewer)
library(readxl)
library(org.Hs.eg.db)
library(devtools)
library(GO.db)
library(AnnotationDbi)
library(GOstats)
library(dplyr)
library(tidyr)
library(sigora)
library(biomaRt)

#set working directory
setwd("C:/Aatmika/PSP_PICK_WGCNA")


options (stringsAsFactors = FALSE)
allowWGCNAThreads()

#input the data 
p_data <-  read.csv("PSP_PICK_newproteins.csv")

p_data <- p_data[,-c(1,3:5)]

rownames(p_data) <- p_data$Gene.names

at100 <- read_xlsx("AT100_new.xlsx")

#take a quick look to data
dim(p_data)
names(p_data)

#remove the sample data and trnaspose the expression data
datExpr0 = as.data.frame(t(p_data))
names(datExpr0) = p_data$`Gene.names`
rownames(datExpr0) = names(p_data)
datExpr0 <- datExpr0[-c(1),]

datExpr0[] <- lapply(datExpr0, as.numeric)

#check for missing values
gsg = goodSamplesGenes(datExpr0, verbose = 3)

#cluster the samples
sampleTree = hclust(dist(datExpr0), method = "average")

#plot the sample tree

par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detcted outliers", sub = "", xlab = "", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

#plot a line to show the cut
abline(h = 0.045, col = "red")

#Determine the cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 0.045, minSize = 0.045)
table(clust)

#clust 1 contains samples we want to keep 
keepSamples = (clust == c(1,2))
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

##Chosing soft threshold power
#choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from= 12, to =20, by =2))

#call the network topology analysis fucntion
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)

#plot the results
sizeGrWindow(9,5)
par(mfrow = c(1,2))
cex1 = 0.9

#scale free topology fit index as a function of soft thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[, 3])*sft$fitIndices[,2],
     xlab = "Soft threshold (power)", 
     ylab = "Scale free Topology Model Fit, signed R^2", type = "n", 
     main = paste("Scale independence"))
text(sft$fitIndices[, 1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels = powers, cex = cex1, col = "red")

#the line corresponds to using R^2 cut-off of h
abline(h = 0.90, col = "red")

#Mean connectivity as a fucntion of the thresholding power 
plot(sft$fitIndices[, 1], sft$fitIndices[, 5], 
     xlab = "Soft threshold (power)", 
     ylab = "Mean Connectivity", type = "n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, cex = cex1, col = "red")

#We calculate adjancies using soft thresholding power 5 
softPower = 16
adjacency = adjacency(datExpr0, power = softPower)

#topological overlap matrix
#Turn adjacency to topological overlap
TOM = TOMsimilarity(adjacency)
dissTOM = 1-TOM

#clustering using TOM

#Call the hierarichical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average")
#plot the resulting clustering tree
sizeGrWindow(12,9)
plot(geneTree, xlab = "", sub = "", main = "Gene clustering on TOM-based dissimilarity", 
     labels = FALSE, hang = 0.04)

#For large modules we set min module size relatively high
minModuleSize = 30

#Module identification using dynamic tree cut
dynamicMods = cutreeDynamic(dendro= geneTree, distM = dissTOM, 
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize
)
table(dynamicMods)

#cut numerical labels into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

#plot the dendogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", 
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05, 
                    main = "Gene dendrogram and module colors")

#Merging modules whose expression profiles are very similar 

#Calculate eigengenes
MEList = moduleEigengenes(datExpr0, colors = dynamicColors)
MEs = MEList$eigengenes

#Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs)

#Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average")

#plot the results
sizeGrWindow(7,6)
plot(METree, main = "Clustering of module eigengenes", 
     xlab = "", sub = "")

MEDissThres = 0.15

#plot the cut line into dendogram
abline(h = MEDissThres, col = "red")

#call an automatic merging function
merge = mergeCloseModules(datExpr0, dynamicColors, cutHeight = MEDissThres, verbose = 3)

#the merged module colors
mergedColors = merge$colors

#Eigengenes of the new merged modules
mergedMEs = merge$newMEs

sizeGrWindow(12,9)

plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

write.csv(mergedMEs, "Eigengenes.csv")

#Rename to module Colors
moduleColors = mergedColors

#Module-trait relationship

nGenes = ncol(datExpr0)
nSamples = nrow(datExpr0)

#moduleColors = moduleColors[moduleColors !='grey']

turq = moduleColors[moduleColors == 'turquoise']
black = moduleColors[moduleColors == 'black']
green = moduleColors[moduleColors == 'green']
pink = moduleColors[moduleColors == 'pink']

#all gene names
probes = names(datExpr0)
inModule = is.finite(match(moduleColors, pink))

#gene for specific module
modProbes = probes[inModule]

#extract ribaq values for each module

pink_data = datExpr0[modProbes]
write.csv(pink_data, "TCIT_new_gene_ribaq.csv")

#recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr0, moduleColors)$eigengenes
#MEs0 = MEs0[MEs0!='grey']
MEs_1 =orderMEs(MEs0)
moduleTraitCor = cor(MEs_1, at100[,-c(1:3)], use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

sizeGrWindow(10,6)

testMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")

dim(testMatrix) = dim(moduleTraitCor)
par(mar = c(9, 8.5, 3, 3))

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(at100[,-c(1:3)]),
               yLabels = names(MEs_1),
               ySymbols = names(MEs_1),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = testMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))


#Hub proteins with 

#PDGFRB
PDGFRB <- as.data.frame(datExpr0$PDGFRB)

rownames(PDGFRB) = rownames(datExpr0)

PDGFRB_new$ribaq <- as.data.frame(PDGFRB[-c(3,4,10:14,18,19),])
PDGFRB_new <- PDGFRB_new[,-1]
colnames(PDGFRB_new) = "ribaq"

plot(PDGFRB_new$ribaq[-c(1,2)] ~ mito_1$`Ratio 2-6h_D3`[-c(1,2)])
cor(PDGFRB_new$ribaq[-c(1,2)] , mito_1$`Ratio 2-6h_D3`[-c(1,2)])

model_3 <- lm(PDGFRB_new$ribaq[-c(1,2)] ~ mito_1$`Ratio 2-6h_D3`[-c(1,2)])
summary(model_3)


#ECE1
ECE1 <- as.data.frame(datExpr0$ECE1)

rownames(ECE1) = rownames(datExpr0)

ECE1_new <- as.data.frame(ECE1[-c(3,4,10:14,18,19),])
#ECE1_new <- ECE1_new[,-1]
colnames(ECE1_new) = "ribaq"

plot(ECE1_new$ribaq[-c(1,2)] ~ mito_1$`Ratio 2-6h_D3`[-c(1,2)])
cor(ECE1_new$ribaq[-c(1,2)] , mito_1$`Ratio 2-6h_D3`[-c(1,2)])

model_4 <- lm(ECE1_new$ribaq[-c(1,2)] ~ mito_1$`Ratio 2-6h_D3`[-c(1,2)])
summary(model_4)

#MTCH2
MTCH2 <- as.data.frame(datExpr0$MTCH2)

rownames(MTCH2) = rownames(datExpr0)

MTCH2_new <- as.data.frame(MTCH2[-c(3,4,10:14,18,19),])
#ECE1_new <- ECE1_new[,-1]
colnames(MTCH2_new) = "ribaq"

plot(MTCH2_new$ribaq[-c(1,2)] ~ mito_1$`Ratio 2-6h_D3`[-c(1,2)])
cor(MTCH2_new$ribaq[-c(1,2)] , mito_1$`Ratio 2-6h_D3`[-c(1,2)])

model_5 <- lm(MTCH2_new$ribaq[-c(1,2)] ~ mito_1$`Ratio 2-6h_D3`[-c(1,2)])
summary(model_5)

#SAMM50
SAMM50 <- as.data.frame(datExpr0$SAMM50)

rownames(SAMM50) = rownames(datExpr0)

SAMM50_new <- as.data.frame(SAMM50[-c(3,4,10:14,18,19),])
#ECE1_new <- ECE1_new[,-1]
colnames(SAMM50_new) = "ribaq"

plot(SAMM50_new$ribaq[-c(1,2)] ~ mito_1$`Ratio 2-6h_D3`[-c(1,2)])
cor(SAMM50_new$ribaq[-c(1,2)] , mito_1$`Ratio 2-6h_D3`[-c(1,2)])

model_5 <- lm(SAMM50_new$ribaq[-c(1,2)] ~ mito_1$`Ratio 2-6h_D3`[-c(1,2)])
summary(model_5)



model_7 <- lm(mito_1$`Ratio 2-6h_D3` ~ PDGFRB_new$ribaq + ECE1_new$ribaq + MTCH2_new$ribaq + SAMM50_new$ribaq)
summary(model_7)



#plot eigengene heatmap

par(cex = 1.0)

plotEigengeneNetworks(MEs_1, "Eigengene Network", marHeatmap = c(3,4,2,2),
                      marDendro = c(0,4,1,2), cex.adjacency = 0.3, plotDendrograms = TRUE,
                      xLabelsAngle = 90, heatmapColors = blueWhiteRed(100)[51:100])


datExpr0$Staging = c("CTL", "CTL", "CTL", "CTL", "PSP", "PSP", "PSP", "PSP","PSP", "PSP", "PSP", "PSP", "PSP", "PSP", "PICK", "PICK", "PICK", "PICK", "PICK" )

toplot=t(MEs_1)
cols=substring(colnames(MEs_1),3,20)
par(mfrow=c(3,3))
par(mar=c(3,2,2,2))


for (i in 1:nrow(toplot)) {
  boxplot(toplot[i,]~factor(as.vector(as.factor(datExpr0$Staging)),c('CTL','PSP', 'PICK')),col=cols[i],ylab="ME",main=rownames(toplot)[i],xlab=NULL,las=2)
  # verboseScatterplot(x=as.numeric(targets.Ref.AD_1$Age),y=toplot[i,],xlab="Age",ylab="ME",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col=cols[i],pch=19)
  # boxplot(toplot[i,]~factor(targets.Ref.AD_1$Gender),col=cols[i],ylab="ME",main=rownames(toplot)[i],xlab=NULL,las=2)
}

#Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs

symbols <- p_data$'Gene.names'

ENSEMBL <- mapIds(org.Hs.eg.db, symbols, 'ENSEMBL', 'SYMBOL')

annot <- data.frame(p_data$'Gene.names', ENSEMBL)

probes <- names(datExpr0)
probes2annot <- match(probes, annot$p_data.Gene.names)


allENT = annot$ENSEMBL[probes2annot]

ME_1 <- as.data.frame(cbind(annot$ENSEMBL, probes, moduleColors ))

colnames(ME_1)[1] = "ENSEMBL.Gene.ID"
colnames(ME_1)[2] = "GeneSymbol"
colnames(ME_1)[3] = "Initially.Assigned.Module.Color"

write.csv(ME_1, "geneInfo.csv")

ME_1$SystemCode = rep("En", length = nrow(ME_1))
background = ME_1[, "ENSEMBL.Gene.ID"]
background = as.data.frame(background)

#output files for GO elite
background <- cbind(background, rep("En", length = length(background)))
colnames(background) <- c("Source Identifier", "SystemCode")
write.table(background, "./geneInfo/background/denominator.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

uniquemodcolors=unique(moduleColors)
uniquemodcolors=uniquemodcolors[uniquemodcolors!='grey']


#i = Number of modules
for (i in 1:length(uniquemodcolors)) {
  thismod= uniquemodcolors[i]
  ind = which(colnames(ME_1)==paste("ME", thismod, sep = ""))
  thisInfo=ME_1[ME_1$'Initially.Assigned.Module.Color'==thismod, ]
  colnames(thisInfo) <- c("Source Identifier", "SystemCode", "ME")
  write.table(thisInfo, file = paste("./geneInfo/input/", thismod, "_Module.txt", sep = ""), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
}

#Run GO elite as nohupped shell script


codedir <- "C:/Aatmika/PSP_PICK_WGCNA/GO-Elite_v.1.2.5-Py"
pathname <- "C:/Aatmika/PSP_PICK/geneInfo"

nperm = 10000

#system(paste("nohup python3.10.9", codedir, "/GO_Elite.py  --species Hs  --mod ENTREZ  --permutations", nperm, "  --method, 'z-score' , --zscore 1.96 --pval 0.01  --num 5 --input ",
#pathname, "/input --denom ", pathname, "/background  --output ", pathname, "/output &", sep = "" ))

python_script_path <- paste(codedir, "GO_Elite.py", sep = "/")

python_version <- "C:/Python27/python.exe"

#command <- paste("nohup", python_version, python_script_path)

command <- paste("nohup", python_version, python_script_path, "--species Hs --mod ENSEMBL --permutations", nperm, "  --method \"z-score\" --zscore 1.96 --pval 0.01 --num 5 --input ",pathname,
                 "/input --denom ",pathname,"/background --output ",pathname,"/output &",sep="")

system(command)

pathname <- "C:/Aatmika/PSP_PICK_WGCNA/geneInfo/Output/GO-ELite_results/CompleteResults"

uniquemodcolors=uniquemodcolors[-c(14)] # For some reason sometimes modules are not run correctly, therefore won't be able to be plotted so they are excluded

pdf("GOElite_plot_Modules_MF.pdf",height=8,width=12)
for(i in 1:length(uniquemodcolors)){
  thismod= uniquemodcolors[i]
  tmp=read.csv(file=paste(pathname,"/ORA_pruned/",thismod,"_Module-GO_z-score_elite.txt",sep=""),sep="\t")
  # tmp=subset(tmp,Ontology.Type!= c('cellular_component'))
  tmp=tmp[,c(2,9)] ## Select GO-terms and Z-score
  tmp=tmp[order(tmp$Z.Score,decreasing=T),] #
  if (nrow(tmp)<10){
    tmp1=tmp ## Take top 10 Z-score
    tmp1 = tmp1[order(tmp1$Z.Score),] ##Re-arrange by increasing Z-score
    par(mar=c(5,40,5,2))
    barplot(tmp1$Z.Score,horiz=T,col="blue",names.arg= tmp1$Ontology.Name,cex.names=1.2,las=1,main=paste("Gene Ontology Plot of",thismod,"Module"),xlab="Z-Score")
    abline(v=2,col="red")
  } else {
    tmp1=tmp[c(1:10),] ## Take top 10 Z-score
    tmp1 = tmp1[order(tmp1$Z.Score),] ##Re-arrange by increasing Z-score
    par(mar=c(5,40,5,2))
    barplot(tmp1$Z.Score,horiz=T,col="blue",names.arg= tmp1$Ontology.Name,cex.names=1.2,las=1,main=paste("Gene Ontology Plot of",thismod,"Module"),xlab="Z-Score")
    abline(v=2,col="red")
  }
  
  cat('Done ...',thismod,'\n')
}

dev.off()


#cell type enrichment
OR <- function(q,k,m,t) {
  q #<-  ## Intersection of test list and reference list, aka number of white balls drawn
  m #<-  ## All genes in reference list, aka number of draws
  k #<-  ## All genes in test list, aka number white balls
  t #<-  ## Total number of genes assessed, aka black plus white balls
  
  fisher.out <- fisher.test(matrix(c(q, k-q, m-q, t-m-k+q), 2, 2),conf.int=TRUE)
  OR <- fisher.out$estimate
  pval <- fisher.out$p.value
  upCI <- fisher.out$conf.int[1]
  downCI <- fisher.out$conf.int[2]
  
  output <- c(OR,pval,upCI,downCI)
  names(output) <- c("OR","Fisher p","-95%CI","+95%CI")
  return(output)
}




ORA <- function(testpath,refpath,testbackground,refbackground) {
  q <- length(intersect(testpath,refpath)) ## overlapped pathway size
  k <- length(intersect(refpath,testbackground))  ## input gene set
  m <- length(intersect(testpath,refbackground)) ## input module
  t <- length(intersect(testbackground,refbackground)) ## Total assessed background (intersect reference and test backgrounds)
  
  empvals <- OR(q,k,m,t)
  
  tmpnames <- names(empvals)
  empvals <- as.character(c(empvals,q,k,m,t,100*signif(q/k,3)))
  names(empvals) <- c(tmpnames,"Overlap","Reference List","Input List","Background","% List Overlap")
  return(empvals)
}


geneInfo <- ME_1

datKME <- geneInfo[,c("ENSEMBL.Gene.ID","Initially.Assigned.Module.Color")] ## Get a list of genes to test for enrichment, e.g. genes with modules defined
testbackground <- as.character(geneInfo$ENSEMBL.Gene.ID) # background list
#datKME=subset(datKME,Initially.Assigned.Module.Color!="grey")
namestestlist <- names(table(datKME[,2])) ## module
multiTest <- vector(mode = "list", length = length(namestestlist))
names(multiTest) <- namestestlist  

for (i in 1:length(multiTest))
{
  multiTest[[i]] <- datKME[datKME[,2]==namestestlist[i],1]
}  

cells <- read.csv("C:/Aatmika/PSP_PICK_WGCNA/cell_type_markers.csv")  

datcell <- pivot_longer(cells, cols = c('Astrocytes', 'Microglia', 'Neuron', 'Oligodendrocytes', 'Endothelia'), names_to = 'Cell_type', values_to = 'Gene_names')

datCells <- as.data.frame(datcell %>% select(Gene_names, everything()))

cell_names <- mapIds(org.Hs.eg.db, datCells$Gene_names, 'ENSEMBL', 'SYMBOL')

datCells$Gene_names <- cell_names


## Set up reference lists
namesreflist <- names(table(datCells[,2])) ## category or module color
multiRef <- vector(mode = "list", length = length(namesreflist))
names(multiRef) <- namesreflist
for (i in 1:length(multiRef))
{
  multiRef[[i]] <- datCells[datCells[,2]==namesreflist[i],1]
}

refbackground<- testbackground  

ORA.OR = matrix(NA,nrow=length(multiTest),ncol=length(multiRef))
colnames(ORA.OR) = names(multiRef)
rownames(ORA.OR) = names(multiTest)  

ORA.P = matrix(NA,nrow=length(multiTest),ncol=length(multiRef))
colnames(ORA.P) = names(multiRef)
rownames(ORA.P) = names(multiTest)  

for (i in 1:length(multiRef)) {
  for (j in 1:length(multiTest)) {
    result = ORA(multiTest[[j]],multiRef[[i]],testbackground,refbackground);
    ORA.OR[j,i] = result[1];
    ORA.P[j,i] = result[2];
  }
}

ORA.OR<-apply(ORA.OR,2,as.numeric)
dim(ORA.OR)<-dim(ORA.P)  

FDRmat.Array <- matrix(p.adjust( ORA.P,method="fdr"),nrow=nrow( ORA.P),ncol=ncol( ORA.P))
rownames(  FDRmat.Array)=rownames(ORA.P)
colnames(  FDRmat.Array)=colnames(ORA.P)

ORA.P=matrix(as.numeric(ORA.P),nrow=nrow( ORA.P),ncol=ncol( ORA.P))
ORA.OR=matrix(as.numeric(ORA.OR),nrow=nrow( ORA.OR),ncol=ncol( ORA.OR))
rownames(ORA.P) <- rownames(ORA.OR) <- rownames(  FDRmat.Array)
colnames(ORA.P) <- colnames(ORA.OR) <- colnames(  FDRmat.Array)  

dispMat <- ORA.OR

txtMat <-  ORA.OR
#txtMat <- round(-log(ORA.P), 3)
txtMat <- signif(ORA.OR, 2)

txtMat[FDRmat.Array >0.05] <- ""
txtMat[FDRmat.Array <0.05&FDRmat.Array >0.01] <- "*"
txtMat[FDRmat.Array <0.01&FDRmat.Array >0.005] <- "**"
txtMat[FDRmat.Array <0.005] <- "***"

txtMat1 <- signif(ORA.P)
txtMat1[txtMat1<2] <- ""  

textMatrix1 = paste( txtMat1, '\n', txtMat , sep = '')
textMatrix1= matrix(textMatrix1,ncol=ncol(ORA.P),nrow=nrow(ORA.P))  

pdf("CellTypeEnrich_WGCNAMods_110618.pdf", width=6,height=10)
labeledHeatmap(Matrix=dispMat,
               yLabels=rownames(dispMat),
               yColorLabels=TRUE,
               xLabels= colnames(dispMat),
               colors=blueWhiteRed(40),
               textMatrix = textMatrix1,
               cex.lab.x=1.0,
               zlim=c(-0.1,3),
               main="Cell-type enrichment Heatmap")
dev.off()  


#call the function for Cell type erncriment
net <- merge
cleanDat <- p_data
numericMeta <- at100

heatmapScale = "minusLogFDR"
heatmapTitle = "My Network Module Overlaps with 5 Brain Cell Type Marker Reference Lists"

paletteColors = "YlGnBu"

FileBaseName = "MyNetworkModules_FET_to_5brainCellTypeMarkerLists"
refDataDescription = "5BrainCellTypes"

refDataFiles <- c("Mc_Kenzie_v1.csv")


geneListFET(FileBaseName="PSP_PICKnetModules_FET_to_5brainCellTypes.barChart", barOption=TRUE,
            heatmapTitle="PSP_PICK Network Module Overlaps with 5 Brain Cell Type Marker Reference Lists",
            modulesInMemory=TRUE,categorySpeciesCode="hsapiens",  #use network in memory; what species code are the symbols in cleanDat rownames? In case symbol interconversion across species is needed...
            refDataFiles=refDataFiles,speciesCode=c("hsapiens"),refDataDescription="5brainCellTypes")

datExpr0[] <- lapply(datExpr0, as.numeric)
p_data[] <- lapply(p_data, as.numeric)
p_data <- p_data[,-1]

