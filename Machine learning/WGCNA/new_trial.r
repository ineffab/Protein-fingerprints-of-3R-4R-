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

at100 <- read_xlsx("AT100.xlsx")

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
#abline(h = 0.045, col = "red")

#Determine the cluster under the line
#clust = cutreeStatic(sampleTree, cutHeight = 0.045, minSize = 0.045)
#table(clust)

#clust 1 contains samples we want to keep 
#keepSamples = (clust == c(1,2))
#datExpr = datExpr0[keepSamples, ]
#nGenes = ncol(datExpr)
#nSamples = nrow(datExpr)

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

#Rename to module Colors
moduleColors = mergedColors

#Module-trait relationship

nGenes = ncol(datExpr0)
nSamples = nrow(datExpr0)

#moduleColors = moduleColors[moduleColors !='grey']

#recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr0, moduleColors)$eigengenes
#MEs0 = MEs0[MEs0!='grey']
MEs_1 =orderMEs(MEs0)
moduleTraitCor = cor(MEs_1, at100[,-c(1:4)], use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

sizeGrWindow(10,6)

testMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")

dim(testMatrix) = dim(moduleTraitCor)
par(mar = c(9, 8.5, 3, 3))

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(at100[,-c(1:4)]),
               yLabels = names(MEs_1),
               ySymbols = names(MEs_1),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = testMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

#plot correlation with mitotimer 

mito <- read_xlsx("mitotimer_PSP_PICK.xlsx")

mito_1 <- mito[,-1]
rownames(mito_1) <- mito$ribaq

MEs_1_new <- MEs_1[-c(3,4,10:14,18,19),]

mito_cor = cor(MEs_1_new, mito_1, use = "p")
mito_cor_P_value = corPvalueStudent(mito_cor, nrow(MEs_1_new))

sizeGrWindow(10,6)

testMatrix_1 = paste(signif(mito_cor, 2), "\n(",
                   signif(mito_cor_P_value, 1), ")", sep = "")

dim(testMatrix_1) = dim(mito_cor)
par(mar = c(9, 8.5, 3, 3))

labeledHeatmap(Matrix = mito_cor,
               xLabels = names(mito_1),
               yLabels = names(MEs_1_new),
               ySymbols = names(MEs_1_new),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = testMatrix_1,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-mitotimer relationships"))

#linear regression

model_1 <- lm(MEs_1_new$MEturquoise[-c(1,2)] ~ #mito_1$`Ratio 2-6h_D1` 
              #+ mito_1$`Ratio 2-6h_ D2` 
              + mito_1$`Ratio 2-6h_D3`[-c(1,2)])
              #+ mito_1$`Ratio 2-6h_D4`)
summary(model_1)

plot(MEs_1_new$MEturquoise[-c(1,2)], mito_1$`Ratio 2-6h_D3`[-c(1,2)], col = "red")
abline(model_1, col = "red", lwd =2)

  
at100_new <- at100[-c(3,4,10:14,18,19), 11]

model_2 <- lm(at100_new$`Coverage %  AT100`[-9] ~ mito_1$`Ratio 2-6h_D4`[-9])
summary(model_2)

plot(at100_new$`Coverage %  AT100`[-9] ~ mito_1$`Ratio 2-6h_D4`[-9])
abline(model_2)


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







###################################################################################################################
# Cross-species FET modified to optionally adjust for symbol lookup inefficiency/loss
# by Eric Dammer & Divya Nandakumar
#-------------------------------------------------------------------------#
# +2/10/19 improved duplicate removal within and across reference lists
# +2/10/19 added toggles and speciescode for biomaRt lookup as parameters
# +8/11/20 fixed calculations to match divya's (swapped moduleList and categories vars, and removed unique() for totProteomeLength
# +5/08/21 added barOption - Divya style barplots
# +1/26/23 Converted to geneListFET function for Levi Wood ALS/FTD network collaboration
# +4/17/23 Added RColorBrewer palette specification to parameter/variable paletteColors (vector of length=# of marker list file inputs),
#          vector character strings must be one of the sequential (first group) or qualitative (third group) palettes shown by:
#          RColorBrewer::display.brewer.all()
#-------------------------------------------------------------------------#
# revisited to define fly cell types in Seurat 87 lists 2/10/2019
# Analysis for Laura Volpicelli, mouse a-Syn Bilaterally Injected Brain Regions 2/15/2019
# LFQ-MEGA Cell Type analysis performed with this code, with grey proteins added back in to totProteome, allGenes  4/5/2019 #***## (2 lines)
#=========================================================================#
geneListFET <- function(modulesInMemory=TRUE,categoriesFile=NA,categorySpeciesCode=NA,resortListsDecreasingSize=FALSE,barOption=FALSE,adjustFETforLookupEfficiency=FALSE,allowDuplicates=TRUE,
                        refDataFiles=NA,speciesCode=NA,refDataDescription="RefList(s)_not_described",FileBaseName="geneListFET_to_RefList(s)",paletteColors="YlGnBu",
                        heatmapTitle="Heatmap Title (not specified)", heatmapScale="minusLogFDR", verticalCompression=3, rootdir="./", reproduceHistoricCalc=FALSE, env=.GlobalEnv) {
  
  require(WGCNA,quietly=TRUE)
  require(RColorBrewer,quietly=TRUE)
  require(biomaRt,quietly=TRUE)
  
  refDataDir<-outputfigs<-outputtabs<-rootdir
  
  
  if(!modulesInMemory) {  # Read in Categories as list 
    # old format: 2 column .csv with Symbol and "ClusterID" columns
    #  enumeratedListsDF<-read.csv(file=paste0(refDataDir,"/",categoriesFile),header=TRUE)
    #  enumeratedLists<-list()
    #  for(eachList in unique(enumeratedListsDF[,"ClusterID"])) { enumeratedLists[[as.character(eachList)]] <- enumeratedListsDF[which(enumeratedListsDF$ClusterID==eachList),"GeneSymbol"] }
    
    # new format: multicolumn .csv with Symbols or UniqueIDs and each cluster's symbols in a separate column with clusterID as column name (in row 1)
    enumeratedLists <- as.list(read.csv(paste(refDataDir,categoriesFile, sep=""),sep=",", stringsAsFactors = FALSE,header=T,check.names=FALSE)) 
    names(enumeratedLists)
    
    #number of entries with no blanks
    length(unlist(lapply(enumeratedLists,function(x) x[!x==''] )))
    #take out blanks from list
    enumeratedLists <- lapply(enumeratedLists,function(x) x[!x==''] )
    #are there symbols duplicated? (yes, if below result is less than above)
    length(unique(unlist(enumeratedLists)))
    
    
    enumeratedLists<-lapply(enumeratedLists,function(x) as.data.frame(do.call("rbind",strsplit(x,"[|]")))[,1] )
    # leave duplicated symbols within each module/category list -
    #enumeratedLists<-lapply(enumeratedLists,unique)
    
    # leave duplicates in the clusters or modules being checked - they should contribute to overlap/enrichmend with reference lists more than once if duplicated.
    #if(!allowDuplicates) {
    #  while( length(unique(unlist(enumeratedLists))) < length(unlist(enumeratedLists)) ) {
    #    duplicatedvec<-unique(unlist(enumeratedLists)[which(duplicated(unlist(enumeratedLists)))])
    #    #remove duplicates from any marker list
    #    enumeratedLists<-lapply(enumeratedLists,function(x) { remIndices=as.vector(na.omit(match(duplicatedvec,x))); if (length(remIndices)>0) { x[-remIndices] } else { x }; } )
    #  }
    #}
    
  } else {  # use WGCNA modules in memory
    
    # Module lookup table
    nModules<-length(table(net$colors))-1
    modules<-cbind(colnames(as.matrix(table(net$colors))),table(net$colors))
    orderedModules<-cbind(Mnum=paste("M",seq(1:nModules),sep=""),Color=WGCNA::labels2colors(c(1:nModules)))
    modules<-modules[match(as.character(orderedModules[,2]),rownames(modules)),]
    #as.data.frame(cbind(orderedModules,Size=modules))
    
    # Recalculate Consensus Cohort Eigengenes, i.e. eigenproteins and their relatedness order
    MEs<-data.frame()
    MEList = WGCNA::moduleEigengenes(t(cleanDat), colors = net$colors)
    MEs = orderMEs(MEList$eigengenes)
    net$MEs <- MEs
    colnames(MEs)<-gsub("ME","",colnames(MEs)) #let's be consistent in case prefix was added, remove it.
    rownames(MEs)<-rownames(numericMeta)
    if("grey" %in% colnames(MEs)) MEs[,"grey"] <- NULL
    
    # Make list of module member gene product official symbols
    enumeratedLists<-sapply( colnames(MEs),function(x) as.vector(data.frame(do.call("rbind",strsplit(rownames(cleanDat),"[|]")))[,1])[which(net$colors==x)] )
    greyToAddToTotProteome<- as.vector( data.frame(do.call("rbind",strsplit(rownames(cleanDat),"[|]")))[,1])[which(net$colors=="grey")] #***##
  }
  moduleList=enumeratedLists
  
  
  
  pdf(file=paste0(outputfigs,"/",FileBaseName,".Overlap.in.",refDataDescription,".pdf"),height=15,width=24) 
  ############
  
  #***iterating through multiple files (each one a page of output PDF):  
  iter=0
  for (refDataFile in refDataFiles) {
    iter=iter+1
    this.heatmapScale<-heatmapScale
    
    
    refData <- as.list(read.csv(paste(refDataDir,refDataFile, sep=""),sep=",", stringsAsFactors = FALSE,header=T)) 
    names(refData)
    
    
    #number of entries with no blanks
    length(unlist(lapply(refData,function(x) x[!x==''] )))
    #take out blanks from list
    refData <- lapply(refData,function(x) x[!x==''] )
    #are there symbols duplicated? (yes, if below result is less than above)
    length(unique(unlist(refData)))
    
    ##Remove duplicates from all lists if allowDuplicates=FALSE
    #remove duplicated exact symbols within each list regardless:
    if (reproduceHistoricCalc) refData<-lapply(refData,unique)
    refData<-lapply(refData,function(x) as.data.frame(do.call("rbind",strsplit(x,"[|]")))[,1] )
    if (!reproduceHistoricCalc) refData<-lapply(refData,unique)
    
    if(!allowDuplicates) {
      while( length(unique(unlist(refData))) < length(unlist(refData)) ) {
        duplicatedvec<-unique(unlist(refData)[which(duplicated(unlist(refData)))])
        #remove duplicates from any marker list
        refData<-lapply(refData,function(x) { remIndices=as.vector(na.omit(match(duplicatedvec,x))); if (length(remIndices)>0) { x[-remIndices] } else { x }; } )
      }
      duplicateHandling="DuplicatesREMOVED"
    } else {
      duplicateHandling="DuplicatesALLOWED"
    }
    length(unlist(refData))
    unlist(refData)[which(duplicated(unlist(refData)))]
    refDataMouse<-refData
    
    groupvec<-placeholders<-vector()
    for (i in 1:length(names(refData))) {
      placeholders=c(placeholders,rep(i,length(refData[[i]])))
      groupvec=c(groupvec,rep(names(refData)[i],length(refData[[i]])))
    }
    categoriesData<-data.frame(UniqueID=unlist(refData),Color=labels2colors(placeholders), Annot=groupvec,Mnum=paste0("M",placeholders)) #,row.names=unlist(refData)) #will not work if allowDuplicates==TRUE
    #^holding refData all lists' items, not modules gene symbols
    
    categoriesNameMatcher<-unique(categoriesData[,2:4])
    rownames(categoriesNameMatcher)<-NULL
    
    
    if(!categorySpeciesCode==speciesCode[iter]) {
      cat(paste0("Converting ",speciesCode[iter]," to ",categorySpeciesCode," for lists in ",refDataFile," ... "))
      
      this.heatmapTitle=paste0(heatmapTitle," in ",categorySpeciesCode," homologs")
      #library(biomaRt)
      
      #human = useEnsembl("genes", dataset = "hsapiens_gene_ensembl", host="https://dec2021.archive.ensembl.org")  #ver=105 equivalent to dec2021
      #mouse = useEnsembl("genes", dataset = "mmusculus_gene_ensembl", host="https://dec2021.archive.ensembl.org") 
      
      category.species = useEnsembl("genes", dataset=paste0(categorySpeciesCode,"_gene_ensembl"), host="https://dec2021.archive.ensembl.org")  #ver=105 equivalent to dec2021  ; old code: #useMart("ensembl",dataset=paste0(categorySpeciesCode,"_gene_ensembl"))
      other = useEnsembl("genes", dataset=paste0(speciesCode[iter],"_gene_ensembl"), host="https://dec2021.archive.ensembl.org")   # old code: #useMart("ensembl",dataset=paste0(speciesCode[iter],"_gene_ensembl"))
      
      #category species to other species conversion (first column is other species)
      if(speciesCode[iter]=="hsapiens") {
        genelist.mouseConv<-getLDS(attributes=c("hgnc_symbol"), filters="hgnc_symbol", values=categoriesData$UniqueID, mart=other, attributesL="external_gene_name",martL = category.species)
        categoriesData$BiomartFlySymbol <- genelist.mouseConv[match(categoriesData$UniqueID,genelist.mouseConv$HGNC.symbol),"Gene.name"]
      } else {
        if(categorySpeciesCode=="hsapiens") {
          genelist.mouseConv<-getLDS(attributes=c("external_gene_name"), filters="external_gene_name", values=categoriesData$UniqueID, mart=other, attributesL="hgnc_symbol",martL = category.species)
          categoriesData$BiomartFlySymbol <- genelist.mouseConv[match(categoriesData$UniqueID,genelist.mouseConv$Gene.name),"HGNC.symbol"]
        } else {
          genelist.mouseConv<-getLDS(attributes=c("external_gene_name"), filters="external_gene_name", values=categoriesData$UniqueID, mart=other, attributesL="external_gene_name",martL = category.species)
          categoriesData$BiomartFlySymbol <- genelist.mouseConv[match(categoriesData$UniqueID,genelist.mouseConv$Gene.name),"Gene.name.1"]
        }
      }
    } else { this.heatmapTitle=heatmapTitle }
    
    categoriesData.original<-categoriesData
    
    categoriesData.reducedFly<-na.omit(categoriesData)
    categoriesData.reducedFly$MouseID.original<-categoriesData.reducedFly$UniqueID
    if(!categorySpeciesCode==speciesCode[iter]) categoriesData.reducedFly$UniqueID<-categoriesData.reducedFly$BiomartFlySymbol
    categoriesData.reducedFly<-categoriesData.reducedFly[,-5]
    
    refDataFly<-list()
    for (i in unique(categoriesData.reducedFly$Annot)) {
      refDataFly[[i]]<-unique(categoriesData.reducedFly$UniqueID[which(categoriesData.reducedFly$Annot==i)])
    }
    refDataFly.original<-refDataFly
    
    #remove within-list duplicates from any marker list (some homologs map to multiple reference species unique genes)
    refDataFly<-lapply(refDataFly,unique)
    #final check
    length(unlist(refDataFly))
    unlist(refDataFly)[which(duplicated(unlist(refDataFly)))] #any duplicates across lists allowed
    length(unlist(refDataFly))==length(unique(unlist(refDataFly))) #false if duplicates across lists.
    
    refData<-refDataFly
    
    #make data frame of all markers
    mouseSymbolVec<-groupvec<-placeholders<-vector()
    for (i in 1:length(names(refData))) {
      placeholders=c(placeholders,rep(i,length(refData[[i]])))
      groupvec=c(groupvec,rep(names(refData)[i],length(refData[[i]])))
      mouseSymbolVec=c(groupvec,rep(names(refData)[i],length(refData[[i]])))
    }
    categoriesData<-data.frame(UniqueID=unlist(refData),Color=labels2colors(placeholders), Annot=groupvec,Mnum=paste0("M",placeholders))  #,row.names=unlist(refData))
    categoriesData$MouseSymbol=categoriesData.reducedFly$MouseID.original[match(categoriesData$UniqueID,categoriesData.reducedFly$UniqueID)]
    
    categoriesNameMatcher<-unique(categoriesData[,2:4])
    rownames(categoriesNameMatcher)<-NULL
    
    
    
    if(modulesInMemory) {
      allGenes<- c(unlist(moduleList), greyToAddToTotProteome)  #unique() here decreases significance, totProteomeLength; Not in original code.
    } else {
      allGenes<- unlist(enumeratedLists)  #ANOVAout$Symbol  #here the background is all measured proteins, #categoriesData$BiomartMouseSymbol
    }
    allGenesNetwork <- as.matrix(allGenes,stringsAsFactors = FALSE) 
    
    categories <- list()
    categoryNames=names(refData) #reference list names
    for (i in 1:length(categoryNames)) {
      element<-categoryNames[i]
      categories[[element]] <- categoriesData$UniqueID[which(categoriesData$Annot==categoryNames[i])]  #categoriesData$BiomartMouseSymbol[which(categoriesData$colors==modcolors[i])]
    }
    
    ##+#+#+#+#+#+#+#+#+#+#+#+#+
    # Final Data Cleaning
    
    nModules <- length(names(moduleList))
    nCategories <- length(names(categories))
    
    for (a in 1:nCategories) {
      categories[[a]] <- unique(categories[[a]][categories[[a]] != ""])
      categories[[a]] <- categories[[a]][!is.na(categories[[a]])]
    }
    for (b in 1:nModules) {
      moduleList[[b]] <- unique(moduleList[[b]][moduleList[[b]] != ""])
    }
    
    if(resortListsDecreasingSize) {
      categories <- categories[order(sapply(categories,length),decreasing=T)]
      #only sort lists if 'moduleList' is not a list of WGCNA modules (keep them in relatedness order if they are)
      if (!modulesInMemory) moduleList <- moduleList[order(sapply(moduleList,length),decreasing=T)]
    } #if FALSE, do not resort lists -- we have them in a precise order already
    
    
    allGenes_cleaned <- na.omit(allGenesNetwork)
    totProteomeLength <- length(allGenes_cleaned)
    if(max(sapply(refData,length))>totProteomeLength) {
      cat (paste0("One of your reference data lists is larger than the background from the categories (WGCNA or specified categories file--all symbols)!\nUsing the bigger number for Fisher Exact would change all stats. Skipping ",refDataFile,".\n\n"))
      next
    }
    
    ### Fisher's Exact Test
    
    cat(paste0("Performing FET for lists [",iter,"] now.\n"))
    
    
    ###swap cell type lists to categories and module members to moduleList (they are backwards up till here)
    #categories.bak<-categories
    #categories<-moduleList
    #moduleList<-categories.bak
    #nModules <- length(names(moduleList))
    #nCategories <- length(names(categories))
    
    FTpVal <- matrix(nrow = nModules, ncol = nCategories)
    categoryOverlap <- matrix(nrow = nModules, ncol = nCategories) 
    numCategoryHitsInDataset <- numCategoryHitsInDataset.UNADJ <- matrix(nrow = nModules, ncol = nCategories) 
    CategoryHitsInDataset <- list()
    hitLists<-matrix(NA,nrow=nModules,ncol=nCategories) #use a matrix of collapsed (";") gene list strings
    ADJRedundancyAfterLookup=1 #length(unlist(refDataFly.original))/length(unlist(refDataFly)) #bigger than 1
    ADJforCrossSpeciesLookupFailure=nrow(categoriesData)/nrow(categoriesData.original) #less than 1
    totProteomeLength.ADJ <- as.integer(totProteomeLength*ADJforCrossSpeciesLookupFailure*ADJRedundancyAfterLookup)
    RefDataElements<-Categories1<-vector()
    
    for (i in 1:nModules){
      sampleSize <- length(moduleList[[i]])
      RefDataElements=c(RefDataElements,names(moduleList)[i])
      for (j in 1:nCategories){
        if(i==1) { Categories1=c(Categories1,names(categories)[j]) }
        #CategoryHitsInProteome <- categories[[j]] ## If using all of the markers and not just markers in proteome
        CategoryHitsInProteome <- intersect(categories[[j]],allGenesNetwork[,1])
        if (!adjustFETforLookupEfficiency) {
          ##Unadjusted calculations:
          numCategoryHitsInProteome <- length(CategoryHitsInProteome) 
          numNonCategoryHitsInProteome <- totProteomeLength - numCategoryHitsInProteome
          overlapGenes <- intersect(moduleList[[i]],CategoryHitsInProteome)
          numOverlap <- length(overlapGenes)
          otherCategories <- sampleSize - numOverlap
          notInModule <- numCategoryHitsInProteome - numOverlap
          notInMod_otherCategories <- totProteomeLength - numCategoryHitsInProteome - otherCategories
        } else {  ##allGenesNetwork has different species Symbols, and categories are also from that full Symbol List, so adjust for comparison to interconverted list overlap
          #&& adjustments noted
          numCategoryHitsInProteome <- as.integer(length(CategoryHitsInProteome)*(ADJforCrossSpeciesLookupFailure*ADJRedundancyAfterLookup)) #&& adjusted down for lookup inefficiency
          numCategoryHitsInProteome.UNADJ <- length(CategoryHitsInProteome)
          numNonCategoryHitsInProteome <- totProteomeLength.ADJ - numCategoryHitsInProteome #first term is adjusted
          overlapGenes <- intersect(moduleList[[i]],CategoryHitsInProteome) #does not need adjustment, both subject to lower lookup efficiency
          
          numOverlap <- length(overlapGenes)
          otherCategories <- sampleSize - numOverlap
          notInModule <- numCategoryHitsInProteome - numOverlap #&&using down-adjusted number numCategoryHitsInProteome for lookup efficiency
          notInMod_otherCategories <- totProteomeLength - numCategoryHitsInProteome - otherCategories #&&first term not adjusted down because this is non-overlap so lookup inefficiency does not apply
          #&& but second term is adjusted because it it the hits subject to lookup efficiency
        }
        hitLists[i,j]<-paste(overlapGenes,collapse=";")
        contingency <- matrix(c(numOverlap,otherCategories,notInModule,notInMod_otherCategories),nrow=2,ncol=2,dimnames=list(c("GenesHit","GenesNotHit"),c("withinCategory","inProteome")))
        #debugging:     if(i==6 & j==3) cat(contingency)
        FT <- fisher.test(contingency,alternative="greater") #variable with presumed explanatory effect should be the row definitions, if known. (can transpose, but no effect on outcome p values)
        FTpVal[i,j] <- FT$p.value
        categoryOverlap[i,j] <- numOverlap
        numCategoryHitsInDataset[i,j] <- numCategoryHitsInProteome
        numCategoryHitsInDataset.UNADJ[i,j] <- if(adjustFETforLookupEfficiency) { numCategoryHitsInProteome.UNADJ } else { numCategoryHitsInProteome }
        if (i==1){
          CategoryHitsInDataset[[j]] <- array(CategoryHitsInProteome)
        }		
      }
    }
    
    
    #moduleList<-categories
    #categories<-categories.bak
    #nModules <- length(names(moduleList))
    #nCategories <- length(names(categories))
    
    
    rownames(FTpVal) <- RefDataElements
    colnames(FTpVal) <- Categories1
    rownames(categoryOverlap) <- RefDataElements
    colnames(categoryOverlap) <- Categories1
    colnames(numCategoryHitsInDataset) <- Categories1
    rownames(numCategoryHitsInDataset) <- RefDataElements
    names(CategoryHitsInDataset) <- Categories1
    rownames(hitLists) <- RefDataElements
    colnames(hitLists) <- Categories1
    
    
    #### Format Data for Plotting ########
    
    NegLogUncorr <- -log10(FTpVal)
    rownames(NegLogUncorr) <- rownames(FTpVal)
    colnames(NegLogUncorr) <- colnames(FTpVal)
    NegLogUncorr <- as.matrix(NegLogUncorr)
    
    nCategories = ncol(FTpVal)
    nModules = nrow(FTpVal)
    
    FisherspVal <- unlist(FTpVal)
    adjustedPVal <- p.adjust(FisherspVal, method = "fdr", n=length(FisherspVal))
    adjustedPval <- matrix(adjustedPVal,nrow=nModules,ncol=nCategories)
    rownames(adjustedPval) <- rownames(FTpVal)
    colnames(adjustedPval) <- colnames(FTpVal)
    NegLogCorr <- -log10(adjustedPval)
    
    ## Transpose above stats and hits matrices
    categoryOverlap<-t(categoryOverlap)
    numCategoryHitsInDataset<-t(numCategoryHitsInDataset)
    numCategoryHitsInDataset.UNADJ<-t(numCategoryHitsInDataset.UNADJ)
    CategoryHitsInDataset<-t(CategoryHitsInDataset)
    hitLists<-t(hitLists)
    NegLogUncorr<-t(NegLogUncorr)
    NegLogCorr<-t(NegLogCorr)
    adjustedPval<-t(adjustedPval)
    FTpVal<-t(FTpVal)
    
    
    ##Make sure colors are in correct (WGCNA) order before changing to numbered modules!
    if(modulesInMemory) {
      orderedLabels<-cbind(paste("M",seq(1:nCategories),sep=""),labels2colors(c(1:nCategories)))
    } else {
      orderedLabels<- cbind(paste("M",seq(1:nModules),sep=""),labels2colors(c(1:nModules))) #these go from M1 to M(# of reference lists)
    }
    
    #if you want the modules in order of relatedness from the module relatedness dendrogram:
    if(!modulesInMemory) {
      orderedLabelsByRelatedness<- orderedLabels #(this is chronol. order)
      if (!length(na.omit(match(orderedLabelsByRelatedness[,2],RefDataElements)))==nrow(orderedLabelsByRelatedness)) orderedLabelsByRelatedness[,2]<- RefDataElements; dummyColors=orderedLabelsByRelatedness[,2]; # our category/cluster names on categoriesFile row 1 are not WGCNA colors.
      NegLogUncorr<-NegLogUncorr[,match(orderedLabelsByRelatedness[,2],colnames(NegLogUncorr))]
      NegLogCorr<-NegLogCorr[,match(orderedLabelsByRelatedness[,2],colnames(NegLogCorr))]
      adjustedPval<-adjustedPval[,match(orderedLabelsByRelatedness[,2],colnames(adjustedPval))]
      FTpVal<-FTpVal[,match(orderedLabelsByRelatedness[,2],colnames(FTpVal))]
    } else {
      orderedLabelsByRelatedness<- cbind( orderedLabels[ match(gsub("ME","",colnames(MEs)),orderedLabels[,2]) ,1] ,gsub("ME","",colnames(MEs)) )
      
      NegLogUncorr<-NegLogUncorr[,match(orderedLabelsByRelatedness[,2],colnames(NegLogUncorr))]
      NegLogCorr<-NegLogCorr[,match(orderedLabelsByRelatedness[,2],colnames(NegLogCorr))]
      adjustedPval<-adjustedPval[,match(orderedLabelsByRelatedness[,2],colnames(adjustedPval))]
      FTpVal<-FTpVal[,match(orderedLabelsByRelatedness[,2],colnames(FTpVal))]
    }
    xlabels <- orderedLabelsByRelatedness[,1]
    
    
    
    
    
    ### Write p Values to a table/file
    #rownames(hitLists)<-categoriesNameMatcher$Annot[match(rownames(categoryOverlap),categoriesNameMatcher$Annot)]
    #rownames(FTpVal)<-categoriesNameMatcher$Annot[match(rownames(FTpVal),categoriesNameMatcher$Annot)]
    #rownames(adjustedPval)<-categoriesNameMatcher$Annot[match(rownames(adjustedPval),categoriesNameMatcher$Annot)]
    #rownames(categoryOverlap)<-categoriesNameMatcher$Annot[match(rownames(categoryOverlap),categoriesNameMatcher$Annot)]
    
    outputData <- rbind("FET pValue", FTpVal,"FDR corrected",adjustedPval,"Overlap",categoryOverlap,"CategoryHitsInDataSet(ADJ)",numCategoryHitsInDataset,"CategoryHitsInDataSet(Unadj)",numCategoryHitsInDataset.UNADJ,"OverlappedGeneLists",hitLists)
    write.csv(outputData,file = paste0(outputtabs,"/",FileBaseName,".Overlap.in.",refDataFile,"-",duplicateHandling,"-hitListStats.csv"))
    
    #auto-check if all FET (BH) calculations = 1, then switch to p value visualization
    if(mean(rowMeans(adjustedPval,na.rm=T),na.rm=T)==1) { this.heatmapScale<-"p.unadj"; addText="-No FDR values lower than 100%"; } else { this.heatmapScale<-heatmapScale; addText=""; }
    
    ## Use the text function with the FDR filter in labeledHeatmap to add asterisks, e.g. * 
    txtMat <- adjustedPval
    txtMat[adjustedPval>=0.05] <- ""
    txtMat[adjustedPval <0.05&adjustedPval >0.01] <- "*"
    txtMat[adjustedPval <0.01&adjustedPval >0.005] <- "**"
    txtMat[adjustedPval <0.005] <- "***"
    
    txtMat1 <- signif(adjustedPval,2)
    txtMat1[adjustedPval>0.25] <- ""
    
    
    textMatrix1 = paste( txtMat1, txtMat , sep = ' ');
    textMatrix1= matrix(textMatrix1,ncol=ncol(adjustedPval),nrow=nrow(adjustedPval))
    
    #for textMatrix of p.unadj
    txtMat <- FTpVal
    txtMat[FTpVal>=0.05] <- ""
    txtMat[FTpVal <0.05&FTpVal >0.01] <- "*"
    txtMat[FTpVal <0.01&FTpVal >0.005] <- "**"
    txtMat[FTpVal <0.005] <- "***"
    
    txtMat.p.unadj <- signif(FTpVal,2)
    txtMat.p.unadj[FTpVal>0.25] <- ""
    
    textMatrix.p.unadj = paste( txtMat.p.unadj, txtMat , sep = ' ');
    textMatrix.p.unadj= matrix(textMatrix.p.unadj,ncol=ncol(FTpVal),nrow=nrow(FTpVal))
    
    
    ## Plotting
    if(!barOption) {
      par(mfrow=c(verticalCompression,1))
      par( mar = c(9.5, 10, 4.5, 2) ) #bottom, left, top, right #text lines
      
      if(exists("colvec")) suppressWarnings(rm(colvec))
      
      #RColorBrewer::display.brewer.all()
      if(iter>length(paletteColors)) {
        cat(paste0("  - paletteColors specified being recycled for additional heatmaps for inputs after #",length(paletteColors),"\n"))
        paletteColors<-c(paletteColors,rep(paletteColors,ceiling(length(refDataFiles)/length(paletteColors))))
      }
      if(!paletteColors[iter] %in% c("YlOrRd","YlOrBr","YlGnBu","YlGn","Reds","RdPu","Purples","PuRd","PuBuGn","PuBu","OrRd","Oranges","Greys","Greens","GnBu","BuPu","BnGn","Blues",
                                     "Spectral","RdYlGn","RdYlBu","RdGy","RdBu","PuOr","PRGn","PiYG","BrBG")) {
        cat(paste0("  - paletteColors specified as '",paletteColors[iter],"' is not in RColorBrewer::display.brewer.all() groups 1 or 3.\n    Using palette 'YlGnBu' (yellow, green, blue)...\n"))
        paletteColors[iter]="YlGnBu"
      }
      if(paletteColors[iter] %in% c("YlOrRd","YlOrBr","YlGnBu","YlGn","Reds","RdPu","Purples","PuRd","PuBuGn","PuBu","OrRd","Oranges","Greys","Greens","GnBu","BuPu","BnGn","Blues")) {
        paletteLength=9
        outOfParkColor=brewer.pal(paletteLength,paletteColors[iter])[paletteLength]
        colvec<- brewer.pal(paletteLength,paletteColors[iter])[1:6]
      } else {
        paletteLength=11
        # pure purple for outOfPark maximum scale color, deprecated.
        # if(paletteColors[iter] %in% c("Spectral","RdYlGn","RdYlBu","RdGy","RdBu","PuOr","BrBG")) { outOfParkColor="#A020F0" } else { outOfParkColor="darkviolet" }
        outOfParkColor=brewer.pal(paletteLength,paletteColors[iter])[paletteLength]
        #		   if(as.boolean(revPalette[iter])) {
        colvec<- rev(brewer.pal(paletteLength,paletteColors[iter])[1:6])  #rev so we take the left side of the palette color swatches
        #		   } else {
        #		      colvec<- brewer.pal(paletteLength,paletteColors[iter])[6:11]  #no rev if we want to take the right half of the palette color swatches
        #		   }
      }
      
      colvecRamped1<- vector()
      for (k in 1:(length(colvec)-1)) {
        gradations <- if (k<4) { 6 } else { 25 }
        temp<-colorRampPalette(c(colvec[k],colvec[k+1]))
        colvecRamped1<-c(colvecRamped1, temp(gradations))
      }
      
      temp2<-colorRampPalette(c(colvecRamped1[length(colvecRamped1)], outOfParkColor)) ## grade to outOfParkColor at top of scale
      colvecRamped1<-c(colvecRamped1, temp2(gradations))
      
      colvecRamped1<-c("#FFFFFF",colvecRamped1)  ## grade to white at bottom of scale
      
      
      if (modulesInMemory) { categoryColorSymbols=paste0("ME",names(moduleList)) } else { if(!length(na.omit(match(orderedLabelsByRelatedness[,2],rownames(NegLogUncorr))))==nrow(orderedLabelsByRelatedness)) { categoryColorSymbols=dummyColors } else { categoryColorSymbols=names(moduleList) } }
      xSymbolsText= ifelse ( rep(modulesInMemory,length(names(moduleList))), paste0(names(moduleList)," ",orderedModules[match(names(moduleList),orderedModules[,2]),1]), names(moduleList) )
      if (this.heatmapScale=="p.unadj") {
        labeledHeatmap(Matrix = FTpVal,
                       yLabels = names(categories), #refData list elements, ordered by size if that option was on
                       xLabels = categoryColorSymbols,
                       xLabelsAngle = 90,
                       xSymbols = xSymbolsText,
                       xColorLabels=FALSE,
                       colors = rev(colvecRamped1),
                       textMatrix =  textMatrix.p.unadj,
                       setStdMargins = FALSE,
                       cex.text = 0.6,
                       cex.lab.y = 0.7,
                       verticalSeparator.x=c(rep(c(1:length(names(moduleList))),nrow(orderedLabelsByRelatedness))),
                       verticalSeparator.col = 1,
                       verticalSeparator.lty = 1,
                       verticalSeparator.lwd = 1,
                       verticalSeparator.ext = 0,
                       horizontalSeparator.y=c(rep(c(1:length(names(categories))),nrow(orderedLabelsByRelatedness))),
                       horizontalSeparator.col = 1,
                       horizontalSeparator.lty = 1,
                       horizontalSeparator.lwd = 1,
                       horizontalSeparator.ext = 0,
                       zlim = c(min(FTpVal),1),
                       main = paste0("Enrichment of ",this.heatmapTitle,"\nof ",refDataFile," Marker Lists by Gene Symbol (",duplicateHandling,")\nHeatmap: Fisher Exact p value, Uncorrected\n (p-values shown",addText,")"),
                       cex.main=0.8)
      }
      
      if (this.heatmapScale=="minusLogFDR") {
        labeledHeatmap(Matrix = NegLogCorr,
                       yLabels = names(categories), #refData list elements, ordered by size if that option was on
                       xLabels = categoryColorSymbols,
                       xLabelsAngle = 90,
                       xSymbols = xSymbolsText,
                       xColorLabels=FALSE,
                       colors = colvecRamped1,
                       textMatrix = textMatrix1, #signif(adjustedPval, 2),
                       setStdMargins = FALSE,
                       cex.text = 0.6,
                       cex.lab.y = 0.7,
                       verticalSeparator.x=c(rep(c(1:length(names(moduleList))),nrow(orderedLabelsByRelatedness))),
                       verticalSeparator.col = 1,
                       verticalSeparator.lty = 1,
                       verticalSeparator.lwd = 1,
                       verticalSeparator.ext = 0,
                       horizontalSeparator.y=c(rep(c(1:length(names(categories))),nrow(orderedLabelsByRelatedness))),
                       horizontalSeparator.col = 1,
                       horizontalSeparator.lty = 1,
                       horizontalSeparator.lwd = 1,
                       horizontalSeparator.ext = 0,
                       zlim = c(0,max(NegLogCorr,na.rm=TRUE)),
                       main = paste0("Enrichment of ",this.heatmapTitle,"\nof ",refDataFile," Marker Lists by Gene Symbol (",duplicateHandling,")\nHeatmap: -log(p), BH Corrected\n (Corrected p-values, FDR, shown)"), #*** Uncorrected\n (p-values shown)"),
                       cex.main=0.8)
      }
    } else {  #if barOption==TRUE:  PLOT BAR PLOTS FOR EACH REFERENCE LIST
      par(mfrow=c(verticalCompression,2))
      par(mar=c(15,7,4,1))
      
      moduleColors= if (modulesInMemory) { names(moduleList) } else { "bisque4" }  # if (modulesInMemory), expect colors for names(moduleList)
      xSymbolsText= ifelse ( rep(modulesInMemory,length(names(moduleList))), paste0(names(moduleList)," ",orderedModules[match(names(moduleList),orderedModules[,2]),1]), names(moduleList) )
      if (this.heatmapScale=="p.unadj") {
        
        for( i in 1:nrow(NegLogUncorr)) {
          plotting <- NegLogUncorr[i,]
          cellType <- rownames(NegLogUncorr)[i]
          barplot(plotting,main = cellType, ylab="",cex.names=1.1, width=1.5,las=2,cex.main=2, legend.text=F,col=moduleColors,names.arg=xSymbolsText)
          mtext(side=2, line=3, "-log(pValue)\n(Uncorrected)", col="black", font=1, cex=1.5)
          abline(h=1.3,col="red")
        }
      }
      
      if (this.heatmapScale=="minusLogFDR") {
        for( i in 1:nrow(NegLogCorr)) {
          plotting <- NegLogCorr[i,]
          cellType <- rownames(NegLogCorr)[i]
          barplot(plotting,main = cellType, ylab="",cex.names=1.1, width=1.5,las=2,cex.main=2, legend.text=F,col=moduleColors,names.arg=xSymbolsText)
          mtext(side=2, line=3, "-log(FDR)\n(Benjamini-Hochberg Correction)", col="black", font=1, cex=1.5)
          abline(h=1.3,col="red")
        }
      }
      
    } #ends if(!barOption)
    
    
    #+#+#+#+#+#+#+#+#+#+#+#+#+
  } #end for(refDataFile ...
  dev.off()
}
 

net <- merge
cleanDat <- p_data
numericMeta <- at100

heatmapScale = "minusLogFDR"
heatmapTitle = "My Network Module Overlaps with 5 Brain Cell Type Marker Reference Lists"

paletteColors = "YlGnBu"

FileBaseName = "MyNetworkModules_FET_to_5brainCellTypeMarkerLists"
refDataDescription = "5BrainCellTypes"

refDataFiles <- c("MyGene-Human-SharmaZhangUnion.csv", "MyGene-Mouse-SharmaZhangUnion.csv")


geneListFET(FileBaseName="PSP_PICKnetModules_FET_to_5brainCellTypes.barChart", barOption=TRUE,
            #heatmapTitle="PSP_PICK Network Module Overlaps with 5 Brain Cell Type Marker Reference Lists",
            modulesInMemory=TRUE,categorySpeciesCode="hsapiens",  #use network in memory; what species code are the symbols in cleanDat rownames? In case symbol interconversion across species is needed...
            refDataFiles=refDataFiles,speciesCode=c("hsapiens","mmusculus"),refDataDescription="5brainCellTypes")
  
datExpr0[] <- lapply(datExpr0, as.numeric)
p_data[] <- lapply(p_data, as.numeric)
p_data <- p_data[,-1]
