library(biomaRt)
library(edgeR)
library(KEGGREST)
library(WGCNA)
library(Boruta)
library(e1071)
library(limma)
library(caret)
library(glmnet)
library(Rcpp)
library(STRINGdb)
library(rstudioapi)
library(hexbin)
library(stringr)

setwd(dirname(getActiveDocumentContext()$path))  

hgPwayGenes = read.csv("../common_files/kegg_genes.csv", header = FALSE, stringsAsFactors = FALSE)
hgPwayGenes2 = hgPwayGenes[,2]
names(hgPwayGenes2) = hgPwayGenes[,1]

cellLineArray =
    read.csv("../common_files/sanger1018_brainarray_ensemblgene_rma.txt", sep="\t")

cellLineWES = read.csv("WES_variants.csv")

cellLineWES_BRAF = cellLineWES[cellLineWES[,"Gene"]=="BRAF", ]

csmBRAFCls = unique(cellLineWES_BRAF[grep("p.V600E", cellLineWES_BRAF[,"AA"]), c("SAMPLE", "COSMIC_ID") ])
csmBRAFWES = cellLineWES[cellLineWES[,1] %in% csmBRAFCls[,1], ]

drugDat = read.csv("v17.3_fitted_dose_response.csv")
drugDatNonBraf = drugDat[!(drugDat[,"COSMIC_ID"] %in% csmBRAFCls[,2]), ]

##
wantBrafDs = c("PLX-4720", "Dabrafenib")

plxthresh =
  quantile(drugDat[grep("PLX-4720", drugDat[,"DRUG_NAME"]), "AUC"], 0.05)

dabrafthresh =
  quantile(drugDat[grep("Dabrafenib", drugDat[,"DRUG_NAME"]), "AUC"], 0.05)

details = read.csv("../common_files/Cell_Lines_Details.csv")
detailsMat = details[match(csmBRAFCls[,2], details[,2]), c("Sample.Name", "COSMIC.identifier", "GDSC.Tissue.descriptor.2")]

melCls = detailsMat[detailsMat[,3] == "melanoma", ]
wantDrugs = drugDat[drugDat[, "COSMIC_ID"] %in% melCls[,2], ]

melCls = detailsMat[detailsMat[,3] == "melanoma", ]

wantVEM = drugDat[grep("PLX-4720", drugDat[,"DRUG_NAME"]), ]

wantVEM = wantVEM[wantVEM[,"CELL_LINE_NAME"] %in% colnames(expr),]

wantVEM2 = wantVEM[!duplicated(wantVEM[,"CELL_LINE_NAME"]) ,]

resist = wantVEM2[wantVEM2[,"COSMIC_ID"] %in% resistCL,]

sensCL = c()
resistCL = c() 

for (cl in melCls[,2]){
  currDrugs = wantDrugs[wantDrugs[, "COSMIC_ID"] == cl, ]

  if( "PLX-4720" %in% currDrugs[,"DRUG_NAME"] ){
    if(mean(currDrugs[currDrugs[,"DRUG_NAME"] == "PLX-4720", "AUC"]) < plxthresh){            
      sensCL = c(sensCL, cl)
    } else {
      if (length(currDrugs[currDrugs[, "DRUG_NAME"] == "Dabrafenib", "AUC"])>0 ) {
        if( !(mean(currDrugs[currDrugs[,"DRUG_NAME"] == "Dabrafenib", "AUC"])==Inf) ){
          if(mean(currDrugs[currDrugs[,"DRUG_NAME"] == "Dabrafenib", "AUC"]) < dabrafthresh){
            sensCL = c(sensCL, cl)
            next
          }
        }
        
        
        resistCL = c(resistCL, cl)
        next
      }
      resistCL = c(resistCL, cl)
      
    }
  }
  
}

allCLs = c(sensCL, resistCL)
allCLscms = paste("X", allCLs, sep="")

noExprCl = allCLscms[!(allCLscms %in% colnames(cellLineArray))]
allCLscms = allCLscms[!(allCLscms %in% noExprCl)]

rownames(cellLineArray) = cellLineArray[,1] 
cellLineArray = cellLineArray[, -1]

melCls[,2] = paste("X", melCls[,2], sep="")

cellLineArray = cellLineArray[, match(allCLscms, colnames(cellLineArray)) ]
colnames(cellLineArray) = melCls[match(colnames(cellLineArray), melCls[,2]), 1]

expr = cellLineArray

## change ensmble rownames to HG
rownames(expr) = sub("\\..*", "", rownames(expr)) 
ensembl=useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")

ensmblHG =
  getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"),
        filters = "ensembl_gene_id", values = rownames(expr), mart = ensembl)

if(any(ensmblHG[,1]=="")){
  ensmblHG = ensmblHG[-which(ensmblHG[,1]==""), ] 
}

ensmblHGunq = unique(ensmblHG[,2])
unqENSmat = ensmblHG[match(ensmblHGunq, ensmblHG[,2]), ]

hgncExpr = expr[match(unqENSmat[,2], rownames(expr)), ]
rownames(hgncExpr) = unqENSmat[,1]

expr = hgncExpr
exprMatName = rownames(expr)

sensCL2 = as.character(melCls[melCls[,2] %in% paste("X",sensCL,sep=""), 1])

sensCLInd = colnames(expr) %in% sensCL2
restCLInd = !(sensCLInd)

## remove genes in pways that have no expression data
hgPwayGenes2 = hgPwayGenes2[hgPwayGenes2 %in% exprMatName]
pwayNames = names(hgPwayGenes2)
unqPways = unique(pwayNames)

######################## A1. test pathway significance ##########################
################################################################################
    
genes = c()

activityMat = matrix(0, nrow = length(unqPways), ncol = ncol(expr))
rownames(activityMat) = unqPways

for (currPway in unqPways ){
  
    currGenes = hgPwayGenes2[ which(pwayNames==currPway) ]
    if(length(currGenes)<10) {next}
    
    currExpr = expr[match(currGenes, exprMatName), ]
    currNormExpr = t(apply(currExpr, 1, scale))

    tResult = apply(currNormExpr, 1, function(x) t.test(x[sensCLInd], x[restCLInd]))
    tPval = unlist(lapply(tResult, function(x) x$p.value) )
    tStat = unlist(lapply(tResult, function(x) x$statistic) )
    
    tResult_simple = unlist(lapply(tResult, function(x) x$statistic) )
    
    if ( mean(tResult_simple)<0 ){
        inds = sort(tResult_simple, index.return=T)$ix
        tResult_simple = tResult_simple[inds]
        tPval = tPval[inds]
    } else {
        inds = sort(tResult_simple, index.return=T, decreasing=T)$ix
        tResult_simple = tResult_simple[inds]
        tPval = tPval[inds]
    }
    
    if(length(tResult_simple)==0){next}
    
    orderedNames = sub("\\..*","",names(tResult_simple))
    currNormExpr = currNormExpr[orderedNames,]

    actScore = tResult_simple[1]
    actScoreNew = actScore

    goodInds = 1
    
    for(ind in 2:length(tResult_simple) ) {

        if(length(tResult_simple)==1){break}
            currActVec = colSums(currNormExpr[c(goodInds, ind), ])/sqrt( length(goodInds)+1 )    
            actScoreNew = t.test(currActVec[sensCLInd], currActVec[restCLInd])$statistic

        if (tResult_simple[1] < 0){
            actScoreNewComp = -1*actScoreNew
            actScoreComp = -1*actScore
        } else {
            actScoreNewComp = actScoreNew
            actScoreComp = actScore
        }
                    
        if (actScoreNewComp > actScoreComp){
            goodInds = c(goodInds, ind)
            actScore = actScoreNew
        } else {
            break
        }
    }    

    ## need to prevent summing columns
    if (length(goodInds) == 1){
        next
    } else {
        activity = colSums(currNormExpr[goodInds,] )/sqrt(length(goodInds))
    }
    
    activityMat[currPway, ] = activity
    genes = c(genes, rownames(currNormExpr[goodInds,]) )
    
}

activityMat = activityMat[-which(rowSums(activityMat)==0),]

pathPvals = apply(activityMat, 1, function(x) t.test(x[sensCLInd], x[restCLInd]))

pathPvals_simple = unlist(lapply(pathPvals, function(x) x$p.value) )
pathPval_adj = p.adjust(pathPvals_simple, method="fdr", n=length(pathPvals_simple))

pathSTAT_simple = unlist(lapply(pathPvals, function(x) x$statistic) )
    
sortedPathPvals = sort(pathPval_adj)
sortedPathTscores = sort(pathSTAT_simple)

names(sortedPathTscores) = sub("\\..*", "", names(sortedPathTscores))  

randomTscores = get(load("random_t_scores_Final.Rdata"))
upperThresh = quantile(randomTscores, .9)
lowerThresh = quantile(randomTscores, .1) 

wantPaths =
  c(names(which(sortedPathTscores > upperThresh)),
    names(which(sortedPathTscores < lowerThresh))) 
      
finalPaths = sortedPathTscores[names(sortedPathTscores) %in% wantPaths]

descripL = length(finalPaths)

dscrip_Tscores = c()
for( i in 1:descripL ){
  dscrip_Tscores = c(dscrip_Tscores, keggGet( names(finalPaths)[i] )[[1]]$NAME )
}

ppp = cbind(dscrip_Tscores, finalPaths)

############################## A2. find gene modules ############################
################################################################################

allGenes = c()

for (pway in rownames(ppp)){
  currGenes = unique(hgPwayGenes2[names(hgPwayGenes2)==pway])
  allGenes = unique(c(allGenes, currGenes))
}

string_db = STRINGdb$new(version = "10", species=9606)
df = data.frame("genes" = allGenes)
mapped = string_db$map(df, "genes")
mapped = mapped[!is.na(mapped[,2]), ]

strAdjMat = matrix(0, nrow(mapped), nrow(mapped))
rownames(strAdjMat) = mapped[,2]
colnames(strAdjMat) = mapped[,2]

for ( ind in 1:nrow(mapped) ){
  
  flag = FALSE
  strID = mapped[ind,2]
  
  tryCatch( string_db$get_neighbors(strID), error=function(e) flag<<-TRUE)
  
  if(flag){ next }
  
  interactors = string_db$get_neighbors(strID)
  interactions = string_db$get_interactions(c(strID, interactors))
  
  targetInts = interactions[ interactions[, "from"] == strID, c("to", "combined_score") ]
  targetInts2 = interactions[ interactions[, "to"] == strID, c("from", "combined_score") ]
  
  ## remove genes that are not in expression dataset
  targetInts = targetInts[targetInts[,1] %in% rownames(strAdjMat), ]
  targetInts2 = targetInts2[targetInts2[,1] %in% rownames(strAdjMat), ]
  
  qq = targetInts[,1] %in% targetInts[,2]
  
  if (any(qq)){print("duplicates detected")}
  
  strAdjMat[ strID, targetInts[,1] ] = targetInts[, 2]
  strAdjMat[ strID, targetInts2[,1] ] = targetInts2[, 2]
  
}

dupGene = which( duplicated(colnames(strAdjMat)) )

if(length(dupGene)>0){strAdjMat = strAdjMat[-dupGene, -dupGene]}

## this file contains the strAdjMat and importGenes used in the paper's analysis
load("VEM.Rdata")
################################################################################

strAdjMat = strAdjMat/1000
diag(strAdjMat) = rep(1, ncol(strAdjMat))
strAdjMat[ strAdjMat>=0.4 ] = 1
strAdjMat[ strAdjMat<0.4 ] = 0

## compute TOM matrix
dissTOM = TOMdist(strAdjMat, TOMType = "unsigned", TOMDenom = "min", verbose = 1, indent = 0)

geneTree = hclust(as.dist(dissTOM), method = 'average')

#module identification using dynamic tree cut algorithm
modules = cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 4, pamRespectsDendro = FALSE,
                        minClusterSize = 30)
#assign module colours
module.colours = labels2colors(modules)

gene_Names = colnames(strAdjMat)
gene_Names_hg = mapped[match(gene_Names, mapped[,2]) ,1]

kIM = intramodularConnectivity(strAdjMat, module.colours, scaleByMax=T)
rownames(kIM) = gene_Names_hg
kIM = kIM[!is.nan(kIM[,2]), ]

module_names = as.character(modules)
names(gene_Names_hg) = module_names 

###################### A3. Identify important genes #############################
################################################################################

unqMods = unique(modules)
y = rep("res", ncol(expr))
y[sensCLInd] = "sens"
y = factor(y)

set.seed(54328)
featList = c()

# i of ten was used in paper
for (iter in 1:1){
  print(iter)
  ## set.seed(sample(1:1000,1))
  importFeatures = c()
  for(ind in 1:length(unqMods)){
    
    currMod = unqMods[ind]
    
    currMod = unqMods[1]
    
    currGenes = gene_Names_hg[names(gene_Names_hg)==currMod]
    currGenes = currGenes[currGenes %in% rownames(expr)]
    currEx = expr[currGenes,]
    currEx = t(currEx)
    bortOut = Boruta(currEx, y)
    importGenes = names(bortOut$finalDecision[bortOut$finalDecision=="Confirmed"]) 
    if (length(importGenes) == 0){next}
    importFeatures = rbind(importFeatures, cbind(importGenes, as.character(currMod) ))
    print(ind)
  }
  
  featList = c(featList, list(importFeatures[,1]))
}

importGenes = Reduce(intersect, featList)
importGenes = as.character( importFeatures[,1] )

## this file contains the strAdjMat and importGenes used for paper's analysis
load("VEM.Rdata")
################################################################################

############################## A4. SVM-RFE #####################################
################################################################################

source("../common_files/msvmRFE_LOO.R")

sensCLInd = !restCLInd

y = rep("res", ncol(expr))
y[sensCLInd] = "sens"
y = factor(y)

set.seed(455)
featList = c()

## revert symbols
rownames(expr)[rownames(expr)=="HLA-DQA1"] = "HLA.DQA1"
rownames(expr)[rownames(expr)=="HLA-DQB1"] = "HLA.DQB1"
importGenes = importGenes[importGenes %in% rownames(expr)]

trainMat = data.frame(t(expr[importGenes, ] ))
trainMat = cbind(y, trainMat)

nfold = 10
nrows = nrow(trainMat)
folds = rep(1:nfold, len=nrows)[sample(nrows)]
folds = lapply(1:nfold, function(x) which(folds == x))

results = lapply(folds, svmRFE.wrap, trainMat, k=1, halve.above=100)

### uncomment code below to estimate error for a given number of features
# featsweep = lapply(1:25, FeatSweep.wrap, results, trainMat)
# save(featsweep, file="feature_error.Rdata")
# errors = sapply(featsweep, function(x) ifelse(is.null(x), NA, x$error))

top.features = WriteFeatures(results, trainMat, save=F)

## find pathways that are enriched
wantMods = names(gene_Names_hg[gene_Names_hg %in% top.features[1:4,1]])
wantGenes = gene_Names_hg[names(gene_Names_hg) %in% wantMods]

ensmblENTZ = getBM(attributes = c("hgnc_symbol", "entrezgene_id"),
                   filters = "hgnc_symbol", values = wantGenes, mart = ensembl)

finalPaths = c()

for (mod in unique(wantMods)){
  
  currGenes = wantGenes[names(wantGenes) == mod]
  currEns = ensmblENTZ[ensmblENTZ[,1] %in% currGenes, 2]
  deKeggPways = kegga(currEns, species = "Hs")
  
  tpkg = list(topKEGG(deKeggPways, n=10))
  
  finalPaths = c(finalPaths, tpkg)
}

names(finalPaths) = unique(wantMods) 

########################## B1. Correlated genes ################################
################################################################################ 

currNormExpr = expr

tResult = apply(currNormExpr, 1, function(x) t.test(x[sensCLInd], x[restCLInd])) 
tResult_simple = unlist(lapply(tResult, function(x) x$p.value) ) 

tResult_simple_cor = p.adjust(tResult_simple, method = "holm")
length(tResult_simple_cor[tResult_simple_cor<0.1])

sigGenes = names(tResult_simple_cor[tResult_simple_cor<0.1])

ensmblENTZ = getBM(attributes = c("hgnc_symbol", "entrezgene_id"), 
                   filters = "hgnc_symbol", values = sigGenes, mart = ensembl)

paths = goana(ensmblENTZ[,2])

paths = paths[order(paths[,"P.DE"]),]
paths2 = paths[paths[,"Ont"]=="BP",]

# padj = p.adjust(paths2[,"P.DE"],method = "holm")

# paths2[,"P.DE"] = padj
head(paths2[order(paths2[,5]),])
paths2Want = paths2[paths2[,"P.DE"]<0.1,]
qq = head(paths2Want[order(paths2Want[,5]),], 20)
outDat = data.frame(GOID=rownames(qq), TERM=qq[,"Term"], pVal=qq[,"P.DE"] )
colnames(outDat)[3] = "p-value"

write.csv(outDat, quote = FALSE, row.names = FALSE, file="ttVEM.csv")


######################### B2. Linear model genes ###############################
################################################################################ 

load("braf_auc.Rdata")

wantCCLEnames = sub("_SKIN", "", colnames(expr))
wantCCLEnames = gsub("-", "", wantCCLEnames)
wantCCLEnames = toupper(wantCCLEnames)

colnames(expr) = wantCCLEnames

BRAFVec = read.csv("VEM.csv")

aucMatch = BRAFVec[match(wantCCLEnames, BRAFVec[,1]), ]
aucNoNa = aucMatch[!is.na(aucMatch[,1]),]

aucNoNa = aucNoNa[!duplicated(aucNoNa[,1]), ]
expr = expr[, aucNoNa[,1]]

control <- trainControl(method = "cv",
                        number = 10,
                        search = "random",
                        verboseIter = TRUE)

mads = apply(expr,1, mad)
wantMads = names(tail(sort(mads), 5000))

trainMat = cbind(AUC=aucNoNa[,2], t(expr[wantMads,]))

elastic_model <- train(AUC ~ .,
                       data = trainMat,
                       method = "glmnet",
                       preProcess = c("center", "scale"),
                       tuneLength = 20,
                       trControl = control)

model2 = elastic_model

var <- data.frame(as.matrix(coef(model2$finalModel, model2$bestTune$lambda)))
sigGenes = rownames(var)[var$s1!=0]

ensembl=useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")

ensmblENTZ = getBM(attributes = c("hgnc_symbol", "entrezgene_id"), 
                   filters = "hgnc_symbol", values = sigGenes, mart = ensembl)

paths = goana(ensmblENTZ[,2])

paths = paths[order(paths[,"P.DE"]),]
paths2 = paths[paths[,"Ont"]=="BP",]

padj = p.adjust(paths2[,"P.DE"],method = "holm")

paths2[,"P.DE"] = padj

outVEM = head(paths2[order(paths2[,5]),],100)

write.csv(data.frame(GO_ID=rownames(outVEM), Term=outVEM[,1], P.value=outVEM[,5]), 
          quote = FALSE, row.names = FALSE, file="linear_model.csv")

############################## C. Visualization ################################
################################################################################

topGenes = as.character(top.features[1:4,1])
processList = list()

for(topGene in topGenes){
  wantMod = names(gene_Names_hg)[gene_Names_hg == topGene]
  modgeneNames = gene_Names_hg[names(gene_Names_hg)==wantMod]
  hgncGenes = getBM(filters= "hgnc_symbol", attributes= c( "entrezgene_id", "hgnc_symbol"), values=modgeneNames, mart=ensembl)
  erPathways = goana(hgncGenes[,1], FDR = 0.05)
  erPathways = erPathways[order(erPathways[,"P.DE"]),]
  erPathways = erPathways[erPathways[,"Ont"]=="BP", ]
  processList = c(processList, list(erPathways))
}

processList[[1]][,"P.DE"] = p.adjust(processList[[1]][,"P.DE"])
processList[[2]][,"P.DE"] = p.adjust(processList[[2]][,"P.DE"])
processList[[3]][,"P.DE"] = p.adjust(processList[[3]][,"P.DE"])
processList[[4]][,"P.DE"] = p.adjust(processList[[4]][,"P.DE"])

processListFinal = processList

allTerms = lapply(processListFinal, function(x) cbind(rownames(x),  x[,"Term"]) )

p1 = data.frame(processListFinal[[1]][1:5, ], mod="25")
p2 = data.frame(processListFinal[[2]][1:5, ], mod="12")
p3 = data.frame(processListFinal[[3]][1:5, ], mod="31")
p4 = data.frame(processListFinal[[4]][c(1:4, 12), ], mod="28")

## pathways enrichment plot
allPaths = rbind(p1, p2, p3, p4)
allPaths[,"P.DE"] = -1*log10(allPaths[,"P.DE"])
allPaths[,"Term"] = str_wrap(allPaths[,"Term"], width = 40)
allPaths[,"Term"] = factor(allPaths[,"Term"], levels=rev(allPaths[,"Term"]) )
allPaths[,"mod"] = factor(allPaths[,"mod"], levels=c("25", "12", "31", "28") )

pdf(file="Pathways.pdf", width=11, height=15)

  ggplot(allPaths, aes(mod, Term)) + geom_point(aes(size = DE, fill = P.DE), shape = 21) +
    scale_fill_viridis_c() +
    scale_size_continuous(range=c(4, 14), breaks=c(10,15,30,45)) + 
    labs( x ="Module", y = "Pathway", fill="-log10(p)", size="# of enriched genes") + theme_bw() +
    theme(axis.text.x = element_text(size=20),
          axis.text.y = element_text(size=15, lineheight = 0.78),
          axis.title.x=element_text(size=26), 
          axis.title.y=element_text(size=26),
          legend.text=element_text(size=20),
          legend.title=element_text(size=22))

dev.off()

modules = unique(names(gene_Names_hg))
totModuleGenes = table(names(gene_Names_hg))
totImportGenes = table(names(gene_Names_hg[gene_Names_hg %in% importGenes]))
totModuleGenes = totModuleGenes[names(totImportGenes)]

df <- data.frame(totGenes = as.numeric(totImportGenes), totModules = as.numeric(totModuleGenes))

pdf(file="modules_hex.pdf", width=4.2, height=3)
ggplot(df, aes(totModules, totGenes)) +
  geom_hex(bins = 6) + scale_fill_gradient(low="palegreen2",high="darkgreen", breaks = c(2,5,8)) + labs( x ="Total # of Genes", y = "Relevant # of Genes", fill = "# of Modules") + theme_bw() +
  theme(axis.text = element_text(size=10), axis.title.x=element_text(size=14), axis.title.y=element_text(size=14))
dev.off()                          

###
importGenesMods = top.features[1:4, "FeatureName"]
importMods = names(gene_Names_hg[gene_Names_hg %in% importGenesMods])

moduleTab = gene_Names_hg[order(names(gene_Names_hg))]
moduleTab2 = cbind(moduleTab, names(moduleTab))

moduleTab2 = moduleTab2[moduleTab2[,2] %in% importMods,]

write.csv(moduleTab2, file = "genes_modules.csv", quote = FALSE, row.names = FALSE)

## rank plot
topFeatMods = gene_Names_hg[match(top.features[,"FeatureName"], gene_Names_hg)]  ##top.features[,"FeatureName"]
topfeaturesOut = cbind(top.features, mods = factor(names(topFeatMods), levels = as.character(c(55:0)) ))

theme_dotplot <- theme_bw(12) +
  theme(axis.text.y = element_text(size = rel(.5)),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = rel(1)),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(size = 0.5),
        panel.grid.minor.x = element_blank())


# pdf(file="rank.pdf", width=3.5, height=3)
myData = topfeaturesOut[,c("FeatureName", "AvgRank", "mods")]

minRanks = c()

for (mod in unique(myData[,"mods"])){
  curr = myData[myData[,"mods"]==mod, ]
  minRank = which.min(curr[,2])
  minRanks = rbind(minRanks, curr[minRank,] )
}

colnames(minRanks) = c("gene", "AvgRank", "mods")
minRanks = data.frame(minRanks)
minRanks[,"gene"] = factor(minRanks[,"gene"], levels = rev(minRanks[,1]))
minRanks[,"mods"] = factor(minRanks[,"mods"], levels = rev(minRanks[,3]) )

pdf(file="dotRanks_mods.pdf", width=7, height=9)
# Horizontal version
ggplot(minRanks, aes(x=mods, y=AvgRank)) +
  geom_segment( aes(x=mods, xend=mods, y=0, yend=AvgRank), color="black",  size = 1.4)+
  geom_point( color="black", size=4) +
  theme_bw() + labs( y ="Average Rank", x = "Module") +
  theme(axis.text = element_text(size=13), axis.text.x = element_text(size=18),
        axis.title.x=element_text(size=24), axis.title.y=element_text(size=24)) + 
  coord_flip()

dev.off()

pdf(file="dotRanks_genes.pdf", width=7, height=9)
ggplot(minRanks, aes(x=gene, y=AvgRank)) +
  geom_segment( aes(x=gene, xend=gene, y=0, yend=AvgRank), color="black",  size = 1.4)+
  geom_point( color="black", size=4) +
  theme_bw() + labs( y ="Average Rank", x = "Module") +
  theme(axis.text = element_text(size=13), axis.text.x = element_text(size=18),
        axis.title.x=element_text(size=24), axis.title.y=element_text(size=24)) + 
  coord_flip()
dev.off()

## rfe plot
load("feature_errorLOO.Rdata")
errors = sapply(featsweep, function(x) ifelse(is.null(x), NA, x$error))

d = data.frame(x = 1:length(errors) , y = errors) 

pdf(file="svm_rfe.pdf", width=3, height=3)
  ggplot(d, aes(x=x, y=y)) +
    geom_line() + geom_point() + geom_hline(yintercept=0.23, linetype = 2) +
    labs( x ="Number of Features", y = "10x Cross Validation") + theme_bw() +
    theme(axis.text = element_text(size=10), axis.title.x=element_text(size=14), 
          axis.title.y=element_text(size=14)) + ylab("10x Cross Validation Error")
dev.off()

