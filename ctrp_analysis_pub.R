setwd("~/Desktop/Projects/Machine_learning")

##     0%     25%     50%     75%    100% 
## 0.0691 11.6200 13.4490 14.5850 29.3500 

cellLines = read.csv("v20.meta.per_cell_line.txt", sep="\t")
cellLines[cellLines[,"ccle_primary_hist"]=="malignant_melanoma", ]
cmpd = read.csv("v20.meta.per_compound.txt", sep="\t")
meta = read.csv("v20.meta.per_experiment.txt", sep="\t")
auc = read.csv("v20.data.curves_post_qc.txt", sep="\t")

cmpds = as.character(unique(cmpd[,"master_cpd_id"]))

varofi = function(i){
    mad(auc[auc[,"master_cpd_id"]==cmpds[i], "area_under_curve"])
}

varss = lapply(1:nrow(cmpd), varofi)
qq = unlist(varss)

names(qq) = cmpd[,"cpd_name"]


drugMat  = matrix(0,length(cmpds), length(cellLines[,"master_ccl_id"]))
rownames(drugMat) = cmpds
colnames(drugMat) = unique(cellLines[, "ccl_name"])

for (i in 1:nrow(cmpd)){
    currAuc = auc[auc[,"master_cpd_id"]==cmpds[i], ]
    cellIDs = meta[match(currAuc[,1], meta[,1]), "master_ccl_id"]
    ccleNames = as.character(cellLines[match(cellIDs, cellLines[,"master_ccl_id"]), "ccl_name"] )
    drugMat[cmpds[i],ccleNames] = currAuc[,"area_under_curve"]
    print(i)
}

cor(drugMat, use="complete.obs")

pheatmap(drugMat)
## figure 4
lol = cbind(1:nrow(cmpd), qq)
qq = qq[order(-qq)[1:100]]
bardata = data.frame(name=factor(names(qq), levels=names(qq[order(-qq)])), mad = qq)

pdf(file="drug_mads.pdf", width=4, height=5)

ggplot(data=bardata, aes(x=name, y=mad)) +
    geom_bar(stat="identity") + theme_classic() + 
    theme(axis.title.x = element_text(size = 15), axis.text.x = element_blank(), 
          axis.text = element_text(size = 15),
          axis.title.y = element_text(size = 15),
          legend.text = element_text(size=18),
          legend.title = element_text(size=20), axis.ticks.x = element_blank() ) + 
    xlab("Drug") + ylab("Mean Absolute Deviation")

dev.off()

# aucs = auc[auc[,"master_cpd_id"]==cmpds[465], "area_under_curve"]
# hist(aucs, 30)

aucs = auc[auc[,"master_cpd_id"]=="660217", "area_under_curve"]
hist(aucs, 30)

wantCMPDs = cmpd[cmpd[,1] %in% cmpds[tail(order(qq),6)], ]

brafID2 = cmpd[grep("BRAF",cmpd[,"gene_symbol_of_protein_target"]),]
brafID = cmpd[grep("^PLX-4032", cmpd[,"cpd_name"]), "master_cpd_id"]
erasId = cmpd[grep("RLS3",cmpd[,"cpd_name"]),"master_cpd_id"]
ML210ID = cmpd[grep("^ML210", cmpd[,"cpd_name"]), "master_cpd_id"]

TO = cmpd[grep("^MK-2206", cmpd[,"cpd_name"]), "master_cpd_id"]
TGID = "377381"
gemID = cmpd[grep("gemcitabine", cmpd[,"cpd_name"]), "master_cpd_id"]
at = cmpd[grep("^AT7867", cmpd[,"cpd_name"]), "master_cpd_id"]
atID = "660409"

# paclitaxelID = cmpd[grep("^paclitaxel", cmpd[,"cpd_name"]), "master_cpd_id"]
# vincristineID = cmpd[grep("^vincristine", cmpd[,"cpd_name"]), "master_cpd_id"]
# docetaxelID = cmpd[grep("^docetaxel", cmpd[,"cpd_name"]), "master_cpd_id"]

paclitaxelID = "26956"
vincristineID = "62602"
docetaxelID = "660364"

erasAUC = auc[auc[,"master_cpd_id"]==erasId, ]
BRAFAUC = auc[auc[,"master_cpd_id"]==brafID, ]
GPX4AUC = auc[auc[,"master_cpd_id"]==ML210ID, ]

paxAUC = auc[auc[,"master_cpd_id"]==paclitaxelID, ]
vinAUC = auc[auc[,"master_cpd_id"]==vincristineID, ]s
docAUC = auc[auc[,"master_cpd_id"]==docetaxelID, ]

TGAUC = auc[auc[,"master_cpd_id"]==TGID, ]
atAUC = auc[auc[,"master_cpd_id"]==atID, ]

intCls = intersect(rownames(tgDatFrame), rownames(scoredf_ov))

intCls = intersect(rownames(atDatFrame), rownames(atDatFrame))

plot(aktSig[intCls,1], atDatFrame[intCls,2])

plot(tgDatFrame[intCls,2], atDatFrame[intCls,2])

plot(tgDatFrame[intCls,2], ovExpr["MECOM", intCls])
plot(tgDatFrame[intCls,2], scoredf_ov[intCls,1])

plot(tgDatFrame[intCls, 2], paxDatFrame[intCls, 2])
plot(atDatFrame[intCls, 2], paxDatFrame[intCls, 2])

# GPX4AUC = BRAFAUC
# cellIDs = meta[match(GPX4AUC[,1], meta[,1]), "master_ccl_id"]
# ccleNames = as.character(cellLines[match(cellIDs, cellLines[,"master_ccl_id"]), "ccl_name"])
# outGPX4data = data.frame(cl=ccleNames, auc=GPX4AUC[,"area_under_curve"])
# outGPX4data = outGPX4data[!duplicated(outGPX4data[,1]),]
# BRAFVec = outGPX4data

# write.csv(BRAFVec, file="VEM.csv", row.names = FALSE)

## write to csv
# cellLinesBK = cellLines
cellLines = cellLinesBK
# cellLines = cellLines[cellLines[,"ccle_primary_site"]!="haematopoietic_and_lymphoid_tissue",]
GPX4AUC = paxAUC
cellIDs = meta[match(GPX4AUC[,1], meta[,1]), "master_ccl_id"]
ccleNames = as.character(cellLines[match(cellIDs, cellLines[,"master_ccl_id"]), "ccl_name"])
outGPX4data = data.frame(cl=ccleNames, auc=GPX4AUC[,"area_under_curve"])
outGPX4data = outGPX4data[!duplicated(outGPX4data[,1]),]
paxVec = outGPX4data
# paxVec_noblood = paxVec
# write.csv(outGPX4data, file="paxAuc.csv", row.names = FALSE)

write.csv(outGPX4data, file="paxAuc_no_blood.csv", row.names = FALSE)

GPX4AUC = vinAUC
cellIDs = meta[match(GPX4AUC[,1], meta[,1]), "master_ccl_id"]
ccleNames = as.character(cellLines[match(cellIDs, cellLines[,"master_ccl_id"]), "ccl_name"])
outGPX4data = data.frame(cl=ccleNames, auc=GPX4AUC[,"area_under_curve"])
outGPX4data = outGPX4data[!duplicated(outGPX4data[,1]),]
vinVec = outGPX4data
write.csv(outGPX4data, file="vinAuc.csv", row.names = FALSE)

GPX4AUC = docAUC
cellIDs = meta[match(GPX4AUC[,1], meta[,1]), "master_ccl_id"]
ccleNames = as.character(cellLines[match(cellIDs, cellLines[,"master_ccl_id"]), "ccl_name"])
outGPX4data = data.frame(cl=ccleNames, auc=GPX4AUC[,"area_under_curve"])
outGPX4data = outGPX4data[!duplicated(outGPX4data[,1]),]
docVec = outGPX4data
write.csv(outGPX4data, file="docAuc.csv", row.names = FALSE)

GPX4AUC = TGAUC
cellIDs = meta[match(GPX4AUC[,1], meta[,1]), "master_ccl_id"]
ccleNames = as.character(cellLines[match(cellIDs, cellLines[,"master_ccl_id"]), "ccl_name"])
outGPX4data = data.frame(cl=ccleNames, auc=GPX4AUC[,"area_under_curve"])
outGPX4data = outGPX4data[!duplicated(outGPX4data[,1]),]
tgVec = outGPX4data
write.csv(outGPX4data, file="tgAuc.csv", row.names = FALSE)

GPX4AUC = atAUC
cellIDs = meta[match(GPX4AUC[,1], meta[,1]), "master_ccl_id"]
ccleNames = as.character(cellLines[match(cellIDs, cellLines[,"master_ccl_id"]), "ccl_name"])
outGPX4data = data.frame(cl=ccleNames, auc=GPX4AUC[,"area_under_curve"])
outGPX4data = outGPX4data[!duplicated(outGPX4data[,1]),]
atVec = outGPX4data
write.csv(outGPX4data, file="atAuc.csv", row.names = FALSE)

allNames = unique(c(as.character(paxVec[,1]), as.character(vinVec[,1]), as.character(docVec[,1])))
    
taxolsMat = matrix(-5, 3, length(allNames))

colnames(taxolsMat) = allNames
taxolsMat[1, as.character(paxVec[,1])] = paxVec[,2]
taxolsMat[2, as.character(vinVec[,1])] = vinVec[,2]
taxolsMat[3, as.character(docVec[,1])] = docVec[,2]

pheatmap(taxolsMat)

clNames = as.character(cellLines[,"ccl_name"])

paxDatFrame = 
    data.frame(cl=as.character(paxVec[,1]), 
               auc=paxVec[,2], 
               type = as.character(cellLines[match(as.character(paxVec[,1]),  clNames), "ccle_primary_site"])) 

levels(paxDatFrame[,3])[1] = "unknown"

vinDatFrame = 
    data.frame(cl=as.character(vinVec[,1]), 
               auc=vinVec[,2], 
               type = as.character(cellLines[match(as.character(vinVec[,1]),  clNames), "ccle_primary_site"])) 

docDatFrame = 
    data.frame(cl=as.character(docVec[,1]), 
               auc=docVec[,2], 
               type = as.character(cellLines[match(as.character(docVec[,1]),  clNames), "ccle_primary_site"])) 

tgDatFrame = 
    data.frame(cl=as.character(tgVec[,1]), 
               auc=tgVec[,2], 
               type = as.character(cellLines[match(as.character(tgVec[,1]),  clNames), "ccle_primary_site"])) 


atDatFrame = 
    data.frame(cl=as.character(atVec[,1]), 
               auc=atVec[,2], 
               type = as.character(cellLines[match(as.character(atVec[,1]),  clNames), "ccle_primary_site"])) 

# paxDatFrame[paxDatFrame[,"type"]=="", "type"] = factor( rep("unknown", length(paxDatFrame[paxDatFrame[,"type"]=="", "type"]) ) )

rownames(paxDatFrame) = paxDatFrame[,1]
rownames(vinDatFrame) = vinDatFrame[,1]
rownames(docDatFrame) = docDatFrame[,1]
rownames(tgDatFrame) = tgDatFrame[,1]
rownames(atDatFrame) = atDatFrame[,1]

intCls = intersect(paxDatFrame[,1],atDatFrame[,1])
cor(paxDatFrame[intCls,2], vinDatFrame[intCls,2])

intCls = intersect(paxDatFrame[,1],docDatFrame[,1])

cancerType = as.character(cellLines[match(intCls,  clNames), "ccle_primary_site"])
cancerTypess = rep("Haematopoietic/Lymphoid", length(cancerType))
cancerTypess[cancerType != "haematopoietic_and_lymphoid_tissue"] = "Solid"
# plot(paxDatFrame[intCls,2], vinDatFrame[intCls,2])
tubbData = data.frame(pax = paxDatFrame[intCls,2], 
                      vin = vinDatFrame[intCls,2], 
                      Type = cancerTypess)

corrVal = cor(paxDatFrame[intCls,2], vinDatFrame[intCls,2], use="complete.obs", method="pearson")

## PTX_vin correlation plot
pdf(file="pax_vin_corr.pdf", width=12, height=8)
ggplot(tubbData, aes(x=pax, y=vin)) +
    geom_point(aes(color=Type)) + 
    geom_smooth(method=lm, color="black") + 
    scale_color_brewer(palette="Dark2") + 
    stat_cor(method = "pearson", label.x = 3, label.y = 16, size=8) + 
    theme_classic() + xlab("PTX Response (AUC)") + ylab("Vincristine Response (AUC)") + 
    ylim(0,17) + xlim(0,17) +
    theme(axis.title.x = element_text(size = 20),
          axis.text = element_text(size = 18),
          axis.title.y = element_text(size = 20),
          legend.text = element_text(size=18),
          legend.title = element_text(size=20))
dev.off()
############################# 4 panel AKT/Tubbulin Plot
## PTX_AT correlation plot

## 1 
intCls = intersect(paxDatFrame[,1], atDatFrame[,1])

cancerType = as.character(cellLines[match(intCls,  clNames), "ccle_primary_site"])
cancerTypess = rep("Haematopoietic/Lymphoid", length(cancerType))
cancerTypess[cancerType != "haematopoietic_and_lymphoid_tissue"] = "Solid"
# plot(paxDatFrame[intCls,2], vinDatFrame[intCls,2])
tubbData = data.frame(pax = paxDatFrame[intCls,2], 
                      vin = atDatFrame[intCls,2], 
                      Type = cancerTypess)

A = ggplot(tubbData, aes(x=pax, y=vin)) +
    geom_point(aes(color=Type)) + 
    geom_smooth(method=lm, color="black") + 
    scale_color_brewer(palette="Dark2") + 
    stat_cor(method = "pearson", label.x = 6, label.y = 7.5, size=8) + 
    theme_classic() + xlab("PTX Response (AUC)") + ylab("AT7867 Response (AUC)") + 
    ylim(5.5,17) + xlim(0,17) +
    theme(axis.title.x = element_text(size = 20),
          axis.text = element_text(size = 18),
          axis.title.y = element_text(size = 20),
          legend.text = element_text(size=30),
          legend.title = element_text(size=30))


## 2
intCls = intersect(paxDatFrame[,1], tgDatFrame[,1])

cancerType = as.character(cellLines[match(intCls,  clNames), "ccle_primary_site"])
cancerTypess = rep("Haematopoietic/Lymphoid", length(cancerType))
cancerTypess[cancerType != "haematopoietic_and_lymphoid_tissue"] = "Solid"
# plot(paxDatFrame[intCls,2], vinDatFrame[intCls,2])
tubbData = data.frame(pax = paxDatFrame[intCls,2], 
                      vin = tgDatFrame[intCls,2], 
                      Type = cancerTypess)

B = ggplot(tubbData, aes(x=pax, y=vin)) +
    geom_point(aes(color=Type)) + 
    geom_smooth(method=lm, color="black") + 
    scale_color_brewer(palette="Dark2") + 
    stat_cor(method = "pearson", label.x = 6, label.y = 7.5, size=8) + 
    theme_classic() + xlab("PTX Response (AUC)") + ylab("MK-2206 Response (AUC)") + 
    ylim(5.5,17) + xlim(0,17) +
    theme(axis.title.x = element_text(size = 20),
          axis.text = element_text(size = 18),
          axis.title.y = element_text(size = 20),
          legend.text = element_text(size=30),
          legend.title = element_text(size=30))

## 3
intCls = intersect(vinDatFrame[,1], atDatFrame[,1])

cancerType = as.character(cellLines[match(intCls,  clNames), "ccle_primary_site"])
cancerTypess = rep("Haematopoietic/Lymphoid", length(cancerType))
cancerTypess[cancerType != "haematopoietic_and_lymphoid_tissue"] = "Solid"
# plot(paxDatFrame[intCls,2], vinDatFrame[intCls,2])
tubbData = data.frame(pax = vinDatFrame[intCls,2], 
                      vin = atDatFrame[intCls,2], 
                      Type = cancerTypess)

C = ggplot(tubbData, aes(x=pax, y=vin)) +
    geom_point(aes(color=Type)) + 
    geom_smooth(method=lm, color="black") + 
    scale_color_brewer(palette="Dark2") + 
    stat_cor(method = "pearson", label.x = 6, label.y = 7.5, size=8) + 
    theme_classic() + xlab("Vincristine Response (AUC)") + ylab("AT7867 Response (AUC)") + 
    ylim(5.5,17) + xlim(0,17) +
    theme(axis.title.x = element_text(size = 20),
          axis.text = element_text(size = 18),
          axis.title.y = element_text(size = 20),
          legend.text = element_text(size=30),
          legend.title = element_text(size=30))

## 4
intCls = intersect(vinDatFrame[,1], tgDatFrame[,1])

cancerType = as.character(cellLines[match(intCls,  clNames), "ccle_primary_site"])
cancerTypess = rep("Haematopoietic/Lymphoid", length(cancerType))
cancerTypess[cancerType != "haematopoietic_and_lymphoid_tissue"] = "Solid"
# plot(paxDatFrame[intCls,2], vinDatFrame[intCls,2])
tubbData = data.frame(pax = vinDatFrame[intCls,2], 
                      vin = tgDatFrame[intCls,2], 
                      Type = cancerTypess)

D = ggplot(tubbData, aes(x=pax, y=vin)) +
    geom_point(aes(color=Type)) + 
    geom_smooth(method=lm, color="black") + 
    scale_color_brewer(palette="Dark2") + 
    stat_cor(method = "pearson", label.x = 6, label.y = 7.5, size=8) + 
    theme_classic() + xlab("Vincristine Response (AUC)") + ylab("MK-2206 Response (AUC)") + 
    ylim(5.5,17) + xlim(0,17) +
    theme(axis.title.x = element_text(size = 20),
          axis.text = element_text(size = 18),
          axis.title.y = element_text(size = 20),
          legend.text = element_text(size=30),
          legend.title = element_text(size=30))

library(ggpubr)

pdf(file="AUC_corrs.pdf", width=14, height=14, onefile = FALSE)
ggarrange(A, B, C, D, ncol=2, nrow=2, common.legend = TRUE, legend="bottom")
dev.off()

# plot(vinDatFrame[intCls,2], atDatFrame[intCls,2])

## PTX_MK correlation plot
intCls = intersect(vinDatFrame[,1],docDatFrame[,1])
cor(vinDatFrame[intCls,2], docDatFrame[intCls,2])
plot(vinDatFrame[intCls,2], docDatFrame[intCls,2])

## VIN_AT correlation plot


## VIN_MK correlation plot

intCls = intersect(paxDatFrame[paxDatFrame[,1] %in% colnames(ovExpr), 1]  , 
                   tgDatFrame[tgDatFrame[,1] %in% colnames(ovExpr), 1])

rownames(paxDatFrame) = paxDatFrame[,1]
rownames(tgDatFrame) = tgDatFrame[,1]

cor(paxDatFrame[intCls,2], tgDatFrame[intCls,2])
plot(paxDatFrame[intCls,2], tgDatFrame[intCls,2])


### hippo predicted scores analysis
# hippoScores = read.csv("hippo_scores.csv", stringsAsFactors = FALSE)
# 
# hippoScores[,1] = gsub("-", "", hippoScores[,1] )
# hippoScores[,1] = gsub(" ", "", hippoScores[,1] )
# hippoScores[,1] = toupper(hippoScores[,1])
# 
# highThresh = quantile(hippoScores[,"Predicted"], 0.75)
# lowThresh = quantile(hippoScores[,"Predicted"], 0.25)
# 
# highCLs = hippoScores[hippoScores[,"Predicted"]>highThresh, 1]
# lowCLs = hippoScores[hippoScores[,"Predicted"]<lowThresh, 1]
# 
# ovarCLs = cellLines[cellLines[, "ccle_primary_site"]=="ovary", "ccl_name"]
# # ovarCLs = cellLines[cellLines[, "ccle_primary_site"]=="lung", "ccl_name"]
# hippoScores2 = hippoScores[hippoScores[,"INVENTORY_SAMPLE_NAME"] %in% ovarCLs,]
# 
# highCLs = highCLs[highCLs %in% ovarCLs]
# lowCLs = lowCLs[lowCLs %in% ovarCLs]
# 
# paxHigh = paxDatFrame[paxDatFrame[,1] %in% highCLs, 2]
# paxLow = paxDatFrame[paxDatFrame[,1] %in% lowCLs, 2]
# 
# hippoPax = data.frame(aucs=c(paxHigh, paxLow), class=factor(c(rep("high", length(paxHigh)), 
#                                                             rep("low", length(paxLow)))) )
# 
# p <- ggplot(hippoPax, aes(x=class, y=aucs)) + geom_violin() + 
#     geom_dotplot(binaxis='y', stackdir='center', dotsize=0.2)

## plotting hippo cluster genes


paxDatFrame[,3] = as.character(paxDatFrame[,3])
paxDatFrame[paxDatFrame[,3]=="",3] = "unknown"
paxDatFrame[paxDatFrame[,3] == "haematopoietic_and_lymphoid_tissue",3] = "haematopoietic_and_lymphoid"
paxDatFrame[,3] = factor(paxDatFrame[,3], levels = rev(unique(paxDatFrame[,3])) )

## fig 4
pdf(file="pax_cancer_types.pdf", width=14, height=8)
ggplot(paxDatFrame, aes(x=type, y=auc)) +
        geom_violin() + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.25) + theme_classic() +
        theme(axis.title.x = element_text(size = 20),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          axis.text = element_text(size = 18),
          axis.title.y = element_text(size = 20)) +
    xlab("Primary Site") + ylab("Area under the curve")
dev.off()

# ## Fig 5 Akt
aktSig = read.csv("akt_scores.csv", header=FALSE, stringsAsFactors = FALSE)
aktSig[,4] = gsub("-", "", aktSig[,4])
rownames(aktSig) = aktSig[,4]

## fig 4
aktGenes = read.csv("cmap_akt.csv", skip = 1, stringsAsFactors = FALSE)
upGenesAkt = aktGenes[grep("up", aktGenes[,1]),4]
downGenesAkt = aktGenes[grep("down", aktGenes[,1]),4]

refinedAkt = read.csv("akt_genes_refined.csv", stringsAsFactors = FALSE, header=FALSE)

upGenesAkt = upGenesAkt[upGenesAkt %in% unlist(refinedAkt)]
downGenesAkt = downGenesAkt[downGenesAkt %in% unlist(refinedAkt)]

normExpressionCCLE_final_Old_ALL_filt = 
    normExpressionCCLE_final_Old_ALL[rowSums(normExpressionCCLE_final_Old_ALL>1)>100,]

normExpressionCCLE_final_Old_ALL_filt = 
    normExpressionCCLE_final_Old_ALL_filt[,!duplicated(colnames(normExpressionCCLE_final_Old_ALL_filt))]

maddss = apply(normExpressionCCLE_final_Old_ALL_filt, 1, mad)
topGenes = head(maddss[order(-maddss)],5000)

normExpressionCCLE_final_Old_ALL_filt_all = 
    normExpressionCCLE_final_Old_ALL_filt[names(topGenes),]

rankData <- rankGenes(normExpressionCCLE_final_Old_ALL_filt_all)

# aktLincs = read.csv("akt_genes_lncs.csv", stringsAsFactors = FALSE)
# rownames(aktLincs) = aktLincs[,1]
# aktLincs = aktLincs[,-1]
# 
# ## compute correlation of 
# corss = rep(0, nrow(aktLincs))
# for(i in 1:nrow(aktLincs)){
#     corss[i] = cor(unlist(aktLincs[i,]), aktSig[,1], use="complete.obs")
# }

# upp = upSet = unlist(refinedAkt)[1:72]
# downn  = unlist(refinedAkt)[73:146]
# 
# upp[upp %in% downGenesAkt]
# 
# downn[downn %in% downGenesAkt]
# downn[downn %in% upGenesAkt]

scoredf <- simpleScore(rankData, downSet = unique(downGenesAkt), 
                       upSet = unique(upGenesAkt), knownDirection = TRUE)

scoredf_yap_cluster <- simpleScore(rankData, upSet = yapClusterGenes, knownDirection = TRUE)

# scoredf <- simpleScore(rankData, upSet = upGenesAkt, downSet=downGenesAkt, knownDirection = TRUE)
hist(scoredf[,1])

lower = quantile(scoredf[,1], 0.2)
upper = quantile(scoredf[,1], 0.8)

clLow = rownames(scoredf[scoredf[,1] < lower,])
clHi = rownames(scoredf[scoredf[,1] > upper,])

# lowAkt = rownames(aktSig[aktSig[,1] < -0.25,])
# hiAkt = rownames(aktSig[aktSig[,1] > 0.25,])

###
paxAktLo = paxDatFrame[paxDatFrame[,1] %in% clLow,2]
paxAktHi = paxDatFrame[paxDatFrame[,1] %in% clHi,2]

aktDatFrame = data.frame(AUC = c(paxAktLo, paxAktHi),
                         AktDependnceScore =  c(rep("Low", length(paxAktLo)), rep("High", length(paxAktHi)))  )

pdf(file="akt_PTX_auc.pdf", width=7, height=7)

ggplot(aktDatFrame, aes(x=AktDependnceScore, y=AUC))  +
    geom_violin() +
    geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) +
    scale_fill_manual(values=c("#999999", "#E69F00")) + theme_classic() +
    theme(axis.text = element_text(size = 22), axis.title.x = element_text(size=26),
          axis.title.y = element_text(size=26)) +
    stat_compare_means(method = "t.test", tip.length = 0.01, size=10,
                       symnum.args=list(cutpoints=c(0,0.05,1), symbols=c("*","ns")), label.y = 19,
                       label.x = 1.2) + ylab("PTX Response (AUC)") + xlab("PI3K/Akt Signature")
dev.off()

highThresh = quantile(scoredf_yap_cluster[,1], 0.75)
lowThresh = quantile(scoredf_yap_cluster[,1], 0.25)

highCLs = rownames(scoredf_yap_cluster[scoredf_yap_cluster[,1]>highThresh, ])
lowCLs = rownames(scoredf_yap_cluster[scoredf_yap_cluster[,1]<lowThresh,])

# ovarCLs = cellLines[cellLines[, "ccle_primary_site"]=="ovary", "ccl_name"]
# ovarCLs = cellLines[cellLines[, "ccle_primary_site"]=="lung", "ccl_name"]
# hippoScores2 = hippoScores[hippoScores[,"INVENTORY_SAMPLE_NAME"] %in% ovarCLs,]

# highCLs = highCLs[highCLs %in% ovarCLs]
# lowCLs = lowCLs[lowCLs %in% ovarCLs]

paxHigh = paxDatFrame[paxDatFrame[,1] %in% highCLs, 2]
paxLow = paxDatFrame[paxDatFrame[,1] %in% lowCLs, 2]

hippoPax = data.frame(aucs=c(paxHigh, paxLow), class=factor(c(rep("High", length(paxHigh)), 
                                                              rep("Low", length(paxLow)))) )

pdf(file="hippo_PTX_auc.pdf", width=7, height=7)
ggplot(hippoPax, aes(x=class, y=aucs)) +
    geom_violin() +
    geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) +
    scale_fill_manual(values=c("#999999", "#E69F00")) + theme_classic() +
    theme(axis.text = element_text(size = 22), axis.title.x = element_text(size=26),
          axis.title.y = element_text(size=26)) +
    stat_compare_means(method = "t.test", tip.length = 0.01, size=10,
                       symnum.args=list(cutpoints=c(0,0.05,1), symbols=c("*","ns")), label.y = 19,
                       label.x = 1.2) + ylab("PTX Response (AUC)") + xlab("Yap/Adhesion Signature")

dev.off()

## akt-yap correlation
intCls = intersect(rownames(tgDatFrame), rownames(scoredf_yap_cluster))

aktSigs = scoredf[intCls,1]
hippoScores = hippoScores[!duplicated(hippoScores[,1]),]
rownames(hippoScores) = hippoScores[,1]
hipSigs = hippoScores[intCls,"Predicted"]

plot(scoredf[intCls, 1], aktSig[intCls,1])

plot(scoredf_yap_cluster[intCls, 1], paxDatFrame[intCls,2])

intCls = intersect(rownames(scoredf), colnames(ovExpr))

plot(tgDatFrame[intCls, 2], ovExpr["MECOM",intCls])

cor(tgDatFrame[intCls, 2], scoredf[intCls,1])

plot(tgDatFrame[intCls, 2], paxDatFrame[intCls,2])
plot(atDatFrame[, 2], paxDatFrame[rownames(atDatFrame),2])

yapClusterGenes = unlist(read.csv("yap_PC1_genes.csv", stringsAsFactors = FALSE, header=FALSE))

teadmrtf = unlist(read.csv("tead_mrtf_genes.csv", stringsAsFactors = FALSE, header=FALSE))

ovarian_clinical = read.csv("nationwidechildrens.org_clinical_drug_ov.txt", 
                            stringsAsFactors = FALSE, sep="\t", skip=1)

## exploratory analyssis RPPA
rppa = read.csv("CCLE_RPPA_20181003.csv")
rownames(rppa) = rppa[,1]
qq = sub("_.*", "", rownames(rppa))
rppa = rppa[!duplicated(qq),]
qq = qq[!duplicated(qq)]
rownames(rppa) = qq
rppa = rppa[,-1]
rppa = data.frame(rppa)
rppaProteins = colnames(rppa)

wantRppa = t(rppa[rownames(rppa) %in% ovarCLs, ])
mads = apply(wantRppa, 1, mad)
wantGenes = names(tail(sort(mads),20))

ovPax = paxDatFrame[paxDatFrame[,1] %in% colnames(wantRppa),]
ovPax = paxDatFrame[which(paxDatFrame[,3] == "ovary"),]
ovPaxNames = ovPax[order(ovPax[,2]),]

ovExpr = 
    normExpressionCCLE_final_Old_ALL[,
                                     colnames(normExpressionCCLE_final_Old_ALL) %in% as.character(ovPaxNames[,1]) ]

lol = rowSums(ovExpr>1)
ovExpr = ovExpr[lol>5,]

qq = apply(ovExpr, 1, mad)==
head(sort(qq,decreasing = TRUE), 10000)

pheatmap(ovExpr[c("CAV1", "CDH1"), as.character(ovPaxNames[,1])], cluster_cols = FALSE ) 
pheatmap(wantRppa[wantGenes,as.character(ovPaxNames[,1])], cluster_cols = FALSE)
pheatmap(ovExpr[grep("^ITGA", rownames(ovExpr)), 
                as.character(ovPaxNames[,1])], cluster_cols = FALSE)

lolwat = rbind(ovExpr[,as.character(ovPax[,1]) ], lol = ovPax[,2])
lolwat = lolwat[rownames(lolwat) %in% names(head(sort(qq,decreasing = TRUE),10000)), ]

lolwat = rbind(lolwat, lol = ovPax[,2])
corlol = cor(t(lolwat))

qq = c("MECOM", "FZD4", "HES1", "PSEN2", "JAG2", "PPARG", "HEY1", "CDC16", "MFNG", "EP300")

corlol["lol",qq[qq %in% rownames(corlol)]]

qq = c("MECOM", "DLL3", "JAG1", "MAML3")
corlol["lol",qq[qq %in% rownames(corlol)]]

corrPax = corlol["lol", ]
corrPax = corrPax[corrPax!=1]

barplot( corrPax[order(-corrPax)] )
tail(corrPax[order(-corrPax)],20)

genes = rownames(lolwat)
corVal = rep(0, length(genes))
pvals = rep(0, length(genes))

for(i in 1:nrow(lolwat)){
    corVal[i] = cor.test(lolwat[i,], ovPax[,2], method=c("pearson"))[[4]]
    pvals[i] = cor.test(lolwat[i,], ovPax[,2], method=c("pearson"))[[3]]
}

names(corVal) = genes

upPax8Genes = c("CPA4", 
                "PCDH10", 
                "LIPH", 
                "NT5E", 
                "DUSP10", 
                "SEMA3C", 
                "NAV3",
                "SLC2A3",
                "FAM198B",
                "INPP4B",
                "DUSP1",
                "AHR",
                "SEMA3A",
                "BMF",
                "AMOT",
                "ID3",
                "NRP1",
                "THBS1",
                "STX11",
                "PDK4",
                "ST3GAL5",
                "LYPD1",
                "KRT18",
                "DCBLD2")

downPax8Genes = c("MECOM",
                  "GPRC5B",
                  "TSPAN1",
                  "EHF",
                  "SLITRK5",
                  "UPK3B",
                  "KLHL13",
                  "SLC18B1",
                  "SORCS1",
                  "LRRC8C",
                  "ADAMTS15",
                  "AKAP6",
                  "BMP7",
                  "PAX8",
                  "MAMDC2",
                  "AKAP12",
                  "PTGS1",
                  "GLDC",
                  "HIST1H1D",
                  "BRI3BP",
                  "ITGB3",
                  "COL4A2",
                  "SNTB1",
                  "PSAT1",
                  "RGS10",
                  "NMU",
                  "FGF18",
                  "CLDN16")

ovExpr_filt = ovExpr[names(head(sort(qq,decreasing = TRUE),20000)),]2

rankData <- rankGenes(ovExpr_filt)

mecomGeneSet = read.csv("mecom_change.csv", stringsAsFactors = FALSE)
mecomGeneSet = mecomGeneSet[!duplicated(mecomGeneSet[,1]),]
vivo = mecomGeneSet[,grep("invivo", colnames(mecomGeneSet))]

vivo = vivo[,grep("PAX8", colnames(vivo))]
rownames(vivo) = as.character(mecomGeneSet[,1])

# vivo[which(rowSums(vivo)>1),]

downPax8Genes2 = rownames(vivo)[ any(vivo < -1) ]
upPax8Genes2 = rownames(vivo)[ any(vivo > 1) ]

downPax8Genes3 = downPax8Genes[downPax8Genes %in% downPax8Genes2]
upPax8Genes3 = upPax8Genes[upPax8Genes %in% upPax8Genes2]

scoredf_ov <- simpleScore(rankData, upSet = downPax8Genes3, 
                       downSet = upPax8Genes3, knownDirection = TRUE)

scoredf_ov_yap <- simpleScore(rankData, upSet = teadmrtf, knownDirection = TRUE)

# scoredf <- simpleScore(rankData, upSet = upGenesAkt, downSet=downGenesAkt, knownDirection = TRUE)
# hist(scoredf_ov_yap[,1])

wantNames = rownames(scoredf_ov_yap)[rownames(scoredf_ov_yap) %in% rownames(hippoScores)]
wantNames = rownames(scoredf_ov_yap)[rownames(scoredf_ov_yap) %in% colnames(ovExpr)]

plot(hippoScores[wantNames,"Predicted"], ovPaxNames[wantNames, 2] )

plot(scoredf[rownames(ovPaxNames),1], ovPaxNames[, 2])

plot(scoredf_ov_yap[wantNames,1], hippoScores[wantNames, "Predicted"] )

plot(hippoScores[wantNames,"Predicted"], scoredf_ov[wantNames,1])

plot(hippoScores[colnames(ovExpr),"Predicted"], ovExpr["MECOM",])

plot(scoredf_ov[wantNames,1], ovExpr["THBS1",] )

plot(ovExpr["MECOM",], ovExpr["CAV2",])

pheatmap(ovExpr[c("MECOM", "ITGB3","COL4A2", "TGFBI"),])
pheatmap(ovExpr[c("MECOM", "ITGB3","COL4A2", "TGFBI"),])

library(GEOquery)
library(limma)
library(umap)

# load series and platform data from GEO

gset <- getGEO("GSE15622", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL8414", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

ex <- exprs(gset)
# log2 transform
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
    (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
ex <- log2(ex) }

rownames(ex) = sub("_at", "", rownames(ex))

wantGenes = rownames(ex)[grep("ENSG", rownames(ex))]
wantGenes = sub("_at", "", wantGenes)

ensembl=useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")  

hgncGenes <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id",
                                                             "hgnc_symbol"), values=wantGenes, mart= ensembl)

swapGenes = hgncGenes[match(wantGenes, hgncGenes[,1]) , 2]
hgncGenes = hgncGenes[-which(hgncGenes[,2]==""),]
hgncGenes = hgncGenes[ !duplicated(hgncGenes[,1]),]

exFinal = ex[rownames(ex) %in% hgncGenes[,1], ]
rownames(exFinal) = hgncGenes[match(rownames(exFinal), hgncGenes[,1]) , 2]

exFinal_scale = t(as.data.frame(apply(exFinal, 1, normalize)))

clinical_dat = read.csv("GSE15622_pt_info.txt", sep="\t", stringsAsFactors = FALSE)
clinica_dat_filt = data.frame(samp = unlist(clinical_dat[1,]), 
                              status = unlist(clinical_dat[7,]), 
                              rx = unlist(clinical_dat[9,]), 
                              response = unlist(clinical_dat[10,]) )

clinica_dat_filt2 = clinica_dat_filt[grep("pre", rownames(clinica_dat_filt)),]

pacPts = clinica_dat_filt2[grep("Pac", clinica_dat_filt2[,"rx"]),]

sens = as.character(pacPts[pacPts[,"response"]=="response: sensitive","samp"])
resist = as.character(pacPts[pacPts[,"response"]=="response: resistant","samp"])

dPax8 = downPax8Genes[downPax8Genes %in% rownames(exFinal)]
uPax8 = upPax8Genes[upPax8Genes %in% rownames(exFinal)]

madds = apply(exFinal,1,mad)

lol = tail(sort(madds), 500)

# lol2 = tail(sort(madds), 12000)
# sort(madds[c(dPax8, uPax8)])

lol2 = lol[names(lol) %in% rownames(vivo)]

yapGenes = yapClusterGenes[yapClusterGenes %in% names(lol)]

lol3 = tail(sort(madds), 5000)

# exFinal2 = exFinal[rownames(exFinal) %in% names(lol2), ]

pheatmap(exFinal_scale[c("MECOM", "PAX8", "NMU", yapGenes), c(sens,resist)], cluster_cols = FALSE)

pheatmap(exFinal_scale[rownames(exFinal) %in% names(lol2), c(sens,resist)], cluster_cols = FALSE)

pheatmap(exFinal_scale[dSet[dSet %in% names(lol)], c(sens)], cluster_cols = TRUE)

pheatmap(exFinal_scale[dSet[dSet %in% names(lol)], c(resist)], cluster_cols = TRUE)

pheatmap(exFinal_scale[c("MECOM", "PAX8"), c(resist)], cluster_cols = TRUE)

pheatmap(exFinal_scale[c(names(lol)), c(sens, resist)], cluster_cols = FALSE)

pheatmap(exFinal_scale[c(dPax8, uPax8), c(sens, resist)], cluster_cols = TRUE)

pheatmap(exFinal_scale[c(dPax8, uPax8), c(resist)], cluster_cols = TRUE)

pheatmap(exFinal_scale[yapClusterGenes[yapClusterGenes %in% rownames(exFinal_scale)], c(sens, resist)], cluster_cols = TRUE)

rankData <- rankGenes(exFinal[names(lol3),])

uSet = names(which(apply(vivo[names(lol2), ]> 2,1,any)))
dSet = names(which(apply(vivo[names(lol2), ] < -2,1,any)))

scoredf_ov_pax <- simpleScore(rankData, downSet = unique(uSet), 
                       upSet = unique(dSet), knownDirection = TRUE)

scoredf_ov_pt_yap <- simpleScore(rankData, upSet = unique(yapClusterGenes), knownDirection = TRUE)


scoredf_ov_pt_yap <- simpleScore(rankData, upSet = unique(upGenesAkt), 
                                 downSet = unique(downGenesAkt), knownDirection = TRUE)

ptData = data.frame(Score=c(scoredf_ov_pax[rownames(scoredf_ov_pax) %in% sens, 1], 
                            scoredf_ov_pax[rownames(scoredf_ov_pax) %in% resist, 1]),
                    Group=c(rep("Sensitive", length(sens)) , rep("Resistant", length(resist)) ))

ptData = data.frame(Score=c(scoredf_ov_pt_yap[rownames(scoredf_ov_pt_yap) %in% sens, 1], 
                            scoredf_ov_pt_yap[rownames(scoredf_ov_pt_yap) %in% resist, 1]),
                    Group=c(rep("Sensitive", length(sens)) , rep("Resistant", length(resist)) ))

ggplot(ptData, aes(x=Group, y=Score)) +
    geom_boxplot() + geom_dotplot(binaxis='y', stackdir='center') + theme_classic() +
    theme(axis.text = element_text(size = 10), axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    stat_compare_means(method = "t.test", label.x = 1.07, label.y = 0.2, size=4)

# library(Rtsne)
# # library(hrbrthemes)
# #### 
# 
# kmeanRes = kmeans(t(normExpressionTCGA_ALL[clustGenes,]), 2, iter.max = 20, nstart = 1 )
# clusts2 = kmeanRes$cluster
# 
# lol = c(uSet, dSet)
# 
# qq = Rtsne(t(exFinal[lol[lol %in% rownames(exFinal)], c(sens,resist)]), labels = as.factor(clusts2), perplexity = 5, epoch=1000)
# 
# clinica_dat_filt = clinica_dat_filt[-1,  ]
# clinica_dat_filt2 = clinica_dat_filt[clinica_dat_filt[,1] %in% c(sens,resist),]
# 
# clusterDatFrame = data.frame(tsne1 = qq$Y[,1], tsne2 = qq$Y[,2], Subtype = clinica_dat_filt2[,"response"] )
# 
# ggplot(clusterDatFrame, aes(x=tsne1, y=tsne2, fill=Subtype)) + geom_point(size=3, alpha=0.5, shape=21) + theme_bw() +
#     theme(axis.text = element_text(size = 16), axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18),
#           panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.text = element_text(size=14),
#           legend.title = element_text(size=16)) + scale_fill_viridis(discrete=T, name="")
# 
# # hist(scoredf_ov_pax[rownames(scoredf_ov_pax) %in% sens,1])
# hist(scoredf_ov_pax[rownames(scoredf_ov_pax) %in% resist,1])

# pheatmap(normExpression_ovarian_TCGA[c("MECOM", "GLDC", "ITGB3","COL4A2"), c(CR, Progressive)], cluster_cols = TRUE)
# plot(normExpression_ovarian_TCGA["MECOM",], normExpression_ovarian_TCGA["COL4A2",])
# responseDat = responseInfo[responseInfo[,"bcr_patient_barcode"] %in% TaxolPts, ]
# responseDat[,]

### pax8 reg genes analysis

metab = data.frame(read.csv("CCLE_metabolomics_20190502.csv", stringsAsFactors = FALSE))
rownames(metab) = metab[,1]
metab = metab[,-c(1:3)]
# metab = metab[colnames(normExpression_filt_CCLE),-c(1:3)]
metab = t(metab)

colnames(metab) = sub("_.*","", colnames(metab))
pheatmap(metab[,as.character(ovPaxNames[,1]) ], cluster_cols = FALSE )
madds = apply(metab[,as.character(ovPaxNames[,1]) ], 1, mad)

# notchgenes = read.csv("notchGenes.txt", stringsAsFactors = FALSE, skip=1)
notchgenes = read.csv("ovarian_notch_genes_full.csv", stringsAsFactors = FALSE)
# notchgenes = read.csv("ovarian_notch_genes.csv", stringsAsFactors = FALSE)[,2]

rownames(notchgenes) = notchgenes[,2]
rownames(notchgenes) = toupper(rownames(notchgenes))

notchUpGenes = rownames(notchgenes[notchgenes[,"log2FoldChange"]>1,])
notchDownGenes = rownames(notchgenes[notchgenes[,"log2FoldChange"]<= 1,])

ovExprScale = t(as.data.frame(apply(ovExpr, 1, normalize)))

pheatmap(ovExprScale[notchgenes[notchgenes %in% rownames(ovExprScale)], 
                as.character(ovPaxNames[,1])], cluster_cols = TRUE)

pheatmap(ovExprScale[c("MECOM", "CPA4", "PCDH10","CLDN16","BMP7","NOS1"), 
                     as.character(ovPaxNames[,1])], , cluster_cols = FALSE)

lol = cor(t(ovExprScale))
lol2 = lol["MECOM",]

hist(scoredf[colnames(ovExprScale),1])

hippoScores[hippoScores[,1] %in% c("MCAS", "OC314", "A2780", " OVSAHO", "OVISE", "OV90"),]

plot(tgDatFrame[colnames(ovExpr_filt),2],
     ovExpr_filt["MECOM",])

intCls = intersect(rownames(scoredf), rownames(tgDatFrame))
cor(tgDatFrame[intCls,2],
     scoredf[intCls,1])

plot(atDatFrame[intCls,2], scoredf[intCls,1])

aktSig = aktSig[!is.na(aktSig[,1]),]
intCls = intersect(rownames(aktSig), rownames(scoredf))

### 
intCls = intersect(rownames(scoredf_ov), rownames(paxDatFrame))

plot(tgDatFrame[intCls,2],
    scoredf_ov[intCls,1])

aktSig = aktSig[!is.na(aktSig[,1]),]

intCls = intersect(rownames(aktSig), rownames(paxDatFrame))

plot(atDatFrame[intCls,2], ovExpr["MECOM",intCls])

cor(tgDatFrame[intCls,2], paxDatFrame[intCls,2])

# (paxDatFrame[intCls,2], aktSig[intCls,1])
## validate in patients
## pan-drug analysis

# i = 1

# drugData2 = read.csv("GDSC2_fitted_dose_response_25Feb20.csv")

# erkData = drugData2[drugData2[,"DRUG_NAME"] == "ML210", ]
# conc10 = erkData[erkData[,"MAX_CONC"]==2, ]

hist(conc10[,"AUC"],20)

GPX4AUC = BRAFAUC
cellIDs = meta[match(GPX4AUC[,1], meta[,1]), "master_ccl_id"]
ccleNames = as.character(cellLines[match(cellIDs, cellLines[,"master_ccl_id"]), "ccl_name"])
outGPX4data = data.frame(ccleNames, GPX4AUC[,"area_under_curve"])
dupRows = outGPX4data
outGPX4data = outGPX4data[!duplicated(outGPX4data[,1]),]

pdf(file="ML210_AUC.pdf", width=6, height=5)

ggplot(gpx4Aucs, aes(x=auc)) + geom_histogram(binwidth=0.5) + geom_vline(xintercept = 9, linetype="dotted") +
    theme_bw() + xlab("AUC") + ylab("Counts") + 
    theme(axis.text.x = element_text(size=15),
          axis.text.y = element_text(size=15),
          axis.title.x=element_text(size=20), 
          axis.title.y=element_text(size=20))

dev.off()
# cancerTypes = cellLines[match(cellIDs, cellLines[,1]), "ccle_primary_hist"]
# badInds = which(cancerTypes=="")
# # 
# # # erasAUCFinal = cbind(erasAUC[-badInds,], cancerTypes = cancerTypes[-badInds])
# # # erasAUCFinal = cbind(BRAFAUC[-badInds,], cancerTypes = cancerTypes[-badInds])
# BRAFAUCFinal = cbind(BRAFAUC[-badInds,], cancerTypes = cancerTypes[-badInds])
# # # erasAUCFinalPlot = erasAUCFinal[, c("cancerTypes","area_under_curve")]
# # 
# wantCells = sub("_SKIN", "", sampleLable)
# wantCellsIDs = cellLines[cellLines[,"ccl_name"] %in% wantCells, 1]
# wantExps = meta[meta[,"master_ccl_id"] %in% wantCellsIDs, "experiment_id"]
# 
# eras_auc = erasAUCFinal[match(wantExps, erasAUCFinal[,1]), "area_under_curve"]
# braf_auc = BRAFAUCFinal[match(wantExps, BRAFAUCFinal[,1]), "area_under_curve"]
BRAFaucF = BRAFAUCFinal[grep("melanoma", BRAFAUCFinal[,"cancerTypes"]), "area_under_curve"]

gpx4AUCs = GPX4AUC[,"area_under_curve"]

plot(braf_auc, eras_auc, xlim=c(8,18), ylim=c(8,18))

e = ggplot(erasAUCFinal, aes(x = cancerTypes, y = area_under_curve))

e + geom_violin(trim = FALSE) +
    geom_dotplot(
        binaxis='y', stackdir='center', binwidth=0.1,
        color = "black", fill = "#999999"
    ) +
    stat_summary(
        fun.data="mean_sdl",  fun.args = list(mult=1),
        geom = "pointrange", color = "#FC4E07", size = 0.4
    ) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), plot.title = element_text(hjust = 0.5)) +
    ggtitle("VEM") +
    xlab("Cancer Type") + ylab("Area under the curve")

### src data
brafID = cmpd[grep("SRC",cmpd[,"gene_symbol_of_protein_target"]), ]

## sara_ID;; need to get rid of 
# BRAFAUC = auc[auc[,"master_cpd_id"]=="660346", ]
# 
# cellIDs = meta[match(auc[,1], meta[,1]), "master_ccl_id"]
# 
# cancerTypes = cellLines[match(cellIDs, cellLines[,1]), "ccle_primary_hist"]
# badInds = which(cancerTypes=="")

# erasAUCFinal = cbind(erasAUC[-badInds,], cancerTypes = cancerTypes[-badInds])
# erasAUCFinal = cbind(BRAFAUC[-badInds,], cancerTypes = cancerTypes[-badInds])

# BRAFAUCFinal = cbind(BRAFAUC[-badInds,], cancerTypes = cancerTypes[-badInds])
# erasAUCFinalPlot = erasAUCFinal[, c("cancerTypes","area_under_curve")]

wantCells = sub("_SKIN", "", rownames(lol2))
wantCellsIDs = cellLines[cellLines[,"ccl_name"] %in% wantCells, ]

badCells = c("COLO679", "COLO783", "HS695T", "G361", "SKMEL24", "SKMEL1", "COLO741")

wantCellsIDs = wantCellsIDs[-which(as.character(wantCellsIDs[,"ccl_name"]) %in% badCells), ]

durgsList = cmpd[,1]
drugMat = c()

for(wantCellsID in wantCellsIDs[,1]){
    # wantCellsID = "472"
    wantExp = meta[meta[,"master_ccl_id"]==wantCellsID, 1][1]
    wantDrugs = auc[auc[,"experiment_id"]==wantExp, ]
    wantDrugsOrdered = wantDrugs[match(durgsList, wantDrugs[,"master_cpd_id"]), "area_under_curve"]
    drugMat = cbind(drugMat, wantDrugsOrdered)
}

rownames(drugMat) = cmpd[,2]
colnames(drugMat) = wantCellsIDs[,2]

wantCPDnames = cmpd[grep("SRC", cmpd[,"gene_symbol_of_protein_target"]), "cpd_name"]
# wantCPDnames = c(as.character(wantCPDnames), "PLX-4032")
# wantCPDnames = cmpd[grep("BRAF", cmpd[,"gene_symbol_of_protein_target"]), "cpd_name"]
# wantCPDnames = c()
# wantRows = which(rowSums(apply(drugMat, 2, is.na))==0)
# drugMatWant = drugMat[wantRows, ]

# drugMat[is.na(drugMat)] = 3
# pheatmap( drugMat[rownames(drugMat) %in% c(wantCPDnames, "erastin"), wantCells[wantCells %in% wantCellsIDs[,2]] ])
# pheatmap( drugMatWant[rownames(drugMatWant) %in% c(wantCPDnames, "erastin"), wantCells[wantCells %in% wantCellsIDs[,2]] ])

resistCells = c("LOXIMVI", "RPMI7951", "MDAMB435S", "WM793", "A2058", "SKMEL5", "K029AX") 

resistCells2 = colnames(normExpression_filt_CCLE)[which(sensCLInd=="RES")]
resistCells2 = sub("_SKIN", "", resistCells2)

pheatmap( 
    drugMatWant[rownames(drugMatWant) %in% c("saracatinib", "dasatinib"), colnames(drugMatWant) %in% resistCells2], cluster_cols = FALSE )

pheatmap( drugMatWant[lol, resistCells] )

# drugMat = t(drugMat)
##
# durgsList = 

# cmpd[grep("BRAF",cmpd[,"gene_symbol_of_protein_target"]), ]
#     
# wantCellsIDs = cellLines[match(wantCells, cellLines[,"ccl_name"]), "master_ccl_id"]
# 
# wantExps = meta[   match(wantCellsIDs, meta[,"master_ccl_id"])   , "experiment_id"]
# 
# braf_auc = BRAFAUCFinal[match(wantExps, BRAFAUCFinal[,1]), "area_under_curve"]
# 
# names(braf_auc) = rownames(lol2)
# 
# plot(braf_auc, lol2[,"TotalScore"]*10)

# cor(braf_auc, lol2[,"TotalScore"], na.rm=TRUE)
# BRAFAUCFinal[match(wantExps, BRAFAUCFinal[,1]), c("experiment_id", "area_under_curve")]

# e + geom_violin(trim = FALSE) + 
#     stat_summary(
#         fun.data = "mean_sdl",  fun.args = list(mult = 1), 
#         geom = "pointrange", color = "black"
#     ) +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# 

# ggplot(erasAUCFinal, aes(cancerTypes, area_under_curve)) +
#     geom_sina(aes(color = cancerTypes), size = 0.7) + 
#     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#     # scale_color_manual(values =  c("#00AFBB", "#E7B800", "#FC4E07", "#00CCFF", "#FF9999"))

# ggplot(erasAUCFinal, aes(cancerTypes, area_under_curve)) +
#     geom_jitter(
#         position = position_jitter(0.2), color = "darkgray"
#     ) + 
#     geom_pointrange(
#         aes(ymin = len-sd, ymax = len+sd),
#         data = df.summary
#     ) + 
#     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

