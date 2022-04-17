library(limma)
library(biomaRt)

load("workdata.Rdata")
load("GP4_working.Rdata")

###
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# topGenes = c("RAC1", "GRM7", "CNGA3", "FGF17")
topGenes = as.character(top.features[1:4,1])
names(gene_Names_hg)[gene_Names_hg == topGenes]

# processList_GoFunc = list()

processList = list()

for(topGene in topGenes){
    wantMod = names(gene_Names_hg)[gene_Names_hg == topGene]
    modgeneNames = gene_Names_hg[names(gene_Names_hg)==wantMod]
    hgncGenes = getBM(filters= "hgnc_symbol", attributes= c( "entrezgene_id", "hgnc_symbol"), values=modgeneNames, mart=human)
    erPathways = goana(hgncGenes[,1], FDR = 0.05)
    erPathways = erPathways[order(erPathways[,"P.DE"]),]
    erPathways = erPathways[erPathways[,"Ont"]=="BP", ]
    processList = c(processList, list(erPathways))
}

# processListFinal = processList

processList[[1]][,"P.DE"] = p.adjust(processList[[1]][,"P.DE"])
processList[[2]][,"P.DE"] = p.adjust(processList[[2]][,"P.DE"])
processList[[3]][,"P.DE"] = p.adjust(processList[[3]][,"P.DE"])
processList[[4]][,"P.DE"] = p.adjust(processList[[4]][,"P.DE"])

processListFinal = processList

allTerms = lapply(processListFinal, function(x) cbind(rownames(x),  x[,"Term"]) )
lol = Reduce(merge,allTerms)

# processListFinal[[1]][,"P.DE"] = p.adjust(processList[[1]][,"P.DE"])
# processListFinal[[2]][,"P.DE"] = p.adjust(processList[[2]][,"P.DE"])
# processListFinal[[3]][,"P.DE"] = p.adjust(processList[[3]][,"P.DE"])
# processListFinal[[4]][,"P.DE"] = p.adjust(processList[[4]][,"P.DE"])

# library(stringr)

p1 = data.frame(processListFinal[[1]][1:5, ], mod="25")
p2 = data.frame(processListFinal[[2]][1:5, ], mod="36")
p3 = data.frame(processListFinal[[3]][1:5, ], mod="16")
p4 = data.frame(processListFinal[[4]][c(1:4, 12), ], mod="8")

## pathways enrichment plot
allPaths = rbind(p1, p2, p3, p4)
allPaths[,"P.DE"] = -1*log10(allPaths[,"P.DE"])
allPaths[,"Term"] = str_wrap(allPaths[,"Term"], width = 40)
allPaths[,"Term"] = factor(allPaths[,"Term"], levels=rev(allPaths[,"Term"]) )

pdf(file="Pathways_gpx4.pdf", width=11, height=15)

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

#### relevant features in important modules Fig 2
modules = unique(names(gene_Names_hg))
totModuleGenes = table(names(gene_Names_hg))
totImportGenes = table(names(gene_Names_hg[gene_Names_hg %in% importGenes]))
totModuleGenes = totModuleGenes[names(totImportGenes)]

df <- data.frame(totGenes = as.numeric(totImportGenes), totModules = as.numeric(totModuleGenes))

pdf(file="modules_hex_gpx4.pdf", width=4.2, height=3)
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

write.csv(moduleTab2, file = "genes_modules_GPX4.csv", quote = FALSE, row.names = FALSE)

## rfe plot

load("feature_error_GPX4.Rdata")
errors = sapply(featsweep, function(x) ifelse(is.null(x), NA, x$error))

source("./SVM-RFE/msvmRFE.R")
no.info = min(prop.table(table(drugStats)))

PlotErrors(errors, no.info=0.39)

d = data.frame(x = 1:length(errors) , y = errors) 

pdf(file="svm_rfe_ML210.pdf", width=3, height=3)

ggplot(d, aes(x=x, y=y)) +
    geom_line() + geom_point() + geom_hline(yintercept=0.39, linetype = 2) +
    labs( x ="Number of Features", y = "10x Cross Validation") + theme_bw() +
    theme(axis.text = element_text(size=10), axis.title.x=element_text(size=14), 
          axis.title.y=element_text(size=14)) + ylab("10x Cross Validation Error")

dev.off()

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

# ggplot(topfeaturesOut, aes(x = AvgRank, y = mods)) + geom_point(size=1.5) + labs( x ="Average Rank", y = "Module") +
#     theme(axis.text = element_text(size=12), axis.text.x = element_text(size=10),
#           axis.title.x=element_text(size=14), axis.title.y=element_text(size=14)) + theme_classic()

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
# qq = myData[which(myData[,"AvgRank"] %in% minRanks[,2]), ]
# qq[,"mods"] = factor(qq[,"mods"], levels = unique(rev(qq[,"mods"])) )
# 
# qq = data.frame(qq, gene=top.features[match(qq[,1], top.features[,"AvgRank"]),1])
# qq[,"gene"] = factor(qq[,"gene"], levels = rev(top.features[,"FeatureName"]))

pdf(file="dotRanks_GPX4_mods2.pdf", width=10, height=9)

# Horizontal version
ggplot(minRanks, aes(x=mods, y=AvgRank)) +
    geom_segment( aes(x=mods, xend=mods, y=0, yend=AvgRank), color="black",  size = 1.4)+
    geom_point( color="black", size=4) +
    theme_bw() + labs( y ="Average Rank", x = "Module") +
    theme(axis.text = element_text(size=12), axis.text.x = element_text(size=18),
          axis.title.x=element_text(size=24), axis.title.y=element_text(size=24)) + 
    coord_flip()

dev.off()


###
load("pax_all_final.Rdata")

topGenes = as.character(top.features[1:4,1])

processList = list()

for(topGene in topGenes){
    wantMod = names(gene_Names_hg)[gene_Names_hg == topGene]
    modgeneNames = gene_Names_hg[names(gene_Names_hg)==wantMod]
    hgncGenes = getBM(filters= "hgnc_symbol", attributes= c( "entrezgene_id", "hgnc_symbol"), values=modgeneNames, mart=human)
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

p1 = data.frame(processListFinal[[1]][1:5, ], mod="7")
p2 = data.frame(processListFinal[[2]][1:5, ], mod="39")
p3 = data.frame(processListFinal[[3]][1:5, ], mod="37")
p4 = data.frame(processListFinal[[4]][c(1:4, 12), ], mod="19")

## pathways enrichment plot
allPaths = rbind(p1, p2, p3, p4)
allPaths[,"P.DE"] = -1*log10(allPaths[,"P.DE"])
allPaths[,"Term"] = str_wrap(allPaths[,"Term"], width = 40)
allPaths[,"Term"] = factor(allPaths[,"Term"], levels=rev(allPaths[,"Term"]) )

pdf(file="Pathways_pax.pdf", width=11, height=15)

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


####

load("pax_final.Rdata")

processList = list()

for(topGene in topGenes){
    wantMod = names(gene_Names_hg)[gene_Names_hg == topGene]
    modgeneNames = gene_Names_hg[names(gene_Names_hg)==wantMod]
    hgncGenes = getBM(filters= "hgnc_symbol", attributes= c( "entrezgene_id", "hgnc_symbol"), values=modgeneNames, mart=human)
    erPathways = goana(hgncGenes[,1], FDR = 0.05)
    erPathways = erPathways[order(erPathways[,"P.DE"]),]
    erPathways = erPathways[erPathways[,"Ont"]=="BP", ]
    processList = c(processList, list(erPathways))
}

# processListFinal = processList

processList[[1]][,"P.DE"] = p.adjust(processList[[1]][,"P.DE"])
processList[[2]][,"P.DE"] = p.adjust(processList[[2]][,"P.DE"])
processList[[3]][,"P.DE"] = p.adjust(processList[[3]][,"P.DE"])
processList[[4]][,"P.DE"] = p.adjust(processList[[4]][,"P.DE"])

processListFinal = processList

# allTerms = lapply(processListFinal, function(x) cbind(rownames(x),  x[,"Term"]) )
# lol = Reduce(merge,allTerms)

topGenes = as.character(top.features[1:3,1])
gene_Names_hg[which(gene_Names_hg %in% topGenes)]

p1 = data.frame(processListFinal[[1]][1:5, ], mod="10")
p2 = data.frame(processListFinal[[2]][1:5, ], mod="18")
p3 = data.frame(processListFinal[[3]][1:5, ], mod="36")
# p4 = data.frame(processListFinal[[4]][c(1:4, 12), ], mod="19")

## pathways enrichment plot
allPaths = rbind(p1, p2, p3)
allPaths[,"P.DE"] = -1*log10(allPaths[,"P.DE"])
allPaths[,"Term"] = str_wrap(allPaths[,"Term"], width = 40)
allPaths[,"Term"] = factor(allPaths[,"Term"], levels=rev(allPaths[,"Term"]) )

pdf(file="Pathways_pax_solid.pdf", width=11, height=10)

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

####
load("./revisions_analysis/GPX4_allPaths_final2.Rdata")

human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
top.features = WriteFeatures(results, trainMat, save=F)  
topGenes = as.character(top.features[1:4,1])
# names(gene_Names_hg)[gene_Names_hg == topGenes]

processList = list()

for(topGene in topGenes){
    wantMod = names(gene_Names_hg)[gene_Names_hg == topGene]
    modgeneNames = gene_Names_hg[names(gene_Names_hg)==wantMod]
    hgncGenes = getBM(filters= "hgnc_symbol", attributes= c( "entrezgene_id", "hgnc_symbol"), values=modgeneNames, mart=human)
    erPathways = goana(hgncGenes[,1], FDR = 0.05)
    erPathways = erPathways[order(erPathways[,"P.DE"]),]
    erPathways = erPathways[erPathways[,"Ont"]=="BP", ]
    processList = c(processList, list(erPathways))
}

# processListFinal = processList

processList[[1]][,"P.DE"] = p.adjust(processList[[1]][,"P.DE"])
processList[[2]][,"P.DE"] = p.adjust(processList[[2]][,"P.DE"])
processList[[3]][,"P.DE"] = p.adjust(processList[[3]][,"P.DE"])
processList[[4]][,"P.DE"] = p.adjust(processList[[4]][,"P.DE"])

processListFinal = processList

# allTerms = lapply(processListFinal, function(x) cbind(rownames(x),  x[,"Term"]) )
# lol = Reduce(merge,allTerms)

topGenes = as.character(top.features[1:4,1])
gene_Names_hg[which(gene_Names_hg %in% topGenes)]

p1 = data.frame(processListFinal[[1]][1:5, ], mod="38")
p2 = data.frame(processListFinal[[2]][1:5, ], mod="30")
p3 = data.frame(processListFinal[[3]][1:5, ], mod="11")
p4 = data.frame(processListFinal[[4]][c(1:5), ], mod="20")

## pathways enrichment plot
allPaths = rbind(p1, p2, p3, p4)
allPaths[,"P.DE"] = -1*log10(allPaths[,"P.DE"])
allPaths[,"Term"] = str_wrap(allPaths[,"Term"], width = 40)
allPaths[,"Term"] = factor(allPaths[,"Term"], levels=rev(allPaths[,"Term"]) )
allPaths[,"mod"] = factor(allPaths[,"mod"], levels=c("38", "30", "11", "20") )

pdf(file="Pathways_ML210_all.pdf", width=9, height=10)

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

# ggplot(topfeaturesOut, aes(x = AvgRank, y = mods)) + geom_point(size=1.5) + labs( x ="Average Rank", y = "Module") +
#     theme(axis.text = element_text(size=12), axis.text.x = element_text(size=10),
#           axis.title.x=element_text(size=14), axis.title.y=element_text(size=14)) + theme_classic()

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
# qq = myData[which(myData[,"AvgRank"] %in% minRanks[,2]), ]
# qq[,"mods"] = factor(qq[,"mods"], levels = unique(rev(qq[,"mods"])) )
# 
# qq = data.frame(qq, gene=top.features[match(qq[,1], top.features[,"AvgRank"]),1])
# qq[,"gene"] = factor(qq[,"gene"], levels = rev(top.features[,"FeatureName"]))

pdf(file="dotRanks_GPX4_all_paths.pdf", width=7, height=9)

# Horizontal version
ggplot(minRanks, aes(x=mods, y=AvgRank)) +
    geom_segment( aes(x=mods, xend=mods, y=0, yend=AvgRank), color="black",  size = 1.4)+
    geom_point( color="black", size=4) +
    theme_bw() + labs( y ="Average Rank", x = "Module") +
    theme(axis.text = element_text(size=13), axis.text.x = element_text(size=18),
          axis.title.x=element_text(size=24), axis.title.y=element_text(size=24)) + 
    coord_flip()

dev.off()

pdf(file="dotRanks_GPX4_all_paths_genes.pdf", width=7, height=9)
ggplot(minRanks, aes(x=gene, y=AvgRank)) +
    geom_segment( aes(x=gene, xend=gene, y=0, yend=AvgRank), color="black",  size = 1.4)+
    geom_point( color="black", size=4) +
    theme_bw() + labs( y ="Average Rank", x = "Module") +
    theme(axis.text = element_text(size=13), axis.text.x = element_text(size=18),
          axis.title.x=element_text(size=24), axis.title.y=element_text(size=24)) + 
    coord_flip()
dev.off()






# ## make AUC plots
# vemAUC = read.csv("vemAUCs.txt", header=TRUE)
# vemAUC = as.numeric(unlist(vemAUC))
# hist(vemAUC,20)


## pan-cancer PTX



## solid-cancer PTX

# DabrafAUC = 

# ggplot(topfeaturesOut, aes(x = AvgRank, y = mods)) + geom_point(size=1.5) + labs( x ="Average Rank", y = "Module") +
#     theme(axis.text = element_text(size=12), axis.text.x = element_text(size=10),
#           axis.title.x=element_text(size=14), axis.title.y=element_text(size=14)) + theme_classic()




