### File with TCGA barcodes ###
setwd("/Users/rbueno/Google Drive/ALL Files_TCGA-OV_FPKM/") # change directory to User's directory #
ids = read.table("TCGA-OV_FPKM.txt", header = F, sep = "", stringsAsFactors=F, skip = 1)

### Filenames for the TCGA-OV FPKM-UQ data ###
setwd("/Users/rbueno/Google Drive/ALL Files_TCGA-OV_FPKM/TCGA-OV_FPKM/") # change directory to the Folder that contains the RNA FPKM dataset downloaded from the TCGA for ovarian cancer #

filenames = list.files(".")
file.id = substr(filenames,1,36)

## Getting TCGA barcodes for expression dat ##

barcodes = sapply(file.id, function(x) grep(x, ids[,1]))

table(names(barcodes) == file.id) # sanity check 

tcga.codes = ids[barcodes,2] # TCGA barcodes for the TCGA Ovarian FPKM data, in the same order as in the folder TCGA-OV-Counts

### Reading in all expression dat files ###
all.patients = NULL
for(i in 1:length(filenames)){
	a = read.table(filenames[i], header = T, sep = "\t", row.names = 1)
	all.patients = cbind(all.patients, a[,1])
	
	print(i)
	
}

rownames(all.patients) = rownames(a)
colnames(all.patients) = tcga.codes

dim(all.patients) 
####################################################################################################
#Getting the gene symbols 

library("biomaRt") # package needed for getting gene names froms ensembl gene ids 

ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl") # Getting the ensembl info 

ens.ids = rownames(all.patients) # getting the ensembl IDs from the OV dataset
id.splt = strsplit(ens.ids, "[.]") #  removing the decimal point 
ID = sapply(id.splt, function(x) x[1]) #keeping only the numbers to the left of the decimal pt
rownames(all.patients) = ID

genes.id = getBM(attributes= c('ensembl_gene_id', 'hgnc_symbol'), filters = 'ensembl_gene_id', values = ID, mart = ensembl) # getting the ensembl gene ids and corresponding gene name 

dim(genes.id)
# 57086    2

genes = genes.id[!genes.id[,2] == "",] # removing ensembl gene ids with no gene symbols 
dim(genes)
# 36200 2

####################################################################################################
# Changing the rownames of the TCGA ovavarian datasets to gene names 

id.genes = sapply(genes[,1], function(x) grep(x, rownames(all.patients)))

allpatsID.dat = all.patients[id.genes,] # RNA seq data of all OV patients, ensemble IDs are rows 

table(rownames(allpatsID.dat) == genes[,1]) # sanity check to see if IDs match 


 
row.names(allpatsID.dat) <- genes[,2]

table(rownames(allpatsID.dat) == genes[,2]) # sanity check to see if genes match 

tcgaOV.dat = log(allpatsID.dat+1, 2)
save.image("FPKM-RNAprocessing.RData")

###### Grabbing BRCA2- and BRCA1- patients 
setwd("/Users/rbueno/Google Drive/TCGA ovarian RNA seq/TCGA_OvarianCancer_BRCAmuts/")

brca1 = read.table(file="mutation_table_BRCA1.tsv", header = F, sep = "\t", fill = T, stringsAsFactors = F, skip = 1)

brca1.id = brca1[,1] 
length(brca1.id) # 37

brca2 = read.table(file="mutation_table_BRCA2.tsv", header = F, sep = "\t", fill = T, stringsAsFactors = F, skip = 1)

brca2.id = brca2[,1]

length(brca2.id) # 35

brca.ids = c(brca1.id, brca2.id)
length(brca.ids) # 72

uniq.brca = unique(brca.ids) # unique BRCA patients (there were 3 patients with both BRCA1- and BRCA2-)

length(uniq.brca) # 69
## making BRCA expression matrix 
brca.index = sapply(uniq.brca, function(x) grep(x, colnames(tcgaOV.dat))) # index of all patients with BRCA1- and BRCA2- genes
length(unlist(brca.index)) #52
brca.index = unlist(brca.index)
name = names(brca.index)

brca1.name = which(name %in% brca1.id == T)

brca2.name = which(name %in% brca2.id == T)

brca.both = which(brca1.name %in% brca2.name == T) # 2 patients with both BRCA1 and BRCA2

BRCApatgenes.dat = tcgaOV.dat[,brca.index]

BRCApatgenes.dat = BRCApatgenes.dat[,-brca.both]#taking out patients with both BRCA1- and BRCA2- 

dim(BRCApatgenes.dat)


### BRCA1- samples 

brca1.dat = BRCApatgenes.dat[,unlist(sapply(brca1.id, function(x) grep(x, colnames(BRCApatgenes.dat))))]
dim(brca1.dat)

### BRCA2- samples 

brca2.dat = BRCApatgenes.dat[,unlist(sapply(brca2.id, function(x) grep(x, colnames(BRCApatgenes.dat))))]
dim(brca2.dat)

### Wild- type samples 

colnames.brca = colnames(BRCApatgenes.dat)
no.brca = sapply(colnames.brca, function(x) grep(x, colnames(all.patients)))

wt.dat = tcgaOV.dat[,-unlist(no.brca)]
dim(wt.dat)

### Sanity Check 
colnames(brca1.dat) %in% colnames(wt.dat)
colnames(brca2.dat) %in% colnames(wt.dat)

######################################################################################
### pathVar 
brca2.wt.log = cbind(brca2.dat, wt.dat)
dim(brca2.wt.log)

# filtering 
keep.genes = apply(brca2.wt.log, 1, function(x){sum(x>log(1,2))>=300})

brca2.wt.log = brca2.wt.log[keep.genes,]

hist(brca2.wt.log, main = "Distribution of all genes in the TCGA ovarian RNA-sequencing (FPKM) dataset", xlab = "Log2 + 1 transformed FPKM", cex.lab = 1.5, cex.axis = 1.25, cex.main = 2)

brca2.cat = colnames(brca2.wt.log)
brca2.cat[1:21] = "brca2"
brca2.cat[22:350] = "wt tumor"
grp.brca2 = c(rep(1, sum(brca2.cat == "brca2")), rep(2, sum(brca2.cat == "wt tumor"))) 
table(grp.brca2)
#grp.brca2
#  1   2 
# 21 329

#### PathVar functions
library(pathVar)
library(reshape)

diagnosticsVarPlotsTwoSample(brca2.wt.log, groups = as.factor(grp.brca2))

res.brca2wt.sd = pathVarTwoSamplesCont(brca2.wt.log, pways.kegg, groups = as.factor(grp.brca2), varStat = "sd", boot = 1000)
sig.brca2wt.sd = sigPway(res.brca2wt.sd, 0.05) 

path.names = names(sig.brca2wt.sd@genesInSigPways1)
length(path.names)

# table of P-values 
pways.cont.pval = cbind(res.brca2wt.sd@tablePway$PwayName, res.brca2wt.sd@tablePway$APval)

DNAdamResp=c("Base excision repair","Mismatch repair","Nucleotide excision repair", "Homologous recombination","Non-homologous end-joining","Cell cycle")

DNAdam.tableCont = pways.cont.pval[sapply(DNAdamResp, function(x) grep(x, pways.cont.pval[,1])),]


### Discrete Gene Count differences in pathways between BRCA2- and Wt-tumor ###
disc.brca2wt.sd = pathVarTwoSamplesDisc(brca2.wt.log, pways.kegg, groups = as.factor(grp.brca2), varStat = "sd", test = "exact")

sig.brca2wt.disc = sigPway(disc.brca2wt.sd, 0.05) 
# Names of the significant pathways
disc.pathnames = names(sig.brca2wt.disc@genesInSigPways1)
length(disc.pathnames)

pways.disc.pval = cbind(disc.brca2wt.sd@tablePway$PwayName, disc.brca2wt.sd@tablePway$APval)


DNAdam.tableDisc = pways.disc.pval[sapply(DNAdamResp, function(x) grep(x, pways.disc.pval[,1])),]





# Mean as the statistic #

res.brca2wt.mean = pathVarTwoSamplesCont(brca2.wt.log, pways.kegg, groups = as.factor(grp.brca2), varStat = "mean", boot = 1000)
sig.brca2wt.mean = sigPway(res.brca2wt.mean, 0.05) # No pathways 
disc.brca2wt.mean = pathVarTwoSamplesDisc(brca2.wt.log, pways.kegg, groups = as.factor(grp.brca2), varStat = "mean", test = "exact")
sig.brca2wt.MeanDisc = sigPway(disc.brca2wt.mean, 0.05) # No significant pathways 



##############################################################################
## Modified pathVar plot functions ##
plotTwoSamplesDisc <- function(pvalue_results, pathway, sig, sizeXtick=NULL,sizeYtick=NULL, sizeLab = NULL) {
  mp <- pathway
  # If the name of the pathway is two long it will cut it into two lines in the plot.
  if (nchar(mp) > 39) {
    mp <- unlist(strsplit(mp, " "))
    if (is.na((which(cumsum(nchar(mp)) > 34)[1]))) {
      l <- length(mp) - 1
    } else {
      l <- which(cumsum(nchar(mp)) > 34)[1] - 1
    }
    mp_1 <- paste(mp[1:l], collapse = " ")
    mp_2 <- paste(mp[(l + 1):length(mp)], collapse = " ")
    mp <- paste(mp_1, "\n", mp_2)
  }
  path1 <- pvalue_results@pwayCounts1[[pathway]]
  path2 <- pvalue_results@pwayCounts2[[pathway]]
  # if we have the result from the sigOneSample, it will be included in the plot
  if (!is.null(sig)) {
    category <- sig@sigCatPerPway
  } else {
    category <- NULL
  }
  # data frame for the reference counts
  path1 <- as.data.frame(path1)
  colnames(path1) <- c("Cluster", "Number_of_genes")
  # data frame for the pathway distribution
  path2 <- as.data.frame(path2)
  colnames(path2) <- c("Cluster", "Number_of_genes")
  yLimMax <- max(path1[, 2], path2[, 2])
  # plot of the reference counts
  plotPath1 <- ggplot(path1, aes(x = Cluster, y = Number_of_genes, fill = Cluster)) + geom_bar(stat = "identity",
                                                                                               color = "black", fill = "darkorchid1") + ylab("Number of genes") + theme(legend.position = "none",axis.text.x = element_text(size=sizeXtick),axis.text.y = element_text(size=sizeYtick), axis.title = element_text(size = sizeLab)) +
    ggtitle(paste(mp, "\n Group 1")) + xlab("") + ylim(0, yLimMax)
  # Plot of the pwathway counts
  plotPath2 <- ggplot(path2, aes(x = Cluster, y = Number_of_genes, fill = Cluster)) + geom_bar(stat = "identity",
                                                                                               color = "black", fill = "steelblue1") + ylab("Number of genes") + theme(legend.position = "none",axis.text.x = element_text(size=sizeXtick),axis.text.y = element_text(size=sizeYtick), axis.title = element_text(size = sizeLab)) +
    ggtitle(paste(mp, "\n Group 2")) + xlab("") + ylim(0, yLimMax)
  # If the pathway is one of the significant ones, the title will be in red. and the categories, if any, we be highlighted with a red star.
  if (pathway %in% names(category)) {
    sigCat <- category[[pathway]]
    plotPathway1 <- plotPath1 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                      panel.background = element_blank(), axis.line = element_line(colour = "black"),
                                      plot.title = element_text(colour = "red"))
    plotPathway2 <- plotPath2 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                      panel.background = element_blank(), axis.line = element_line(colour = "black"),
                                      plot.title = element_text(colour = "red"))
    if (sum(sigCat) > 0) {
      plotPathway1 <- plotPathway1 + annotate("text", x = sigCat, y = path1[sigCat +
                                                                              0.1, 2], label = rep("*", length(sigCat)), color = "red", size = 15)
      plotPathway2 <- plotPathway2 + annotate("text", x = sigCat, y = path2[sigCat +
                                                                              0.1, 2], label = rep("*", length(sigCat)), color = "red", size = 15)
    }
  } else {
    plotPathway1 <- plotPath1 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                      panel.background = element_blank(), axis.line = element_line(colour = "black"))
    plotPathway2 <- plotPath2 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                      panel.background = element_blank(), axis.line = element_line(colour = "black"))
  }
  # plot the reference and pathway counts side by side
  grid.arrange(arrangeGrob(plotPathway1, plotPathway2, nrow = 1))
}

############################################################################################
plotTwoSamplesCont <- function(pvalue_results, pathway, sig, sizeXtick=NULL,sizeYtick=NULL, sizeLab = NULL) {
  mp <- pathway
  # If the name of the pathway is two long it will cut it into two lines in the plot.
  if (nchar(mp) > 39) {
    mp <- unlist(strsplit(mp, " "))
    if (is.na((which(cumsum(nchar(mp)) > 34)[1]))) {
      l <- length(mp) - 1
    } else {
      l <- which(cumsum(nchar(mp)) > 34)[1] - 1
    }
    mp_1 <- paste(mp[1:l], collapse = " ")
    mp_2 <- paste(mp[(l + 1):length(mp)], collapse = " ")
    mp <- paste(mp_1, "\n", mp_2)
  }
  xtab <- pvalue_results@tablePway
  # If the number of genes of the pathway is less than 3, it is not possible to draw a density and it will return an empty plot with this message.
  if (xtab[PwayName == pathway, NumOfGenesFromDataSetInPathway] < 3) {
    df <- data.frame()
    plotPway <- ggplot(df) + geom_point() + xlim(0, 3) + ylim(0, 1) + theme(panel.grid.major = element_blank(), 
                                                                            panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
      ylab("Density") + theme(legend.position = "none", axis.text.x = element_text(size=sizeXtick),axis.text.y = element_text(size=sizeYtick), axis.title = element_text(size = sizeLab)) + xlab(varStat) + ggtitle(mp) + 
      geom_text(aes(1.5, 0.5, label = "Less than 3 genes, not possible to draw a density.", 
                    col = "red"))
  } else {
    var_1 <- pvalue_results@var1
    var_2 <- pvalue_results@var2
    varStat <- pvalue_results@varStat
    genes <- pvalue_results@genesInPway[[pathway]]
    grp_1 <- as.vector(var_1[genes])
    grp_2 <- as.vector(var_2[genes])
    df <- data.frame(cbind(c(grp_1, grp_2), c(rep("Group 1", length(grp_1)), rep("Group 2", 
                                                                                 length(grp_2)))))
    df.m <- melt(df)
    df.m[, 1] <- as.numeric(as.character(df.m[, 1]))
    colnames(df.m) <- c("value", "group")
    color <- c("steelblue1", "grey")
    # Plot of the two densities (one for each group) of the variability of the genes inside the pathway.
    d <- ggplot(na.omit(df.m)) + geom_density(aes(x = value, y = ..density.., colour = group, 
                                                  fill = group), alpha = 0.3) + scale_colour_manual(values = color) + scale_fill_manual(values = color) + 
      theme(legend.position = "bottom", legend.title = element_blank(), axis.text.x = element_text(size=sizeXtick),axis.text.y = element_text(size=sizeYtick), axis.title = element_text(size = sizeLab)) + ylab("Density") + 
      xlab(varStat) + ggtitle(mp)
    if (!is.null(sig)) {
      significant <- names(sig@genesInSigPways1)
    } else {
      significant <- NULL
    }
    # If we included the results of sigTwoSamplesCont, it will verify if the pathway is one of them and if yes the title will be printed in red.
    if (pathway %in% significant) {
      plotPathway <- d + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                               panel.background = element_blank(), axis.line = element_line(colour = "black"), 
                               plot.title = element_text(colour = "red"))
    } else {
      plotPathway <- d + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                               panel.background = element_blank(), axis.line = element_line(colour = "black"))
    }
  }
  plot(plotPathway)
}


##############################################################################
## Plots of pathways involved in the DNA damage response ##
plotTwoSamplesCont(res.brca2wt.sd, "Base excision repair", sig.brca2wt.sd, sizeXtick = 22, sizeYtick = 22, sizeLab = 24)
plotTwoSamplesCont(res.brca2wt.sd, "Mismatch repair", sig.brca2wt.sd, sizeXtick = 22, sizeYtick = 22, sizeLab = 24)
plotTwoSamplesCont(res.brca2wt.sd, "Nucleotide excision repair", sig.brca2wt.sd, sizeXtick = 22, sizeYtick = 22, sizeLab = 24)
plotTwoSamplesCont(res.brca2wt.sd, "Homologous recombination", sig.brca2wt.sd, sizeXtick = 22, sizeYtick = 22, sizeLab = 24)
plotTwoSamplesCont(res.brca2wt.sd, "Non-homologous end-joining", sig.brca2wt.sd, sizeXtick = 22, sizeYtick = 22, sizeLab = 24)

# other pathway involved in DNA damage response that's not a repair pathway 
plotTwoSamplesCont(res.brca2wt.sd, "Cell cycle", sig.brca2wt.sd, sizeXtick = 22, sizeYtick = 22, sizeLab = 24)

## Plots of repair pathways with significant differences in gene counts ##
plotTwoSamplesDisc(disc.brca2wt.sd, "Base excision repair", sig.brca2wt.disc, sizeXtick = 26, sizeYtick = 26, sizeLab = 30)
plotTwoSamplesDisc(disc.brca2wt.sd, "Mismatch repair", sig.brca2wt.disc, sizeXtick = 26, sizeYtick = 26, sizeLab = 30)
plotTwoSamplesDisc(disc.brca2wt.sd, "Nucleotide excision repair", sig.brca2wt.disc, sizeXtick = 26, sizeYtick = 26, sizeLab = 30)
plotTwoSamplesDisc(disc.brca2wt.sd, "Homologous recombination", sig.brca2wt.disc,sizeXtick = 26, sizeYtick = 26, sizeLab = 30)
plotTwoSamplesDisc(disc.brca2wt.sd, "Non-homologous end-joining", sig.brca2wt.disc,sizeXtick = 26, sizeYtick = 26, sizeLab = 30) # No significant difference between clusters 

# other pathway involved in DNA damage response that's not a repair pathway 
plotTwoSamplesDisc(disc.brca2wt.sd, "Cell cycle", sig.brca2wt.disc, sizeXtick = 26, sizeYtick = 26, sizeLab = 30) # focus on cell cycle pathway  

###################

## Grabbing the genes associated with the repair pathways ##

# Base excision repair #
base.repair = cbind(disc.brca2wt.sd@genesInPway1$"Base excision repair", disc.brca2wt.sd@genesInPway2$"Base excision repair")
colnames(base.repair) = c("brca2-", "wt-tumor")

baserep.genes = rownames(base.repair)
length(baserep.genes) # 33 genes 

table(base.repair[,1]) # Genes in BRCA2- (group1)


table(base.repair[,2]) # Genes in Wt-tumor (group2)


# Mismatch repair #
mismatch.repair = cbind(disc.brca2wt.sd@genesInPway1$"Mismatch repair", disc.brca2wt.sd@genesInPway2$"Mismatch repair")
colnames(mismatch.repair) = c("brca2-", "wt-tumor")

mismatch.genes = rownames(mismatch.repair)
length(mismatch.genes) # 23 genes 

table(mismatch.repair[,1]) # Genes in BRCA2- (group1)

table(mismatch.repair[,2]) # Genes in Wt-tumor (group2)

# Nucleotide excision repair #
nucl.repair = cbind(disc.brca2wt.sd@genesInPway1$"Nucleotide excision repair", disc.brca2wt.sd@genesInPway2$"Nucleotide excision repair")
colnames(nucl.repair) = c("brca2-", "wt-tumor")

nucl.genes = rownames(nucl.repair)
length(nucl.genes) # 45 genes 

table(nucl.repair[,1]) # Genes in BRCA2- (group1)

table(nucl.repair[,2]) # Genes in Wt-tumor (group2)

# Homologous Recombination #
hr.repair = cbind(disc.brca2wt.sd@genesInPway1$"Homologous recombination", disc.brca2wt.sd@genesInPway2$"Homologous recombination")
colnames(hr.repair) = c("brca2-", "wt-tumor")

hr.genes = rownames(hr.repair)
length(hr.genes) # 26 genes 

table(hr.repair[,1]) # Genes in BRCA2- (group1)

table(hr.repair[,2]) # Genes in Wt-tumor (group2)

# Non-homologous end joining #
nhej.repair = cbind(disc.brca2wt.sd@genesInPway1$"Non-homologous end-joining", disc.brca2wt.sd@genesInPway2$"Non-homologous end-joining") 
colnames(nhej.repair) = c("brca2-", "wt-tumor")

nhej.genes = rownames(nhej.repair)
length(nhej.genes) # 11 genes 

table(nhej.repair[,1])

table(nhej.repair[,2])

# cell cycle #  
cellcyc = cbind(disc.brca2wt.sd@genesInPway1$"Cell cycle", disc.brca2wt.sd@genesInPway2$"Cell cycle")
colnames(cellcyc) = c("brca2-", "wt-tumor")
cellcyc.genes = rownames(cellcyc)
length(cellcyc.genes) # 123 genes 

### Grabbing all genes from the five DNA repair pathways ###

all.genes = c(baserep.genes, mismatch.genes, nucl.genes, hr.genes, nhej.genes, cellcyc.genes)
length(all.genes) # 261 genes 
all.genes = sort(all.genes)

uniq.allgenes = unique(all.genes)
length(uniq.allgenes) # 219 geenes 
uniq.allgenes = sort(uniq.allgenes)


###############################################################################
# Grabbing genes that transition from High/Medium Variability to Low Variability #
# Non-homologous End Joining #
dim(nhej.repair)
medVar.nhej = nhej.repair[which(nhej.repair[,2]!=1),]
dim(medVar.nhej)

lowVar.nhej = medVar.nhej[which(medVar.nhej[,1] == 1),]
dim(lowVar.nhej)

sliNhej.brca2 = rownames(lowVar.nhej)
length(sliNhej.brca2)

# Nucleotide Excision Repair #
dim(nucl.repair)

medVar.nucl = nucl.repair[which(nucl.repair[,2]!=1),]
dim(medVar.nucl)

lowVar.nuclrepair = medVar.nucl[medVar.nucl[,1]==1,] 
dim(lowVar.nuclrepair)

sliNuclrepair.brca2 = rownames(lowVar.nuclrepair)
length(sliNuclrepair.brca2)

# Mismatch Repair pathway #
dim(mismatch.repair)  

medVar.mismatch = mismatch.repair[which(mismatch.repair[,2]!=1),]
dim(medVar.mismatch)

lowVar.mismatch = medVar.mismatch[which(medVar.mismatch[,1] == 1),]
dim(lowVar.mismatch)

sliMismatch.brca2 = rownames(lowVar.mismatch)
length(sliMismatch.brca2)

# Homologous Recombination Pathway # 
dim(hr.repair)

medVar.hr = hr.repair[which(hr.repair[,2]!=1),]
dim(medVar.hr)

lowVar.hr = medVar.hr[which(medVar.hr[,1] == 1),]
dim(lowVar.hr)

sliHr.brca2 = rownames(lowVar.hr)
length(sliHr.brca2)

# Base Repair Excision #
dim(base.repair)

medVar.base = base.repair[which(base.repair[,2]!=1),]
dim(medVar.base)

lowVar.base = medVar.base[which(medVar.base[,1] == 1),]
dim(lowVar.base)

sliBase.brca2 = rownames(lowVar.base)
length(sliBase.brca2)

# Cell Cycle #
dim(cellcyc)

medVar.cellcyc = cellcyc[which(cellcyc[,2]!=1),]
dim(medVar.cellcyc)

lowVar.cellcyc = medVar.cellcyc[which(medVar.cellcyc[,1] == 1),]
dim(lowVar.cellcyc)

sliCellcyc.brca2 = rownames(lowVar.cellcyc)
length(sliCellcyc.brca2)

sli.genes = unique(c(sliCellcyc.brca2, sliBase.brca2, sliHr.brca2, sliMismatch.brca2, sliNuclrepair.brca2, sliNhej.brca2))
length(sli.genes)

######################################################################################
# HR.repair genes overlap 
hr.base = intersect(hr.genes, baserep.genes)
hr.mismatch = intersect(hr.genes, mismatch.genes)
hr.nhej = intersect(hr.genes, nhej.genes)
hr.nucl = intersect(hr.genes, nucl.genes)
hr.cellcyc = intersect(hr.genes, cellcyc.genes) # 0 genes

all.hr.overlap = c(hr.base,hr.mismatch,hr.nhej,hr.nucl)
all.hr.overlap = unique(all.hr.overlap)

##############################################################################
### Intersecting unique genes with SLI gene from SynLethDB database ###
setwd("/Users/rbueno/Google Drive/TCGA ovarian RNA seq/")


synlethdb = read.table("synLethDB-BRCA2.txt", header = F, sep = "\t", fill = T, stringsAsFactors = F)
names = synlethdb[1,]
synlethdb = synlethdb[-1,]
dim(synlethdb) # 42 9 

db.genes = unique(c(synlethdb[,3], synlethdb[,1]))

db.genes = db.genes[-which(db.genes == " ")]
db.genes = db.genes[-which(db.genes == "")]
db.genes = db.genes[-which(db.genes%in%"BRCA2" == T)]
length(db.genes)

lowvar.slgenes = intersect(sli.genes, synlethdb[,3])
lowvar.slgenes2 = intersect(sli.genes, synlethdb[,1])

sli.db.pathVar = c(lowvar.slgenes, lowvar.slgenes2) 
length(sli.db.pathVar)# 
sli.db.pathVar

slipathVar.notdb = sli.genes[which(sli.genes %in% sli.db.pathVar==F)]
length(slipathVar.notdb)

sliDB.NotpathVar = setdiff(db.genes, sli.db.pathVar)
length(sliDB.NotpathVar)

NotSLI.pathVar.db = setdiff(uniq.allgenes, sli.genes)
NotSLI.pathVar.db = setdiff(NotSLI.pathVar.db, db.genes)
length(NotSLI.pathVar.db)

# Fisher's Exact Test #
fisher.matrix = matrix(nrow = 2, ncol = 2)
colnames(fisher.matrix) = c("pathVar+", "pathVar-")
rownames(fisher.matrix) = c("SynLethDB+", "SynLethDB-")
fisher.matrix[1,] = c(length(sli.db.pathVar), length(sliDB.NotpathVar))
fisher.matrix[2,] = c(length(slipathVar.notdb), length(NotSLI.pathVar.db))

fisher.test(fisher.matrix, alternative = "greater")
############################################################
## KS analysis for SD values for SLI genes in SynLethDB and pathVar
allgenesSDvals.brca2 = disc.brca2wt.sd@var1
allgenesSDvals.wt = disc.brca2wt.sd@var2

synlethdb.brca2 = disc.brca2wt.sd@var1[which(names(disc.brca2wt.sd@var1)%in%db.genes==T)]
synlethdb.wt = disc.brca2wt.sd@var2[which(names(disc.brca2wt.sd@var2)%in%db.genes == T)]
plot(density(synlethdb.brca2), col = "red", xlab = "sd values", main = "Density plot of the sd values for SynLethDB SLI genes")
lines(density(synlethdb.wt))
legend("topright", "P-value = 0.007", cex = 1.25)


pathVar.brca2 = disc.brca2wt.sd@var1[which(names(disc.brca2wt.sd@var1)%in%sli.genes==T)]
pathVar.wt = disc.brca2wt.sd@var2[which(names(disc.brca2wt.sd@var2)%in%sli.genes==T)]
plot(density(pathVar.brca2), col = "blue", xlab = "sd values", main = "Density plot of the sd values for pathVar SLI genes", xlim = c(0.1,0.9))
lines(density(pathVar.wt))
legend("topright", "P-value = 2.2e-16", cex = 1.5)


plot(density(synlethdb.brca2), col = "red", main = "Density plot of sd values of pathVar vs. SynLethDB SLI genes", ylim = c(0,9), xlab = "sd values", cex.axis = 1.25, cex.main = 1.25, cex.lab = 1.25)
lines(density(pathVar.wt), col = "blue")
legend("topright", "P-value = 0.42", cex = 1.5)

par(mfrow=c(3,1))
# SynLethDB 
plot(density(synlethdb.brca2), col = "red", xlab = "sd values", main = "Density plot of the sd values for SynLethDB SLI genes", cex.axis = 1.25, cex.main = 2, cex.lab = 1.25, lwd = 3)
lines(density(synlethdb.wt), lwd = 3)
legend("topright", "P-value = 0.01", cex = 1.75)
# pathVar 
plot(density(pathVar.brca2), col = "blue", xlab = "sd values", main = "Density plot of the sd values for pathVar SLI genes", xlim = c(0.1,0.9), cex.axis = 1.25, cex.main = 2, cex.lab = 1.25, lwd = 3)
lines(density(pathVar.wt), lwd = 3)
legend("topright", "P-value = <2.2e-16", cex = 1.75)
# SynLethDB vs. pathVar
plot(density(synlethdb.brca2), col = "red", main = "Density plot of sd values of pathVar vs. SynLethDB SLI genes", ylim = c(0,9), xlab = "sd values", cex.axis = 1.25, cex.main = 2, cex.lab = 1.25, lwd = 3)
lines(density(pathVar.brca2), col = "blue", lwd = 3)
legend("topright", "P-value = 7.385e-09", cex = 1.75)

plot(ecdf(synlethdb.brca2))
plot(ecdf(pathVar.), add = T)

plot(ecdf(synlethdb.brca2))
plot(ecdf(synlethdb.wt), add = T, col =  "red")

##
ks.boot(pathVar.brca2,synlethdb.brca2) 
ks.test(synlethdb.brca2,synlethdb.wt, exact = F)
ks.test(pathVar.brca2,pathVar.wt, exact = F)
############################################################
# T-test of SLI genes #
sli.exprs = which(rownames(brca2.wt.log)%in%sli.genes == T)
sli.exDat = brca2.wt.log[sli.exprs,]
dim(sli.exDat)

tTest = function(x){
	test = t.test(x[1:21],x[22:350])
	pval = test$p.value
	return(pval)
}


mean.brca2 = apply(sli.exDat[,1:21],1,mean)
mean.wt = apply(sli.exDat[,22:350], 1, mean)
fchg = mean.brca2 - mean.wt

sum(fchg>=0)
sum(fchg<0)

tpvals = apply(sli.exDat, 1, tTest)
length(tpvals)
sum(tpvals<0.05)
tpvals[which(tpvals<0.05)]

adj.tpvals = p.adjust(tpvals, "BH")
sum(adj.tpvals<0.05)
adj.tpvals[which(adj.tpvals<0.05)]

log.adjP = -log(adj.tpvals,10)
plot(fchg,log.adjP, xlab = "Log2(BRCA2-/BRCA2+)", ylab = "-Log10(FDR adj. P-value)", main = "Volcano Plot")
abline(h = 1.3, col = "red")
text(fchg,log.adjP, labels = names(fchg), cex = 0.6, pos = 1)

#############################################################
load("FPKM-pathVar.RData")
## Network Analyses using iGraph ##
library(igraph)
setwd("/Users/rbueno/Google Drive/TCGA ovarian RNA seq/")

length(sli.genes) # ALL genes that are defined to be SLI using pathVar

ddr.genes = uniq.allgenes[which(uniq.allgenes %in% sli.genes == F)] # genes that are not considered to be SLI but function in the DNA damage response 

# Protein-Protein interactions from the Human Reference Protein database #
hrpd = read.table("BINARY_PROTEIN_PROTEIN_INTERACTIONS.txt", header = F, sep = "\t", stringsAsFactors = F)

hrpd.edgeNode = hrpd[,c(1,2,4,5)]
dim(hrpd.edgeNode)
head(hrpd.edgeNode)

sli.Node = hrpd.edgeNode[hrpd[,1] %in% sli.genes,]
dim(sli.Node)

sli.Node2 = hrpd.edgeNode[hrpd[,4] %in% sli.genes,]
dim(sli.Node2)
sli.Nodes2 = sli.Node2[,c(3,4,1,2)]


sli.Nodes = rbind(sli.Node,sli.Node2)
dim(sli.Nodes)

sli.ddr.net = sli.Nodes[sli.Nodes[,3]%in%ddr.genes,]
sli.ddr.net = sli.ddr.net[order(sli.ddr.net[,1]),]
dim(sli.ddr.net)

## Setting up igraph ##
# Links #
LINK = sli.ddr.net[,c(2,4)]
dim(LINK)

# Nodes #
NODE = unique(c(sli.ddr.net[,1], sli.ddr.net[,3])) # grabbing the unique nodes 
node.id1 = sapply(NODE, function(x) grep(x, sli.ddr.net[,1])) 
node.id2 = sapply(NODE, function(x) grep(x, sli.ddr.net[,3]))
index.id1 = lapply(node.id1, function(x){if(length(x)!= 1){return(x[1])} else{return(x)}}) # finding index of each unique gene 

index.id1.na = index.id1[!is.na(unlist(index.id1))] #removing NAs

index.id2 =  lapply(node.id2, function(x){if(length(x)!= 1){return(x[1])}else{return(x)}}) # finding index of each unique gene

index.id2.na = index.id2[!is.na(unlist(index.id2))] # removing NAs


id1 = sli.ddr.net[unlist(index.id1.na),c(1,2)]
id2 = sli.ddr.net[unlist(index.id2.na),c(3,4)]
colnames(id1) = colnames(id2)

NODE.id = rbind(id1, id2)
NODE.id = NODE.id[,c(2,1)]

dim(NODE.id) # data of all unique Nodes and the string ID associated that node 

# Creating shapes - Circle = SLI, Square = gene in DDR #
NODE.sli = cbind(NODE.id, NODE.id[,2]%in%sli.genes)
shape.func = function(x){
	if(x == T){
		return("circle")
	}else{
		return("square")
	}
}

shape = sapply(NODE.sli[,3], shape.func )

NODE.shape = NODE.sli[,c(1,2)]
NODE.shape = cbind(NODE.shape, shape)

# Creating size vector to represent number of edges - large node represent large number of edges #

SZ1 = table(sli.ddr.net[,1])
SZ2 = table(sli.ddr.net[,3])
sizes = c(SZ1, SZ2)
sizes = as.matrix(sizes)
sizes.index = sapply(NODE.shape[,2], function(x) grep(x, rownames(sizes)))

sz.func = function(x){
	if(x[1]!=1){
		return(x[1])
	}else{
		return(x)
	}
}

sz.ind = lapply(sizes.index, sz.func)
Sizes = as.matrix(sizes[unlist(sz.ind),])
NODE.shape.sz = cbind(NODE.shape, Sizes)
dim(NODE.shape.sz)

# Functional Annotations #
setwd("/Users/rbueno/Google Drive/ALL Files_TCGA-OV_FPKM/")

ipa = read.table("SLI_annotation.txt", header = T, sep = "\t", stringsAsFactors = F, fill = T)

ipa.func = ipa[,c(1,5)]

func.index = sapply(NODE.shape.sz[,2], function(x) grep(x, ipa.func[,1]))
func = function(x){
	if(x[1]!=1){
		return(x[1])
	}else{
		return(x)
	}
}
func.ind = lapply(func.index, func)

type.func = ipa.func[unlist(func.ind),]

all.equal(type.func[,1], NODE.shape.sz[,2])

NODE.shapeSzFunc = cbind(NODE.shape.sz, type.func[,2])

genefunc = function(x){ if(x == "kinase"){return("orange")}else{if(x == "enzyme"){return("light blue")}else{if(x == "transcription regulator"){return("purple")}else{if(x == "other"){return("white")}else{return("green")}}}}}

func.color = sapply(NODE.shapeSzFunc[,5], genefunc)
NODE.shapeSzFuncCol = cbind(NODE.shapeSzFunc,func.color)

# Identifying validated SLIs (red) and Non-validated SLIs (pink) determined by pathVar #
head(NODE.shapeSzFuncCol) 
NODE.shapeSzFuncCol[,6] = as.character(NODE.shapeSzFuncCol[,6])


NODE.shapeSzFuncCol[which(NODE.shapeSzFuncCol[,2]%in%sli.db.pathVar == T),6] <- "red" 

NODE.shapeSzFuncCol[which(NODE.shapeSzFuncCol[,2]%in%slipathVar.notdb == T),6] <- "pink"


# Network figure using igraph #
NET = graph_from_data_frame(d = LINK, vertices = NODE.shape.sz, directed = F)

set.seed(2)
plot(NET,vertex.label = NODE.shape.sz[,2], vertex.size = NODE.shape.sz[,4]*3, vertex.label.cex = 0.8,edge.arrow.size = .5, layout = layout_nicely(NET), vertex.shape = as.character(NODE.shape.sz[,3]), vertex.color = NODE.shapeSzFuncCol[,6])

legend("topright", c("Validated BRCA2-SLI genes", "Non-validated BRCA2-SLI genes", "Enzyme", "other", "Transcription Regulator", "Phosphatase", "Kinase"), pch=21, pt.bg=unique(NODE.shapeSzFuncCol[,6]), pt.cex= 2, cex= 1, bty="n", ncol=1)

# SLI nodes #
SLI.nodes = NODE.shapeSzFuncCol[which(NODE.shapeSzFuncCol[,3] == "circle"),c(2,4,5,6)]
colnames(SLI.nodes) = c("Protein"," Number of Connections", "Function", "Graph Color")
write.table(SLI.nodes, file = "SLItable.txt", quote = F, col.names = T, sep = "\t", row.names = F)

# Graphing EP300 connection only #
# Node and Edges #
edgelist <- get.edgelist(NET)

ep300.node = NODE.shapeSzFuncCol[grep("EP300", NODE.shapeSzFuncCol[,2]),]

ep300link = edgelist[grep(ep300.node[,1],edgelist[,1]),]

DNAdam.nodes = NODE.shapeSzFuncCol[sapply(ep300link[,2], function(x) grep(x, NODE.shapeSzFuncCol[,1])),]

ep300DNA.nodes = rbind(ep300.node, DNAdam.nodes)

# Network figure using igraph #
ep300.NET = graph_from_data_frame(d = ep300link, vertices = ep300DNA.nodes[,1:2], directed = F)

plot(ep300.NET, vertex.label = ep300DNA.nodes[,2], vertex.size = ep300DNA.nodes[,4]*5, vertex.label.cex = 1, edge.arrow.size = .5, vertex.shape = as.character(ep300DNA.nodes[,3]), vertex.color =  ep300DNA.nodes[,6], layout = layout_nicely(ep300.NET))
legend("topleft", c("Non-validated BRCA2-SLI genes","Transcription Regulator","Enzyme", "Kinase"), pch=21, pt.bg=unique(ep300DNA.nodes[,6]), pt.cex= 2, cex= 1, bty="n", ncol=1)

# p53 #
NODE.shapeSzFuncCol[grep("TP53", NODE.shapeSzFuncCol[,2]),]
p53.edge = edgelist[grep(1859,edgelist[,2]),]
NODE.shapeSzFuncCol[sapply(p53.edge[,1], function(x) grep(x, NODE.shapeSzFuncCol[,1])),]

#############################################################
#setwd("/Users/rbueno/Google Drive/ALL Files_TCGA-OV_FPKM/")
#save.image("FPKM-pathVar.RData")
#############################################################
