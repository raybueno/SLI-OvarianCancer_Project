load("exDat.wt-BRCA.RData")

### Pre-processing/filtertering the data using functions from edgeR ###
library("edgeR")

x = cbind(brca1.dat, brca2.dat, wt.dat)

ov.cat = colnames(x)
ov.cat[1:29] = "brca1"
ov.cat[30:50] = "brca2"
ov.cat[51:379] = "wt"

table(ov.cat)
#brca1 brca2    wt 
#   29    21   329 

group = as.factor(ov.cat)

y = DGEList(counts = x, group = group)

### Filtering lowly expressed genes 
keep = rowSums(cpm(y)> 15) >= 15 # filtering out lowly expressed genes

y = y[keep, , keep.lib.sizes=FALSE]# recalculating libary sizes

### Normalization 

y = calcNormFactors(y)

### estimating dispersions 

y = estimateDisp(y)

### organizing the data in a dataframe ###

## BRCA1- and wt in a single dataframe ##
brca1.wt = y$counts[,c(1:29,51:379)]
dim(brca1.wt)
#[1] 12743   358

brca1.wt.log = log10(brca1.wt+1) # log10 transformed + 1

# sanity check to see what the data looks like 
hist(brca1.wt.log)
plot(density(brca1.wt.log))


## designing groups for pathVar analyses ##
brca1.cat = colnames(brca1.wt)
brca1.cat[1:29] = "brca1"
brca1.cat[30:358] = "wt tumor"
grp.brca1 = c(rep(1, sum(brca1.cat == "brca1")), rep(2, sum(brca1.cat == "wt tumor")))
table(grp.brca1)
#grp.brca1
#  1   2 
# 29 329

## BRCA2- and wt in a single dataframe ##
brca2.wt = y$counts[,c(30:50,51:379)] 
dim(brca2.wt)
#[1] 12743   350

brca2.wt.log = log10(brca2.wt+1) # log10 tranformed + 1

# sanity check to see what the data looks like 
hist(brca2.wt.log)
plot(density(brca2.wt.log), main = "Distribution of all genes in the TCGA ovarian cancer RNA-seq Data", xlab = "Read counts of all genes (log10 transformed + 1)")


#############################################################################################################################
## pathVar analyses to determine differences in variability between BRCA2- ovarian cancer patients and BRCA2+ ovarian cancer patients 

library(pathVar)
## designing groups for pathVar analyses ##

brca1.cat = colnames(brca1.wt)
brca1.cat[1:29] = "brca1"
brca1.cat[30:358] = "wt tumor"
grp.brca1 = c(rep(1, sum(brca1.cat == "brca1")), rep(2, sum(brca1.cat == "wt tumor")))
table(grp.brca1)
#grp.brca1
#  1   2 
# 29 329

brca2.cat = colnames(brca2.wt.log)
brca2.cat[1:21] = "brca2"
brca2.cat[22:350] = "wt tumor"
grp.brca2 = c(rep(1, sum(brca2.cat == "brca2")), rep(2, sum(brca2.cat == "wt tumor")))
table(grp.brca2)
#grp.brca2
#  1   2 
# 21 329

# BRCA1
diagnosticsVarPlotsTwoSample(brca1.wt.log , groups = as.factor(grp.brca1))
# SD(R = -0.568), MAD(R = -0.531), CV(R = -0.68)

# BRCA2
diagnosticsVarPlotsTwoSample(brca2.wt.log , groups = as.factor(grp.brca2))
# SD(R = -0.566), MAD(R = -0.529), CV(R = -0.683)

## Distribution Difference in Mean between BRCA1- vs. Wt-tumor ##
res.brca1wt.mean = pathVarTwoSamplesCont(brca1.wt.log, pways.kegg, groups = as.factor(grp.brca1), varStat = "mean", boot = 1000)
sig.brca1wt.mean = sigPway(res.brca1wt.mean, 0.05) # No significant pathways  

## Distribution Difference in SD between BRCA1- vs. Wt-tumor
res.brca1wt.sd = pathVarTwoSamplesCont(brca1.wt.log, pways.kegg, groups = as.factor(grp.brca1), varStat = "sd", boot = 1000)
sig.brca1wt.sd = sigPway(res.brca1wt.sd, 0.05) # No significant pathways 

## Distribution Difference in Mean between BRCA2- vs. Wt-tumor ##
res.brca2.wt.mean = pathVarTwoSamplesCont(brca2.wt.log, pways.kegg, groups = as.factor(grp.brca2), varStat = "mean", boot = 1000)
sig.brca2wt.mean = sigPway(res.brca2.wt.mean, 0.05) # No significant pathways

## Distribution Difference in SD between BRCA2= vs. Wt-tumor ##
res.brca2wt.sd = pathVarTwoSamplesCont(brca2.wt.log, pways.kegg, groups = as.factor(grp.brca2), varStat = "sd", boot = 1000)

sig.brca2wt.sd = sigPway(res.brca2wt.sd, 0.05) # 194 pathways have significant differences in expression variability 

# Names of the significant pathways 
path.names = names(sig.brca2wt.sd@genesInSigPways1)


## Plots of pathways involved in the DNA damage response ##
plotPway(res.brca2wt.sd, "Base excision repair", sig.brca2wt.sd)
plotPway(res.brca2wt.sd, "Mismatch repair", sig.brca2wt.sd, sizeXtick = 22)
plotPway(res.brca2wt.sd, "Nucleotide excision repair", sig.brca2wt.sd)
plotPway(res.brca2wt.sd, "Homologous recombination", sig.brca2wt.sd)
plotPway(res.brca2wt.sd, "Non-homologous end-joining", sig.brca2wt.sd)
plotPway(res.brca2wt.sd, "Cell cycle", sig.brca2wt.sd)

### Discrete Gene Count differences in pathways between BRCA2- and BRCA2+ ovarian cancers ###
disc.brca2wt.sd = pathVarTwoSamplesDisc(brca2.wt.log, pways.kegg, groups = as.factor(grp.brca2), varStat = "sd", test = "exact")

sig.brca2wt.disc = sigPway(disc.brca2wt.sd, 0.05) # 204 Significant Pathways

# Names of the significant pathways
disc.pathnames = names(sig.brca2wt.disc@genesInSigPways1) # 204 significant pathways 

## Plots of repair pathways with significant differences in gene counts ##
plotPway(disc.brca2wt.sd, "Base excision repair", sig.brca2wt.disc)
plotPway(disc.brca2wt.sd, "Mismatch repair", sig.brca2wt.disc)
plotPway(disc.brca2wt.sd, "Nucleotide excision repair", sig.brca2wt.disc)
plotPway(disc.brca2wt.sd, "Homologous recombination", sig.brca2wt.disc)
plotPway(disc.brca2wt.sd, "Non-homologous end-joining", sig.brca2wt.disc) # No significant difference between clusters 
plotPway(disc.brca2wt.sd, "Cell cycle", sig.brca2wt.disc) # focus on cell cycle pathway 

#####################################################################################################################################

## Grabbing the genes associated with the repair pathways ##

# Base excision repair #
base.repair = cbind(disc.brca2wt.sd@genesInPway1$"Base excision repair", disc.brca2wt.sd@genesInPway2$"Base excision repair")
colnames(base.repair) = c("brca2-", "wt-tumor")

baserep.genes = rownames(base.repair)
length(baserep.genes) # 32 genes 

table(base.repair[,1]) # Genes in BRCA2- (group1)
# 1  2  3 # var clusters
#27  4  1 

table(base.repair[,2]) # Genes in Wt-tumor (group2)
# 1  2  3 # var clusters
#14 17  1 

# Mismatch repair #
mismatch.repair = cbind(disc.brca2wt.sd@genesInPway1$"Mismatch repair", disc.brca2wt.sd@genesInPway2$"Mismatch repair")
colnames(mismatch.repair) = c("brca2-", "wt-tumor")

mismatch.genes = rownames(mismatch.repair)
length(mismatch.genes) # 22 genes 

table(mismatch.repair[,1]) # Genes in BRCA2- (group1)
# 1  2 # var clusters
#20  2

table(mismatch.repair[,2]) # Genes in Wt-tumor (group2)
# 1  2  3 # var clusters 
# 9 12  1 

# Nucleotide excision repair #
nucl.repair = cbind(disc.brca2wt.sd@genesInPway1$"Nucleotide excision repair", disc.brca2wt.sd@genesInPway2$"Nucleotide excision repair")
colnames(nucl.repair) = c("brca2-", "wt-tumor")

nucl.genes = rownames(nucl.repair)
length(nucl.genes) # 39 genes 

table(nucl.repair[,1]) # Genes in BRCA2- (group1)
# 1  2 # var clusters
#36  3 

table(nucl.repair[,2]) # Genes in Wt-tumor (group2)
# 1  2 # var clusters 
#22 17 

# Homologous Recombination #
hr.repair = cbind(disc.brca2wt.sd@genesInPway1$"Homologous recombination", disc.brca2wt.sd@genesInPway2$"Homologous recombination")
colnames(hr.repair) = c("brca2-", "wt-tumor")


hr.genes = rownames(hr.repair)
length(hr.genes) # 22 genes 

table(hr.repair[,1]) # Genes in BRCA2- (group1)
# 1 2 # var clusters 
#19 3

table(hr.repair[,2]) # Genes in Wt-tumor (group2)
# 1 2 3 # var clusters 
#9 10 3

# Non-homologous end joining #
nhej.repair = cbind(disc.brca2wt.sd@genesInPway1$"Non-homologous end-joining", disc.brca2wt.sd@genesInPway2$"Non-homologous end-joining") 
colnames(nhej.repair) = c("brca2-", "wt-tumor")


nhej.genes = rownames(nhej.repair)
length(nhej.genes) # 9 genes 

table(nhej.repair[,1])
# 1 2 # var clusters 
# 8 1

table(nhej.repair[,2])
# 1 2 # var clusters 
# 5 4

# cell cycle #  
cellcyc = cbind(disc.brca2wt.sd@genesInPway1$"Cell cycle", disc.brca2wt.sd@genesInPway2$"Cell cycle")
colnames(cellcyc) = c("brca2-", "wt-tumor")
cellcyc.genes = rownames(cellcyc)

### Grabbing all genes from the five DNA repair pathways ###

all.genes = c(baserep.genes, mismatch.genes, nucl.genes, hr.genes, nhej.genes, cellcyc.genes)
length(all.genes) # 242 genes 

all.genes = sort(all.genes) # sorting genes in alphabetical order

uniq.allgenes = unique(all.genes) # getting only unique genes 

length(uniq.allgenes) # 202 genes 

uniq.allgenes = sort(uniq.allgenes) # sorting unique genes in alphabetical order


#########################################################################################################################
### Searching for genes that transition from High/Medium variability (coded 3 and 2, respectively) to Low Variability (coded 1) in the BRCA2- ovarian cancers 

## Non-Homologous End Joining Pathway 
medVar.Wt = nhej.repair[which(nhej.repair[,2]!=1),]
lowVar.brca2 = medVar.Wt[which(medVar.Wt[,1] == 1),]
#
lowVar.nhej = medVar.Wt[which(medVar.Wt[,1] == 1),]
#
lowVarNhej.brca2 = rownames(lowVar.brca2)
length(lowVarNhej.brca2) # 3

## Nucleotide Excision Repair 
medVar.Wt = nucl.repair[which(nucl.repair[,2]!=1),]
lowVar.brca2 = medVar.Wt[medVar.Wt[,1]==1,] 
#
lowVar.nucl = medVar.Wt[medVar.Wt[,1]==1,] 
#
lowVarNuclrep.brca2 = rownames(lowVar.brca2)
length(lowVarNuclrep.brca2) # 14

## Mismatch Repair Pathway 
medVar.Wt = mismatch.repair[which(mismatch.repair[,2]!=1),]
lowVar.brca2 = medVar.Wt[which(medVar.Wt[,1] == 1),]
#
lowVar.mismatch = medVar.Wt[which(medVar.Wt[,1] == 1),]
#
lowVarMisrep.brca2 = rownames(lowVar.brca2)

length(lowVarMisrep.brca2) # 11

## Homologous Recombination Pathway 
medVar.Wt = hr.repair[which(hr.repair[,2]!=1),]
lowVar.brca2 = medVar.Wt[which(medVar.Wt[,1] == 1),]
#
lowVar.hr = medVar.Wt[which(medVar.Wt[,1] == 1),]
#
lowVarHR.brca2 = rownames(lowVar.brca2)

length(lowVarHR.brca2) # 10

## Base Excision Repair Pathway 
medVar.Wt = base.repair[which(base.repair[,2]!=1),]
lowVar.brca2 = medVar.Wt[which(medVar.Wt[,1] == 1),]
#
lowVar.base = medVar.Wt[which(medVar.Wt[,1] == 1),]
#
lowVarBaserep.brca2 = rownames(lowVar.brca2)

length(lowVarBaserep.brca2) #15

## Cell Cycle Pathway 
medVar.Wt = cellcyc[which(cellcyc[,2]!=1),]
lowVar.brca2 = medVar.Wt[which(medVar.Wt[,1] == 1),]
#
lowVar.cellcyc = medVar.Wt[which(medVar.Wt[,1] == 1),]
#
lowVarCellcyc.brca2 = rownames(lowVar.brca2)
length(lowVarCellcyc.brca2) #60

## Getting all Low Var genes 
alluniq.lowVar.genes = unique(c(lowVarBaserep.brca2, lowVarHR.brca2, lowVarMisrep.brca2, lowVarNuclrep.brca2, lowVarNhej.brca2, lowVarCellcyc.brca2))

length(alluniq.lowVar.genes) # 97 genes 

alluniq.lowVar.genes = sort(alluniq.lowVar.genes) # sorting genes alphabetically 

########################################################################################################
### Intersecting unique genes with SLI gene from SynLethDB database ###

synlethdb = read.table("synLethDB-BRCA2.txt", header = F, sep = "\t", fill = T, stringsAsFactors = F) # this is the text file that has all genes that have evidence for genes having synthetic lethal interaction with BRCA2
names = synlethdb[1,]
synlethdb = synlethdb[-1,]
dim(synlethdb) # 42 9 

### Finding overlap between genes found with Low variability using pathVar and the genes from the SynLethDB
lowvar.slgenes = intersect(alluniq.lowVar.genes, synlethdb[,3])
lowvar.slgenes2 = intersect(alluniq.lowVar.genes, synlethdb[,1])

ALLlowvar.slgenes = c(lowvar.slgenes, lowvar.slgenes2) 
length(ALLlowvar.slgenes)# 12 genes 

index = alluniq.lowVar.genes %in% ALLlowvar.slgenes
other.sligenes = alluniq.lowVar.genes[which(index == F)]
length(other.sligenes) # 85 

index2 = synlethdb[,3] %in% ALLlowvar.slgenes
index3 = synlethdb[,1] %in% ALLlowvar.slgenes
slgenes.validate1 = synlethdb[index2,]
dim(slgenes.validate1) # 10 9 
slgenes.validate2 = synlethdb[index3,]
dim(slgenes.validate2) # 2 9 
slgenes.validate = rbind(slgenes.validate1, slgenes.validate2)
dim(slgenes.validate) # 12 9


order.sl = order(slgenes.validate[,3])

slgenes.Validate = slgenes.validate[order.sl,]
colnames(slgenes.Validate) = names
dim(slgenes.Validate) # 12 9 
###################################################################################################






















































