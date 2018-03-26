pruneTree = function(tree, subData) {
	species = unique(as.character(tree$species[!is.na(tree$species)]))
	curTree = tree
	##edit tree##
	for (curSpecies in species) {
		if (is.na(subData[curSpecies])) {
			curNode = curTree[curTree$species==curSpecies & !is.na(curTree$species),1]
			curParent = curTree[curTree$species==curSpecies & !is.na(curTree$species),3]
			curGrandparent = curTree[curTree$node==curParent & !is.na(curTree$node),3]
			curSib = curTree[curTree$ancestor==curParent & !is.na(curTree$ancestor) & curTree$node != curNode,1]     
			if (is.na(curTree[curTree$node == curParent,3])) {
				#remove opossum
				curTree = curTree[-which(curTree$node ==curNode),]
				#remove 1
				curTree = curTree[-which(curTree$node ==curParent),]
				#make opossum sib terminal node
				curTree[curTree$node == curSib,3] = NA
				curTree$time = curTree$time - curTree[curTree$node == curSib,4]
			} else {
				#remove terminal node
				curTree = curTree[-which(curTree$node ==curNode),]
				#connect sib to grandparent
				curTree[curTree$node == curSib,3] = curGrandparent
				#remove parent node
				curTree = curTree[-which(curTree$node ==curParent),]
			}
		}
	}
	tree.tr = seq(1, nrow(curTree))
	names(tree.tr) = curTree$node
	tree.tr = as.data.frame(tree.tr)
	curTree$node = tree.tr[as.character(curTree$node),]
	curTree$ancestor = tree.tr[as.character(curTree$ancestor),]
	rownames(curTree) = curTree$node
	return(curTree)
}

args <- commandArgs(trailingOnly=TRUE)

if (length(args) != 4) {
	cat("usage: Rscript ./fitOUModel.R exp_matrix index tree out\n")
	cat ("example: Rscript ./fitOUModel.R data/liver.exp_matrix.txt data/liver.index.txt data/mammals.tree.txt liver.ou_stats.txt\n")
	
	quit()
}

library(edgeR)
library(ouch)

exp <- read.delim(args[1], header=T, comment.char="", row.names=1, stringsAsFactors = F)
index <- read.delim(args[2], header=F, comment.char="")
ouch <- read.delim(args[3])
outPath <- args[4]

##TODO: check data matches index##

##get average of species
data=exp[2:ncol(exp)]
data_mean = c()
for (i in unique(index[,1])) {
	curSpecies = data[,index[,1] == i]
	if (sum(index[,1] == i) == 1) {
		data_mean = cbind(data_mean, curSpecies)
	} else {
		data_mean = cbind(data_mean, apply(curSpecies, 1, function(x) mean(x, na.rm=T)))
	}
}
colnames(data_mean) = unique(index[,1])

mammals.tree = ouchtree(ouch$node, ouch$ancestor, ouch$time, labels = as.character(ouch$species))
#plot(mammals.tree)

##normalize data##
data_mean.noNA = data_mean
data_mean.noNA[is.na(data_mean)] = 0
scale = calcNormFactors(data_mean.noNA, refColumn = 1, method="TMM")
scale = scale / scale[1]
data_mean.norm = sapply(seq(1, length(scale)), function(x) data_mean[,x] / scale[x])
data_mean.norm = log10(data_mean.norm+0.01)
colnames(data_mean.norm) = colnames(data_mean)
#print(head(data_mean.norm))

sigmas = rep(NA, nrow(data_mean.norm))
alphas = rep(NA, nrow(data_mean.norm))
logLik = rep(NA, nrow(data_mean.norm))
pvalues = rep(NA, nrow(data_mean.norm))
thetas = rep(NA, nrow(data_mean.norm))
brownSigmas = rep(NA, nrow(data_mean.norm))

#check that exp_matrix names matches tree names
numIntersect = length(intersect(colnames(data_mean.norm), as.character(ouch$species[!is.na(ouch$species)])))
if (numIntersect == 0) {
	cat("\nNone of the names from your index file matched a species name from your tree! Exiting.\n")
	quit()
}

cat(paste0("\n", numIntersect, " species found in index and tree.\n"))
print(intersect(colnames(data_mean.norm), as.character(ouch$species[!is.na(ouch$species)])))
cat("\n")

cat("\nFitting...\n")
for (i in seq(1, nrow(data_mean.norm))) {
	if (i %% 1000 == 0) cat(paste0(i, "\tgenes\r"))
	subData = data_mean.norm[i,intersect(colnames(data_mean.norm), as.character(ouch$species[!is.na(ouch$species)]))]

	if (sum(!is.na(subData)) < 3) next
	if (sd(subData, na.rm=T) == 0) next
	curTree = pruneTree(ouch, subData)
	curOuchTree = ouchtree(curTree$node, curTree$ancestor, curTree$time, labels = as.character(curTree$species))

#rearrange subData
	subData = subData[as.character(curTree$species[!is.na(curTree$species)])]
#add ancestral nodes
	curData = c(rep(NA, sum(is.na(curTree$species))), subData[!is.na(subData)])
	names(curData) = curTree$node
	curRegime = as.factor(rep("ns", nrow(curTree)))
	names(curRegime) = curTree$node

	h = hansen(data=curData, tree=curOuchTree, regimes=curRegime, sqrt.alpha = 1, sigma = 1, reltol=1e-5)
	b = brown(curData, curOuchTree)
	brownSigmas[i] = coef(b)$sigma.sq.matrix[,1]
	sigmas[i] = coef(h)$sigma.sq.matrix[,1]
	alphas[i] = coef(h)$alpha.matrix[,1]
	thetas[i] = coef(h)$theta$curData
	logLik[i] = logLik(h)
	pvalues[i] =  1 - pchisq((logLik(h) - logLik(b)) * 2, 1)
}
#print(head(sigmas))   
names(sigmas) = rownames(data_mean.norm)
names(alphas) = rownames(data_mean.norm)
names(logLik) = rownames(data_mean.norm)
names(pvalues) = rownames(data_mean.norm)  
names(thetas) = rownames(data_mean.norm)
names(brownSigmas) = rownames(data_mean.norm)
qvalues = p.adjust(pvalues, method="fdr")
var = sigmas / (alphas*2)
gene_name = exp[names(thetas),1]

stats = cbind(gene_name, qvalues, thetas, var, brownSigmas)

write.table(stats, file=outPath, sep="\t", quote=F, row.names=T, col.names=T) 


