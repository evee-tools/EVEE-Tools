pruneTree = function(tree, subData) {
	species = unique(as.character(tree$species[!is.na(tree$species)]))
	curTree = tree
	##edit tree##
	for (curSpecies in species) {
		if (!(curSpecies %in% names(subData)) | is.na(subData[curSpecies])) {
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

if (length(args) != 5) {
	cat("usage: Rscript ouRegimes.R exp_matrix index tree regimes outPrefix\n")
	cat("example: Rscript ouRegimes.R data/liver.exp_matrix.txt data/liver.index.txt data/mammals.tree.txt data/regimes.all.txt liver.regimes\n")
	quit()
}

library(edgeR)
library(ouch)

exp <- read.delim(args[1], header=T, comment.char="", row.names=1, stringsAsFactors = F)
index <- read.delim(args[2], header=F, comment.char="")
ouch <- read.delim(args[3])
regimes <- read.delim(args[4], header=T, comment.char="")
outPrefix <- args[5]

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

pvalues = matrix(NA, nrow = nrow(data_mean.norm), ncol = ncol(regimes))
thetas = matrix(NA, nrow = nrow(data_mean.norm), ncol = sum(apply(regimes, 2, function(x) length(unique(x)))))
aic = matrix(NA, nrow = nrow(data_mean.norm), ncol = ncol(regimes))
bic = matrix(NA, nrow = nrow(data_mean.norm), ncol = ncol(regimes))

#check that exp_matrix names matches tree names
numIntersect = length(intersect(colnames(data_mean.norm), as.character(ouch$species[!is.na(ouch$species)])))
if (numIntersect == 0) {
	cat("\nNone of the names from your index file matched a species name from your tree! Exiting.\n")
	quit()
}

cat(paste0("\n", numIntersect, " species found in index and tree.\n"))
print(intersect(colnames(data_mean.norm), as.character(ouch$species[!is.na(ouch$species)])))
cat("\n")

ouch.species = as.character(ouch$species[!is.na(ouch$species)])
data.species = colnames(data_mean.norm)
intersect.species = ouch.species[ouch.species %in% data.species]

cat("\nFitting...\n")
for (i in seq(1, nrow(data_mean.norm))) {
	if (i %% 1000 == 0) cat(paste0(i, "\tgenes\r"))
	
	subData = data_mean.norm[i,intersect.species]
	
	#subData = data_mean.norm[i,ouch[!is.na(ouch$species),]$species]
	if (sum(!is.na(subData)) < 3) next
	if (sd(subData, na.rm=T) == 0) next
						    
	curNullTree = pruneTree(ouch, subData)
	curOuchNullTree = ouchtree(curNullTree$node, curNullTree$ancestor, curNullTree$time, labels = as.character(curNullTree$species))
	curData = c(rep(NA, sum(is.na(curNullTree$species))), subData[!is.na(subData)])
	names(curData) = curNullTree$node
	b = brown(curData, curOuchNullTree)
	
	#print(subData)
	#print(curData)
	#print(curNullTree)

	#column indices for thetas
	thetaStart = 1
	
	thetaNames = c()

	for (x in 1:ncol(regimes)) {
		numThetas = length(unique(regimes[,x]))
		thetaNames = c(thetaNames, rep(colnames(regimes)[x], numThetas))

		ouch.curRegime = cbind(ouch, regimes[,x])
		curTree = pruneTree(ouch.curRegime, subData)
		colnames(curTree) = c("node", "species", "ancestor", "time", "regimes")
		curOuchTree = ouchtree(curTree$node, curTree$ancestor, curTree$time, labels = as.character(curTree$species))
		curData = c(rep(NA, sum(is.na(curTree$species))), subData[!is.na(subData)])
		names(curData) = curTree$node
		curRegime = as.factor(curTree$regimes)
		names(curRegime) = curTree$node
		
		#print(curData)
		h = hansen(data=curData, tree=curOuchTree, regimes=curRegime, sqrt.alpha = 1, sigma = 1, reltol=1e-5)
		pvalues[i,x] = 1 - pchisq((logLik(h) - logLik(b)) * 2, summary(h)$dof-1)
		thetas[i,(thetaStart:(thetaStart+numThetas-1))] = coef(h)$theta$curData
		thetaStart = thetaStart + numThetas
		aic[i,x] = summary(h)$aic
		bic[i,x] = summary(h)$sic
	}
}

rownames(pvalues) = rownames(data_mean.noNA)
rownames(thetas) = rownames(data_mean.noNA)
rownames(aic) = rownames(data_mean.noNA)
rownames(bic) = rownames(data_mean.noNA)

colnames(thetas) = thetaNames
colnames(pvalues) = colnames(regimes)
colnames(aic) = colnames(regimes)
colnames(bic) = colnames(regimes)

write.table(pvalues, file = paste0(outPrefix, ".pvalues.txt"), sep="\t", quote=F)
write.table(thetas, file = paste0(outPrefix, ".thetas.txt"), sep="\t", quote=F)
write.table(aic, file = paste0(outPrefix, ".aic.txt"), sep="\t", quote=F)
write.table(bic, file = paste0(outPrefix, ".bic.txt"), sep="\t", quote=F)

