args <- commandArgs(trailingOnly=TRUE)

if (length(args) != 4) {
	cat("usage: Rscript ./scoreGenes.R ou_stats ref_exp test_exp outPrefix\n")
	cat("example: Rscript ./scoreGenes.R data/liver.ou_stats.txt data/ref_exp/liver.ref_exp.txt data/patient.test.txt patient\n")	
	quit()
}

library(edgeR)

cat("\nreading ou_stats...\n")
stats <- read.delim(args[1], header=T, comment.char="", stringsAsFactors = F, row.names=1)
#get only unique genes
cat(paste0("      removed ", sum(duplicated(stats[,1])), " genes with duplicated names\n"))
stats = stats[!duplicated(stats[,1]),]
rownames(stats) = stats[,1]

cat("reading ref_exp...\n")
ref_exp.in <- read.delim(args[2], header=T, comment.char="", stringsAsFactors = F)
ref_exp.in = ref_exp.in[!duplicated(ref_exp.in[,1]),]
ref_exp = as.data.frame(ref_exp.in[,2])
rownames(ref_exp) = ref_exp.in[,1]

cat("reading test_exp...\n")
exp.in <- read.delim(args[3], header=T, comment.char="", stringsAsFactors = F)
exp.in = exp.in[!duplicated(exp.in[,1]),]
exp = as.data.frame(exp.in[,2])
rownames(exp) = exp.in[,1]
colnames(exp) = colnames(exp.in)[2]
#print(head(exp))
outPath <- args[4]

gene_set = intersect(rownames(stats), rownames(exp))
#print(head(gene_set))

stats.filt = stats[gene_set,]
ref_exp.filt = ref_exp[gene_set,]
exp.filt = as.data.frame(exp[gene_set,])
rownames(exp.filt) = gene_set
colnames(exp.filt) = colnames(exp)
#print(head(exp.filt))

##normalize gene expression
cat("\nnormalizing gene expression...\n")
cur_mean = cbind(ref_exp.filt, exp.filt)
cur_mean.noNA = cur_mean
cur_mean.noNA[is.na(cur_mean)] = 0
scale = calcNormFactors(cur_mean.noNA, refColumn = 1, method="TMM")
scale = scale / scale[1]

exp.norm = sapply(seq(2, length(scale)), function(x) cur_mean[,x] / scale[x])
rownames(exp.norm) = rownames(exp.filt)
colnames(exp.norm) = colnames(exp.filt)

cat(paste0("\nscoring ", length(gene_set), " genes...\n"))
zscores = c()
for (i in 1:nrow(exp.norm)) {
	curRPKM = log10(as.vector(as.matrix(exp.norm[i,]))+0.01)
	curMean = stats.filt[i,]$thetas
	curVar = stats.filt[i,]$var
	curQvalue = stats.filt[i,]$qvalues
	if (!is.na(curVar) & curMean > .69 & curQvalue<0.05) {
		curZ = (curRPKM - curMean) / sqrt(curVar)
	} else {
		curZ = rep(NA, length(curRPKM))
	}
	zscores = rbind(zscores, curZ)
}
rownames(zscores) = rownames(exp.norm)
colnames(zscores) = colnames(exp.norm)

pvalues = apply(zscores, 2, function(x) 2*pnorm(-abs(x)))
qvalues = apply(pvalues, 2, function(x) p.adjust(x, method="fdr"))

write.table(qvalues, file = paste0(outPath, ".qvalues.txt"), sep="\t", quote=F)
write.table(zscores, file = paste0(outPath, ".zscores.txt"), sep="\t", quote=F)
