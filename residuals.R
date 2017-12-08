args <- commandArgs(trailingOnly=TRUE)

if (length(args) != 4) {
	cat("usage: Rscript residuals.R data_matrix index reference out\n")
	cat("example: Rscript residuals.R data/brain.exp_matrix.one2one.txt data/brain.index.txt human brain.residuals.human.txt\n")
	quit()
}

exp <- read.delim(args[1], header=T, comment.char="", row.names=1, stringsAsFactors = F)
index <- read.delim(args[2], header=F, comment.char="")
refSpecies <- args[3]
outPath <- args[4]

data = exp[2:ncol(exp)]
#print(head(data))

residuals = c()
ref = as.matrix(data[,index[,1] == refSpecies])
if (sum(index[,1] == refSpecies) == 1) {
	ref.mean = ref
	rownames(ref.mean) = rownames(data)
} else {
	ref.mean = as.matrix(apply(ref, 1, function(x) mean(x, na.rm=T)))
}

cat("Processing...\n")
for (i in 1:ncol(data)) {
	cat("\t")
	cat(colnames(data)[i])
	cat("\n")
	other = as.matrix(data[,i])
	ref.data = ref.mean[ref.mean > 0 & other > 0,]
	other.data = other[ref.mean > 0 & other > 0,]
	t = cbind(log10(ref.data+0.01), log10(other.data+0.01))
	t = t[apply(t, 1, function(x) sum(is.na(x)))  == 0,]
	m = prcomp(t)
	res = m$x[,2]
	r = res[rownames(data)]
	if (m$rotation[2,2]<0) r = r*-1
	residuals = cbind(residuals, r)
}
colnames(residuals) = colnames(data)
rownames(residuals) = rownames(data)
write.table(cbind(exp[,1], residuals), file=outPath, sep="\t", quote=F)
