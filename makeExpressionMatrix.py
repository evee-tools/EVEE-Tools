#!/usr/bin/env python
import os
import sys
import argparse

def main():
	parser = argparse.ArgumentParser(description='makeExpressionMatrix - Given a file of ortholog groups and rsem expression files across species, creates an expression matrix of gene x species and an index file for easily accessing a species of interest.')
	parser.add_argument('ortholog_groups', type=file, help='File of ortholog groups. First column is ortholog group name. Remaining columns are genes in each ortholog group, prefixed by species ID such that it matches the first column in rsem_list (e.g. speciesID|geneID).')
	parser.add_argument('rsem_list', type=file, help='File pointing to path of rsem output files (e.g. *.genes.results). First column is sample name (used for column header), second column is species ID (used for ortholog_groups lookup), and third column is path to rsem file. If you are using --gene_name flag, the gene ID from the --gene_name file should belong to the first species listed.')
	parser.add_argument('out_prefix', type=str, help="prefix for output files")
	parser.add_argument('--gene_name', type=file, help="File of species gene ID (first column) to gene name (second column). Gene ID should be from the first species listed in rsem_list.")
	args = parser.parse_args()
	
	outfile = open(args.out_prefix + ".exp_matrix.txt", 'w')
	outindex = open(args.out_prefix + ".index.txt", 'w')
	
	#species to dictionary
	speciesToSamples = {}
	samplesToTPM = {}

	species = []
	samples = []
	paths = []

	print "reading..."
	for line in args.rsem_list.readlines():
		line = line.split()
		
		curSample = line[0].strip()
		curSpecies = line[1].strip()
		curPath = line[2].strip()

		print ("\t%s\t%s" % (curSample, curPath))
		species.append(curSpecies)
		paths.append(curPath)
		samples.append(curSample)
		rsemFile = open(curPath, 'r')
		
		if curSpecies in speciesToSamples:
			curSamplesSet = speciesToSamples[curSpecies]
		else:
			curSamplesSet = set()
		curSamplesSet.add(curSample)
		speciesToSamples[curSpecies] = curSamplesSet
		
		if curSample in samplesToTPM:
			curDict = samplesToTPM[curSample]
		else:
			curDict = {}
		
		for line in rsemFile.readlines():
			line = line.split()
			gene = line[0].strip()
			tpm = line[5].strip()
			if tpm == "TPM": continue
			if gene in curDict:
				expArray = curDict[gene]
				expArray.append(float(tpm))
			else:
				expArray = [float(tpm)]
			curDict[gene] = expArray
		
		samplesToTPM[curSample] = curDict

	#speciesUniq 
	speciesUniq = []
	for curSpecies in species:
		if curSpecies not in speciesUniq: speciesUniq.append(curSpecies)
	
	#header/index
	outfile.write("#OG\tGENE_NAME")
	index = {}
	for curSpecies in speciesUniq:
		for curSample in speciesToSamples[curSpecies]:
			index[curSample] = samples.count(curSample)
			for x in range(samples.count(curSample)):
				outindex.write("%s\t%s\n" % (curSample, curSpecies))
				outfile.write("\t%s" % curSample)
	outfile.write("\n")

	#print species
	#print speciesUniq
	#print index

	#####GENE NAMES#####
	geneNames = {}
	if args.gene_name is not None:
		print "\nUsing %s as reference species for finding gene names" % species[0]
		for line in args.gene_name.readlines():
			line = line.split()
			geneNames[line[0].strip()] = line[1].strip()


	for line in args.ortholog_groups.readlines():
		line = line.split()
		og = line[0]
		genes = line[1:]
		
		curDict = {}
		
		geneName = ""
		for gene in genes:
			gene = gene.split("|")
			curSpecies = gene[0].strip()
			geneID = gene[1].strip()
			if curSpecies not in speciesToSamples: continue

			if curSpecies == speciesUniq[0]:
				if geneID in geneNames:
					geneName = geneName + geneNames[geneID] + "_"
					#print geneName

			for curSample in speciesToSamples[curSpecies]:
				rsemDict = samplesToTPM[curSample]
				if geneID in rsemDict:
					if curSpecies in curDict:
						expArray = curDict[curSample]
						for i in range(len(expArray)):
							expArray[i] = expArray[i] + rsemDict[geneID][i]
					else:
						curDict[curSample] = rsemDict[geneID]
		
		#reformat geneName
		if geneName == "": geneName = "Unannotated"
		else: geneName = geneName[:-1]
		outfile.write("%s\t%s" % (og, geneName))
		
		for curSpecies in speciesUniq:
			for curSample in speciesToSamples[curSpecies]:
				if curSample in curDict:
					curExp = curDict[curSample]
					for entry in curExp:
						outfile.write("\t%s" % entry)
				else:
					#print curSpecies
					for i in range(index[curSample]):
						outfile.write("\tNA")
		outfile.write("\n")

if __name__ == "__main__":
	main()
