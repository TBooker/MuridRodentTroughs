#### Script to generate a table of positions and chromosomes with CpG status and alleles for an ingroup and two outgroups
import argparse, pysam
from Bio import SeqIO
import pandas as pd

def main():

	parser = argparse.ArgumentParser(description="")

	parser.add_argument("--alleles", 
		required = True,
		dest = "alleles",
		type =str, 
		help = "The name(s) of the FASTA files with Alleles in them",
		nargs = '+')
	parser.add_argument("--CpG", 
		required = True,
		dest = "cpg",
		type =str, 
		help = "The name(s) of the FASTA files with CpG status in them",
		nargs = '+')
	parser.add_argument("--names", 
		required = True,
		dest = "names",
		type =str, 
		help = "The name of the sequences that you are making a table for ",
		nargs = '+')
	parser.add_argument("--chrom", 
		required = True,
		dest = "chrom",
		type =str, 
		help = "The chromosome you want ot make a table for")
	parser.add_argument("--output", 
		required = True,
		dest = "output",
		type =str, 
		help = "The name for the output file")
	args = parser.parse_args()


	seqs = []
	for allele_file_i in args.alleles:
		fasta = pysam.FastaFile(allele_file_i)	
		seqs.append([i.upper() for i in fasta.fetch(args.chrom)])	
#		allele_record_dict = SeqIO.to_dict(SeqIO.parse(allele_file_i, "fasta"))
#		seqs.append( allele_record_dict[args.chrom].seq )
	
		
	cpgs = []
	for cpg_file_i in args.cpg:
		cpg_fasta = pysam.FastaFile(cpg_file_i)	
		cpgs.append([i.upper() for i in cpg_fasta.fetch(args.chrom)])	

#		cpg_record_dict = SeqIO.to_dict(SeqIO.parse(cpg_file_i, "fasta"))
#		cpgs.append( cpg_record_dict[args.chrom].seq )
	dfs = []
	names = []
	for i in range(len(seqs)):
		temp = pd.DataFrame({args.names[i]: list(seqs[i]),
				 args.names[i]+'_cpg':list(cpgs[i])})
		dfs.append(temp)

		names.append(args.names[i])
		names.append(args.names[i]+'_cpg')

	full = pd.concat(dfs, axis=1)
	pos = len(full)
	chrom = args.chrom
	full['POS'] = [i+1 for i in range(pos)]
	full['CHROM'] = chrom


	full = full[['CHROM','POS'] + names]
	full.to_csv(args.output, header = False, index = False, sep = '\t')

if '__name__':
	main()
