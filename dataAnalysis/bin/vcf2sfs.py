
import pysam,argparse

def main():

	parser = argparse.ArgumentParser(description="")

	parser.add_argument("--vcf", 
		required = True,
		dest = "vcf",
		type =str, 
		help = "the name of the VCF file")
	parser.add_argument("--table", 
		required = True,
		dest = "table",
		type =str, 
		help = "the name of the TSV file that has the diergence and CpG status stuff")
	args = parser.parse_args()


	vcf = pysam.Tabixfile(args.vcf)
	table = pysam.Tabixfile(args.table)
if '_name__':
	main()
