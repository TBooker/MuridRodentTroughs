
import pysam, argparse, vcf, random
from collections import Counter

def SFS_from_all_frequencies(frequencies,N):
	SFS = [0]*(N+1)
	length = len(frequencies)
	for i in frequencies:
		if type(i) != int:
			print 'Why are you treating strings like integers you weirdo?'
			return
		if i > N:
			print "SFS_from_frequencies: Error in your frequencies vector: One of the values is greater than the number of individuals\nThe offending value is: " + str(i) +" and the sample is "+str(N)
			return
		if i > 0 and i <= N:
			SFS[i] += 1
			
	SFS[0] = length - sum(SFS)
	if sum(SFS) < length:
		print"SFS_from_frequencies: Error in your frequencies vector: Fewer items in the SFS than the length of the region"
		return
	if sum(SFS) > length:
		print"SFS_from_frequencies: Error in your frequencies vector: More items in the SFS than the length of the region"
		return
	return SFS

def getGenotypes(genos, allele_dict):
## This is a really horrible nested list comprehension within a list comprehension that is used to flatten the list
	return [g for sublist in [[allele_dict[i.split('/')[0]] , allele_dict[i.split('/')[1]]] for i in genos] for g in sublist]		
	
def genos2freqs(genos):
## get the MAF for the alleles, or return 0 if the site is invariant 
	freqs = Counter(genos)
	if len(freqs) >2:
		return None
	if len(freqs) == 1:
	#	print freqs, len(freqs)
		return 0
	else:
	#	print freqs, min(freqs[freqs.keys()[0]], freqs[freqs.keys()[1]])
		return min(freqs[freqs.keys()[0]], freqs[freqs.keys()[1]])

def indel_quest(vcf_line):

	if vcf_line.ALT == '[None]':
		return False
	elif len(vcf_line.REF) > 1 or len(vcf_line.ALT) >1:
		return True
	elif len(vcf_line.REF) == 1 and len(vcf_line.ALT)  == 1:
		return False
	else:
		print 'special case'
		return False

def divergent(genos, outgroup_base):
	if len(set(genos)) == 1:
		if genos[0] != outgroup_base:
			return 1
		else:
			return 0
	elif len(set(genos)) > 1:

		if random.choice(genos) != outgroup_base:
			return 1
		else:
			return 0

def testCpG( cpg_status , index):
	if set( [cpg_status[0][index], cpg_status[1][index], cpg_status[2][index] ]) == set(['0']):

		return True
	else:
		return False

def analyse_chunk(vcf_chunk, start, divergence, cpg_status, min_QUAL = 20, max_Depth = 700, min_GQ = 15, N = 10):#,  outgroup_1, outgroup_2, cpg_outgroup_1, cpg_outgroup_2, exclude_cpg = False,):
	all_freqs = []
	ncpg_freqs = []
	index = 0
        previous_POS = start-1
        print 'starting at base:',previous_POS
	diverged_all_1 = 0 ## Is outgroup 1 divergent (all sites)?
	diverged_all_2 = 0 ## Is outgroup 2 divergent (all sites)?
	diverged_ncpg_1 = 0 ## Is outgroup 1 divergent (ncpg sites)?
	diverged_ncpg_2 = 0 ## Is outgroup 2 divergent (ncpg sites)?
        current_position = 0
	for i in vcf_chunk:
                
#                print previous_POS,i.POS
                while i.POS != start + current_position:
                        current_position +=1
  #              print i.POS, start + current_position
                
                fasta_pos = current_position
#                index +=1

                if divergence[0][fasta_pos] =='.' or divergence[1][fasta_pos] =='.':
                        continue
                #print i.REF, divergence[0][fasta_pos], divergence[1][fasta_pos]

#		index +=1
		allele_dict = {'0' : i.REF}
#		if len(i.ALT)> 1:
#			print 'more than 1 alternate allele'
#			print 'I only analyse biallelic sites'
#			continue
		if i.ALT != [None]:

			variant = True
			for q in range(len(i.ALT)):
				allele_dict[str(q+1)] = str(i.ALT[q])		
		else:
			variant = False
	 			### FILTER BED
		if i.QUAL < min_QUAL: continue
		## Is the site an IN/DEL?
		if indel_quest(i): continue
		## Does the number of called individuals match the requirement
		if i.num_called < N: continue

		if i.INFO['DP'] > max_Depth: continue

		if variant:
			if float(i.INFO['ExcessHet']) > 30: continue
			## If the site is a variant AND the min depth matches the min DP and GQ
			genos = getGenotypes( [j['GT'] for j in i.samples if j['GQ'] >= min_GQ ], allele_dict ) 
			if len(genos) < N: continue
		elif not variant:
			genos = [str(i.REF) , str(i.REF)] * N
		
		all_frequency = genos2freqs(genos)

		if all_frequency != None:
			all_freqs.append( all_frequency )

		diverged_all_1 += divergent(genos, divergence[0][fasta_pos].upper())
		diverged_all_2 += divergent(genos, divergence[1][fasta_pos].upper())
		## Here you add a CpG-prone filter 

		if testCpG( cpg_status , fasta_pos):
                        continue

		ncpg_freq = all_frequency

		if ncpg_freq != None:
			ncpg_freq =  all_frequency
                        #print all_frequency
			ncpg_freqs.append(all_frequency)
			diverged_ncpg_1 += divergent( genos, divergence[0][fasta_pos].upper() )
			diverged_ncpg_2 += divergent( genos, divergence[1][fasta_pos].upper() )

		

	all_sfs = SFS_from_all_frequencies( all_freqs, N )
	ncpg_sfs = SFS_from_all_frequencies( ncpg_freqs, N )
	
	print all_sfs, float(diverged_all_1)/sum(all_sfs), float(diverged_all_2)/sum(all_sfs)
        print ncpg_sfs, float(diverged_ncpg_1)/sum(ncpg_sfs), float(diverged_ncpg_2)/sum(ncpg_sfs)
	return

def concatBed(bed_line):
	return bed_line[0] + ':' + str(bed_line[1]) + '-' + str(bed_line[2])

def main():

	parser = argparse.ArgumentParser(description="")

	parser.add_argument("--vcf", 
		required = True,
		dest = "vcf",
		type =str, 
		help = "the name of the VCF file")

	parser.add_argument("--cpg", 
		required = True,
		dest = "ingroup_cpg",
		type =str, 
		help = "A fasta file with the CpG status of the ingroup")

	parser.add_argument("--outgroup_1", 
		required = True,
		dest = "outgroup_1",
		type =str, 
		help = "The name of the FASTA files with alleles in them for outgroup_1")

	parser.add_argument("--outgroup_2", 
		required = True,
		dest = "outgroup_2",
		type =str, 
		help = "The name of the FASTA files with alleles in them for outgroup_2")

	parser.add_argument("--outgroup_cpg_1", 
		required = True,
		dest = "outgroup_cpg_1",
		type =str, 
		help = "The name of the FASTA files with CpG status in them for outroup1")

	parser.add_argument("--outgroup_cpg_2", 
		required = True,
		dest = "outgroup_cpg_2",
		type =str, 
		help = "The name of the FASTA files with CpG status in them for outgroup2")

	parser.add_argument("--min_individuals", 
		required = True,
		dest = "min_individuals",
		type = int, 
		help = "What is the minimum number of individuals ")

	args = parser.parse_args()
 
	testChunks = [['chr19',20000000,20020000] ,
			['chr19',10000000,10020000]]
	
	vcf_object = vcf.Reader(open(args.vcf, 'r'))
	outgroup_1_fasta = pysam.FastaFile(args.outgroup_1)	
	outgroup_2_fasta = pysam.FastaFile(args.outgroup_2)	
	ingroup_cpg_fasta = pysam.FastaFile(args.ingroup_cpg)	
	outgroup_1_cpg_fasta = pysam.FastaFile(args.outgroup_cpg_1)	
	outgroup_2_cpg_fasta = pysam.FastaFile(args.outgroup_cpg_2)	


	for i in testChunks:
		print concatBed(i)


		divergence = [	outgroup_1_fasta.fetch(i[0],start = i[1]-1, end = i[2]),
				outgroup_2_fasta.fetch(i[0],start = i[1]-1, end = i[2]) ]

		cpg_status = [	ingroup_cpg_fasta.fetch(i[0],start = i[1]-1, end = i[2]),
				outgroup_1_cpg_fasta.fetch(i[0],start = i[1]-1, end = i[2]),
				outgroup_2_cpg_fasta.fetch(i[0],start = i[1]-1, end = i[2]) ]

## PyVCF is 'half-open' so the first position is not included in the bit that is fetched from the file
		analyse_chunk( vcf_object.fetch(i[0],start = i[1]-1, end = i[2])  , i[1], divergence, cpg_status)



if '_name__':
	main()
