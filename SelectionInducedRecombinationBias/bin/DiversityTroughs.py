from tom import brace
import argparse, pandas as pd, sys, random, gzip, pylab, glob, scipy
import site_frequency_spectrum as SFS
import matplotlib
import numpy as np


def gammaPDF(Nes,lamb):
	return (lamb) * np.exp(-1.*lamb * Nes)


def sweep( Nes, Ner4):
	return math.pow(Nes ,-1.*Ner4/(Nes/2)) ## Rho is in terms of 4Ner

def combined(x_Nes, Ner4, lamb, pa, mut_rate): #lamb here is the mean of the exponential midribution
	Va = mut_rate * pa * x_Nes # the rate of sweeps


	scipy.special.gamma()

	return Va * sweep(x_Nes, Ner4) * exponentialPDF(x_Nes, 1/lamb)


def BGS(mid):
	"""Use the Nordborg approx to get the diersity reduction around the exon"""

	t_mean = 0.025 # Mean sh for 
	shape = 0.2
	pa = 0.75
	sites = 1500

	Ne = 10000
	mut_rate = 0.01/(4.*Ne)

	r = mid / (4*Ne)

	s = (1.*Ya) / (2*Ne)

	Va = mut_rate * pa * Ya * 2

	S =  2 * Ne * sites * Va *  (Ya ** (-4.*r/s))
	
	model = 1./ (1. + S)

	return model 


def combinedSel(mid):
	"""model the trough in diverstiy around Exons"""
	"""Provide recombination midances as the 'x'"""
	"""Provide the reduction in diversity as the 'data'"""

	Ya = 700
	pa = 0.0003124
	sites = 1600

	Ne = 10000
	mut_rate = 0.01/(4.*Ne)

	r = mid / (4*Ne)

	s = (1.*Ya) / (2*Ne)

	Va = mut_rate * pa * Ya * 2

	S =  2 * Ne * sites * Va *  (Ya ** (-4.*r/s))
	
	model = 1./ (1. + S)

	return model 


def mergeManySFS(spectra):
	count = 0
	for i in spectra:
		if count == 0:
			sfs = i
		else:
			sfs = SFS.merge_SFS(sfs, i)
		count +=1
	return sfs 
## Quick script to turn a slim output file into a fasta

def processSLiM(input_file):
	intro = []
	introTrig = True
	mutations = []
	mutationsTrig = False
	genomes = []
	genomesTrig = False
	for i in gzip.open(input_file):
		if i == '// Starting run at generation <start>:\n':
			introTrig = False
			continue
		if i == 'Mutations:\n':
			mutationsTrig = True
			continue
		if i == 'Genomes:\n':
			mutationsTrig = False
			genomesTrig = True
			continue
		if introTrig == True:
			intro.append(i.strip())
		if mutationsTrig == True:
			mutations.append(i.strip())
		if genomesTrig == True:
			genomes.append(i.strip())

	
	

	mutations_raw = [i.split(' ') for i in mutations if i != 'Mutations:\n']	
	mutations_raw = [[int(i[0]),int(i[1]),i[2],int(i[3]),float(i[4]),float(i[5]),i[6],int(i[7]),int(i[8])] for i in mutations_raw if i != 'Mutations:\n']
	mutations = pd.DataFrame(mutations_raw, columns = ['ID','Num','Type','POS','s','h','Pop','genInit','freq'])

	architechture_raw = [i.strip(');').split('(')[1].split(', ') for i in intro if i.split('(')[0] == 'initializeGenomicElement']
	architechture_raw = [[i[0],int(i[1]),int(i[2])] for i in architechture_raw]
	architechture = pd.DataFrame(architechture_raw, columns = ['element','start','stop'])


	genomesDict = {}
	for i in genomes:
		g = i.split(' ')
		genomesDict[g[0]] = map(int,g[2:])
	if architechture.start.min() == 0:
		length = architechture.stop.max() + 1
	elif architechture.start.min() == 1:
		length = architechture.stop.max()
	else:
		print "we don't deal with your weird chromosomes here, bucko"
		sys.exit() 
	## Make the default DNA sequence for the region
	
	## Assume functional elements are named: g1
	focals = architechture[architechture.element == 'g1']

	focal_mid = list((focals.stop + focals.start)/2)[0] +1
	## Assume a single recombination rate ... VERY BREAKABLE
	r_rate_raw = [i for i in intro if i.startswith('initializeRecombinationRate')]
	r_rate = float(r_rate_raw[0].strip(');').split('(')[1])
	## Assume that the chromosome is 100,000bp long with the 1Kbp exon in the exact centre
	bins = [i for i in range(0,49000,100)] + [i for i in range(51000,100000, 100)]

	analysis = []
	for b in bins:
		start = b
		end = b + 99
		muts_in_range = list( mutations[(mutations.POS >= start) & (mutations.POS <= end)].freq )
		sfs = SFS.SFS_from_frequencies(muts_in_range, 100, 20)
		if focal_mid > end: 
			mid = focal_mid - ((1.*start+end)/2)
		elif focal_mid < end: 
			mid = ((1.*start+end)/2) - focal_mid

		r_dist_true = mid*r_rate * 40000		
		analysis.append( [start, end, mid, r_dist_true, sfs] )


	data = pd.DataFrame(analysis, columns = ['start','end','mid','r_dist_true','sfs'])
	return data	

def main():
	parser = argparse.ArgumentParser(description="Give a directory and I'll analyse the patterns of diversity around a simulated exon")

	parser.add_argument("-i","--input", 
			required = True,
			dest = "input",
			type =str, 
			help = "The name of the file that contains the SLiM output")
	args = parser.parse_args()
	count = 0
	r_bins = [float(i)/10 for i in range(1,500, 1)] + [1.*i for i in range(50,500)]

	approx = combinedSel(np.array(r_bins))
	pylab.plot(r_bins, approx*0.01, 'r')

	for i in glob.glob(args.input + '*.out.gz'):

		if count == 0:
			data = processSLiM(i)
		else:
			temp = processSLiM(i)
			data = pd.concat([data, temp])
		count += 1
		if count == 50:
			break
	top = max(range(len(r_bins)))
	analysis = []
	for r in range(len(r_bins)):
		bin = r_bins[r]

		if r == top:
			next = 1e8
		else:
			next = r_bins[r+1]

		in_range = list( data[(data.r_dist_true >= bin) & (data.r_dist_true < next)].sfs )
		if len(in_range) == 0: continue
		analysis.append([ r_bins[r], SFS.pi(mergeManySFS(in_range)) ])
	true_r = pd.DataFrame(analysis, columns = ['dist','pi'])
	pylab.plot(true_r.dist, true_r.pi, 'b')
	pylab.show()
if '__name__':
	main()

