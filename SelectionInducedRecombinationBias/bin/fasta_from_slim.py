from tom import brace
import argparse, pandas as pd, sys, random, copy, gzip
## Quick script to turn a slim output file into a fasta

def initSequence(length):
	DNA = ['A','T','G','C']
	return [random.choice(DNA) for i in range(length)]

def mutate(base):
	DNA = ['A','T','G','C']
	mutation = random.choice(DNA)
	while mutation == base:
		mutation = random.choice(DNA)
		
	return mutation


def main():
	parser = argparse.ArgumentParser(description="From a SLiM utput file, generate a FASTA file")

	parser.add_argument("-i","--input", 
			required = True,
			dest = "input",
			type =str, 
			help = "The name of the file that contains the SLiM output")
	parser.add_argument("-o","--output", 
			required = True, 
			metavar ="output", 
			type = str, 
			help = "What name do you want to give the output file?")
	args = parser.parse_args()
	intro = []
	introTrig = True
	mutations = []
	mutationsTrig = False
	genomes = []
	genomesTrig = False
	for i in gzip.open(args.input):
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
	TemplateSequence = initSequence(length)
	changesDict = {}
	outputFASTA = open(args.output,'w')
	for i in genomesDict.keys():
		mySeq = copy.copy(TemplateSequence)
		for m in genomesDict[i]:
			mPOS = list(mutations[mutations.ID == m].POS)
			if len(mPOS) > 1:
				print "There should be one mutation per ID"
				return
			try:
				mySeq[mPOS[0]] = changesDict[mPOS[0]]
			except KeyError:				
				changesDict[mPOS[0]] = mutate(mySeq[mPOS[0]])
				mySeq[mPOS[0]] = changesDict[mPOS[0]]
		outputFASTA.write('>'+i+'\n'+''.join(mySeq)+'\n')
	outputFASTA.close()

if '__name__':
	main()
