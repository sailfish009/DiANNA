#!/usr/bin/python
# F.Ferre
# Jan 2006: Stand-alone version 1.0
# Cysteine connectivity prediction with neural networks
# 1 - PSI-Blast profile and secondary structure prediction using psiPred
# 2 - cysteines state (half or free) prediction using NN (using evolutionary information)
# 3 - cysteine connectivity prediction using NN (using evolutionary and secondary structure information)
# 4 - weighted graph matching using the Gabow algorithm


import cgi, re, sys, os, string



########################
### Global Variables ###
########################

aminoAcids        = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
AvgAaFreq         =  '   8   5   5   6   2   4   7   7   2   6   9   6   2   4   4   6   6   1   3   7'
emptyProfileRow   =  '   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0'
secStructEncoding = {'H':'1.000000 0.000000 0.000000', 'C':'0.000000 1.000000 0.000000', 'E':'0.000000 0.000000 1.000000'}

# The name of the BLAST data bank
dbname  = '/cluster/home/clotelab/blast/swissprot_red.nt'
# Where the NCBI programs have been installed
ncbidir = '/cluster/home/clotelab/blast'
# Where the PSIPRED V2 programs have been installed
execdir = '/cluster/home/clotelab/psipred/bin'
# Where the PSIPRED V2 data files have been installed
datadir = '/cluster/home/clotelab/psipred/data'
# Psipred folder
rootdir = '/cluster/home/clotelab/psipred'



##################################################
### Check the cysteine number and the seq name ###
##################################################
def checkCys():
	sequenceFile = open("seqFile.faa","r")
	sequence     = sequenceFile.readlines()

	if sequence[0][0] != '>':    ###Check; this case should be impossible
		print "The submitted sequence must be in FASTA format"
		sys.exit(1)
	else:
		seq     = ""
		seqName = sequence[0][1:-1]
		for i in range(1,len(sequence)):
			seq += re.sub("\n","",sequence[i])
		seq.replace("\r","")
		print "Sequence %s Length %d residues" % (seqName,len(seq))
		print "Cysteines in this sequence: %d\n" % seq.count('C')
		
		if seq.count('C') == 0 or seq.count('C') == 1:
			print "No disulfide bonds are possible for this sequence, the prediction will not be attempted"
			sys.exit(1)

	return(seqName)


	
####################################
### Compute secondary structure ####
####################################
def runPsipred(seqName,sel):
	os.system('cp seqFile.faa psitmp.fasta')

	print "Step 1: Running PSI-BLAST with input sequence"

	comm = "%s/blastpgp -b 0 -j 3 -h 0.001 -d %s -i psitmp.fasta -C psitmp.chk -Q out.ckp > temp.blast" % (ncbidir, dbname)
	os.system(comm)

	if sel != 'ModuleA':
		print "Step 2: Predicting secondary structure using PSIPRED"

	os.system('echo psitmp.chk > psitmp.pn')
	os.system('echo psitmp.fasta > psitmp.sn')

	comm = '%s/makemat -P psitmp' % ncbidir
	os.system(comm)

	comm = '%s/psipred psitmp.mtx %s/weights.dat %s/weights.dat2 %s/weights.dat3 %s/weights.dat4 > out.ss' % (execdir,datadir,datadir,datadir,datadir)
	os.system(comm)

	comm ='%s/psipass2 %s/weights_p2.dat 1 0.98 1.09 out.ss2 out.ss > out.horiz' % (execdir, datadir)
	os.system(comm)
### Remove temporary files
#	os.system('rm -f '+rootdir+'/psitmp.*')
#	os.system('rm -f psitmp.*')
#	print "\nSecondary Structure Predicted.\n"



############################
### read psiblast output ###
############################
def readProfile(seq):
	#### Profiles directory ####
	profDir = 'profileOut/'
	prof    = ''

###append fake profle rows before the beginning of the sequence profile rows (avg residue frequencies)
	for i in range(10):
		prof += AvgAaFreq

###read blast output
	profile = open('out.ckp','r').readlines()
	for i in range(3,len(profile)-6):
		if profile[i][70:150] == emptyProfileRow:
			prof += writeProfileRow(seq[i-13])
		else:
			prof += profile[i][70:150]

###append fake profle rows at the end of the sequence profile rows
	for i in range(10):
		prof += AvgAaFreq
	return(prof)



############################################################################################
### handle empty profile rows (happens when input sequence have no homoogs in swissprot) ###
############################################################################################
def writeProfileRow(r):
	row = ''
	for i in range(20):
		if r == aminoAcids[i]:
			row += ' 100'
		else:
			row += '   0'
	return(row)



###########################
### read psipred output ###
###########################
def readSecStruct():
	#### Secondary Structure predictions directory ####

	secStruct    = 'CCCCCCCCCC'
	secStructOut = open('out.ss2','r').readlines()
	for i in range(2,len(secStructOut)):
		secStruct += secStructOut[i][7]
	secStruct += 'CCCCCCCCCC'
	return(secStruct)



###############################
### read seq from input file ###
###############################
def readSeq():
	seq      = 'XXXXXXXXXX'   ###append X at the beginning of the sequence (to handle C too close to the edges)
	sequence = open('seqFile.faa','r').readlines()
	for i in range(1,len(sequence)):
		seq += re.sub("\n","",sequence[i])
	seq =  string.replace(seq,'\r','')
	seq += 'XXXXXXXXXX'       ###append X at the end of he sequence
	return(seq)



#######################################################################
### write input file for Module A and return cysteines for ModuleBC ###
#######################################################################
def createDataset(seq,prof,secStruct,window):
	if len(seq) == len(prof)/80 == len(secStruct):
		outFile   = open('outFileA','w')
		cysteines = []
		for i in range(len(seq)):
			if seq[i] == 'C':
				if i-window >= 0 and i+window+1 <= len(seq):
					neighbor          = seq[i-window:i+window+1]
					secStructNeighbor = secStruct[i-window:i+window+1]
					profNeighbor      = prof[i*80-window*80:i*80+window*80+80]
					outFile.write(neighbor+' '+secStructNeighbor+' '+profNeighbor+'\n')
					cysteines.append([i+1-10,neighbor,secStructNeighbor,profNeighbor])
		outFile.close()
		return(cysteines)
	else:
		print "An error occurred! Please check the input sequence"
		sys.exit(1)



###################################
### Oxidation State prediction ####
###################################
def disulfideStatePrediction(cysteines,window,Print):
	outFile = open('outFileB','w')
	if Print:
		print "Step 2: Disulfide Oxidation State Prediction using a trained Neural Network"

	### create foo.pat for state prediction
	createFooPatA('outFileA',window)

	### test foo.pat
	os.system('runNetOxState foo.test.pat foo.test.res')

	### read scores
	results = open('foo.test.res','r').readlines()
	scores  = []
	halfC   = 0

	for line in results:
		NNOutput = float(line)
		if NNOutput >= 0.5:
			halfC += 1
		scores.append(NNOutput)

	### write results in a table
	positives = []
	i         = 0
	print
	print 'Sequence     Cysteine         Score'
	print 'Position     neighborhood     -----'
	print '--------     ------------'
	for i in range(0,len(scores)):
		print '   %d%s     %s     %f' % (cysteines[i][0],' '*(5-len(repr(cysteines[i][0]))),cysteines[i][1],scores[i])
		if scores[i] > 0.5:
			positives.append([cysteines[i][1],cysteines[i][2],cysteines[i][3]])	

	for j in range(len(positives)):
		for k in range(j+1,len(positives)):
			outFile.write(positives[j][0]+positives[k][0]+' '+positives[j][1]+positives[k][1]+' '+positives[j][2]+positives[k][2]+'\n')
	outFile.close()

	return(halfC)



#######################
### Create foo.patA ###
#######################
def createFooPatA(posFileName,window):
	num    = 0
	out    = open('foo.test.pat','w')
	file   = open(posFileName,"r")
	line   = file.readline()
	seqLen = len(string.split(line)[0])
	while line:
		num += 1
		line = file.readline()
	file.close()
  #num is now total number of sequences

	num      = 0
	file     = open(posFileName,"r")
	line     = file.readline()
	while line:
		num += 1
		word = line[:-1]
		sec  = window*2+2
		ind  = (window*2+2)*2
		while ind+80 < (len(line)+1):
			stretch = ''
			for i in range(ind,ind+80,4):
				value    = string.strip(line[i:i+4])
				stretch += "%f " % (float(value)/100.0)
			ind += 80
			out.write(stretch)
		out.write('\n')
		line = file.readline()
	file.close()
	out.close()



#################################
### create ModuleB input file ###
#################################
def createBonds(cysteines):
	outFile   = open('outFileB','w')
	i         = 0
	positives = []
	bonds     = []
	for i in range(0,len(cysteines)):
		positives.append([cysteines[i][1],cysteines[i][2],cysteines[i][3],cysteines[i][0]])
	for j in range(len(positives)):
		for k in range(j+1,len(positives)):
			outFile.write(positives[j][0]+positives[k][0]+' '+positives[j][1]+positives[k][1]+' '+positives[j][2]+positives[k][2]+'\n')
			distance = int(positives[k][3])-int(positives[j][3])
			ID       = '%s - %s' % (positives[j][3],positives[k][3])
			bonds.append([positives[j][0]+'-'+positives[k][0],distance,ID])
	outFile.close()
	return(bonds)



###############################
### Connectivity Prediction ###
###############################
def disulfideBondPrediction(bonds,window,PRINT):
	outFile = open('bondPrediction','w')
	print
	print "Step 4: Disulfide Bonds Prediction using a trained Neural Network"

	### set foo.net for connectivity prediction
#	os.system('cp foo.trained.net_bond'+repr(window)+' foo.trained.net')
	os.system('cp foo.trained.net_2.18.05 foo.trained.net')
#	os.system('cp foo.trained.net.3_stored foo.trained.net')

	### create foo.pat
	createFooPatB('outFileB',window)

	### test foo.pat
	os.system('runNetConnectivity foo.test.pat foo.test.res')

	### read results
	results = open('foo.test.res','r').readlines()
	scores  = []
	for line in results:
		NNOutput = float(line)
		scores.append(NNOutput)

	### print results in a table
	i         = 0
	positives = []
	print
	print 'Sequence   Sequence              Bond               Score'
	print 'position   Distance              ----               -----'
	print '--------   --------'
	for i in range(len(scores)):
		print '%s%s    %s     %s%s     ' % (bonds[i][2],' '*(12-len(repr(bonds[i][2]))),bonds[i][1],' '*(5-len(repr(bonds[i][1]))),bonds[i][0]),
		outFile.write(bonds[i][0]+' '+repr(scores[i])+'\n')
		if scores[i] >= 0.5:
			positives.append(bonds[i][0][:11]+bonds[i][0][12:])
		print '%s' % (scores[i])	
	outFile.close()	



#####################
### createFooPatB ###
#####################
def createFooPatB(posFileName,window):
	num    = 0
	file   = open(posFileName,"r")
	out    = open('foo.test.pat','w')
	line   = file.readline()
	seqLen = len(string.split(line)[0])
	while line:
		num += 1
		line = file.readline()
	file.close()
  #num is now total number of sequences

	num    = 0
	file   = open(posFileName,"r")
	line   = file.readline()
	while line:
		num += 1
		word = line[:-1]
		sec  = (window*2+2)*2-1
		ind  = (window*2+1)*4+2
		while ind+80 < len(line)+1:
			str = ''
			for i in range(ind,ind+80,4):
				str += '%f ' % (float(string.strip(line[i:i+4]))/100.0)
			ind += 80
			str += '%s ' % (secStructEncoding[word[sec]])
			sec += 1
			out.write(str)
		out.write('\n')
		line = file.readline()
	file.close()
	out.close()



#########################################
### Gabow's maximum weighted matching ###
#########################################
def weightedMatch(cysteines):
	print
	print "Step 5: Weighted matching"

	### read ModuleB output
	inFile   = open('bondPrediction','r').readlines()

	### create wmatch output file
	outFile  = open('bondPrediction.in','w')

	### initialize stuff
	num      = len(cysteines)
	tot      = 0
	output   = []
	alphabet = []
	bond2num = {}
	num2bond = {}
	count    = {}

	for i in range(1,num+1):
		alphabet.append(i)

	### write wmatch input file
	ind = 0
	for i in inFile:
		bond  = string.split(i)[0]
		bondL = bond[:len(bond)/2]
		bondR = bond[len(bond)/2+1:]
		if bondL not in bond2num.keys():
			bond2num[bondL]=alphabet[ind]
			num2bond[alphabet[ind]]=bondL
			count[alphabet[ind]]=1
			ind = ind+1
		else:
			count[bond2num[bondL]] += 1
		if bondR not in bond2num.keys():
			bond2num[bondR]         = alphabet[ind]
			num2bond[alphabet[ind]] = bondR
			count[alphabet[ind]]    = 0
			ind += 1

	for j in count.keys():
		tot += count[j]

	outFile.write(repr(len(count))+' '+repr(tot)+' U\n')

	for i in inFile:
		bond  = string.split(i)[0]
		bondL = bond[:len(bond)/2]
		bondR = bond[len(bond)/2+1:]
		score = string.split(i)[1]
		if bond2num[bondL] in alphabet:
			del alphabet[alphabet.index(bond2num[bondL])]
			outFile.write(repr(count[bond2num[bondL]])+' '+bondL+' 0 0\n')
		outFile.write(repr(bond2num[bondR])+' '+repr(int(float(score)*100000))+'\n')

	outFile.write('0 1 0 0\n')
	outFile.close()

	### run wmatch
	os.system('/cluster/home/clotelab/weightedMatch/weighted-match/wmatch bondPrediction.in > bondPrediction.out')

	### read wmatch output and print results
	inFile       = open('bondPrediction.out','r').readlines()
	connectivity = ''
	print '\nPredicted bonds:'
	for line in inFile:
		left  = int(string.split(line)[0])
		right = int(string.split(line)[1])
		if left < right and left != 0 and right != 0:
			for k in cysteines:
				if k[1] == num2bond[int(left)]:
					L = repr(k[0])
				elif k[1] == num2bond[int(right)]:
					R = repr(k[0])
			connectivity += '%s-%s, ' % (left,right)
		
			output.append(num2bond[int(left)]+num2bond[int(right)])
			print '%s-%s (%s,%s)' % (L,R,num2bond[int(left)],num2bond[int(right)])

	print '\nPredicted connectivity:'
	print '%s' % (connectivity[:-2])

	return(output)				



############
### main ###
############
def main(window,sel):

	print
	if sel == 'ModuleA':	
		print '###########################################'
		print '### Cysteine Oxidation state prediction ###'
		print '###########################################\n'
	else:
		print '#########################################'
		print '### Disulfide Connectivity prediction ###'
		print '#########################################\n'

	seqName   = checkCys()        #check if the input contains at least 2 cys, and returns the seq name
	runPsipred(seqName,sel)	      #compute the secondary structure, stored in output files
	seq       = readSeq()         #open and store the input sequence
	prof      = readProfile(seq)  #read and store the evolutionary information
	secStruct = readSecStruct()   #read and store secondary structure information

	cysteines = createDataset(seq,prof,secStruct,window)


	if sel == 'ModuleA':
		halfC = disulfideStatePrediction(cysteines,window,1)

	else:
		print "Step 3: Disulfide Oxidation State Prediction"
		halfC = disulfideStatePrediction(cysteines,window,0)
		bonds = createBonds(cysteines)
		if halfC < 2:
			print 'Warning! The number of predicted half-cystines is lower than 2'
		if bonds:
			graph  = disulfideBondPrediction(bonds,window,0)
			output = weightedMatch(cysteines)



##################################################
### Here it reads the input and calls the main ###
##################################################
if __name__ == "__main__":
	if len(sys.argv) != 3:
		print "Usage:  %s  inputSequenceFile  1/2[1:Ox State Prediction; 2:Connectivity prediction]" % (sys.argv[0])
		sys.exit(1)
																 
	### The input sequence is stored in a file
	inputSeq = open('seqFile.faa','w')
	if os.path.exists(sys.argv[1]):
		inputFile = open(sys.argv[1],'r').readlines()
		if len(inputFile) > 1:
			if inputFile[0][0] != '>':
				inputSeq.write('>inputSeq\n')
				startSeq = 0
			else:
				inputSeq.write(inputFile[0])
				startSeq = 1
			for i in range(startSeq,len(inputFile)):
				inputSeq.write(inputFile[i][:-1].upper())
		else:
			inputSeq.write('>inputSeq\n')
			inputSeq.write(inputFile[0].upper())
		inputSeq.close()
	else:
		print "File %f not found" % (sys.argv[1])
		sys.exit(1)	

	### call the main function
	selections = {1:'ModuleA',2:'ModuleBC'}
	sel        = selections.get(int(sys.argv[2]),2)
	window     = 5
	main(window,sel)
	print

	