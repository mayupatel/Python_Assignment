
Name: Mayuri Patel


#importing the modules for use
import random
#importing regular expression module
import re 


#dictionary of amino acids as a reference.
aa_dict = {'Met':['ATG'], 'Phe':['TTT', 'TTC'], 'Leu':['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'], 'Cys':['TGT', 'TGC'], 'Tyr':['TAC', 'TAT'], 'Trp':['TGG'], 'Pro':['CCT', 'CCC', 'CCA', 'CCG'], 'His':['CAT', 'CAC'], 'Gln':['CAA', 'CAG'], 'Arg':['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'], 'Ile':['ATT', 'ATC', 'ATA'], 'Thr':['ACT', 'ACC', 'ACA', 'ACG'], 'Asn':['AAT', 'AAC'], 'Lys':['AAA', 'AAG'], 'Ser':['AGT', 'AGC', 'TCT', 'TCC', 'TCA', 'TCG'], 'Val':['GTT', 'GTC', 'GTA', 'GTG'], 'Ala':['GCT', 'GCC', 'GCA', 'GCG'], 'Asp':['GAT', 'GAC'], 'Glu':['GAA', 'GAG'], 'Gly':['GGT', 'GGC', 'GGA', 'GGG'], '*':['TAA','TAG','TGA']}


#this function generates a new sequence.
def newSeq(average_seq,freq_AA,finaldict,userchoice):
	
	#takes average value and check for divisible by 3 or not.
	#Done for the codon count as each codon contain three characters
	if average_seq % 3 == 0:
		codon_lt = int(average) #if divisible, then codon length is average
	else:
		#we have to create a new average by taking +-10 values from given average.
		while True:
			a= average_seq-10
			b= average_seq+10
			avgcheck = random.randint(int(a),int(b))
			#checking again for divisible by 3...
			if avgcheck % 3 == 0:
				codon_lt = avgcheck
				break

	#Creating mutliple Seq based on user choice			
	dnaBaselist =[]
	dnaMapList={}
	usercount =0
	#while loop until the condition is met. Userchoice to generate sequence is taken
	while(usercount!=userchoice):
		usercount =usercount+1
		#within a function calling another function for sequence generation.
		baseseq = creatnew_seq(codon_lt,average_seq,freq_AA,finaldict)
		#each sequence is counted and stored as a dictionary
		dnaMapList[usercount] =baseseq
		
	#End of the creating multiple sequence
	return dnaMapList


#generating a sequence code.
def creatnew_seq(codon_lt,average_seq,freq_AA,finaldict):

	#taking list of stop codons
	lstStopCodon = ["TAA","TGA","TAG"]

	codon_lt = int(average_seq/3)
	#substracting the number with 2
	finalLt = codon_lt-2
	
	seq_read = []
	final_seq_read = []
	dnaStr = ""
	for n in range(1,finalLt):
			
			number1 = random.randint(1,10) # selected random Amino acid freq
			for key,value in freq_AA.items():# find the amino acid
				if key == "*": # if stop codon, then continue
					continue
				if number1 ==value:
					for fkey,fvalue in finaldict.items():# for AA we are looking for codon freq
						if key==fkey:
							for fskey,fsvalue in fvalue.items():
								while True:
									number2 =random.randint(0,100)
									if number2 == fsvalue:
										seq_read.append(fskey)
										break
	#appending start codon							
	final_seq_read.append("ATG")
	ctSeq =0
	for m in seq_read:
		ctSeq =ctSeq+1
		final_seq_read.append(m)
		if ctSeq==finalLt:
			#opting random stop codon from given list
			stopcodon = random.choice(lstStopCodon)
			final_seq_read.append(stopcodon)# appending stop codon
			break

	#new sequence is converted into a string.
	dnaBase = dnaStr.join(final_seq_read)
	
	#exit the function by returing dna sequence
	return dnaBase


#this function mutates the sequence
def subsitute_base(dnaBase):

	newNucleotide = ""
	lstStopCodon = "TAA" or "TGA" or"TAG"
	mutate_count = 0
	seqlen = len(dnaBase)

	final_stopcodon = False

	while True:
		substitute_dna = ""
		mutate_count = mutate_count + 1
		
		seqlen = len(dnaBase)
		lst = []

		#random position is selected from the sequence
		position1 = random.randint(0, len(dnaBase)-1)
		firstRandom = dnaBase[position1]
		
		count = 0
		#dna sequence is converted to a list
		nucleotideList = list(dnaBase)# string converted to list
		position2 = random.randint(0, len(nucleotideList)-1)# generated second random number
		#take the base from the second random position
		secondRandom = nucleotideList[position2]
		#replace the position 1 to second random dna base
		nucleotideList[position1] = secondRandom.lower() 
		newNucleotide = "".join(nucleotideList)# list is converted to string
		
		mutatedSeq = newNucleotide.upper()# new mutated sequence 

		#there are three conditions to met when a sequence is mutated
		# First condition
		if position1 ==0 or position1 ==1 or position1 ==2:
			print("FIRST CONDITION - START CODON DISRUPTION MATCH")
			#print("Number OF Iteration took to match The Condition: " + str(mutate_count))
			break

		elif position1 ==(seqlen-1)or position1 ==(seqlen-2) or position1 ==(seqlen-3):
			print("SECOND CONDITION - STOP CODON DISRUPTION MATCH")
			#print("Number OF Iteration took to match The Condition: " + str(mutate_count))
			break

		else:
			#finding the premature stop codon in the sequence
			codon = ""
			count = 0
			mutated_codon = []
			mutate_lt = len(mutatedSeq)
			
			#each time looking for the 3 bases as codon 
			for b in mutatedSeq:
				codon = codon + b
				count = count + 1
				if count == 3:
					mutated_codon.append(codon)
					count = 0
					codon = ""

			#using the regular expression to find stop codon
			mutated_codonLt = len(mutated_codon)
			for i in range(0,mutated_codonLt):

				if i != mutated_codonLt-1:
					#finding the codon using finall function
					if bool(re.findall("[A-Z]",mutated_codon[i])):
						if mutated_codon[i] in lstStopCodon:
							print("THIRD CONDITION - PREMATURE STOP CODON FOUND")
							#print("Number OF Iteration took to match The Condition: " + str(mutate_count))
							final_stopcodon = True
							break
			if final_stopcodon:
				break

	#Exit the code with return of count and sequence
	return mutate_count,mutatedSeq
		


#this function executes amino acid counts and its frequencies.
def aminoA_freq(fastaDict):
	#assigning empty variables
	lis = []
	seq_read = ""
	sumLt = 0
	count = 0

	#looping dictionary of parsed file containing header and sequence
	for tup in fastaDict.values():
		#join all sequences.	
		lis.append(tup)
	
	#joining all the sequences and taking the length of it.
	final_seq = seq_read.join(lis)
	numberBase = len(final_seq)
	

	#count of Amino acids. 
	AA_dict ={}
	
	for b in range(0, numberBase):
		codon = final_seq[ b:b + 3]#sliding window is used
		#reference amino acid dictionary is looped here to find the codon in the given dict
		for key, value in aa_dict.items():
			if codon in value and key in AA_dict.keys():
				AA_dict[key] += 1
			elif codon in value:
				AA_dict[key] = 1
	
	#average/freq of Amino acids.	
	vsum = 0
	Each_val = []
	keyAA = []	
	valueAA = []

	#taking new dictionary created from previous part
	for key, valuelst in AA_dict.items():
		Each_val.append(valuelst) # count is appended into a list
		vsum = vsum + valuelst
		keyAA.append(key)
		
	for val in Each_val:
		average_Aa = int((val/vsum)*100) #calculating average for amino acid frequencies
		valueAA.append(average_Aa)
	
	#zip function zips the key and value into a dictionary
	freq_AA = dict(zip(keyAA,valueAA))
	
	#exit the function by returning the AA frequencies.
	return freq_AA



#this function executes codon counts and frequencies.
def codon_freq(fastaDict):
	lis = []
	seq_read = ""

	#joins all the sequences into one.
	for tup in fastaDict.values():
		lis.append(tup)
		#joining all the sequences and taking the length of it.
		final_seq = seq_read.join(lis)
		numberBase = len(final_seq)
	
	#codon count
	newdict = {}	
	for b in range(0, numberBase):
		codon = final_seq[ b:b + 3]#sliding window is used
		
		for key, value in aa_dict.items():
			if codon in value and codon in newdict.keys():
				#reference amino acid dictionary is looped here to find the codon in the given dict
				newdict[codon] += 1
			elif codon in value:
				newdict[codon] = 1


	#frequencies of codons.
	finaldict = {}
	#using the reference amino acids
	for key, valuelst in aa_dict.items():
		freqDict = {}
		count = 0
		vsum = 0
		#taking values from the dictionary as looped them again
		for value in valuelst:
			count +=1
			#comparing the value to key in the dictionary 
			if value in newdict.keys():
				vsum = vsum + newdict[value] #summing the number count.
		if count == len(valuelst):
			for v in valuelst:
				if v in newdict.keys():
					#calculating the frequency of codon for same Amino acid.
					freqDict[v] = int((newdict[v]/vsum)*100)
		finaldict[key] = freqDict
	
	#exit the function by returing the codons and their frequencies.
	return finaldict




#this function parses fasta file and
#concatenate the sequence reads.
def read_fasta(appleGen):
	#this uses yield as a generator.

	#assigning empty variables
	header,sequence = None,[]
	
	#looping each line of the file
	for line in appleGen:
		#replacing the new line character
		line = line.rstrip().replace('\n','')
		
		#parsing for header and then joining its particular sequence
		if line.startswith(">"):
			if header:
				#yielding header and sequence
				yield (header,"".join(sequence))
				
			header,sequence =line,[]
		else:
			#all the sequence line is appended into one list for each header
			sequence.append(line)
 
	if header:
		yield (header,"".join(sequence))
	#exit the function by yielding header and sequence into the main function




#The main function executes the code.
def main():

	#Open and read the fasta file
	#whole of the dataset file is renamed to a shorter version
	with open("Mdomestica_491_v1.1.cds_primaryTranscriptOnly.fa","r") as appleGen:
		
		#assigning the empty variables
		count = 0
		seq_list =[]

		#generator is used to parse fasta file.
		fastaDict = {}	

		#each time yielding the header and sequence from a file	
		for header,sequence in read_fasta(appleGen):
			#building a dictionary of header and sequence
			fastaDict[header] = sequence

		#calculating the average score of the coding sequences
			#calculate the length and append into a list to sum it up.
			seq_lt = len(sequence)
			seq_list.append(seq_lt)
			count = count + 1

		#average is sum of the sequence length divided by total number of sequences
		average_seq = int(sum(seq_list)/count)


		#function takes parsed fasta file and returns the output of AA frequencies.
		freq_AA = aminoA_freq(fastaDict)
		

		#function takes parsed fasta file and returns the output of codon frequencies.
		finaldict = codon_freq(fastaDict)
		

		
		#this function takes all the arguments,
		#and return new sequence into the main function.
		#providing the userinput to generate new sequence
		userchoice = int(input("How Many Sequences You Want To Generate(in number): "))
		#function takes called here
		dnaMapList = newSeq(average_seq,freq_AA,finaldict,userchoice)
		#for loop the new sequences generated for key and value
		for key, value in dnaMapList.items():
			#format the output of it to better understanding for the user
			print(str(key)  + "-" +  str(value))
		
		# ask the user which seq you want to mutate
		userchoiceSeq = int(input("Which One of the Sequences You Want To Mutate, Please Select Number As Shown: " ))
		
		#loop until the condition is satisfied
		while True:
			#proper user input is provided then sequence is mutated and breaked the code
			if userchoiceSeq <= userchoice: 
				selectUserSeq = dnaMapList.get(userchoiceSeq)
				break
			else:
				#if user won't provide the proper number as it is shown in the output, error is thrown
				userchoiceSeq = int(input("Invalid Choice, Please Select The Number From The Output: " ))
				selectUserSeq = dnaMapList.get(userchoiceSeq)
				continue

		# Pass User selected Sequence here to use it in substitute method
		# mutating a sequence
		#returns the count and mutated sequence for a function
		mutate_count, mutatedSeq = substitute_base(selectUserSeq)
		print("The Mutated Sequence: " + str(mutatedSeq))
		print("Number OF Iteration took to match The Condition: " + str(mutate_count))

		
if __name__ == "__main__":
	main()