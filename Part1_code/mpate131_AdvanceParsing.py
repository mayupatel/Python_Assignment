
# Name: Mayuri Patel


#amino acids as reference.
aa_dict = {'Met':['ATG'], 'Phe':['TTT', 'TTC'], 'Leu':['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'], 'Cys':['TGT', 'TGC'], 'Tyr':['TAC', 'TAT'], 'Trp':['TGG'], 'Pro':['CCT', 'CCC', 'CCA', 'CCG'], 'His':['CAT', 'CAC'], 'Gln':['CAA', 'CAG'], 'Arg':['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'], 'Ile':['ATT', 'ATC', 'ATA'], 'Thr':['ACT', 'ACC', 'ACA', 'ACG'], 'Asn':['AAT', 'AAC'], 'Lys':['AAA', 'AAG'], 'Ser':['AGT', 'AGC', 'TCT', 'TCC', 'TCA', 'TCG'], 'Val':['GTT', 'GTC', 'GTA', 'GTG'], 'Ala':['GCT', 'GCC', 'GCA', 'GCG'], 'Asp':['GAT', 'GAC'], 'Glu':['GAA', 'GAG'], 'Gly':['GGT', 'GGC', 'GGA', 'GGG'], '*':['TAA','TAG','TGA']}


#this function parse fasta file and
#concatenate the sequence read.
def read_fasta(appleGen):

	#assigning empty variables
	header,sequence = None,[]
	
	#looping each line of the file.
	for line in appleGen:	
		#removing the next line character.	
		line = line.rstrip().replace('\n','')

		# parsing each line of the file
		if line.startswith(">"):
			#storing the header into a variable
			if header:
				# yielding each header and sequence
				yield (header,"".join(sequence))
				
			header,sequence =line,[]

		else:
			#appending each sequence into one part
			#sequence += line.rstrip()
			sequence.append(line)

	#yield header, sequence each time into the main function
	if header:
		yield (header,"".join(sequence))



#main function executes the code and handles the flow.
def main():

	# file is open and read
	# this file of apple genes are used
	with open("Mdomestica_491_v1.1.cds_primaryTranscriptOnly.fa","r") as appleGen:
	
		fastaDict = {}
		newdict = {}
		
		# calling the function and yielding HEADER and SEQUENCE.
		for header,sequence in read_fasta(appleGen): 
			#generated a dictionary of the given fasta file
			fastaDict[header] = sequence
		

		lis = []
		seq_read = ""
		#looping the sequence and concatenating
		for tup in fastaDict.values():			
			lis.append(tup)
		final_seq = seq_read.join(lis) # joining the sequences of each header
		
		#length of final concatenating sequence
		numberBase = len(final_seq)
		
		#codon are checked and appended in dictionary
		newdict = {}
		for b in range(0, numberBase):
			codon = final_seq[ b:b + 3]#sliding window is used
			
			#searching the codon within the amino acid dictionary and 
			#building a new codon dictionary	
			for key, value in aa_dict.items():
				if codon in value and codon in newdict.keys():
					newdict[codon] += 1
				elif codon in value:
					newdict[codon] = 1
		
		#amino acid dictionary is created
		finaldict = {}

		#looping amino acid dictionary
		for key, valuelst in aa_dict.items(): 
			freqDict = {}
			count = 0
			vsum = 0
			#again looping the dictionary values
			for value in valuelst:
				count +=1
				if value in newdict.keys():
					vsum = vsum + newdict[value] #summing the number count.
			if count == len(valuelst):
				for v in valuelst:
					if v in newdict.keys():
						#calculating the frequency of codon for same Amino acid.
						freqDict[v] = newdict[v]/vsum 

			finaldict[key] = freqDict
		
		
		# codon usage table is created
		#this gives table format to the output.
		print("AA\t " + "Codon\t " + "Frequency\n----\t------\t----------")
		for fkey,fvalue in finaldict.items():
				for fskey,fsvalue in fvalue.items():
					print(  str(fkey)  +   "    \t "  +  str(fskey)  + "   \t "  + str(fsvalue) )		



if __name__ == "__main__":
	main()