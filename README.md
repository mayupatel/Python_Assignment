# Python_Assignment3

PART 1:  Advanced file parsing lab

For this assignment, you will be reading in this file of apple genes and, based on these coding sequences, generated a codon usage bias table for this species. 

You will probably want to use the dictionary of amino acids and codons for this.

We're looking for the frequency of the occurrence of a codon relative to other codons of the same amino acid.  As such, you will have multiple counts to keep track of.  Think about how you want to keep track of these counts, and what combinations of dictionaries, lists, and tuples are best suited to this task. Plan this BEFORE you start.

Just as important as it is to get these counts and generate a result, you must also consider how best to present these results. Your code should generate a human-readable file that shows your frequencies in a way that is verbose and makes sense.  Your success in presenting the results will be considered just as important as how you arrived to them.


PART 2: 

Using the advanced file parsing lab as a starting point:

 
1. Rewrite your code so that the parsing of the FASTA file is handled.  Update the code that records codon bias frequencies so that is also records overall frequencies of amino acids within all of the sequences, and the average length of each sequence.  At the end of this, you should have a breakdown of the frequency of amino acids, and the frequency of each codon for each amino acid, and average length of coding sequences.

2. Write a function that will, using those frequencies, generate a series of randomized coding sequences.  The regular rules describing a coding sequence should apply (always starts with a start codon, ends in a stop codon) but the composition of each sequence should be randomly generated using the frequencies from step one.  For instance, if Cys had an overall frequency of 12%, with codon frequencies of TGT: 75% and TGC: 25%, then your code should have a 12% chance to pick a Cys codon, and once that happens have a 75% chance for the codon to be TGT. You should build these sequences to have a random length, but center the length around the average you found in step 1.  Keep in mind this length will always be divisible by 3 to maintain frame.

3. Write a function that will take one of your randomized sequences, and introduce a randomized substitution within the sequence.  Select a base position at random along its entire length, and change that base to a different nucleotide. Continue this until 1 of three conditions occurs: the start codon is disrupted, the stop codon is disrupted, or a premature stop codon is introduced.  You should use a regular expression to detect a premature stop codon. Record how many substitutions had to occur before one of these disruptions was detected.  Write code that will allow you to repeat this test a large number of times, and output the results.

 
 
 - For both of the parts, same file is used to write a code. 