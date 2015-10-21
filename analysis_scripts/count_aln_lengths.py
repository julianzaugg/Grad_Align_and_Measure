import sys
sys.path.append("../")
import os

import numpy as np
from collections import Counter

from sequence import *

f1 = "linsi"
f2 = "muscle"
f3 = "t_coffee"

folders = [f1,f2,f3]
# folders = [f2]
mutation_positions = [215,219,244,249,317,318,349,350]
# Set to true to calculate average entropy across alignment and across specific mutation positions for
# a reduced alignment (trimmed gap positions in respect to reference sequence)
do_entropy = True 

if do_entropy:
	print "Method\tSeqs\tLength\tMeanEnt\tMeanMutEnt\tP215_Ent\tP219_Ent\tP244_Ent\tP249_Ent\tP317_Ent\tP318_Ent\tP349_Ent\tP350_Ent\tLratio"
else: 
	print "Method\tSeqs\tLength\tLratio"
for folder_name in folders:
	for aln_file in os.listdir(folder_name):
		if not aln_file.startswith("."):
			try:
				aln = read_clustal_file(folder_name + "/" + aln_file, Protein_Alphabet)
				if do_entropy:
					trimmed_aln = aln.get_ungapped_using_reference("aspni-hyl1")
					col_entropies = [trimmed_aln.get_shannon_entropy(p) for p in xrange(trimmed_aln.alignlen)]
					mut_pos_entropies = [col_entropies[p-1] for p in mutation_positions]
					trimmed_length = len(trimmed_aln)
					print "%s\t%i\t%i\t%.3f\t%.3f\t%s\t%0.3f" % (folder_name, len(aln), aln.alignlen, 
					np.mean(col_entropies), np.mean(mut_pos_entropies), "\t".join(map(str, mut_pos_entropies)), float(len(aln))/trimmed_aln.alignlen)
				else:
					print "%s\t%s\t%i\t%0.3f" % (folder_name, aln_file.split(".")[0].split("_")[-1], aln.alignlen, float(len(aln))/aln.alignlen)
			except:
				# print aln_file + " is incomplete or corrupt"
				continue
