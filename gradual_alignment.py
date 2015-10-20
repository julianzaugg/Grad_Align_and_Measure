"""
Gradual alignment of sequences using alignment methods.
Requires MAFFT, MUSCLE, T-COFFEE installed to path (comment out or change paths 
as necessary)

Will take a fasta file and do a gradual alignment of the sequences, storing the 
length for each alignment and the CLUSTAL formatted alignments.

"""
import argparse
import subprocess

from sequence import *



# Path locations for alignment algorithms, and additional parameters.
# Currently assume mafft, muscle and t_coffee are install to path, but
# comment out or change as desired.
ALN_ALG_PATH_ARGS = { 
    "linsi" : [], #Mafft
    "muscle" : ["-maxiters 50"],
    "t_coffee" : [],
}


def _dir_check(directory):
    if not os.path.exists(directory):
        raise StandardError("The directory %s does not exist or the location \
        	was incorrectly specified" % directory)

def _parse_arguments(my_parser, my_args):
    # Load sequences from fasta file
	input_seqs = read_fasta_file(my_args.input, Protein_Alphabet)


	for s in input_seqs:
		print s.name, len(s)


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Removes gaps from a clustal \
                                            formmated sequence alignment file')
	parser.add_argument('-i', '--input', help='Input alignment file', 
    	required=True)
	parser.add_argument('-o', '--output', help='Output Location', 
    	required=False, default="./")

	parser.add_argument('-alnm', '--alignment_methods', help='Alignment \
		algorithms to use', required=True, choices = ("linsi", "muscle", 
													"t_coffee", "all"))

	parser.add_argument('-r', '--reference', help='Reference sequence name, if \
	 provided will remove gaps in respect to this reference', required=False)

    parser.add_argument('-order_file', help='if an order file is provided, ', required=False)

	# args = parser.parse_args()
	
	input_file = "/Users/julianzaugg/Documents/University/Phd/Projects/" +\
        "Evolutionary_Pathway/Data/ANEH/Sequences/ESTHER/epoxide_fasta.txt"
	args = parser.parse_args(["-i" , input_file, "-alnm", "linsi"])

	_parse_arguments(parser, args)




