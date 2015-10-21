"""
Gradual alignment of sequences using alignment methods.
Currently assumes MAFFT, MUSCLE and T-COFFEE installed to path (comment out or change paths
as necessary; T-Coffee typically comes with with these other alignment tools provided)

Will take a fasta file and do a gradual alignment of the sequences, storing the 
CLUSTAL formatted alignments at each step. By gradual I mean - take first 2 sequences and align and save; take first 3
sequences and align and save...etc.,

"""
import argparse
import subprocess
import os
import sys

from sequence import *
from annotation import *

OUT_LOCATION = "./"
ANNOTATIONS = None
ANNOTATION_DICT = None
ORDER_ANNOTATION_NAME = None

# Path locations for alignment algorithms, and additional parameters.
# Currently assume mafft, muscle and t_coffee are install to path, but
# comment out or change as desired.
# {PATH : PARAMETERS}
ALN_ALG_PATH_ARGS = {
    "linsi": ["linsi", "--quiet", "--clustalout"], #Mafft
    "muscle": ["muscle", "-maxiters 50", "-clw", "-clwstrict", "-quiet"],
    "t_coffee": ["t_coffee", "-quiet", "-n_core=4"],
}



def _isfloat(x):
    try:
        a = float(x)
    except ValueError:
        return False
    else:
        return True

def _isint(x):
    try:
        a = float(x)
        b = int(a)
    except ValueError:
        return False
    else:
        return a == b


def _make_dir(name):
    global OUT_LOCATION
    if not os.path.exists(OUT_LOCATION + name):
        os.makedirs(OUT_LOCATION + name)

def _dir_check(directory):
    if not os.path.exists(directory):
        raise StandardError("The directory %s does not exist or the location \
            was incorrectly specified" % directory)

def _parse_arguments(my_parser, my_args):
    global OUT_LOCATION, ANNOTATIONS, ANNOTATION_DICT, ORDER_ANNOTATION_NAME
    OUT_LOCATION = my_args.output

    # Load sequences from fasta file
    input_seqs = read_fasta_file(my_args.input, Protein_Alphabet)

    # If an order file was provided, we order the input sequences
    if my_args.order_file:
        ANNOTATIONS = Annotation(my_args.order_file[0])
        ORDER_ANNOTATION_NAME = my_args.order_file[1]
        _name_column = ANNOTATIONS.get_column("Name")
        _annotation_column = ANNOTATIONS.get_column(ORDER_ANNOTATION_NAME)
        if _isfloat(_annotation_column[0]):
            _annotation_column = map(float, _annotation_column)
        ANNOTATION_DICT = dict(zip(_name_column, _annotation_column))
        input_seqs = sorted(input_seqs, key=lambda x: ANNOTATION_DICT[x.name])

    cnt = my_args.seqnumber
    if my_args.alignment_methods == "all":
        map(_make_dir, ["linsi", "muscle", "t_coffee"])
    else:
        _make_dir(my_args.alignment_methods)
    while cnt < len(input_seqs):
        # ALIGN THE INPUT WITH METHOD
        cur_seqs = input_seqs[:cnt]
        # Write temporary fasta file for infile
        temp_in_file = OUT_LOCATION + "temp_cur_seqs.txt"
        write_fasta_file(temp_in_file, cur_seqs)
        print "Aligning first %i sequences" % len(cur_seqs)
        if my_args.alignment_methods == "linsi":
            # Call aligner
            print "linsi\t%i" % cnt
            command_args = ALN_ALG_PATH_ARGS["linsi"] + ["%s" % temp_in_file, ">", "%slinsi_%i.txt" % (OUT_LOCATION +
                                                                                                        "linsi/", cnt)]
            subprocess.call(" ".join(command_args),
                            shell=True)
        if my_args.alignment_methods == "muscle":
            # Call aligner
            print "muscle\t%i" % cnt
            command_args = ALN_ALG_PATH_ARGS["muscle"] + ["-in", "%s" % temp_in_file, "-out", "%smuscle_%i.txt" %
                                                           (OUT_LOCATION + "muscle/", cnt)]
            subprocess.call(" ".join(command_args),
                            shell=True)
        if my_args.alignment_methods == "t_coffee":
            print "t_coffee\t%i" % cnt
            # Call aligner
            command_args = ALN_ALG_PATH_ARGS["t_coffee"] + ["-infile=%s" % temp_in_file, "-output=aln",
                                                             "-outfile=%st_coffee_%i.txt" % (OUT_LOCATION +
                                                                                             "t_coffee/", cnt)]
            subprocess.call(" ".join(command_args),
                            shell=True)
        if my_args.alignment_methods == "all":
            # Call aligner
            command_args1 = ALN_ALG_PATH_ARGS["linsi"] + ["%s" % temp_in_file, ">", "%slinsi_%i.txt" % (OUT_LOCATION +
                                                                                                        "linsi/", cnt)]
            command_args2 = ALN_ALG_PATH_ARGS["muscle"] + ["-in", "%s" % temp_in_file, "-out", "%smuscle_%i.txt" %
                                                           (OUT_LOCATION + "muscle/", cnt)]
            command_args3 = ALN_ALG_PATH_ARGS["t_coffee"] + ["-infile=%s" % temp_in_file, "-output=aln",
                                                             "-outfile=%st_coffee_%i.txt" % (OUT_LOCATION +
                                                                                             "t_coffee/", cnt)]
            print "linsi\t%i" % cnt
            subprocess.call(" ".join(command_args1),
                            shell=True)
            print "muscle\t%i" % cnt
            subprocess.call(" ".join(command_args2),
                            shell=True)
            print "t_coffee\t%i" % cnt
            subprocess.call(" ".join(command_args3),
                            shell=True)
        cnt += my_args.skip
    os.remove(temp_in_file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Performs gradual alignment of sequences, saving alignments at each '
                                                 'sequence adding step')
    parser.add_argument('-i', '--input', help='Input FASTA file', required=True)
    parser.add_argument('-o', '--output', help='Output Location', required=False, default="./")

    parser.add_argument('-alnm', '--alignment_methods', help='Alignment algorithms to use', required=True,
                        choices=("linsi", "muscle", "t_coffee", "all"))

    parser.add_argument('-skip', '--skip', help='Skip through input sequences, aligning every Nth set',
                        required=False, type=int, default="1")

    parser.add_argument('--order_file', nargs=2, help='if an order file is provided, sequences will be ordered '
                                                        'first; Format should be a tab delimited text file where at '
                                                        'least one column must have a header of "Name" with content '
                                                        'matching sequence name. The user must also nominate another '
                                                        'column to order sequences by, e.g., '
                                                      'header =  "Name   Length" '
                                                      'args = --order_file order_file_name Length',
                                                        required=False)
    parser.add_argument('-sn', '--seqnumber', help='Number of sequences to start alignment from', type=int, default=2)
    # args = parser.parse_args()

    input_file = "/Users/julianzaugg/Documents/University/Phd/Projects/Evolutionary_Pathway/" \
                 "Data/ANEH/Sequences/ESTHER/epoxide_fasta.txt"
    annotation_data = "/Users/julianzaugg/Documents/University/Phd/Projects/Evolutionary_Pathway/Data/ANEH/Sequences/BLAST_results/epoxide_blast_out_simple.txt"
    # args = parser.parse_args(["-i", input_file, "-alnm", "linsi", "-o", "/Users/julianzaugg/Desktop/my_test/", "--order_file", annotation_data, "EValue", "-sn", "179", "-skip", "20"])
    # args = parser.parse_args(["-i", input_file, "-alnm", "muscle", "-o", "/Users/julianzaugg/Desktop/my_test/", "--order_file", annotation_data, "EValue", "-sn", "398", "-skip", "1"])
    # args = parser.parse_args(["-i", input_file, "-alnm", "t_coffee", "-o", "/Users/julianzaugg/Desktop/my_test/", "--order_file", annotation_data, "EValue", "-sn", "246", "-skip", "1"])
    # args = parser.parse_args(["-i", input_file, "-alnm", "all", "-o", "/Users/julianzaugg/Desktop/my_test/", "--order_file", annotation_data, "EValue", "-sn", "86", "-skip", "20"])

    _parse_arguments(parser, args)




