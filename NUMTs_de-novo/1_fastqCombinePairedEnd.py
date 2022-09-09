#!/usr/bin/env python

# FILE: 1_fastqCombinePairedEnd.py --------------------------------------------
#
# DESCRIPTION: Resynchronize 2 fastq or fastq.gz files (R1 and R2) after they 
#              have been trimmed and cleaned. Resynchronization is needed 
#              because sometimes during trimming only one of the reads from 
#              pair is kept, and we need both of them.
# 
# USAGE: python fastqCombinePairedEnd.py R1.fastq R2.fastq None
# OPTIONS: input1 = LEFT  fastq or fastq.gz file (R1)
#          input2 = RIGHT fastq or fastq.gz file (R2)
#          separator = character that separates the name of the read from the 
#          part that describes if it goes on the left or right, usually with 
#          characters '1' or '2'.  The separator is often a space, but could be
#          another character. A space is used by default. If the sequence names
#          do not contain two parts and you want to use the full name info to 
#          pair your sequences, use 'None' (as text) for the separator. 
#          Eg: python fastqCombinePairedEnd.py R1.fastq R2.fastq None
# 
# REQUIREMENTS: python > 2.7, gzip, sys
#               
# BUGS:   ---
# NOTES:  ---
# AUTHOR:  Scipt is orinally published at https://raw.githubusercontent.com/enormandeau/Scripts/master/fastqCombinePairedEnd.py
# COMPANY:  UCL, London, the UK
# VERSION:  2
# CREATED:  13.07.2019
# REVISION: 04.09.2022

# Importing modules
import gzip
import sys

# Parsing user input
try:
    in1 = sys.argv[1]
    in2 = sys.argv[2]
except:
    print(__doc__)
    sys.exit(1)

try:
    separator = sys.argv[3]
    if separator == "None":
        separator = None
except:
    separator = " "

# Defining classes
class Fastq(object):
    """Fastq object with name and sequence
    """

    def __init__(self, name, seq, name2, qual):
        self.name = name
        self.seq = seq
        self.name2 = name2
        self.qual = qual

    def getShortname(self, separator):
        if separator:
            self.temp = self.name.split(separator)
            del(self.temp[-1])
            return separator.join(self.temp)
        else:
            return self.name

    def write_to_file(self, handle):
        handle.write(self.name + "\n")
        handle.write(self.seq + "\n")
        handle.write(self.name2 + "\n")
        handle.write(self.qual + "\n")

# Defining functions
def myopen(infile, mode="r"):
    if infile.endswith(".gz"):
        return gzip.open(infile, mode=mode)
    else:
        return open(infile, mode=mode)

def fastq_parser(infile):
    """Takes a fastq file infile and returns a fastq object iterator
    """
    
    with myopen(infile) as f:
        while True:
            name = f.readline().strip()
            if not name:
                break

            seq = f.readline().strip()
            name2 = f.readline().strip()
            qual = f.readline().strip()
            yield Fastq(name, seq, name2, qual)

# Main
if __name__ == "__main__":
    seq1_dict = {}
    seq2_dict = {}
    seq1 = fastq_parser(in1)
    seq2 = fastq_parser(in2)
    s1_finished = False
    s2_finished = False

    if in1.endswith('.gz'): 
    	outSuffix='.fastq.gz'
    else:
    	outSuffix='.fastq'
    	
    with myopen(in1 + "_pairs_R1" + outSuffix, "w") as out1:
        with myopen(in2 + "_pairs_R2" + outSuffix, "w") as out2:
            with myopen(in1 + "_singles" + outSuffix, "w") as out3:
                while not (s1_finished and s2_finished):
                    try:
                        s1 = seq1.next()
                    except:
                        s1_finished = True
                    try:
                        s2 = seq2.next()
                    except:
                        s2_finished = True

                    # Add new sequences to hashes
                    if not s1_finished:
                        seq1_dict[s1.getShortname(separator)] = s1
                    if not s2_finished:
                        seq2_dict[s2.getShortname(separator)] = s2

                    if not s1_finished and s1.getShortname(separator) in seq2_dict:
                        seq1_dict[s1.getShortname(separator)].write_to_file(out1)
                        seq1_dict.pop(s1.getShortname(separator))
                        seq2_dict[s1.getShortname(separator)].write_to_file(out2)
                        seq2_dict.pop(s1.getShortname(separator))

                    if not s2_finished and s2.getShortname(separator) in seq1_dict:
                        seq2_dict[s2.getShortname(separator)].write_to_file(out2)
                        seq2_dict.pop(s2.getShortname(separator))
                        seq1_dict[s2.getShortname(separator)].write_to_file(out1)
                        seq1_dict.pop(s2.getShortname(separator))
                        
                # Treat all unpaired reads
                for r in seq1_dict.values():
                    r.write_to_file(out3)

                for r in seq2_dict.values():
                    r.write_to_file(out3)