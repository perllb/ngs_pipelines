#!/usr/bin/env python
# coding: utf-8


# Take reverse complement of index 5 well list provided by iCell8x

# ## Install biopython
# pip install biopython


from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

import numpy as np
import pandas as pd
import csv

import sys, getopt

usage = '>>> Usage: \n> test.py -i <inputfile> -o <outputfile>'

def main(argv):
    # Read WellList file
    wellf = ""
    newName = ""
    
    try:
        opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
    except getopt.GetoptError:
        print(usage)
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print(usage)
            sys.exit()
        elif opt in ("-i", "--ifile"):
            wellf = arg
        elif opt in ("-o", "--ofile"):
            newName = arg
    
    if wellf == "":
        print(">Error: no input file!")
        print(usage)
        sys.exit()
    if newName == "":
        print(">Error: no output file!")
        print(usage)
        sys.exit()
    
            
    print('Input file is: "%s"' % wellf)
    print('Output file is "%s"' % newName)
    
    welldf = pd.read_csv(wellf,sep="\t")

    # Get barcodes
    barcodes = welldf["Barcode"]
    bcSpl = barcodes.str.split("+")

    # Create new barcode list with rev comp i5s
    newBCs = []
    for bc in bcSpl:
        myseq = Seq(bc[1],generic_dna)
        i5_rc = str(myseq.reverse_complement())
        newbc = "%s+%s" % (bc[0],i5_rc)
        newBCs.append(newbc)    

    # Create new wellList df with new barcodes
    newwell = welldf
    newwell["Barcode"] = newBCs

    # Write new well List file
    newwell.to_csv(newName,sep="\t",index=False)

if __name__ == "__main__":
    main(sys.argv[1:])




