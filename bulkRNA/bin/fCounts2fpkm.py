#!/bin/python

import pandas as pd
import numpy as np
import sys, getopt
import os

def main(argv):
    inputfile = ''
    outputfile = ''
    rewrite = False

    usage='> Usage: fCounts2fpkm.py -i <inputfile> [-r (write new count table with modified column names) '
    try:
        opts, args = getopt.getopt(argv,"hi:r",["ifile=", "rewrite="])
    except getopt.GetoptError:
        print(usage)
        sys.exit(2)
    if len(sys.argv) <= 1:
        print("> Error: No input file entered:")
        print(usage)
        sys.exit()
    for opt, arg in opts:
        if opt == '-h':
            print(usage)
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-r", "--rewrite"):
            rewrite=True
        

                  
    outputfile = os.path.splitext(inputfile)[0] + ".fpkm"
    rwoutput =  os.path.splitext(inputfile)[0] + ".csv"
    csvoutput = os.path.splitext(inputfile)[0] + ".mod.csv"
    

    print('Input file is "',inputfile,'"')
    print('Output FPKM file is "',outputfile,'"')
    
    data = pd.read_table(inputfile,header=1,index_col=0)
        
    counts = data.iloc[:,6:]

    lengths = data.iloc[:,4]

    # Remove path and Aligned* suffix for colnames
    colnames = counts.columns
    newcol = [os.path.basename(x) for x in colnames]
    newcol2 = [x.replace('Aligned.sortedByCoord.out.bam','') for x in newcol]
    counts.columns = newcol2

    if rewrite:
        colnames = data.columns
        newcol = [os.path.basename(x) for x in colnames]
        newcol2 = [x.replace('Aligned.sortedByCoord.out.bam','') for x in newcol]
        data.columns = newcol2
        print('Output csv with mod. sample names: "',rwoutput,'"')
        data.to_csv(rwoutput)

        print('Output csv with mod. sample names, without position data: "',csvoutput,'"')
        counts.to_csv(csvoutput)

        

    fpkm = 1000000000*counts.div(counts.sum()).div(lengths,axis="rows")

    fpkm.to_csv(outputfile)

if __name__ == "__main__":
    main(sys.argv[1:])
               

