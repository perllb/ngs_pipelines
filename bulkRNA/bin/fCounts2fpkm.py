

import pandas as pd
import numpy as np
import sys, getopt
import os

def main(argv):
    inputfile = ''
    outputfile = ''
    try:
        opts, args = getopt.getopt(argv,"hi:",["ifile="])
    except getopt.GetoptError:
        print('fCounts2fpkm.py -i <inputfile> ')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('fCounts2fpkm.py -i <inputfile> ')
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
            
    outputfile = os.path.splitext(inputfile)[0] + ".fpkm"
                
    print('Input file is "', inputfile)
    print('Output file is "', outputfile)
    
    data = pd.read_table(inputfile,header=1,index_col=0)
        
    counts = data.iloc[:,6:]

    lengths = data.iloc[:,4]

    # Remove path and Aligned* suffix for colnames
    colnames = counts.columns
    newcol = [os.path.basename(x) for x in colnames]
    newcol2 = [x.replace('Aligned.sortedByCoord.out.bam','') for x in newcol]
    counts.columns = newcol2

    fpkm = 1000000000*counts.div(counts.sum()).div(lengths,axis="rows")

    fpkm.to_csv(outputfile)

if __name__ == "__main__":
    main(sys.argv[1:])
               

