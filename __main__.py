from datetime import datetime
import sys
import os

from sort import *
from parse import parse
from shared import shared
from gnames import gnames
startTime = datetime.now()

print('geneparser v1.62')
print(' ')

# check for numpy and pandas modules
try:
    import pandas as pd
    try:
        import numpy as np
    except ImportError:
        print('Missing numpy module, please install')
        sys.exit()
except ImportError:
    print('Missing pandas module, please install')
    sys.exit()

# get pivar, pcvar and evar from sys.argv, check if they are in the correct range
try:
    pivar = float(sys.argv[1])
    if pivar < 0 or pivar > 100:
        print('Percent identity must be a number between 0 and 100')
        sys.exit()
except ValueError:
    print('Percent identity must be a number between 0 and 100')
try:
    pcvar = float(sys.argv[2])
    if pcvar < 0 or pcvar > 100:
        print('Percent coverage must be a number between 0 and 100')
        sys.exit()
except ValueError:
    print('Percent coverage must be a number between 0 and 100')
    sys.exit()
try:
    evar = float(sys.argv[3])
    if evar < 0:
        print('Expected value must be greater than 0')
        sys.exit()
except ValueError:
    print('Expected value must be a number greater than 0')
    sys.exit()

if '-unw' in sys.argv:
    unwanted = True
else:
    unwanted = False

print('Percent identity cutoff set at ',pivar)
print('Percent coverage cutoff set at ',pcvar)
print('Expected value cutoff set at ',evar)
print('')

#run sort
rd_files, orgnum, gn_file, data_dir, base_dir = sort(startTime)
print(datetime.now()-startTime)
print('')

#find common substrings
orgstring = gene_lcs(rd_files,base_dir)
print(datetime.now()-startTime)
print('')

#parse for shared genes
parse(orgstring,pivar,pcvar,evar,rd_files,orgnum,base_dir,unwanted)
print(datetime.now()-startTime)
print('')

#find shared genes
output = shared(orgnum,base_dir)
print(datetime.now()-startTime)
print('')

#make shared gene name files
gnames(orgnum,orgstring,base_dir,output)
print('Parsing finished in',datetime.now()-startTime)