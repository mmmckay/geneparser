import sys
import os
import glob
import numpy as np
import csv
from datetime import datetime

def sort(startTime):
    base_dir = os.getcwd()+'/'
    gname_dir = base_dir+'genenames/'
    data_dir = base_dir+'rawdata/'
    os.chdir(base_dir)
    file_list = glob.glob('*.csv')

    #create new directories for sorting
    if not os.path.exists(gname_dir):
        os.makedirs(gname_dir)
    if not os.path.exists(data_dir):
        os.makedirs(data_dir)

    #move files into new folders
    for file in file_list:
        if 'genes' in file.lower():
            os.rename(base_dir+file,gname_dir+file)
        else:
            os.rename(base_dir+file,data_dir+file)

    #check file count
    gn_filecount = len(glob.glob(gname_dir+'*.csv'))
    rd_filecount = len(glob.glob(data_dir+'*.csv'))

    if gn_filecount == rd_filecount and len(file_list) == 0:
        print('Files already sorted')
    elif gn_filecount == rd_filecount and len(file_list) > 0:
        print('Files sorted in ',datetime.now() - startTime)
    elif gn_filecount != rd_filecount:
        print('Unequal sets of files')
        print('%s gene name files to %s raw data files'%(gn_filecount,rd_filecount))
        sys.exit()

    #create file lists
    os.chdir(data_dir)
    rd_files = glob.glob('*.csv')
    os.chdir(gname_dir)
    gn_files = glob.glob('*.csv')
    return rd_files, rd_filecount, gn_files, data_dir, base_dir

def gene_lcs(genomes,base_dir):
    #create substring folder, check if substrings have previously been calculated
    if not os.path.exists(base_dir+'substrings/'):
        os.makedirs(base_dir+'substrings/')
    os.chdir(base_dir+'substrings/')

    #import previous substring file if it exists
    substring_file = glob.glob('*.csv')
    orgstring = []
    if len(substring_file) == 1:
        with open('substrings.csv', newline='') as csvfile:
            reader = csv.reader(csvfile, delimiter=',', quotechar='|')
            for row in reader:
                orgstring.append(row[0])

    if len(orgstring) == len(genomes):
        print('Organism substrings already calculated')
    else:
        print('Finding common substrings')
        #define long substring function
        def long_substr(data):
            substr = ''
            if len(data) > 1 and len(data[0]) > 0:
                for i in range(len(data[0])):
                    for j in range(len(data[0])-i+1):
                        if j > len(substr) and all(data[0][i:i+j] in x for x in data):
                            substr = data[0][i:i+j]
            return substr
        orgstring = []
        for genome in genomes:
            os.chdir(base_dir+'rawdata/')
            with open(genome,'rU') as csvfile:
                reader = csv.reader(csvfile,delimiter=',',dialect=csv.excel_tab)
                m=list(reader)
                odata=np.array(m).astype('object')
            orgstring.append(long_substr(odata[1:,0]))
        #write orgstring file
        os.chdir(base_dir+'substrings/')
        with open('substrings.csv', 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, delimiter=',',quotechar='|', quoting=csv.QUOTE_MINIMAL)
            for line in orgstring:
                writer.writerow([line])

    return orgstring