import csv
from datetime import datetime
import os
import glob
import sys
from multiprocessing import Pool as ThreadPool
from itertools import repeat

#Functions
def long_substr(data):
    substr = ''
    if len(data) > 1 and len(data[0]) > 0:
        for i in range(len(data[0])):
            for j in range(len(data[0])-i+1):
                if j > len(substr) and all(data[0][i:i+j] in x for x in data):
                    substr = data[0][i:i+j]
    return substr
def extract(genome,base_dir):
    os.chdir(base_dir+'/rawdata')
    with open(genome,'rU') as csvfile:
        reader = csv.reader(csvfile,delimiter=',',dialect=csv.excel_tab)
        m=list(reader)
        odata=np.array(m).astype('object')
    string = long_substr(odata[1:,0])
    return string
def name_parse(file, shared_file, orgnum, base_dir):
    #fix sequence names
    os.chdir(base_dir+'genenames/')
    name_data = pd.read_csv(file)
    name_data.columns = ['gene', 'sequence', 'amino']
    fixed_names = []
    for name in name_data['sequence'].tolist():
        name = name.replace('SEED:','')
        fixed_names.append(name.replace(' ','_'))
    name_data['sequence'] = fixed_names

    #put together separate amino acid name files
    for n in range(orgnum):
        if shared_file[n].tolist()[0] in name_data['sequence'].tolist():
            column = n

    aminos = []
    names = []
    sequence_names = shared_file[column].tolist()
    name_set = name_data['sequence'].tolist()
    missing_names = 0
    for name in sequence_names:
        if name in name_set:
            aminos.append(name_data[name_data['sequence'] == name]['amino'].tolist()[0])
            names.append(name_data[name_data['sequence'] == name]['gene'].tolist()[0])
        else:
            if verbose:
                print('No data in name file for',name)
            aminos.append('No data')
            names.append('No data')
            missing_names += 1

    os.chdir(base_dir+'output/')
    output = pd.DataFrame({'gene':names,'sequence':sequence_names,'amino':aminos},index=np.arange(len(aminos)))
    output = output[['gene','sequence','amino']]
    output.to_csv(file.split('genes',1)[0]+'_aminos.csv',header=False,index=False)

    return aminos, names, column
def all_names(base_dir, gn_files):
    os.chdir(base_dir+'genenames/')
    namedf_list = []
    for file in gn_files:
        namedf_list.append(pd.read_csv(file, skiprows=1, names=['Name', 'ID', 'Sequence']))
    return namedf_list

#Sort method
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
    if not os.path.exists(base_dir+'shared/'):
        os.makedirs(base_dir+'shared/')
    if not os.path.exists(base_dir+'odata/'):
        os.makedirs(base_dir+'odata/')

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
        print('Sorting files... ')
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

#Common substring method
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
        print('Finding common substrings...')
        pool = ThreadPool(len(genomes))
        orgstring = pool.starmap(extract,zip(genomes,repeat(base_dir)))
        pool.close()
        pool.join()

        #write orgstring file
        os.chdir(base_dir+'substrings/')
        with open('substrings.csv', 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, delimiter=',',quotechar='|', quoting=csv.QUOTE_MINIMAL)
            for line in orgstring:
                writer.writerow([line])

    return orgstring

#Parse method
def parse(orgstring,pivar,pcvar,evar,file,orgnum,base_dir,unwanted):
    #Do loop for each file, apply cutoff and look for potential shared genes

    os.chdir(base_dir+'rawdata/')
    orig = pd.read_csv(file,names=['query','subject','pi','pc','e'],dtype={'pi':float, 'pc':float, 'e':float})

    #apply cutoffs
    new_data = orig[(orig['pi'] >= pivar) & (orig['pc'] >= pcvar) & (orig['e'] <= evar)]

    gene_list = new_data['query'].tolist()
    genes = set(gene_list)
    for gene in genes:
        if gene_list.count(gene) < orgnum:
            new_data = new_data[new_data['query'] != gene]

    os.chdir(base_dir+'odata/')
    new_data.to_csv(file.split('vs',1)[0]+'_occ.csv',header=False,index=False)

    #remove organisms that are unwanted
    if unwanted == True:
        sub_list = new_data['subject'].tolist()
        orgindex = []
        for sub in sub_list:
            for orgname in orgstring:
                if orgname in sub:
                    orgindex.append(new_data[new_data['subject'] == sub].index[0])

        new_data = new_data.ix[orgindex]

    #Create output list with each gene matching the best gene from the other genomes
    out_list=[]
    gene_list = new_data['query'].tolist()
    genes = sorted(list(set(gene_list)))
    for gene in genes:
        gene_data = new_data[new_data['query'] == gene]
        match_list = gene_data['subject'].tolist()
        match_hits = []
        multiples = []
        multiples_names = []
        for orgname in orgstring:
            n = 0
            multi_names = []
            for match in match_list:
                if orgname in match:
                    n += 1
                    multi_names.append(match)
            if n > 1:
                multiples.append(orgname)
                multiples_names.append(multi_names)
            match_hits.append(n)
        #append lists that contain one and only one of each other genome
        if 0 not in match_hits and sum(match_hits) == orgnum:
            for match in sorted(match_list):
                out_list.append(match)
        #select the best of multiples, parse and append
        if 0 not in match_hits and sum(match_hits) > orgnum:
            temp_list = []
            for multiple in multiples:
                sindex = []
                for subject in gene_data['subject'].tolist():
                    if multiple in subject:
                        sindex.append(gene_data[gene_data['subject'] == subject].index[0])
                multi_genes = gene_data.ix[sindex]
                if gene in multi_genes['subject'].tolist():
                    drop_genes = multi_genes[multi_genes['subject'] != gene]['subject'].tolist()
                else:
                    drop_genes = gene_data.ix[sindex].sort(['pi','pc','e','subject'], ascending = [False,False,True,True])['subject'].tolist()[1:]
                for drop_gene in drop_genes:
                    gene_data = gene_data[gene_data['subject'] != drop_gene]
            for match in sorted(gene_data['subject'].tolist()):
                temp_list.append(match)
            if len(temp_list) == orgnum:
                for match in temp_list:
                    out_list.append(match)

    output = np.reshape(out_list,(int(len(out_list)/orgnum),int(orgnum)))
    output = pd.DataFrame(output,columns=[0]*orgnum)

    os.chdir(base_dir+'shared/')
    output.to_csv(file.split('vs',1)[0]+'_shared.csv',header=False,index=False)
    print(file, 'parsed for shared genes')

#shared gene method
def shared(orgnum, base_dir, pivar, pcvar):
    print('Compiling final shared list...')
    if not os.path.exists(base_dir+'output/'):
        os.makedirs(base_dir+'output/')

    os.chdir(base_dir+'shared/')
    shared_files = glob.glob('*.csv')
    shared_cells = []
    shared_count = []

    #find set of all shared genes for each organism across each individual share file
    for number in range(orgnum):
        cells = []
        for file in shared_files:
            temp_data = pd.read_csv(file,names = np.arange(orgnum))
            cells.append(temp_data[number].tolist())

        shared_set = set(cells[0])

        for cell in cells[1:]:
            shared_set.intersection_update(cell)

        shared_set = sorted(list(shared_set))
        shared_count.append(len(shared_set))
        shared_cells.append(shared_set)

    #choose shared set with smallest set of genes to guarentee conservation across all sets
    min_col = shared_count.index(min(shared_count))
    final_shared = shared_cells[min_col]

    #check each shared gene set from each shared file against every other
    #fins shared gene sets with smallest shared set from above
    cells = []
    for file in shared_files:
        temp_data = pd.read_csv(file,names = np.arange(orgnum))
        cells.append(temp_data)

    final_output = []
    for gene in final_shared:
        shared_cells = []
        for cell in cells:
            shared_cells.append(cell[cell[min_col] == gene].values[0])

        shared_set = set(shared_cells[0])
        dif_set = []
        for shared_cell in shared_cells[1:]:
            if list(shared_set.difference(shared_cell)) != []:
                dif_set.append(list(shared_set.difference(shared_cell)))
            shared_set.intersection_update(shared_cell)

        #pass all check shared gene sets to final output
        if len(shared_set) == orgnum:
            final_output.append(shared_set)

    checked_output = pd.DataFrame(columns = np.arange(orgnum),index=np.arange(len(final_output)))
    for n in range(len(final_output)):
        checked_output.loc[n] = [sorted(list(final_output[n]))[i] for i in range(orgnum)]

    print('%s total shared genes'%len(final_output))
    os.chdir(base_dir+'output/')
    if parameter_file == True:
        with open('parameter_file.csv', 'a') as param_file:
            param_file.write('{}, {}, {}\n'.format(pivar, pcvar, len(final_output)))

    checked_output.to_csv('output.csv',header=False,index=False)

    return checked_output

#Gene name matching method
def gnames(orgnum,orgstring,base_dir,input):
    print('Creating amino acid string files...')
    shared_file = input
    os.chdir(base_dir+'genenames/')
    name_files = glob.glob('*.csv')
    all_amino = []
    col_order = []

    pool = ThreadPool(len(name_files))
    out = pool.starmap(name_parse, zip(name_files, repeat(shared_file),repeat(orgnum),repeat(base_dir)))
    for i in range(len(out)):
        all_amino.append(out[i][0])
        col_order.append(out[i][2])
    names = out[0][1]
    pool.close()
    pool.join()

    shared_wnames = pd.DataFrame({'Core Gene Name':names})
    os.chdir(base_dir+'output/')
    shared_wnames.to_csv('core_names.csv',index=False)

    #join aminos into single string
    print('Creating concatenated amino acid string file...')
    amino_strings = []
    for aminos in all_amino:
        string = ''.join(aminos)
        amino_strings.append(string)

    org_order = []
    for column in col_order:
        for string in orgstring:
            if string in shared_file[column].tolist()[0]:
                org_order.append(string)

    output = pd.DataFrame({'organism':org_order,'aminos':amino_strings},index = np.arange(orgnum))
    output = output[['organism','aminos']]
    output.to_csv('all_aminos.csv',header=False,index=False)

#Find unique genes
def unique(orgstring, pivar, pcvar, evar, file, namedf_list):
    if not os.path.exists(base_dir+'output/'):
        os.makedirs(base_dir+'output/')

    os.chdir(base_dir+'rawdata/')
    orig = pd.read_csv(file,names=['query','subject','pi','pc','e'],dtype={'pi':float, 'pc':float, 'e':float})

    #find matching orgstring
    gene_sample = orig.ix[0, 'query']
    for org in orgstring:
        if org in gene_sample:
            gene_org = org

    cut_data = orig[(orig['pi'] >= pivar) & (orig['pc'] >= pcvar) & (orig['e'] <= evar)]
    gene_list = cut_data['query'].tolist()
    uniques_ids = []
    unique_names = []
    duplicates = []

    for gene in gene_list:
        org_count_list = []
        match_list = cut_data[cut_data['query'] == gene]['subject'].tolist()
        if len(match_list) == 1:
            uniques_ids.append(gene)
        else:
            if gene not in duplicates and gene not in uniques_ids:
                count = 0
                for org in orgstring:
                    for match in match_list:
                        if org in match:
                            org_count_list.append(org)
                if org_count_list.count(gene_org) == len(org_count_list):
                    count += 1
                    uniques_ids.append(gene)
                    for match in match_list:
                        if match != gene:
                            duplicates.append(match)

    #find name file
    checked_indices = []
    os.chdir(base_dir+'genenames/')
    for i, name_file in enumerate(namedf_list):
        if i not in checked_indices:
            name_check = name_file.ix[0, 'ID']
            if gene_org in name_check.replace(' ','_'):
                id_list = name_file['ID'].tolist()
                fixed_ids = []
                for id in id_list:
                    fixed_ids.append(id.replace(' ', '_'))
                name_file['Fixed ID'] = pd.Series(fixed_ids, index = name_file.index)
                checked_indices.append(i)
                for id in uniques_ids:
                    if id in fixed_ids:
                        unique_names.append(name_file[name_file['Fixed ID'] == id]['Name'].tolist()[0])
                    else:
                        unique_names.append('No reference name')

    os.chdir(base_dir+'output/')
    uni_df = pd.DataFrame({'Name': unique_names,'ID': uniques_ids},index=list(range(len(uniques_ids))))
    uni_df.to_csv(file+'_uniques.csv',header=False,index=False)
    print(file, 'parsed for unique genes')
    if verbose == True:
        print('{} duplicates'.format(len(duplicates)))
        print('{} total unique genes'.format(len(unique_names)))

startTime = datetime.now()

print('geneparser v1.81')
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

#Add extra commands here
#for designating fewer than full BLAST of files
if '-unw' in sys.argv:
    unwanted = True
else:
    unwanted = False

if '-v' in sys.argv:
    verbose = True
    print('Verbose printout is on')
else:
    verbose = False

if '-u' in sys.argv:
    uniques = True
    print('Finding unique genes is on')
else:
    uniques = False

if '-pof' in sys.argv:
    parameter_file = True
    print('Parameter file output is on')
else:
    parameter_file = False

if '-t' in sys.argv:
    try:
        threads = int(sys.argv[sys.argv.index('-t')+1])
        print('Number of threads set to {}'.format(threads))
    except ValueError:
        print('Number of threads invalid, defaulting to 1')
        threads = 1

print('')
print('Percent identity cutoff set at ',pivar)
print('Percent coverage cutoff set at ',pcvar)
print('Expected value cutoff set at ',evar)
print('')

#run sort
rd_files, orgnum, gn_files, data_dir, base_dir = sort(startTime)
print(datetime.now()-startTime)
print('')

#find common substrings
orgstring = gene_lcs(rd_files,base_dir)
print(datetime.now()-startTime)
print('')

#parse for shared genes
print('Finding shared genes...')
pool = ThreadPool(threads)
pool.starmap(parse, zip(repeat(orgstring),repeat(pivar),repeat(pcvar),repeat(evar),rd_files,repeat(orgnum),repeat(base_dir),repeat(unwanted)))
pool.close()
pool.join()
print(datetime.now()-startTime)
print('')

#Find uniques
if uniques == True:
    print('Finding unique genes...')
    namedf_list = all_names(base_dir,gn_files)
    pool2 = ThreadPool(threads)
    pool2.starmap(unique, zip(repeat(orgstring), repeat(pivar), repeat(pcvar), repeat(evar), rd_files, repeat(namedf_list)))
    pool2.close()
    pool2.join()
    print(datetime.now()-startTime)
    print('')

#find shared genes
output = shared(orgnum,base_dir, pivar, pcvar)
print(datetime.now()-startTime)
print('')

#make shared gene name files
gnames(orgnum,orgstring,base_dir,output)
print('Parsing finished in',datetime.now()-startTime)