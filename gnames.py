import pandas as pd
import glob
import numpy as np
import os

def gnames(orgnum,orgstring,base_dir,input):
    print('Creating amino acid string files')
    shared_file = input
    os.chdir(base_dir+'genenames/')
    name_files = glob.glob('*.csv')
    all_amino = []
    col_order = []

    for file in name_files:
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
                print('No data in name file for',name)
                aminos.append('No data')
                names.append('No data')
                missing_names += 1

        all_amino.append(aminos)
        col_order.append(column)

        os.chdir(base_dir+'output/')
        output = pd.DataFrame({'gene':names,'sequence':sequence_names,'amino':aminos},index=np.arange(len(aminos)))
        output = output[['gene','sequence','amino']]
        output.to_csv(file.split('genes',1)[0]+'_aminos.csv',header=False,index=False)

    shared_wnames = pd.DataFrame({'Core Gene Name':names})
    shared_wnames.to_csv('core_names.csv',index=False)

    #join aminos into single string
    print('Creating concatenated amino acid string file')
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