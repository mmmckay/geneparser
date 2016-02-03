import os
import glob
import pandas as pd
import numpy as np

def shared(orgnum, base_dir):
    print('Compiling final shared list')
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
    checked_output.to_csv('output.csv',header=False,index=False)

    return checked_output