import pandas as pd
import os
import numpy as np

def parse(orgstring,pivar,pcvar,evar,rd_files,orgnum,base_dir,unwanted):
    #Do loop for each file, apply cutoff and look for potential shared genes
    if not os.path.exists(base_dir+'shared/'):
        os.makedirs(base_dir+'shared/')
    if not os.path.exists(base_dir+'odata/'):
        os.makedirs(base_dir+'odata/')

    for file in rd_files:
        os.chdir(base_dir+'rawdata/')
        print('Parsing',file)
        orig = pd.read_csv(file,names=['query','subject','pi','pc','e'],dtype={'pi':float, 'pc':float, 'e':float})

        #apply cutoffs
        print('Applying cutoffs')
        new_data = orig[(orig['pi'] >= pivar) & (orig['pc'] >= pcvar) & (orig['e'] <= evar)]

        gene_list = new_data['query'].tolist()
        genes = set(gene_list)
        for gene in genes:
            if gene_list.count(gene) < orgnum:
                new_data = new_data[new_data['query'] != gene]

        os.chdir(base_dir+'odata/')
        new_data.to_csv(file.split('vs',1)[0]+'_occ.csv',header=False,index=False)

        #remove organism that are unwanted
        if unwanted == True:
            sub_list = new_data['subject'].tolist()
            orgindex = []
            for sub in sub_list:
                for orgname in orgstring:
                    if orgname in sub:
                        orgindex.append(new_data[new_data['subject'] == sub].index[0])

            new_data = new_data.ix[orgindex]

        #Create output list with each gene matching the best gene from the other genomes
        print('Compiling shared genes')
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

        output = np.reshape(out_list,(len(out_list)/orgnum,orgnum))
        output = pd.DataFrame(output,columns=[0]*orgnum)

        os.chdir(base_dir+'shared/')
        output.to_csv(file.split('vs',1)[0]+'_shared.csv',header=False,index=False)
