from Bio import SeqIO
import os
import glob
from multiprocessing import Pool as ThreadPool
from datetime import datetime
import time
import argparse

def parse_flags():
    parser = argparse.ArgumentParser()
    # handle flags
    parser.add_argument("-t", "--threads",      help="Thread count for multiprocessing (don't specify more than cores available)", type=int, dest="threads", default=1)
    parser.add_argument("-o", "--output",       help="Desired output folder location",                                             type=str, default="")
    parser.add_argument("-d", "--directory",    help="Folder containing raw genbank files, output will appear here unless -o is specified",  type=str, default="./", dest="directory")

    return parser.parse_args()

flags = parse_flags()

def prepare_directories():
    folders = ['genbanks', 'database', 'output', 'fastas']
    for folder in folders:
        if not os.path.exists(folder):
            os.mkdir(folder)

    # move genbank files
    gbk_files = glob.glob('*.gb*')
    for file in gbk_files:
        os.rename(file, os.path.join("genbank", file))

    return

def main():
    start_time = time.time()

    # move to directory
    os.chdir(flags.directory)

    # create necessary directories and move raw files
    prepare_directories()



startTime = datetime.now()
base = os.getcwd()+'/'
print(base)
#make folders
folders = ['/genbanks', '/database', '/output', '/fastas']
for folder in folders:
    if not os.path.exists(base+folder):
        os.mkdir(base+folder)

#fix file names and move them
gbk_files = glob.glob('*.gb*')
for file in gbk_files:
    os.rename(base+file,base+'/genbanks/'+file.replace(' ','_'))
os.chdir(base+'/genbanks')
gbk_files = glob.glob('*.gb*')

seqids = []

for file in gbk_files:
    with open(file, 'r') as file_obj:
        # get sequence name
        for line in file_obj:
            if "ORGANISM" in line:
                organism_name = line.split("ORGANISM")[-1].lstrip().rstrip() + "_" + file.split(".")[0]
                break

        file_obj.seek(0)

        input_gbk = SeqIO.parse(file_obj, 'genbank')
        db_output = open(os.path.join(base+'/fastas/'+'all.fasta'), 'a')
        indv_output = open(os.path.join(base+'/fastas/'+'{}.fasta'.format(organism_name.replace(' ','_'))),'w')
        genename_filename = '{}_genes.csv'.format(organism_name.replace(' ','_'))
        genename_file = open(os.path.join(base+'/output/',genename_filename),'w')
        genename_file.write('Name,locus_tag,translation\n')

        count = 0
        for sequence in input_gbk:
            for feature in sequence.features:
                if feature.type == 'CDS':
                    if feature.qualifiers['locus_tag'][0] not in seqids and "translation" in feature.qualifiers:
                        try:
                            db_output.write('>{}\n{}\n'.format(feature.qualifiers['locus_tag'][0],feature.qualifiers['translation'][0]))
                            indv_output.write('>{}\n{}\n'.format(feature.qualifiers['locus_tag'][0],feature.qualifiers['translation'][0]))
                            genename_file.write('{},{},{}\n'.format(feature.qualifiers['product'][0].replace(',','_'),feature.qualifiers['locus_tag'][0],feature.qualifiers['translation'][0]))
                            count += 1
                        except KeyError:
                            continue
                    seqids.append(feature.qualifiers['locus_tag'][0])
        print(count)

os.system('makeblastdb -in {}/all.fasta -out {}/protblastdb -parse_seqids -dbtype prot'.format(base+'/fastas',base+'/database'))

os.chdir(base+'/fastas')
fasta_files = glob.glob('*.fasta')

blastp_list = []
for file in fasta_files:
    if file != 'all.fasta':
        string = 'blastp -db {} -query {} -outfmt "10 qseqid sseqid pident qcovs evalue" -out {}vsall.csv'.format(base+'/database/protblastdb',base+'/fastas/'+file,base+'/output/'+file.rsplit('.')[0])
        blastp_list.append(string)
pool = ThreadPool(4)
pool.map(os.system, blastp_list)
pool.close()
pool.join()

print(datetime.now()-startTime)