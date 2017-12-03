import os
from collections import defaultdict
import glob
import numpy as np
import dill
import pickle
import datetime
import argparse
import hashlib
import sys
import time

VERSION = "2.1.1"

def verbose_print(text, end="\n"):
    if flags.verbose:
        print(text, end=end)
    sys.stdout.flush()

def parse_flags():
    parser = argparse.ArgumentParser()
    # handle flags
    parser.add_argument("-i", "--percent-identity",      help="Percent identity cutoff for core genome",                type=float, dest="pi", default=90)
    parser.add_argument("-c", "--percent-coverage",      help="Percent coverage cutoff for core genome",                type=float, dest="pc", default=90)
    parser.add_argument("-e", "--expected-value",        help="Expected value cutoff for core genome",                  type=float, dest="ev", default=.0001)
    parser.add_argument("-p", "--pan-genome",            help="Construct pan genome as part of output",                 action="store_true", dest="pan")
    parser.add_argument("-v", "--verbose",               help="Increase output verbosity",                              action="count")
    parser.add_argument("-s", "--save-gp",               help="Save gp file for test runs",                             action="store_true", dest="gp")
    parser.add_argument("-g", "--core-pan-progression",  help="Generate csv of core and pan genome progression sizes",  action="store_true", dest="prog")
    parser.add_argument("-u", "--unique",                help="Generate list of uniques",                               action="store_true", dest="unique")
    parser.add_argument("-l", "--list-values",           help="Append all found values w/ cutoffs to csv for graphing", action="store_true", dest="list")
    parser.add_argument("-w", "--overwrite-list",        help="Overwrite a previous list",                              action="store_true", dest="overwrite")
    parser.add_argument("-x", "--core-pan-intersection", help="Calculate intersection of core and pan genome",          action="store_true", dest="intersect")
    parser.add_argument("-r", "--rule-out-similar",      help="Don't rule out similar genes from the same genome",      action="store_true", dest="rule_out")
    parser.add_argument("-o", "--output",                help="Desired output folder location",                         type=str, default="")
    parser.add_argument("-d", "--input-directory",       help="Folder containing raw data, output will be placed here as well if -o is not passed",  type=str, default="./", dest="directory")

    return parser.parse_args()

def find_best_match(matches):
    if len(matches) == 1:
        return matches[0]["name"]
    else:
        maximum_pi = max(arg["values"]["pi"] for arg in matches)
        max_pi_dicts = sorted([arg for arg in matches if arg["values"]["pi"] == maximum_pi], key=lambda k: k["name"])
        if len(max_pi_dicts) > 1:
            maximum_pc = max(arg["values"]["pc"] for arg in max_pi_dicts)
            max_pc_dicts = sorted([arg for arg in max_pi_dicts if arg["values"]["pc"] == maximum_pc], key=lambda k: k["name"])
            if len(max_pc_dicts) > 1:
                minimum_ev = min(arg["values"]["ev"] for arg in max_pc_dicts)
                min_ev_dicts = sorted([arg for arg in max_pc_dicts if arg["values"]["ev"] == minimum_ev], key=lambda k: k["name"])
                return min_ev_dicts[0]["name"]
            else:
                return max_pc_dicts[0]["name"]
        else:
            return max_pi_dicts[0]["name"]

def meets_cutoffs(match_dict):
    """
    Checks for a cutoff match in the rows
    Ordered to prevent unneccesary comparisons
    """
    cutoffs = {
        "pi": flags.pi,
        "pc": flags.pc,
        "ev": flags.ev
    }
    matches = True
    for key in ["pi", "pc", "ev"]:
        if key == "ev":
            if not match_dict[key] <= cutoffs[key]:
                matches = False
                break
        else:
            if not match_dict[key] >= cutoffs[key]:
                matches = False
                break
    return matches

def pickle_checksum(gp_file, gp_hash):
    """
    Check if the pickle is valid to avoid ingesting bad or malicious data
    """
    # get the saved hash first for comparison
    with open(gp_hash, "r") as file:
        for line in file:
            saved_hash = line.rstrip()
            break

    # hash the pickle file for comparison
    BLOCKSIZE = 65536
    hasher = hashlib.sha256()
    with open(gp_file, 'rb') as file:
        buf = file.read(BLOCKSIZE)
        while len(buf) > 0:
            hasher.update(buf)
            buf = file.read(BLOCKSIZE)
    
    return saved_hash == hasher.hexdigest()

def hash_pickle(parse_bytes):
    """
    Hash the pickled data before write for later checksumming
    """
    return hashlib.sha256(parse_bytes).hexdigest()

def loading_bar(done, max_value, max_length=30):
    bar = "["
    length = int(round(done / max_value % max_length, 1) * max_length)
    for n in range(max_length):
        if n < length:
            bar += "="
        else:
            bar += " "
    return bar + "]"

def sort_files():
    verbose_print("Sorting files... ", end="")

    # directory information
    raw_folder = os.path.exists("rawdata")
    name_folder = os.path.exists("genenames")
    unsorted_files = glob.glob("*.csv")

    if raw_folder and name_folder and not len(unsorted_files):
        verbose_print("already sorted.")
        return
    else:
        if not raw_folder:
            os.mkdir("rawdata")
        if not name_folder:
            os.mkdir("genenames")

        raw_count, name_count = 0, 0
        for file in unsorted_files:
            if "_genes" in file:
                os.rename(file, f"genenames/{file}")
                name_count += 1
            else:
                os.rename(file, f"rawdata/{file}")
                raw_count += 1

        if raw_count != name_count:
            exit("Unequal raw data and gene name files, something is missing")

        verbose_print("sorting finished.")

        return

def parse_gp_file():
    """
    Looks for pickled data file and reads it in if found
    """
    raw_parse = None
    hash_file = os.path.exists(".gphash")
    parse_file = os.path.exists("parse.gp")

    verbose_print("Looking for gp file and checksum... ", end="")
    if parse_file and hash_file:
        valid_pickle = pickle_checksum("parse.gp", ".gphash")
        if valid_pickle:
            verbose_print("found! Loading in...")
            with open("parse.gp", "rb") as raw:
                raw_parse = dill.load(raw)
        else:
            verbose_print("previous data found, but SHA1 hash is invalid, possible corruption.")
            preparsed_data = False
    else:
        verbose_print("none found.")
        preparsed_data = False

    return raw_parse

def ingest_raw_data():
    """
    Converts the raw data files into gp format for parsing
    """
    raw_gene_dict = defaultdict(lambda: defaultdict(dict))
    gene_sets = defaultdict(set)
    gene_names = {}

    # parse raw csvs
    verbose_print("Parsing raw data CSVs...")
    raw_file_list = glob.glob("rawdata/*.csv")
    for n, raw_csv in enumerate(raw_file_list):
        verbose_print(f"{loading_bar(n, len(raw_file_list))} {raw_csv.lstrip('rawdata/')}       ", end="\r")
        with open(raw_csv, "r") as file:
            for line in file:
                parts = line.rstrip().split(",")
                temp = {
                    "pi": float(parts[2]),
                    "pc": float(parts[3]),
                    "ev": float(parts[4])
                }
                if meets_cutoffs(temp):
                    raw_gene_dict[raw_csv][parts[0]][parts[1]] = temp
                gene_sets[raw_csv].add(parts[0])

    #new line for next loading bar
    verbose_print("")

    # parse gene name csvs
    verbose_print("Parsing gene name CSVs...")
    gene_name_file_list = glob.glob("genenames/*.csv")
    for n, gene_name_csv in enumerate(gene_name_file_list):
        verbose_print(f"{loading_bar(n, len(gene_name_file_list))} {gene_name_csv.lstrip('genenames/')}       ", end="\r")
        with open(gene_name_csv, "r") as file:
            for line in file:
                parts = line.rstrip().split(",")
                gene_names[parts[1]] = {"name": parts[0], "amino": parts[2]}

    # combine and write output
    raw_parse = {"raw": raw_gene_dict, "gene_sets": gene_sets, "gene_names": gene_names}
    verbose_print("")
    if flags.gp:
        verbose_print("Writing out to parse.gp...")
        with open("parse.gp", "wb") as raw, open(".gphash", "w") as gphash:
            raw_parse_bytes = dill.dumps(raw_parse)
            hashed = hash_pickle(raw_parse_bytes)
            verbose_print(f"Parse signature {hashed}")
            raw.write(raw_parse_bytes)
            gphash.write(f"{hashed}\n")

    return raw_parse

def gene_parse(raw_parse):
    """
    Actually do the gene parsing.
    This function takes in the custom dictionary data structure of genomes and finds the shared (and possibly pan) genes
    """

    # progression is the list of genomes for this particular parse
    # form the parse structures from the progression
    verbose_print("Processing genomes... ")
    genome_list = list(raw_parse["raw"].keys())
    if flags.prog:
        flags.pan = True
        flags.unique = True
        genome_progression = [genome_list[:n] for n in range(1, len(genome_list)+1)]
    else:
        genome_progression = [genome_list]

    if flags.intersect:
        flags.pan = True
        if flags.prog:
            flags.prog = False
            verbose_print("-g is disabled when -i is enabled.")

    # loop through genomes, break after first loop unless looking for pan genome
    core_pan_counts = defaultdict(list)
    for n, progression in enumerate(genome_progression):
        if flags.prog:
            verbose_print(f"{loading_bar(n, len(genome_progression))} {n + 1} genomes       ", end="\r")
        raw_genes = {organism: raw_parse["raw"][organism] for organism in progression}
        temp_gene_sets = {organism: raw_parse["gene_sets"][organism] for organism in progression}
        shared_gene_groups = set()
        any_shared_gene_groups = set()
        unique = defaultdict(set)
        first_loop = True
        for i, organism in enumerate(progression):
            if not flags.prog and (flags.pan or flags.unique):
                verbose_print(f"{loading_bar(i, len(progression))} {organism.lstrip('rawdata/').rstrip('vsall.csv')}      ", end="\r")
            ruled_out_similar_genes = set()
            for m, gene in enumerate(sorted(list(raw_genes[organism]))):
                if not flags.prog and not flags.pan and not flags.unique:
                    verbose_print(f"{loading_bar(m, len(raw_genes[organism]))} {gene}      ", end="\r")
                if gene not in ruled_out_similar_genes:
                    # create set of all comparisons matching the cutoff
                    matching_gene_set = set(raw_genes[organism][gene].keys())
                    if not flags.rule_out:
                        for ruled_out_gene in matching_gene_set & temp_gene_sets[organism]:
                            ruled_out_similar_genes.add(ruled_out_gene)
                    best_matches = [gene]
                    other_organism_match_count = 0
                    for other_organism in [org for org in progression if org != organism]:
                        gene_intersection = matching_gene_set & temp_gene_sets[other_organism]
                        if gene_intersection:
                            # if there is a match between the cutoff genes and the corpus of the other genome
                            # grab the best possible of the matches
                            other_organism_match_count += 1
                            best_match = find_best_match([{"name": match_gene, "values": raw_genes[organism][gene][match_gene]} for match_gene in gene_intersection])
                            # add to list of best valid matches from each other genome
                            best_matches.append(best_match)

                    # we've found a shared gene if it matches the cutoff in all other organisms
                    # aka total in the progression - 1
                    if first_loop:
                        if other_organism_match_count == len(progression) - 1:
                            # make sure it's not close to the same
                            if not any(len(frozenset(best_matches) & group) for group in shared_gene_groups): 
                                shared_gene_groups.add(frozenset(best_matches))
                    if flags.pan:
                        any_shared_gene_groups.add(frozenset(best_matches))
                    if flags.unique and len(best_matches) == 1:
                        unique[organism].add(best_matches[0])

            first_loop = False
            if not flags.pan and not flags.unique:
                break

        overlap = 0
        if flags.intersect:
            verbose_print("")
            verbose_print("Calculating intersection...")
            for n, group in enumerate(shared_gene_groups):
                verbose_print(f"{loading_bar(n, len(shared_gene_groups))} {n + 1}      ", end="\r")
                for gene in group:
                    for other_group in any_shared_gene_groups:
                        if gene in other_group:
                            overlap += 1

        # add to counts, this only matters if the -g flag is passed
        core_pan_counts["pan"].append(len(any_shared_gene_groups))
        core_pan_counts["core"].append(len(shared_gene_groups))

    verbose_print("")
    return shared_gene_groups, any_shared_gene_groups, core_pan_counts, unique, overlap

def create_output(core_genome, pan_genome, core_pan_counts, unique, overlap_count, raw_parse):
    """
    Creates all the output files
    """

    # change directory
    if flags.output:
        os.chdir(flags.output)

    # create output directory if it doesn't exist
    if not os.path.exists("output"):
        os.mkdir("output")

    # main core output
    # also create single set of all shared for other file creation
    all_shared_set = set()
    with open("output/core_genome.csv", "w") as core_output:
        for core_gene_group in core_genome:
            for gene in core_gene_group:
                all_shared_set.add(gene)
            core_output.write(",".join(sorted(list(core_gene_group))))
            core_output.write("\n")

    # pan output
    pan_genome_processed = []
    with open("output/pan_genome.csv", "w") as pan_output, open("output/pan_genome_concat.csv", "w") as pan_concat:
        for pan_gene_group in pan_genome:
            size = 0
            item_list = []
            for gene_set in sorted(list(raw_parse["gene_sets"])):
                found = False
                for item in pan_gene_group:
                    if item in raw_parse["gene_sets"][gene_set]:
                        found = True
                if found:
                    item_list.append("1")
                    size += 1
                else:
                    item_list.append("0")
            pan_genome_processed.append({"size": size, "items": item_list})
        pan_sorted = sorted(pan_genome_processed, key=lambda k: k["size"], reverse=True)
        # make outputs
        transposed = defaultdict(list)
        gene_list = sorted([item.lstrip("rawdata/").rstrip("vsall.csv") for item in list(raw_parse["gene_sets"])])
        pan_output.write(f"{','.join(gene_list)}\n")
        for line in pan_sorted:
            pan_output.write(f"{','.join(line['items'])}\n")
            for n, item in enumerate(line["items"]):
                transposed[gene_list[n]].append(item)
        for genome in gene_list:
            pan_concat.write(f"{genome},{''.join(transposed[genome])}\n")

    # gene name files, name file and all aminos file
    # trim filenames from raw
    gene_name_files = {}
    file_names_fixed = {}
    for file in raw_parse["gene_sets"]:
        file_name_fixed = file.lstrip('rawdata/').rstrip('vsall.csv')
        gene_name_files[file] = open(f"output/{file_name_fixed}_aminos.csv", "w")
        file_names_fixed[file] = file_name_fixed

    # get genome order when sorted
    # we use this to shorten the search time for correct gene name file for writing
    gene_order = []
    for core_gene_group in core_genome:
        for gene in sorted(list(core_gene_group)):
            for genome in raw_parse["gene_sets"]:
                if gene in raw_parse["gene_sets"][genome]:
                    gene_order.append(genome)
        break

    all_aminos = defaultdict(str)
    with open("output/core_names.csv", "w") as core_names:
        for core_gene_group in core_genome:
            for n, gene in enumerate(sorted(list(core_gene_group))):
                # write to individual gene file
                gene_name_files[gene_order[n]].write(",".join([raw_parse["gene_names"][gene]["name"], gene, raw_parse["gene_names"][gene]["amino"]]))
                gene_name_files[gene_order[n]].write("\n")
                # concat to all aminos
                all_aminos[file_names_fixed[gene_order[n]]] += raw_parse["gene_names"][gene]["amino"]
                # write name file
                if n == 0:
                    core_names.write(raw_parse["gene_names"][gene]["name"])
                    core_names.write("\n")

    # output all aminos
    with open("output/all_aminos.csv", "w") as all_aminos_file:
        for name in sorted(list(all_aminos.keys())):
            all_aminos_file.write(f"{name},{all_aminos[name]}\n")

    if flags.prog:
        with open("output/core_pan.csv", "w") as core_pan_file:
            for n in range(len(core_pan_counts["pan"])):
                core_pan_file.write(f"{core_pan_counts['pan'][n]},{core_pan_counts['core'][n]}\n")

    if flags.list:
        io_value = "a"
        if flags.overwrite:
            io_value = "w"
        with open("output/values_list.csv", io_value) as file:
            file.write(f"{flags.pi},{flags.pc},{flags.ev},{len(core_genome)},{len(pan_genome)},{sum([len(unique[group]) for group in unique])},{overlap_count}\n")

    # close the files
    for file in gene_name_files:
        gene_name_files[file].close()

def main():
    logo = """      _____                 _____                         
     / ____|               |  __ \                        
    | |  __  ___ _ __   ___| |__) |_ _ _ __ ___  ___ _ __ 
    | | |_ |/ _ \ '_ \ / _ \  ___/ _` | '__/ __|/ _ \ '__|
    | |__| |  __/ | | |  __/ |  | (_| | |  \__ \  __/ |   
     \_____|\___|_| |_|\___|_|   \__,_|_|  |___/\___|_|   \n"""

    verbose_print(logo)
    verbose_print(f"GeneParser {VERSION}")
    verbose_print(f"Cutoffs specified:\n - {flags.pi}% identity\n - {flags.pc}% coverage\n - {flags.ev} expected value")
    start_time = time.time()

    # check directories
    if not os.path.exists(flags.directory):
        exit(f"Invalid input directory {flags.directory}")
    if flags.output and not os.path.exists(flags.output):
        exit(f"Invalid output directory {flags.output}")
    
    # change directory
    os.chdir(flags.directory)

    # ingest raw data or load in previous data
    raw_parse = parse_gp_file()
    if not raw_parse:
        # sort files
        sort_files()
        # parse the data
        raw_parse = ingest_raw_data()

    # parse
    core_genome, pan_genome, core_pan_counts, unique, overlap_count = gene_parse(raw_parse)

    # read outs
    verbose_print(f"Core genome size: {len(core_genome)}")
    if flags.pan:
        verbose_print(f"Pan genome size: {len(pan_genome)}")
    if flags.unique:
        verbose_print(f"Unique genome average size: {round(np.mean([len(unique[group]) for group in unique]), 2)}")
    if flags.intersect:
        verbose_print(f"Core/Pan genome intersection size: {overlap_count}")

    # generate core output
    verbose_print("Generating output...")
    create_output(core_genome, pan_genome, core_pan_counts, unique, overlap_count, raw_parse)

    verbose_print(f"Elapsed time: {round(time.time() - start_time, 2)}s")

if __name__ == "__main__":
    # make flags global
    flags = parse_flags()
    
    main()
