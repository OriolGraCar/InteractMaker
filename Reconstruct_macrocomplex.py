#!/usr/bin/env python3


from argparse import ArgumentParser
from BioMacromplex.PDBaligner import good_chain_names, reconstruct_macrocomplex, delete_overlapping_chains
from BioMacromplex.PDB import ProteinStructure as PS
from random import shuffle
import os
import sys
import re


# ---Estos funcionan ---
# folder = 'pdb/proteasoma' #from 1pma.pdb
# folder = 'pdb/deconstruct'#from 1pma.pdb
# folder = 'pdb/nucl'#from 3kuy.pdb
# folder = 'pdb/hemo_deconstruct' #from 1a3n.pdb
# folder = 'pdb/hemoglobin' #from 1a3n.pdb
# folder = 'pdb/phosphate'  # from 2f1d.pdb
# folder = 'pdb/actine'#from 6bno.pdb
# folder = 'pdb/capsid2'  # from 5wk1.pdb

# ---Estos no funcionan, demomento ---
# folder = 'pdb/capsid' #from 1gav.pdb
# folder = 'pdb/ATP'  # from 5wk1.pdb

parser = ArgumentParser(description='%s is a script that uses the BioMacromplex module to reconstruct Protein, DNA or RNA complexes through a set of interactions in pdb format.' %sys.argv[0])
parser.add_argument("-i", dest= "folder", action = "store", type = str, help = "Path to a folder where pdb are stored", required = True)
parser.add_argument("-o", dest="output_pdb", action = "store", type = str, default='macrocomplex.pdb', help = "Output file where the macrocomplex pdb will be stored")
parser.add_argument('-v', '--verbose', dest = 'verbose', action = "store_true", default = False, help = "Print the progress of the program and the log")
parser.add_argument('-s', '--steps', dest = 'tmp_steps', action = "store_true", default= False, help = "Save a temporary pdb in tmp/ each time a chain is added to track the process")
opt = parser.parse_args()

pdb_list = list()
if not os.path.isdir(opt.folder):
    sys.stderr.write("You didn't provide a folder as an input, please check the usage with %s -h\n" %sys.argv[0])
    exit(1)
is_pdb = re.compile('(.*/(.*)\.pdb|(.*)\.pdb)')
id_list = os.listdir(opt.folder)
for pdb_file in id_list:
    pdb_match = is_pdb.match(pdb_file)
    if pdb_match:
        input_pdb = PS(pdb_match.groups(0)[2], '%s/%s' % (opt.folder, pdb_file))
        if len(input_pdb) == 2:
            pdb_list.append(input_pdb)
        else:
            sys.stderr.write('The inputed %s file does not have two chains, which is a required condition of the input '
                             'files so we will not use it for the macrocomplex reconstruction' %pdb_file)
for pdb in pdb_list:
    pdb.childs[0].id = 'A'
    pdb.childs[1].id = 'B'
shuffle(pdb_list)
if opt.tmp_steps:
    print('All new states of the macrocomplex will be saved in tmp/ with the name of partX, where X is the ordinal step when a chain was added.')
    if not os.path.exists('tmp'):
        os.mkdir('tmp')
homo_chains = good_chain_names(pdb_list)
new_pdb, chain_id_dict = reconstruct_macrocomplex(pdb_list, homo_chains, verbose = opt.verbose, steps = opt.tmp_steps)
delete_overlapping_chains(new_pdb, verbose = opt.verbose)
new_pdb.save_to_file(opt.output_pdb)
print('The END')

