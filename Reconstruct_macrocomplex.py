import argparse
from BioMacromplex.PDBaligner import good_chain_names, reconstruct_macrocomplex, delete_overlapping_chains
from BioMacromplex.PDB import ProteinStructure as PS
from random import shuffle
import os




pdb_list = list()

# ---Estos funcionan ---
# folder = 'pdb/proteasoma' #from 1pma.pdb
# folder = 'pdb/deconstruct'#from 1pma.pdb
# folder = 'pdb/nucl'#from 3kuy.pdb
# folder = 'pdb/hemo_deconstruct' #from 1a3n.pdb
# folder = 'pdb/hemoglobin' #from 1a3n.pdb
folder = 'pdb/phosphate'  # from 2f1d.pdb
# folder = 'pdb/actine'#from 6bno.pdb
# folder = 'pdb/capsid2'  # from 5wk1.pdb

# ---Estos no funcionan, demomento ---
# folder = 'pdb/capsid' #from 1gav.pdb
# folder = 'pdb/ATP'  # from 5wk1.pdb




id_list = os.listdir(folder)
for pdb_id in id_list:
    pdb_list.append(PS(pdb_id, '%s/%s' % (folder, pdb_id)))
for pdb in pdb_list:
    pdb.childs[0].id = 'A'
    pdb.childs[1].id = 'A'
shuffle(pdb_list)
if not os.path.exists('tmp'):
    os.mkdir('tmp')
homo_chains = good_chain_names(pdb_list)
new_pdb, chain_id_dict = reconstruct_macrocomplex(pdb_list, homo_chains)
delete_overlapping_chains(new_pdb)
new_pdb.save_to_file('pdb/final.pdb')
print('THE END')

