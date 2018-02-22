from PDB4 import ProteinStructure as PS
from Bio.PDB import Superimposer
import copy
import numpy as np
from Bio.pairwise2 import align, format_alignment
from random import shuffle

def find_common_atoms(chain_fix, chain_mov, bb = True):
    seq_fix = chain_fix.get_sequence_str()
    seq_mov = chain_mov.get_sequence_str()
    atoms_fix = list()
    atoms_mov = list()
    '''print('Alignment from the chain %s from the %s pdb and chain %s from %s pdb' %(
         chain_fix.get_id(), chain_fix.parent.id,chain_mov.get_id(), chain_mov.parent.id))'''
    alignment = align.globalxx(seq_fix, seq_mov)
    count_fix = 0
    count_move = 0
    for i in range(len(alignment[0][0])):
        if alignment[0][0][i] == '-':
            count_fix -= 1
        if alignment[0][1][i] == '-':
            count_move -= 1
        elif alignment[0][0][i] != '-' and alignment[0][1][i] != '-':
            atoms_fix.extend(chain_fix.childs[i + count_fix].backbone())
            atoms_mov.extend(chain_mov.childs[i + count_move].backbone())
    #print(format_alignment(*alignment[0]))
    return atoms_fix, atoms_mov


def superimpose_pdb_by_chain(chain_fix, chain_mov):
    sup = Superimposer()
    atoms_fix, atoms_mov = find_common_atoms(chain_fix, chain_mov)
    sup.set_atoms(atoms_fix, atoms_mov)
    chain_mov.parent.transform(rot = sup.rotran[0], tran = sup.rotran[1])



if __name__== '__main__' :
    pdb_list = list()
    id_list = ['AC', 'AB', 'DC', 'AD', 'BC']
    sqs = dict()
    sup = Superimposer()
    new_pdb = None
    for id in id_list:
        pdb_list.append(PS(id, 'pdb/%s.pdb' %id))
    for i in range(10):
        shuffle(pdb_list)
        for pdb1 in pdb_list:
            for pdb2 in pdb_list:
                if pdb1 is pdb2:
                    break
                for chain1 in pdb1:
                    for chain2 in pdb2:
                        if chain1.get_id() == chain2.get_id():
                            superimpose_pdb_by_chain(chain1, chain2)
    for pdb in pdb_list:
        for chain in pdb:
            if new_pdb is None:
                new_pdb = pdb
                break
            else:
                new_pdb.add_chain(chain, chain.get_id())
    new_pdb.save_to_file('pdb/junto.pdb')
    print(new_pdb['A'].interacting_residues(new_pdb['B']))
    print(new_pdb['B'].interacting_residues(new_pdb['A']))
    print('THE END')
