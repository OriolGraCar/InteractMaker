from PDB import ProteinStructure as PS
from Bio.PDB import Superimposer
import copy
from Bio.pairwise2 import align, format_alignment
import pir

def find_common_atoms(chain_fix, chain_mov):
    '''
    :param chain_fix: Chain object
    :param chain_mov: Chain object
    :return: A two lists of Atom objects prepared to be superimposed.
    That means the atoms are coherent in their position with the other list.
    '''
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
    '''
    Superimpose the pdb that owns the chain_mov on the pdb that owns the chain_fix
    :param chain_fix: Chain object
    :param chain_mov: Chain object
    '''
    sup = Superimposer()
    atoms_fix, atoms_mov = find_common_atoms(chain_fix, chain_mov)
    sup.set_atoms(atoms_fix, atoms_mov)
    chain_mov.parent.transform(rot = sup.rotran[0], tran = sup.rotran[1])
def good_chain_names(pdb_list):
    homo_dict = dict()
    seq_dict = dict()
    for pdb in pdb_list:
        chain = pdb.childs[0]
        other_chain = pdb.childs[1]
        if chain is not other_chain:
            chain_seq = chain.get_sequence_str()
            other_seq = other_chain.get_sequence_str()
            if chain_seq == other_seq:  #La interaccion es de un homodimero
                for known_id in seq_dict:
                    if seq_dict[known_id] == chain_seq: #Esta cadena ya existe con otro nombre y habra que trackear la conversion del homomero
                        homo_dict.setdefault(other_chain.get_id(),set()).add(known_id)
                        chain.id = known_id
                        seq_dict[other_chain.get_id()] = other_seq
                        break
                if chain.get_id() not in seq_dict:
                    seq_dict[chain.get_id()] = chain_seq
                    seq_dict[other_chain.get_id()] = chain_seq
                    homo_dict.setdefault(other_chain.get_id(),set()).add(chain.get_id())
            else:#Las cadenas son distintas, asi que chequeamos por separado
                for known_id in seq_dict:
                    if seq_dict[known_id] == chain_seq:
                        chain.id = known_id
                    elif seq_dict[known_id] == other_seq:
                        other_chain.id = known_id
                if chain.get_id() not in seq_dict:
                    seq_dict[chain.get_id()] = chain_seq
                elif other_chain.get_id() not in seq_dict:
                    seq_dict[other_chain.get_id()] = other_seq
        pdb.restablish_dict()
    return homo_dict
def all_interactions(pdb_list):
    '''

    :param pdb_list: A list of pdbs (Protein Structure objects)
    :return: A dictionary with the chain_id as keys, each of those is a dictionary with tupples of
    interacting residues as key and a pointer to the chain responsible for the interactions as value.
    '''
    r = dict()
    for pdb in pdb_list:
        for chain in pdb:
            for other_chain in pdb:
                if chain is not other_chain:
                    r.setdefault(chain.get_id(), dict())[tuple(chain.interacting_residues(other_chain))] = other_chain
    return r
def construct_macrocomplex(PDB_list):
    '''
    It gets the interactions
    :param PDB_list: A list of pdb (Protein Structure objects)
    :return: A pdb object with all the interactions completed
    '''
    interactions_dict = all_interactions(PDB_list)
    new_pdb = copy.deepcopy(PDB_list[0])
    new_pdb.id = 'Macrocomplex'
    runing = True
    while runing:
        completed = 0
        for chain in new_pdb:
            intermedias = list()
            for other_chain in new_pdb:
                if chain is not other_chain:
                    intermedias.append(tuple(chain.interacting_residues(other_chain)))
            for border in interactions_dict[chain.get_id()]:
                if border not in intermedias:
                    for i in range(100):
                        superimpose_pdb_by_chain(chain,interactions_dict[chain.get_id()][border].parent[chain.get_id()])
                    new_pdb.add_chain(interactions_dict[chain.get_id()][border], interactions_dict[chain.get_id()][border].get_id())
                    intermedias.append(border)
                if intermedias.sort() == list(interactions_dict[chain.get_id()]).sort() and \
                                len(intermedias) == len(interactions_dict[chain.get_id()]) and \
                                intermedias.sort(reverse=True) == list(interactions_dict[chain.get_id()]).sort(reverse=True):
                    completed += 1
                    break
            if intermedias.sort() == list(interactions_dict[chain.get_id()]).sort() and \
                            len(intermedias) == len(interactions_dict[chain.get_id()]) and \
                            intermedias.sort(reverse=True) == list(interactions_dict[chain.get_id()]).sort(reverse=True):
                completed += 1
        if completed == len(new_pdb):
            runing = False
    return new_pdb


if __name__== '__main__' :
    pdb_list = list()
    id_list = ['AC', 'AB', 'DC', 'AD', 'BC']
    sqs = dict()
    for id in id_list:
        pdb_list.append(PS(id, 'pdb/%s.pdb' %id))

    new_pdb = construct_macrocomplex(pdb_list)
    #r = good_chain_names(pdb_list)
    pir.superimpose_to_pir(new_pdb, pdb_list, 'kk.pir')


    '''
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
    new_pdb.remove_chain(new_pdb['Z'])
    new_pdb.remove_chain(new_pdb['Y'])
    new_pdb.remove_chain(new_pdb['X'])
    new_pdb.remove_chain(new_pdb['W'])
    new_pdb.remove_chain(new_pdb['V'])
    new_pdb.remove_chain(new_pdb['U'])
    new_pdb.save_to_file('pdb/junto.pdb')
    new_pdb.get_mw()
    '''
    print('THE END')
