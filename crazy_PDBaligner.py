from PDB import ProteinStructure as PS
from Bio.PDB import Superimposer
import copy
from Bio.pairwise2 import align, format_alignment
import pir
from random import shuffle

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
    '''
      Changed the chain id from the chains of the ProteinStructure objects to have concordancy with their sequence.
      If two chains are equal and are in the same pdb both will be named diferently and it will be tracked down.
      :param pdb_list: a list of ProteinStructure objects from each pdb
      :return: a dict with the id equivalence of homo-interactons (two equal chains)
      '''
    homo_dict = dict()
    seq_dict = dict()
    for pdb in pdb_list:
        chain = pdb.childs[0]
        other_chain = pdb.childs[1]
        if chain is not other_chain:
            chain_seq = chain.get_sequence_str()
            other_seq = other_chain.get_sequence_str()
            if chain_seq == other_seq: # It is an homo-interaction (both chains are the same)
                for known_id in seq_dict:#Check if out sequence is found in another chain id from other pdb
                    if seq_dict[known_id] == chain_seq: #Esta cadena ya existe con otro nombre y habra que trackear la conversion del homomero
                        homo_dict.setdefault(other_chain.get_id(),set()).add(known_id)
                        chain.id = known_id
                        seq_dict[other_chain.get_id()] = other_seq
                        break#As we have found the chain sequence in the known sequences we don't need to keep looking
                if chain.get_id() not in seq_dict:#Check if out sequence is found in another chain id from other pdb
                    seq_dict[chain.get_id()] = chain_seq
                    seq_dict[other_chain.get_id()] = chain_seq
                    homo_dict.setdefault(other_chain.get_id(),set()).add(chain.get_id())
            else:#The chains are diferent so we'll check separatedly for their possible identicals in other pdb's
                for known_id in seq_dict:
                    if seq_dict[known_id] == chain_seq:
                        chain.id = known_id
                    elif seq_dict[known_id] == other_seq:
                        other_chain.id = known_id
                if chain.get_id() not in seq_dict:
                    seq_dict[chain.get_id()] = chain_seq
                if other_chain.get_id() not in seq_dict:
                    seq_dict[other_chain.get_id()] = other_seq
        pdb.restablish_dict()
    return homo_dict
def all_interactions(pdb_list, homolog_chains_dict = None):
    '''

    :param pdb_list: A list of pdbs (Protein Structure objects)
    :param homolog_chains_dict: Default = None. A dictionary with id transformation for homolog chains interacting
    :return: A dictionary with the chain_id as keys, each of those is a dictionary with tupples of
    interacting residues as key and a pointer to the chain responsible for the interactions as value.
    '''
    if homolog_chains_dict is not None:
        temporary_dict = dict()
        for homolog_chain in homolog_chains_dict:
            for homolog_other_chain in homolog_chains_dict[homolog_chain]:
                temporary_dict[homolog_other_chain] = homolog_chain
                temporary_dict[homolog_chain] = homolog_other_chain
        homolog_chains_dict = temporary_dict
    r = dict()
    for pdb in pdb_list:
        for chain in pdb:
            for other_chain in pdb:
                if chain is not other_chain:
                    interacting_res_tuple = tuple(chain.interacting_residues(other_chain))
                    r.setdefault(chain.get_id(), dict())[interacting_res_tuple] = other_chain
                    if chain.get_id() in homolog_chains_dict:
                        r.setdefault(homolog_chains_dict[chain.get_id()], dict())[interacting_res_tuple] = other_chain

    return r
def construct_macrocomplex(PDB_list,  homolog_chains_dict = {}):
    '''
    It gets the interactions
    :param PDB_list: A list of pdb (Protein Structure objects)
    :param homolog_chains_dict: Default = None. A dictionary with id transformation for homolog chains interacting
    :return: A pdb object with all the interactions completed
    '''
    interactions_dict = all_interactions(PDB_list, homolog_chains_dict)
    new_pdb = copy.deepcopy(PDB_list[0])
    new_pdb.id = 'Macrocomplex'
    chain_id_dict = {new_pdb.childs[0].get_id(): new_pdb.childs[0].get_id(),
                     new_pdb.childs[1].get_id(): new_pdb.childs[1].get_id()}
    runing = True

    #reversed_homolog_dict = {value : key for key, value in homolog_chains_dict.items()}
    reversed_homolog_dict = dict()
    for homolog_chain in homolog_chains_dict:
        for homolog_other_chain in homolog_chains_dict[homolog_chain]:
            reversed_homolog_dict[homolog_other_chain] = homolog_chain
    while runing:
        completed = 0
        for chain in new_pdb:
            print(chain.get_id())
            intermedias = list()
            for other_chain in new_pdb:
                if chain is not other_chain:
                    tmp_interaction = chain.interacting_residues(other_chain)
                    if tmp_interaction is not None:
                        intermedias.append(tuple(tmp_interaction))
            for border in interactions_dict[chain_id_dict[chain.get_id()]]:
                if border not in intermedias:
                    if chain_id_dict[chain.get_id()] in interactions_dict[chain_id_dict[chain.get_id()]][border].parent and chain_id_dict[chain.get_id()] != interactions_dict[chain_id_dict[chain.get_id()]][border].get_id():
                        for i in range(100):
                            superimpose_pdb_by_chain(chain,interactions_dict[chain_id_dict[chain.get_id()]][border].parent[chain_id_dict[chain.get_id()]])
                    elif reversed_homolog_dict[chain_id_dict[chain.get_id()]] in interactions_dict[chain_id_dict[chain.get_id()]][border].parent:
                        for i in range(100):
                            superimpose_pdb_by_chain(chain,interactions_dict[chain_id_dict[chain.get_id()]][border].parent[reversed_homolog_dict[chain_id_dict[chain.get_id()]]])
                    new_name = new_pdb.add_chain(interactions_dict[chain_id_dict[chain.get_id()].upper()][border],
                                                 interactions_dict[chain_id_dict[chain.get_id()].upper()][
                                                     border].get_id(), track_name=True)
                    chain_id_dict[new_name] = interactions_dict[chain_id_dict[chain.get_id()].upper()][border].get_id()
                    print(chain_id_dict)
                    intermedias.append(border)
                if intermedias.sort() == list(interactions_dict[chain_id_dict[chain.get_id()]]).sort() and \
                                len(intermedias) == len(interactions_dict[chain_id_dict[chain.get_id()]]) and \
                                intermedias.sort(reverse=True) == list(interactions_dict[chain_id_dict[chain.get_id()]]).sort(reverse=True):
                    completed += 1
                    print('Completed %s and completed = %s' %(chain.get_id(), completed))
                    new_pdb.save_to_file('pdb/%s.pdb' %completed)
                    break
        if completed == len(new_pdb):
            runing = False
    return new_pdb, chain_id_dict


if __name__== '__main__' :
    pdb_list = list()
    #id_list = ['AC', 'AB', 'DC', 'AD', 'BC']
    id_list = ['AB', 'AC', 'AD']
    sqs = dict()
    for id in id_list:
        pdb_list.append(PS(id, 'pdb/%s.pdb' %id))
    homo_chains = good_chain_names(pdb_list)
    new_pdb = construct_macrocomplex(pdb_list, homo_chains)
    pir.superimpose_to_pir(new_pdb, pdb_list, 'kk.pir')

    print('THE END')
