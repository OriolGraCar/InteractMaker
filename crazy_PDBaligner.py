from PDB import ProteinStructure as PS
from Bio.PDB import Superimposer
import copy
from Bio.pairwise2 import align, format_alignment
import os
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

def new_unused_id(used_ids):
    '''
    :param used_ids: an iterable (list, dict or set) with the keys you don't want
    :return: a new id not in the used_ids
    '''
    possible_id = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789ñÑçÇ'
    for id in possible_id:
        if id not in used_ids:
            return id
def good_chain_names(pdb_list):
    '''
      Changed the chain id from the chains of the ProteinStructure objects to have concordancy with their sequence.
      If two chains are equal and are in the same pdb both will be named diferently and it will be tracked down.

      :param pdb_list: a list of ProteinStructure objects from each pdb
      :return: a dict with the id equivalence of homo-interactons (two equal chains)
      '''

    #This function is only done to prevent the case of all chains having the same id or all pdb with the same id combinations
    #also it allows to work with consistent chain id
    homo_dict = dict()
    seq_dict = dict()

    for pdb in pdb_list:
        chain = pdb.childs[0]
        other_chain = pdb.childs[1]
        chain_seq = chain.get_sequence_str()
        other_seq = other_chain.get_sequence_str()

        #Check if id aren't used already in other chains from past pdb's
        if chain.get_id() in seq_dict:
            chain.id = new_unused_id(seq_dict)
        if chain.get_id() == other_chain.get_id():
            other_chain.id = new_unused_id(chain.get_id())
        if other_chain.get_id() in seq_dict:
            id_lists = [x for x in seq_dict]
            id_lists.append(chain.get_id())
            other_chain.id = new_unused_id(id_lists)

        if chain_seq == other_seq: # It is an homo-interaction (both chains are the same)
            for known_id in seq_dict:#Check if out sequence is found in another chain id from other pdb
                if seq_dict[known_id] == chain_seq: #This chain already exist with an other name so we'll change it
                    chain.id = known_id
                    seq_dict[other_chain.get_id()] = 'this name is forbiden'
                    homo_dict[other_chain.get_id()] = known_id
                    break#As we have found the chain sequence in the known sequences we don't need to keep looking
            if chain.get_id() not in seq_dict:#Check if out sequence is found in another chain id from other pdb
                seq_dict[chain.get_id()] = chain_seq
                seq_dict[other_chain.get_id()] = 'this name is forbiden'
                homo_dict[other_chain.get_id()] = chain.get_id()
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
        pdb.restablish_dict()#Keep child_dict with good ids
    return homo_dict
def delete_nonsymetryc_interactions_from_dict(interactions_dict, homolog_chains):
    '''
    If there are homolog chains, it removes the pointer to one self of the homolog chains. Used mainly to prevent
    homologs having to fulfill both non-symmetrical interactions in the same pdb.

    :param interactions_dict: A dictionary with the chain_id as keys, each of those is a dictionary with tupples of
    interacting residues as key and a pointer to the chain responsible for the interactions as value.
    :param homolog_chains:  A dictionary with id transformation for homolog chains interacting
    '''
    new_dictionary = copy.deepcopy(interactions_dict)
    if len(homolog_chains) > 0:
        for chain in new_dictionary:
            for interaction in new_dictionary[chain]:
                if new_dictionary[chain][interaction].get_id() == chain and chain not in homo_chains:
                    del interactions_dict[chain][interaction]


def all_interactions(pdb_list, homolog_chains_dict = {}):
    '''
    A function that returns a dictionary with the interactions of each chain with the others from the pdbs in the pdb_list
    :param pdb_list: A list of pdbs (Protein Structure objects)
    :param homolog_chains_dict: Default = Empty dictionary. A dictionary with id transformation for homolog chains interacting
    :return: A dictionary with the chain_id as keys, each of those is a dictionary with tupples of
    interacting residues as key and a pointer to the chain responsible for the interactions as value.
    '''
    mirror_homolog_dict = copy.copy(homolog_chains_dict)
    mirror_homolog_dict.update({chain : homolog for homolog, chain in homolog_chains_dict.items()})

    r = dict()
    for pdb in pdb_list:
        for chain in pdb:
            for other_chain in pdb:
                if chain is not other_chain:
                    interacting_res_tuple = tuple(chain.interacting_residues(other_chain))
                    r.setdefault(chain.get_id(), dict())[interacting_res_tuple] = other_chain
                    if chain.get_id() in mirror_homolog_dict:
                        r.setdefault(mirror_homolog_dict[chain.get_id()], dict())[interacting_res_tuple] = other_chain
    delete_nonsymetryc_interactions_from_dict(r, mirror_homolog_dict)
    return r
def is_residue_interacting(residue, distance):
    '''
    Checks if the residue given is interacting with another residue from a diferent chain of the same ProteinStructure object
    :param residue: a residue object
    :param distance (int): the max distance you allow to consider an interaction
    :return: boolean
    '''
    pdb = residue.parent.parent
    for chain in pdb:
        if chain is not residue.parent:
            for other_residue in chain:
                    for my_atom in residue:
                        for other_atom in other_residue:
                            if my_atom - other_atom < distance:
                                return True
    return False
def construct_macrocomplex(PDB_list,  homolog_chains_dict = {}):
    '''
    Builds a macrocomplex pdb parting from a list of the interactions from a pdb (protein structure object)
    :param PDB_list: A list of pdb (Protein Structure objects)
    :param homolog_chains_dict: Default = Empty dict. A dictionary with id transformation for homolog chains interacting
    :return: A pdb object with all the interactions completed, a dictionary with the change of names done in the process
    '''

    interactions_dict = all_interactions(PDB_list, homolog_chains_dict)#All known interactions that each chain should have
    new_pdb = copy.deepcopy(PDB_list[0]) #Use one of the inputs as the base to build the macrocomplex
    new_pdb.id = 'Macrocomplex'
    chain_id_dict = {new_pdb.childs[0].get_id(): new_pdb.childs[0].get_id(),
                     new_pdb.childs[1].get_id(): new_pdb.childs[1].get_id()}

    '''     
        Dictionaries explanation (because its quite messy)
        
    - interactions_dict = A dictionary with the chain id from the basic unique chains and the homologs as primary key. 
    Each one is a dict with sets of numbers (interacting residues) as key and a pointer to a chain object responsible to 
    that interaction            
    
    - chain_id_dict = when we change the name of a chain in our new_pdb because there's already a chain with that name
    it fills the dictionary with the equivalence. To prevent problems from happening, the chains their id are also in the 
    dictionary but doing nothing (e.g. {'A':'A'}) 
    
    -homolog_chains_dict = a dictionary with the following structure {homolog id: original_chain id} used to find the correct 
    chain in the pdb with homolog interactions
    '''


    runing = True
    while runing: #While at least one chain has a missing interaction /main loop? what if we run out of names?/
        completed_chain = 0

        for chain in new_pdb: #It goes chain by chain of the new pdb looking for missing interactions
            completed_borders = 0
            for border in interactions_dict[chain_id_dict[chain.get_id()]]:
                #Look in the residues known to interact if they actually do (lax is with a broader margin, to have into account litle deviations in the fit
                confirmed_residues, lax_residues = 0, 0
                for residue_number in border:
                    if is_residue_interacting(chain.get_residue_by_num(residue_number), 4):
                        confirmed_residues += 1
                    elif is_residue_interacting(chain.get_residue_by_num(residue_number), 5):
                        lax_residues += 1
                if confirmed_residues+lax_residues < len(border) or lax_residues > confirmed_residues: #If we have un-interacting atoms or the fitting was too bad. /how to go back and redo?/
                    if chain_id_dict[chain.get_id()] in homolog_chains_dict:
                        for i in range(100):#Fitting 100 times (too much? it's too quick and we don't want deviations to accumulate)
                            superimpose_pdb_by_chain(chain,interactions_dict[chain_id_dict[chain.get_id()]][border].parent[homolog_chains_dict[chain_id_dict[chain.get_id()]]])
                        print('He usado %s del pdb %s para obtener %s y tapar %s de %s' % (
                            interactions_dict[chain_id_dict[chain.get_id()]][border].parent[homolog_chains_dict[chain_id_dict[chain.get_id()]]].get_id(),
                            interactions_dict[chain_id_dict[chain.get_id()]][border].parent.id,
                            interactions_dict[chain_id_dict[chain.get_id()]][border].get_id(), border, chain.get_id()))
                    else:
                        for i in range(100):#Fitting 100 times (too much? it's too quick and we don't want deviations to accumulate)
                            superimpose_pdb_by_chain(chain,interactions_dict[chain_id_dict[chain.get_id()]][border].parent[chain_id_dict[chain.get_id()]])
                        print('He usado %s del pdb %s para obtener %s y tapar %s de %s' % (
                            interactions_dict[chain_id_dict[chain.get_id()]][border].parent[chain_id_dict[chain.get_id()]].get_id(),
                            interactions_dict[chain_id_dict[chain.get_id()]][border].parent.id,
                            interactions_dict[chain_id_dict[chain.get_id()]][border].get_id(), border, chain.get_id()))
                    #Add the chain to fill the interaction and track the new name if necessary
                    new_name = new_pdb.add_chain(interactions_dict[chain_id_dict[chain.get_id()]][border],
                                                 interactions_dict[chain_id_dict[chain.get_id()]][border].get_id(), track_name=True)
                    new_pdb.save_to_file('tmp/%s.pdb'%new_name)
                    chain_id_dict[new_name] = interactions_dict[chain_id_dict[chain.get_id()].upper()][border].get_id()
                completed_borders += 1
            #Are all the interactions for the chain fullfilled?
            if completed_borders == len(interactions_dict[chain_id_dict[chain.get_id()]]):
                completed_chain += 1
        if completed_chain == len(new_pdb):
            runing = False
    return new_pdb, chain_id_dict


if __name__== '__main__' :
    pdb_list = list()
    #id_list = ['AC', 'AB'] #hemoglobin
    id_list = ['1U', 'MN', 'LY'] #proteasoma

    for id in id_list:
        pdb_list.append(PS(id, 'pdb/%s.pdb' %id))
    if not os.path.exists('tmp'):
        os.mkdir('tmp')
    homo_chains = good_chain_names(pdb_list)
    new_pdb, chain_id_dict = construct_macrocomplex(pdb_list, homo_chains)
    new_pdb.save_to_file('pdb/proteasoma.pdb')
    pir.superimpose_to_pir(new_pdb, pdb_list, 'kk.pir', chain_id_dict)

    print('THE END')
