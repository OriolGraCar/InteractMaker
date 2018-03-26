# -*- coding: utf-8 -*-

from Bio.PDB import Superimposer, NeighborSearch
import copy
from Bio.pairwise2 import align, format_alignment
import sys


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
            try:
                atoms_fix.extend(chain_fix.childs[i + count_fix].backbone())
                atoms_mov.extend(chain_mov.childs[i + count_move].backbone())
            except:
                atoms_fix.extend(chain_fix.childs[i + count_fix].backbone(True, 'C'))
                atoms_mov.extend(chain_mov.childs[i + count_move].backbone(True, 'C'))
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

    sys.stderr.write("Internal Error: Limit of valid sequence names reached.\n")
    exit(1)

def good_chain_names(pdb_list):
    '''
      Changed the chain id from the chains of the ProteinStructure objects to have concordancy with their sequence.
      If two chains are equal and are in the same pdb both will be named diferently and it will be tracked down.

      :param pdb_list: a list of ProteinStructure objects from each pdb
      :return: a dict with the id equivalence of homo-interactons (two equal chains)
    '''

    # This function is only done to prevent the case of all chains having the same id or all pdb with the same id combinations
    # also it allows to work with consistent chain id
    homo_dict = dict()
    seq_dict = dict()

    for pdb in pdb_list:
        chain = pdb.childs[0]
        other_chain = pdb.childs[1]
        chain_seq = chain.get_sequence_str()
        other_seq = other_chain.get_sequence_str()

        # Check if id aren't used already in other chains from past pdb's
        if chain.get_id() in seq_dict:
            chain.id = new_unused_id(seq_dict)
        if chain.get_id() == other_chain.get_id():
            other_chain.id = new_unused_id(chain.get_id())
        if other_chain.get_id() in seq_dict:
            id_lists = [x for x in seq_dict]
            id_lists.append(chain.get_id())
            other_chain.id = new_unused_id(id_lists)

        if chain_seq == other_seq:  # It is an homo-interaction (both chains are the same)
            for known_id in seq_dict:  # Check if out sequence is found in another chain id from other pdb
                if seq_dict[known_id] == chain_seq:  # This chain already exist with an other name so we'll change it
                    chain.id = known_id
                    seq_dict[other_chain.get_id()] = 'this name is forbiden'
                    homo_dict[other_chain.get_id()] = known_id
                    break  # As we have found the chain sequence in the known sequences we don't need to keep looking
            if chain.get_id() not in seq_dict:  # Check if out sequence is found in another chain id from other pdb
                seq_dict[chain.get_id()] = chain_seq
                seq_dict[other_chain.get_id()] = 'this name is forbiden'
                homo_dict[other_chain.get_id()] = chain.get_id()
        else:  # The chains are diferent so we'll check separatedly for their possible identicals in other pdb's
            for known_id in seq_dict:
                if seq_dict[known_id] == chain_seq:
                    chain.id = known_id
                elif seq_dict[known_id] == other_seq:
                    other_chain.id = known_id
            if chain.get_id() not in seq_dict:
                seq_dict[chain.get_id()] = chain_seq
            if other_chain.get_id() not in seq_dict:
                seq_dict[other_chain.get_id()] = other_seq
        pdb.restablish_dict()  # Keep child_dict with good ids
    return homo_dict

def delete_nonsymetryc_interactions_from_dict(interactions_dict, homo_chains ):
    '''
    If there are homolog chains, it removes the pointer to one self of the homolog chains. Used mainly to prevent
    homologs having to fulfill both non-symmetrical interactions in the same pdb.

    :param interactions_dict: A dictionary with the chain_id as keys, each of those is a dictionary with tupples of
    interacting residues as key and a pointer to the chain responsible for the interactions as value.
    :param homolog_chains:  A dictionary with id transformation for homolog chains interacting
    '''
    new_dictionary = copy.deepcopy(interactions_dict)
    if len(homo_chains) > 0:
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
    r = dict()
    for pdb in pdb_list:
        for chain in pdb:
            for other_chain in pdb:
                if chain is not other_chain:
                    interacting_res = chain.interacting_residues(other_chain)
                    if interacting_res:
                        interacting_res_tuple = tuple(interacting_res)
                        r.setdefault(chain.get_id(), dict())[interacting_res_tuple] = other_chain
                        if chain.get_id() in homolog_chains_dict:
                            r.setdefault(homolog_chains_dict[chain.get_id()], dict())[interacting_res_tuple] = other_chain
    # In our initial thoughts we didn't want to add this two lines because they induce many wrong supperimposition of chains
    # that we'll have to filter through clash checking. Nevertheless we didn't find another way to fix the unsolved interactions
    # that appeared on multiple homolog and asymmetric interactions (very common in, for example, vira campside)
    for homo in homolog_chains_dict:
        r[homo] = r[homolog_chains_dict[homo]]
    #delete_nonsymetryc_interactions_from_dict(r, homolog_chains_dict)
    return r

def fast_is_residue_interacting(residue, distance):
    '''
     Checks if the residue given is interacting with another residue from a diferent chain of the same ProteinStructure object
     :param residue: a residue object
     :param distance (int): the max distance you allow to consider an interaction
     :return: boolean
     '''
    if residue is not None:
        pdb = residue.parent.parent
        other_atoms = []
        for chain in pdb:
            if chain is not residue.parent:
                other_atoms.extend(chain.get_atoms_list())
        ns = NeighborSearch(other_atoms)
        for atom in residue:
            result = ns.search(center= atom.get_coord(), radius=distance)
            if len(result) > 0:
                return True
    return False

def reconstruct_macrocomplex(PDB_list, homolog_chains_dict = {}, verbose = False, steps = False):
    '''
    Builds a macrocomplex pdb parting from a list of the interactions from a pdb (protein structure object)
     :param PDB_list: A list of pdb (Protein Structure objects)
     :param homolog_chains_dict: Default = Empty dict. A dictionary with id transformation for homolog interactions
     :return: A pdb object with all the interactions completed, a dictionary with the change of names done in the process
    '''

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
    interactions_dict = all_interactions(PDB_list, homolog_chains_dict)#All known interactions that each chain should have
    new_pdb = copy.deepcopy(PDB_list[0]) #Use one of the inputs as the base to build the macrocomplex
    new_pdb.id = 'Macrocomplex'
    chain_id_dict = {new_pdb.childs[0].get_id(): new_pdb.childs[0].get_id(),
                        new_pdb.childs[1].get_id(): new_pdb.childs[1].get_id()}
    if verbose:
        print('The %s.pdb is used as seed to construct the Macrocomplex' %PDB_list[0].get_id())
    runing = True
    tmp_count = 0
    while runing:  # While at least one chain has a missing interaction /main loop? what if we run out of names?/
        completed_chain = 0
        for chain in new_pdb:
            completed_borders = 0
            for border in interactions_dict[chain_id_dict[chain.get_id()]]:
                # Look in the residues known to interact if they actually do (lax is with a broader margin, to have into account litle deviations in the fit
                confirmed_residues, lax_residues = check_interactions(chain, border)

                # If we have un-interacting atoms or the fitting was too bad it tries to fill the border with a chain from the input
                ok = True
                if confirmed_residues + lax_residues < len(border) or (lax_residues > confirmed_residues and confirmed_residues < (len(border) / 2)):
                    tmp_count, ok = fill_interaction(chain, homolog_chains_dict, interactions_dict, chain_id_dict, new_pdb, border, tmp_count, verbose = verbose, save_steps= steps)
                if ok: #If we filled the interaction correctly
                    completed_borders += 1
            if completed_borders == len(interactions_dict[chain_id_dict[chain.get_id()]]):#if we filled all necesary borders
                completed_chain += 1
        if completed_chain == len(new_pdb):#if all chains have all interactions filled.
            runing = False
    return new_pdb, chain_id_dict

def check_interactions(chain,border):
    '''
    Look if the @chain has all the aminoacids in @border interacting with another chain
    :param chain: a Chain instance inside a ProteinStructure
    :param border: a tupple of aminoacid numbers (we'll use it as an identificator)
    :return:confirmed_residues: number of @chain residues are interacting correctly,
            lax_residue: number of @chain residues are interacting but too far than they should to (bad superimpossing?)
    '''
    confirmed_residues, lax_residues = 0, 0
    for residue_number in border:
        if fast_is_residue_interacting(chain.get_residue_by_num(residue_number), 4):
            confirmed_residues += 1
        elif fast_is_residue_interacting(chain.get_residue_by_num(residue_number), 5):
            lax_residues += 1
    return confirmed_residues, lax_residues

def fill_interaction(chain, homolog_chains_dict, interactions_dict, chain_id_dict, new_pdb, border, tmp_count, verbose = False, save_steps = False):
    '''
    :param chain: Chain object
    :param homolog_chains_dict:A dictionary with id transformation for homolog interactions
    :param interactions_dict:  A dictionary with the chain id from the basic unique chains and the homologs as primary key.
         Each one is a dict with sets of numbers (interacting residues) as key and a pointer to a chain object responsible to
         that interaction
    :param chain_id_dict: A dictionary with the name transformation for chain_id in the new_pdb and the original ones
    :param new_pdb: ProteinStructure object
    :param border: Tupple of residue numbers that need to be interacting with another residues but they arent
    :param tmp_count: a count for the temporal pdb names that track the process everytime a new chain is added
    :return: tmp_count, boolean ( if the interaction has correctly been filled)
    '''
    if chain_id_dict[chain.get_id()] in homolog_chains_dict: #Use chains that are involved in homointeractions as the original
        for i in range(20):  # Fitting 100 times (too much? it's too quick and we don't want deviations to accumulate)
            superimpose_pdb_by_chain(chain, interactions_dict[chain_id_dict[chain.get_id()]][border].parent[
                homolog_chains_dict[chain_id_dict[chain.get_id()]]])
    else:
        for i in range(20):  # Fitting 100 times (too much? it's too quick and we don't want deviations to accumulate)
            superimpose_pdb_by_chain(chain, interactions_dict[chain_id_dict[chain.get_id()]][border].parent[
                chain_id_dict[chain.get_id()]])

    # Add the chain to fill the interaction and track the new name if necessary. Try different approaches of new chains to attach
    # and check if they clash, if they do, then pop them (not very efficient but we had the problem with the interaction_dict
    # mentioned in the function all_interactions thus we have to do it this way
    new_name = new_pdb.add_chain(interactions_dict[chain_id_dict[chain.get_id()]][border],
                                    interactions_dict[chain_id_dict[chain.get_id()]][border].get_id(), track_name=True, verbose=verbose)
    if is_chain_clashing(new_name, new_pdb):#if this chain causes clashes we try the other one of the input and pop this one
        new_pdb.childs.pop()
        new_pdb.restablish_dict()
        if verbose:
            sys.stderr.write('The chain %s recently added produced clashes with another chain so it has been deleted\n' %new_name)
        new_name = new_pdb.add_chain(interactions_dict[chain_id_dict[chain.get_id()]][border].parent.get_other_chain(interactions_dict[chain_id_dict[chain.get_id()]][border].get_id()),
                                     interactions_dict[chain_id_dict[chain.get_id()]][border].get_id(),
                                     track_name=True, verbose=verbose)
        if is_chain_clashing(new_name, new_pdb):#If this chain also causes clashes mainly means we superimposed the wrong chain? or the fitting wasn't good
            new_pdb.childs.pop()
            new_pdb.restablish_dict()
            if verbose:
                sys.stderr.write('The chain %s recently added produced clashes with another chain so it has been deleted\n' % new_name)
            new_name = 'NADA'
        elif verbose:
            print('The chain %s from the %s has been added to fulfill the interaction of the following residues %s of chain %s' % (
                interactions_dict[chain_id_dict[chain.get_id()]][border].parent.get_other_chain(interactions_dict[chain_id_dict[chain.get_id()]][border].get_id()).get_id(),
                interactions_dict[chain_id_dict[chain.get_id()]][border].parent.id, border, chain.get_id()))
    if new_name is None: #if no names can be given it means we've ran out of names to give to the chain so the program finishes abruptly
        sys.stderr.write("Last obtained structure is part%s.pdb" % tmp_count)
        exit(1)
    elif new_name == 'NADA':#none of the possible chains to add were good (they had clashes)
        tmp_count, False
    else:#everything is fine
        if save_steps:
            new_pdb.save_to_file('tmp/part%s.pdb' % tmp_count)
            tmp_count += 1
        chain_id_dict[new_name] = interactions_dict[chain_id_dict[chain.get_id()]][border].get_id()
    return tmp_count, True

def is_chain_clashing(chain_name, pdb):
    '''
    :param chain_name (str): chain.id of the chain you want to check clashes of
    :param pdb: a ProteinStructure instance
    :return: The position of the chain that clashes with @chain_name on the pdb.childs list, or False (boolean) if no
    clash found
    '''
    for chain_position in range(len(pdb)):
        other_chain = pdb.childs[chain_position]
        if other_chain.get_id() == chain_name:
            continue
        contador, n_of_atoms = 0, 0
        ns = NeighborSearch(other_chain.get_atoms_list())
        for atom in pdb[chain_name].iter_atoms():
            clash = ns.search(atom.get_coord(), radius=1.5)
            if clash:
                contador += 1
            n_of_atoms += 1
        if contador > n_of_atoms * 0.10:
            return chain_position
    return False
def delete_overlapping_chains(pdb, verbose = False):
    """
    Looks for clashes between chains and if the clash is big it deletes the latest chain added to the pdb
    :param pdb: A Protein Structure object
    """
    for chain in pdb:
        position_to_pop = is_chain_clashing(chain.get_id(), pdb)
        if position_to_pop:
            if verbose:
                sys.stderr.write('Chain %s clashes significantly with chain %s so we will delete the latest\n' %(chain.get_id(), pdb.childs[position_to_pop].get_id()))
            pdb.childs.pop(position_to_pop)

if __name__== '__main__' :
    '''Mas pruebas <3'''