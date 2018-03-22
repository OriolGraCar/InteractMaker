from PDB import ProteinStructure as PS
from Bio.PDB import Superimposer, NeighborSearch
import copy
from Bio.pairwise2 import align, format_alignment
import os
import sys
import datetime
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

def reconstruct(PDB_list, homolog_chains_dict = {}):
    interactions_dict = all_interactions(PDB_list, homolog_chains_dict)#All known interactions that each chain should have
    new_pdb = copy.deepcopy(PDB_list[0]) #Use one of the inputs as the base to build the macrocomplex
    new_pdb.id = 'Macrocomplex'
    chain_id_dict = {new_pdb.childs[0].get_id(): new_pdb.childs[0].get_id(),
                        new_pdb.childs[1].get_id(): new_pdb.childs[1].get_id()}
    runing = True
    tmp_count = 0
    while runing:  # While at least one chain has a missing interaction /main loop? what if we run out of names?/
        completed_chain = 0
        for chain in new_pdb:
            completed_borders = 0
            for border in interactions_dict[chain_id_dict[chain.get_id()]]:
                # Look in the residues known to interact if they actually do (lax is with a broader margin, to have into account litle deviations in the fit
                confirmed_residues, lax_residues = check_interactions(chain, border)
                # If we have un-interacting atoms or the fitting was too bad. /how to go back and redo?/
                ok = True
                if confirmed_residues + lax_residues < len(border) or (lax_residues > confirmed_residues and confirmed_residues < (len(border) / 2)):
                    tmp_count, ok = superpose(chain, homolog_chains_dict, interactions_dict, chain_id_dict, new_pdb, border, tmp_count)
                if ok:
                    completed_borders += 1
            if chain_id_dict[chain.get_id()] in homolog_chains_dict:
                if completed_borders == len(interactions_dict[homolog_chains_dict[chain_id_dict[chain.get_id()]]]):
                    completed_chain += 1
            else:
                if completed_borders == len(interactions_dict[chain_id_dict[chain.get_id()]]):
                    completed_chain += 1
        if completed_chain == len(new_pdb):
            runing = False
    return new_pdb, chain_id_dict

def check_interactions(chain,border):
    confirmed_residues, lax_residues = 0, 0
    for residue_number in border:
        if fast_is_residue_interacting(chain.get_residue_by_num(residue_number), 4):
            confirmed_residues += 1
        elif fast_is_residue_interacting(chain.get_residue_by_num(residue_number), 5):
            lax_residues += 1
    return confirmed_residues, lax_residues

def superpose(chain, homolog_chains_dict, interactions_dict, chain_id_dict, new_pdb, border, tmp_count):
    if chain_id_dict[chain.get_id()] in homolog_chains_dict:
        is_homolog = True
        for i in range(200):  # Fitting 100 times (too much? it's too quick and we don't want deviations to accumulate)
            superimpose_pdb_by_chain(chain, interactions_dict[chain_id_dict[chain.get_id()]][border].parent[
                homolog_chains_dict[chain_id_dict[chain.get_id()]]])
        print('He usado %s del pdb %s para obtener %s y tapar %s de %s' % (
            interactions_dict[chain_id_dict[chain.get_id()]][border].parent[
                homolog_chains_dict[chain_id_dict[chain.get_id()]]].get_id(),
            interactions_dict[chain_id_dict[chain.get_id()]][border].parent.id,
            interactions_dict[chain_id_dict[chain.get_id()]][border].get_id(), border, chain.get_id()))

    else:
        for i in range(100):  # Fitting 100 times (too much? it's too quick and we don't want deviations to accumulate)
            superimpose_pdb_by_chain(chain, interactions_dict[chain_id_dict[chain.get_id()]][border].parent[
                chain_id_dict[chain.get_id()]])
        print('He usado %s del pdb %s para obtener %s y tapar %s de %s' % (
            interactions_dict[chain_id_dict[chain.get_id()]][border].parent[chain_id_dict[chain.get_id()]].get_id(),
            interactions_dict[chain_id_dict[chain.get_id()]][border].parent.id,
            interactions_dict[chain_id_dict[chain.get_id()]][border].get_id(), border, chain.get_id()))
    # Add the chain to fill the interaction and track the new name if necessary

    new_name = new_pdb.add_chain(interactions_dict[chain_id_dict[chain.get_id()]][border].parent.get_other_chain(chain_id_dict[chain.get_id()]),
                                    interactions_dict[chain_id_dict[chain.get_id()]][border].get_id(), track_name=True)
    if clash_maker(new_name, new_pdb):
        new_pdb.childs.pop()
        new_pdb.restablish_dict()
        new_name = new_pdb.add_chain(interactions_dict[chain_id_dict[chain.get_id()]][border],
                                     interactions_dict[chain_id_dict[chain.get_id()]][border].get_id(),
                                     track_name=True)
        if clash_maker(new_name, new_pdb):
            new_pdb.childs.pop()
            new_pdb.restablish_dict()
            new_name = 'NADA'

    if new_name is None:
        sys.stderr.write("Last obtained structure is part%s.pdb" % tmp_count)
        exit(1)
    elif new_name == 'NADA':
        tmp_count, False
    else:
        new_pdb.save_to_file('tmp/part%s.pdb' % tmp_count)
        tmp_count += 1
        chain_id_dict[new_name] = interactions_dict[chain_id_dict[chain.get_id()]][border].get_id()
    return tmp_count, True

def clash_maker(chain_name, pdb):
    for other_chain in pdb:
        if other_chain.get_id() != chain_name:
            contador, n_of_atoms = 0, 0
            ns = NeighborSearch(other_chain.get_atoms_list())
            for atom in pdb[chain_name].iter_atoms():
                clash = ns.search(atom.get_coord(), radius=1.5)
                if clash:
                    contador += 1
                n_of_atoms += 1
            if contador > n_of_atoms*0.35:
                return True
    return False

def delete_overlapping_chains(pdb):
    """
    Looks for clashes between chains and if the clash is big it deletes the latest chain added to the pdb
    :param pdb: A Protein Structure object
    """
    for chain in pdb:
        to_pop = list()
        for chain_position in range(len(pdb)):
            other_chain = pdb.childs[chain_position]
            contador, n_of_atoms = 0, 0
            if chain is other_chain:
                continue
            ns = NeighborSearch(other_chain.get_atoms_list())
            for atom in chain.iter_atoms():
                clash = ns.search(atom.get_coord(), radius=1.5)
                if len(clash) > 0:
                    contador += 1
                n_of_atoms += 1
            if contador > n_of_atoms*0.35:
                sys.stderr.write('Chain %s clashes significantly with chain %s so we will delete the latest\n' %(chain.get_id(), other_chain.get_id()))
                to_pop.append(chain_position)
        to_pop.sort(reverse=True)
        for position_to_pop in to_pop:
            pdb.childs.pop(position_to_pop)

if __name__== '__main__' :
    pdb_list = list()

    # ---Estos funcionan ---
    #folder = 'pdb/proteasoma'
    #folder = 'pdb/deconstruct'
    #folder = 'pdb/nucl'
    #folder = 'pdb/hemo_deconstruct'
    #folder = 'pdb/hemoglobin'
    #folder = 'pdb/phosphate'

    # ---Estos no funcionan ---
    folder = 'pdb/capsid'



    id_list = os.listdir(folder)
    for pdb_id in id_list:
        pdb_list.append(PS(pdb_id, '%s/%s' %(folder, pdb_id)))
    for pdb in pdb_list:
        pdb.childs[0].id = 'A'
        pdb.childs[1].id = 'A'
    if not os.path.exists('tmp'):
        os.mkdir('tmp')
    homo_chains = good_chain_names(pdb_list)
    new_pdb, chain_id_dict = reconstruct(pdb_list, homo_chains)
    delete_overlapping_chains(new_pdb)
    new_pdb.save_to_file('pdb/final.pdb')
    print('THE END')
