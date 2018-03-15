from PDB import ProteinStructure as PS
import copy
import numpy as np
import sys
def refine_interactions(interaction_dict):
    refined_interacciones = dict()
    for pair_type in interaction_dict:
        refined_interacciones.setdefault(pair_type, dict())
        for interaccion in interaction_dict[pair_type]:
            valid = True
            for otra_interaccion in refined_interacciones[pair_type]:
                intersect =  set(interaccion).intersection( set(otra_interaccion))
                if len(intersect) >= max([len(interaccion), len(otra_interaccion)])*0.2:
                    if len(interaccion) >= len(otra_interaccion):
                        refined_interacciones.setdefault(pair_type, dict())[interaccion] = interaction_dict[pair_type][interaccion]
                        del refined_interacciones[pair_type][otra_interaccion]
                        break
                    else:
                        valid = False
                        continue
            if interaccion not in refined_interacciones[pair_type] and valid:
                refined_interacciones.setdefault(pair_type, dict())[interaccion] = interaction_dict[pair_type][interaccion]
    return refined_interacciones
def independize_from_mainpdb(interaction_dict):
    independent_dict = dict()
    for pair_type in interaction_dict:
        for interaccion in interaction_dict[pair_type]:
            independent_dict.setdefault(pair_type, dict())[interaccion] = [copy.deepcopy(interaction_dict[pair_type][interaccion][0]), copy.deepcopy(interaction_dict[pair_type][interaccion][1])]
    return  independent_dict
def deconstruct_macrocomplex_by_interactions(pdb_file, output_folder = './', translation = np.array([0,0,0]), rotation= np.array([[1,0,0], [0,1,0], [0,0,1]])):
    '''
    :param pdb_file(str): A path to a pdb file
    :param output_folder (str): A path to a directory where the interactions will be stored
    :param translation (numpy array 3x1): a np array defining the translation of the new pdb's
    :param rotate (numpy array 3x3): a np array defining the rotation of the new pdb's
    '''
    pdb = PS('original', pdb_file)
    pdb_list = list()
    todas_interacciones = dict()
    seqs = list()
    for chain in pdb:
        if chain.get_sequence_str() not in seqs:
            seqs.append(chain.get_sequence_str())
        for other_chain in pdb:
            if other_chain.get_sequence_str() not in seqs:
                seqs.append(other_chain.get_sequence_str())
            interacting_res_tuple = None
            if chain is not other_chain:
                int = chain.interacting_residues(other_chain)
                if int:
                    interacting_res_tuple = tuple(int)
                if interacting_res_tuple:
                    todas_interacciones.setdefault(''.join(sorted(chain.get_sequence_str()+other_chain.get_sequence_str())), dict())[interacting_res_tuple] = chain.get_id()+other_chain.get_id()
    interacciones_unicas = refine_interactions(todas_interacciones)
    for i in range(100):
        interacciones_unicas = refine_interactions(interacciones_unicas)

    for pair_type in interacciones_unicas:
        for interaccion in interacciones_unicas[pair_type]:
            pdb.transform(tran=translation, rot=rotation)
            pdb.save_to_file('%s/%s.pdb' % (output_folder,interacciones_unicas[pair_type][interaccion]), chain_name=interacciones_unicas[pair_type][interaccion])

'''
Rotation matrix for a angle = A is the following:

        1       0       0
Rx(A) = 0       cosA    -sinA
        0       sinA    cosA

        cosA    0       sinA
Ry(A) = 0       1       0
        -sinA   0       cosA

        cosA    -sinA   0
Rz(A) = sinA    cosA    0
        0       0       1
                                                        '''

if __name__ == '__main__':

    #deconstruct_macrocomplex_by_interactions(sys.argv[1], sys.argv[2])
    deconstruct_macrocomplex_by_interactions('pdb/1pma.pdb', 'pdb/deconstruct')

