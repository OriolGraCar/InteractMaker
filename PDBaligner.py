from PDB4 import ProteinStructure as PS
from Bio.PDB import Superimposer
import copy
import numpy as np
from Bio.pairwise2 import align, format_alignment

def find_common_atoms(chain_fix, chain_mov, bb = True):
    seq_fix = chain_fix.get_sequence_str()
    seq_mov = chain_mov.get_sequence_str()
    atoms_fix = list()
    atoms_mov = list()
    print('Alignment from the chain %s from the %s pdb and chain %s from %s pdb' %(
        chain_fix.get_id(), chain_fix.parent.id,chain_mov.get_id(), chain_mov.parent.id))
    alignment = align.globalxx(seq_fix, seq_mov)
    count_fix = 0
    count_move = 0
    for i in range(len(alignment[0][0])):
        if alignment[0][0][i] == '-':
            count_fix -= 1
        if alignment[0][1][i] == '-':
            count_move -= 1
        elif alignment[0][0][i] != '-' and alignment[0][1][i] != '-':
            atoms_fix.append(chain_fix.childs[i + count_fix]['CA'])
            atoms_mov.append(chain_mov.childs[i + count_move]['CA'])
    #print(format_alignment(*alignment[0]))
    return atoms_fix, atoms_mov


def superimpose_pdb_by_chain(chain_fix, chain_mov):
    sup = Superimposer()
    atoms_fix, atoms_mov = find_common_atoms(chain_fix, chain_mov)
    sup.set_atoms(atoms_fix, atoms_mov)
    chain_mov.parent.transform(rot = sup.rotran[0], tran = sup.rotran[1])



if __name__== '__main__' :
    pdb_list = list()
    id_list = ['ABXX', 'AXCX', 'AXXD']
    sqs = dict()
    sup = Superimposer()

    for id in id_list:
        pdb_list.append(PS(id, 'pdb/%s.pdb' %id))
    for i in range(100):
        for pdb1 in pdb_list:
            for pdb2 in pdb_list:
                if pdb1 is pdb2:
                    break
                for chain1 in pdb1:
                    for chain2 in pdb2:
                        if chain1.get_id() == chain2.get_id():
                            superimpose_pdb_by_chain(chain1, chain2)
    for pdb in pdb_list:
        pdb.save_to_file('pdb/%s_fit.pdb'%pdb.get_id())




'''
#rot = np.identity(3).astype('f')
rot = np.array([[  1.00000000e+00,  -5.22362039e-10,  0.000009],
                [  0.00000009,   1.00000000e+00,  -3.29020445e-7],
                [  7.64685491e-7,   5.90978817e-7,   1.00000000e+00]])
tran = np.array((30.0, 20.0, 30.0), 'f')

pdb1 = PS('1a3n_fix', 'pdb/1a3n.pdb')
pdb2 = PS('1a3n_mov', 'pdb/1a3n.pdb')
pdb2.transform(rot=rot, tran=tran)
pdb2.save_to_file('pdb/1a3n_movido.pdb')

a_fix = pdb1.get_atoms()
a_mov = pdb2.get_atoms()

sup = Superimposer()

sup.set_atoms(a_fix, a_mov)

print(sup.rotran)
print(sup.rms)

sup.apply(a_mov)
pdb2.save_to_file('pdb/1a3n_fit.pdb')
'''
print('THE END')
