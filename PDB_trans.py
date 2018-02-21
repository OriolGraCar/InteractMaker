from PDB4 import ProteinStructure as PS
import numpy as np
import math as m

pdb = PS('original', 'pdb/6bb5.pdb')
id_list = ['AC', 'AD', 'BC', 'DC', 'AB']
i = 3
tran = np.array((30.0, i*5.0, 30.0), 'f')
X = 40
RotX = np.array([[1, 0, 0],[0, m.cos(m.radians(X)), -m.sin(m.radians(X))], [0, m.sin(m.radians(X)), m.cos(m.radians(X))]])
for i in range(len(id_list)):
    pdb.transform(tran=tran, rot=RotX)
    pdb.save_to_file('pdb/%s.pdb'%id_list[i], chain_name=id_list[i])


'''Rotation matrix for a angle = A is the following:

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