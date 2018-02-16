#!/bin/env python3

from Sequences import ProteinSequence
from Bio.PDB import protein_letters_3to1
import numpy
import copy

class BASE(object):
    """Base class to build the ProteinStructure on top of it"""
    def __init__(self, id):
        self.id = id
        self.childs = []
        self.__child_dict = self._get_childs_dict(self.childs)
        self.parent = None
    def _get_childs_dict(self, list_of_childs):
        """Private method to generate a dict of pointers to the childs (just to iterate better if needed a dict)"""
        d_child = {}
        for child in list_of_childs:
            d_child[child.id] = child
        return d_child
    def __len__(self):
        """Return the number of childs"""
        return len(self.child)
    def __eq__(self, other):
        if self.__class__ == other.__class__:
            return self.get_id() == other.get_id()
        else:
            raise ValueError("%s and %s can't be compared" %(self.__class__.__name__, other.__class__.__name__))
    def __lt__(self, other):
        if self.__class__ == other.__class__:
            return self.get_id() < other.get_id()
        else:
            raise ValueError("%s and %s can't be compared" %(self.__class__.__name__, other.__class__.__name__))
    def __hash__(self):
        return hash(self.get_id())
    def __iter__(self):
        """Iterate over children."""
        for child in self.childs:
            yield child
    def __getitem__(self, id):
        """Return the child with given id."""
        return self.__child_dict[id]
    def __contains__(self, id):
        """True if there is a child element with the given id."""
        return (id in self.__child_dict)
    def parenting(self):
        """Sets recursively the parents of the childs as self"""
        for child in self.childs:
            child.set_parent(self)
    def set_parent(self, papa):
        """Set own parent and starts parenting if childs"""
        self.parent = papa
        self.parenting()
    def get_id(self):
        """Return the id."""
        return self.id
    def get_parent(self):
        """Return the id."""
        return self.parent
    def transform(self, rot, tran):
        """
        Apply rotation and translation to the atomic coordinates.
        It goes all until Atoms and it transforms them. (this method
        overrided in the ATOM class

        @param rot: A right multiplying rotation matrix
        @type rot: 3x3 Numeric array

        @param tran: the translation vector
        @type tran: size 3 Numeric array
        """
        for o in self:
            o.transform(rot, tran)
class ProteinStructure(BASE):
    """Protein Structure class with the typical hierarchical structure:
            Structure
                ·Chain
                    ·Residue
                        ·Atom
            It inherits the attributes from BASE with some changes: its childs are a list of chain objects"""
    def __init__(self, id ,struct_dict):
        BASE.__init__(self, id)
        self.childs = self._init_chains(struct_dict)
        self.mw = None
        self.parenting()
        self.get_mw()
    def _init_chains(self, struct_dict):
        """Private method to generate and return the child chain objects for initialization"""
        c = []
        for chain in struct_dict:
            c.append(Chain(chain, struct_dict[chain]))
        return c
    def get_mw(self):
        """Returns the molecular weight as the sum of the molecular weight of it's chains"""
        if self.mw is None:
            mweight = 0
            for child in self:
                mweight += child.get_mw()
            self.mw = round(mweight, 4)
        return self.mw
    def get_chains(self):
        """Returns a list of chains"""
        return self.childs
    def remove_chain(self, chain_id):
        """Removes a chain from the structure"""
        self.childs.remove(chain_id)
    def save_to_file(self, outfile, atom_name = None):
        """Saves the structure file in a pdb format to the outfile"""
        if type(atom_name) == str :
            atom_name = [atom_name]
        with open(outfile, "w") as out_pdb:
            for chain in self:
                for residue in chain:
                    for atom in residue:
                        if atom_name is None or atom.get_name() in atom_name:
                            out_pdb.write("%-6s%5s %4s %3s %s%4s    %8s%8s%8s\n" %('ATOM', atom.num, atom.name,
                                                                               residue.name, chain.id, residue.num,
                                                                               atom.coords[0], atom.coords[1], atom.coords[2]))
class Chain(BASE):
    """Chain class in the typical hierarchical structure:
               Structure
                   ·Chain
                       ·Residue
                           ·Atom
        It inherits the attributes from BASE with some changes: its childs are a list of residues objects
        additionally its sequence is a ProteinSequence object"""
    def __init__(self, id, chain_dict):
        BASE.__init__(self, id)
        self.childs = self._init_residues(chain_dict)
        self.sequence = self._obtain_sequenceobj(chain_dict)
        self.mw = self.sequence.get_mw()
    def _init_residues(self, chain_dict):
        """Private method to generate and return the child residue objects for initialization"""
        r = []
        for res in chain_dict:
            r.append(Residue(res, chain_dict[res]))
        return r
    def _obtain_sequenceobj(self, chain_dict):
        """Private method to generate the ProteinSequence object for initialization"""
        seq = ''
        for res_ky in chain_dict:
            seq += protein_letters_3to1[res_ky[1]]
        return ProteinSequence(self.id, seq)
    def get_mw(self):
        """Return molecular weight"""
        return self.mw
    def get_residues(self):
        """Returns a list of residues"""
        return self.childs
class Residue(BASE):
    """Residue class in the typical hierarchical structure:
                   Structure
                       ·Chain
                           ·Residue
                               ·Atom
            It inherits the attributes from BASE with some changes: its childs are a list of Atom object.
            additionally its sequence is a ProteinSequence object"""
    def __init__(self, id_tupple, res_list):
        BASE.__init__(self, id_tupple)
        self.num = id_tupple[0]
        self.name = id_tupple[1]
        self.childs = self._ini_atoms(res_list)
    def _ini_atoms(self, res_list):
        """Private method to generate and return the child atom objects for initialization"""
        a = []
        for aa in res_list:
            a.append(Atom(aa))
        return a
    def get_atoms(self):
        """Returns a list of atoms"""
        return self.childs
    def get_specific_atom(self, atom_name):
        """Returns the atom object with the specified name."""
        for atom in self:
            if atom.get_name() == atom_name:
                return atom

class Atom(BASE):
    """Atom class in the typical hierarchical structure:
                       Structure
                           ·Chain
                               ·Residue
                                   ·Atom
        It inherits the attributes from BASE with some changes:
            · It has no child (it's the bottom of the hierarchy
            · It overrides the transform method with an actuall changing its coordinates."""
    def __init__(self, info):
        BASE.__init__(self, info[1])
        self.num = info[0]
        self.name = info[1]
        self.coords = info[2]
    def __sub__(self, other):
        """
        Calculate distance between two atoms.

        Example:
            · distance = atom1-atom2
        """
        if self.__class__ == other.__class__:
            diff = self.coords - other.coords
            return numpy.sqrt(numpy.dot(diff, diff))
        else:
            raise ArithmeticError('Impossible to calculate the distance from an Atom to a %s' %other.__class__.__name__)
    def get_name(self):
        """Return atom name."""
        return self.name
    def get_coords(self):
        """Returns a list of coords"""
        return self.coords
    def transform(self, rot = numpy.array([[0,0,0],[0,0,0],[0,0,0]]), tran = numpy.array([0,0,0])):
        """Apply rotation and translation to the atomic coordinates.
        @param rot: A right multiplying rotation matrix
        @type rot: 3x3 Numeric array

        @param tran: the translation vector
        @type tran: size 3 Numeric array
        """
        self.coords = tuple(numpy.dot(self.coords, rot) + tran)

def load_pdb(file_name):
    pdb = dict()
    sequence = dict()
    last = None
    with open(file_name, "r") as file:
        for line in file:
            line = line.rstrip()
            if line.startswith('ATOM'):
                cols = line.split()
                if last is None:
                    last = [cols[5], cols[4]]
                    sequence.setdefault(cols[4], []).append(cols[3])
                elif cols[5] > last[0] or cols[4] != last[1]:
                    last = None
                pdb.setdefault(cols[4], dict()).setdefault((cols[5], cols[3]), list()).append(
                    ( cols[1], cols[2], (float(cols[6]), float(cols[7]), float(cols[8]))))
        return ProteinStructure( file_name , pdb)

if __name__ == '__main__':
    my_pdb= load_pdb('pdb/1a3n.pdb')
    my_otherpdb = copy.deepcopy(my_pdb)
    my_otherpdb.transform(rot = numpy.array([[0,-1,0],
                                            [1,0,0],
                                            [0,0,1]]),
                          tran=numpy.array([0,0,0]))
    my_otherpdb.save_to_file('pdb/1a3n_tranformado.pdb', ['CA', 'CB'])
    print(my_pdb > my_otherpdb)
