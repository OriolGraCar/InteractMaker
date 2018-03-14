#!/bin/env python3

from Sequences import *
from Bio.PDB import protein_letters_3to1
import copy
import numpy
from sys import stderr
from Bio.PDB import NeighborSearch

ascii_lowercase = 'abcdefghijklmnopqrstuvwxyz'
ascii_uppercase = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
numbers = '01234567890'
ltr = ascii_uppercase[::-1] + ascii_lowercase + numbers



class BASE(object):
    """Base class to build the ProteinStructure on top of it"""
    def __init__(self, id):
        self.id = id
        self.child_dict = self._get_childs_dict(self.childs)
        self.parent = None
    def _get_childs_dict(self, list_of_childs):
        """Private method to generate a dict of pointers to the childs (just to iterate better if needed a dict)"""
        d_child = {}
        for child in list_of_childs:
            d_child[child.id] = child
        return d_child
    def __len__(self):
        """Return the number of childs"""
        return len(self.childs)
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
        return self.child_dict[id]
    def __contains__(self, id):
        """True if there is a child element with the given id"""
        return (id in self.child_dict)
    def parenting(self):
        """Sets recursively the parents of the childs as self"""
        for child in self.childs:
            child.set_parent(self)
    def set_parent(self, papa):
        """Set own parent and starts parenting if childs"""
        self.parent = papa
        self.parenting()
    def get_id(self):
        """Return the id"""
        return self.id
    def get_parent(self):
        """Return the parent object (one step higher in the hierarchy"""
        return self.parent
    def transform(self, rot = numpy.array([[1,0,0],[0,1,0],[0,0,1]]), tran = numpy.array([0,0,0])):
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
    def __init__(self, id ,pdb_file):
        self.childs = self._init_chains(pdb_file)
        BASE.__init__(self, id)
        self.mw = None
        self.parenting()
        self.get_mw()
        self.find_gaps()
    def _init_chains(self, pdb_file):
        """Private method to generate and return the child chain objects for initialization"""
        c = []
        pdb = self._read_pdb(pdb_file)
        for chain in pdb:
            c.append(Chain(chain, pdb[chain]))
        return c
    def _read_pdb(self, pdb_file):
        cols = list()
        pdb = dict()
        with open(pdb_file, "r") as file:
            for line in file:
                line = line.rstrip()
                if line.startswith('ATOM'):
                    if len(line) > 64:
                        cols = [line[:6].strip(), line[6:11].strip(), line[12:16].strip(), line[17:20].strip(),
                                line[21].strip(), line[22:26].strip(), line[30:38].strip(), line[38:46].strip(),
                                line[46:54].strip(), line[54:60].strip(), line[60:66].strip()]
                    elif len(line) > 53:
                        cols = [line[:6].strip(), line[6:11].strip(), line[12:16].strip(), line[17:20].strip(),
                                line[21].strip(), line[22:26].strip(), line[30:38].strip(), line[38:46].strip(),
                                line[46:54].strip()]
                    if len(cols) > 10: # include Temperature Factor
                        pdb.setdefault(cols[4], dict()).setdefault((cols[5], cols[3]), list()).append(
                        (cols[1], cols[2], (float(cols[6]), float(cols[7]), float(cols[8])), float(cols[9]), float(cols[10])))
                    elif len(cols) > 9: # include atom occupancy
                        pdb.setdefault(cols[4], dict()).setdefault((cols[5], cols[3]), list()).append(
                            (cols[1], cols[2], (float(cols[6]), float(cols[7]), float(cols[8])), float(cols[9])))
                    else: # standart pdb just with coordinates
                        pdb.setdefault(cols[4], dict()).setdefault((cols[5], cols[3]), list()).append(
                            (cols[1], cols[2], (float(cols[6]), float(cols[7]), float(cols[8]))))
        return pdb
    def get_mw(self):
        """Returns the molecular weight as the sum of the molecular weight of it's chains"""
        mweight = 0
        for child in self:
            mweight += child.get_mw()
        self.mw = round(mweight, 4)
        return self.mw
    def get_chains(self):
        """Returns a list of chains"""
        return self.childs
    def get_residues(self):
        """Returns a list with all structure residues"""
        r = []
        for chain in self.get_chains():
            for res in chain:
                r.append(res)
        return r
    def add_chain(self, nw_chain, cid, track_name = False):
        """
            Adds the input chain object to the Structure.

         It does a deep copy not to mess with the original chain object.
         If the given chain id (cid) is already in use by another chain in
          the structure it will look for an alternative name.

        :param nw_chain: Chain object to add to the structure object
        :param cid: chain id you want to give
        :param track_name(Boolean): default False;
        :return: if track_name = True returns the id given to the new chain.
        """
        nw_id = None
        if cid in self.child_dict:
            for letter in ltr:
                if letter not in self.child_dict:
                    nw_id = letter
                    break
            if nw_id is None:
                stderr.write("Internal Error: Limit of valid sequence names reached. Program finished abrubtly.\n")
                return None
            else:
                stderr.write('WARNING!: Tried to add a chain to %s with an already existing Chain id (%s). I will try to change it to %s.\n' %(cid , self.id, nw_id))
                cid = self.add_chain(nw_chain, nw_id, track_name= track_name)
        else:
            my_nw_chain = copy.deepcopy(nw_chain)
            my_nw_chain.id = cid
            self.childs.append(my_nw_chain)
            self.child_dict = self._get_childs_dict(self.childs)
            self.parenting()
        if track_name:
            return cid
    def get_atoms(self):
        """Returns a list with all structure atoms"""
        a = []
        for r in self.get_residues():
            for aa in r:
                a.append(aa)
        return a
    def restablish_dict(self):
        '''To use from outside.. does the same as self.child_dict = _get_childs_dict'''
        self.child_dict = self._get_childs_dict(self.childs)
    def remove_chain(self, chain_id):
        """Removes a chain from the structure."""
        self.childs.remove(chain_id)
        self.child_dict = self._get_childs_dict(self.childs)
    def find_gaps(self):
        """Check if Structure has gaps in its chains. Only looks if the residu num is consecutively."""
        for chain in self:
            current_residue = chain.childs[0].num #current_residue is initialized as the first residue - 1, so the first one is
            if current_residue != str(1):
                stderr.write("WARNING!: Chain %s in pdb %s doesn't start in residue 1\n" %(chain.get_id(), self.get_id()))
            for residue in chain:
                if residue.num != current_residue:
                    stderr.write("WARNING!: Gap found between residue %s and residue %s in the chain %s of the pdb %s\n" %(
                        current_residue, residue.num, chain.get_id(), self.get_id()))
                current_residue = str(int(residue.num) + 1)
    def find_clashes(self):
        """
        Check every pair of backbone atoms of diferent chains to see if they clash.

        A clash is defined when the distance between two atoms is greater than the sum of their VanderWaal radii
        """
        vwr = {'C':1.8, 'O':1.4, 'N':1.7, 'S':2, 'CA': 1.8}
        bb = ('C', 'O', 'N', 'CA')
        for chain in self:
            for other_chain in self:
                if chain is not other_chain:
                    contador_de_clashes = 0
                    ns = NeighborSearch(other_chain.get_atoms_list())
                    for atom in chain.iter_atoms():
                        clashes = ns.search(atom.get_coord, radius=1.5)
                        if len(clashes)> 0:
                            contador_de_clashes += 1
                    if contador_de_clashes > 0:
                        stderr.write("The chain %s has %s atoms clashing with chain %s in the pdb %s" %(chain.get_id(), contador_de_clashes, other_chain.get_id(), self.parent.get_id()))
    def save_fasta(self, outfile):
        """Saves the sequence of all the chains in the pdb """
        with open(outfile, "w") as out_fa:
            n = 80
            for chain in self:
                out_fa.write(">%s:%s\n" %(self.get_id(), chain.get_id()))
                seq = chain.sequence
                seq_list = [seq[i:i+n] for i in range(0, len(seq), n)]
                for line in seq_list:
                    out_fa.write("%s\n"%line)
    def save_to_file(self, outfile, atom_name = None, chain_name = None):
        """
        Saves the structure file in a pdb format to the outfile

        IF @atom_name given : it will select only the especified atoms.
        IF @chain_name given: it will select only the especified chains.
        """
        line_num = 1
        if type(atom_name) == str :
            atom_name = [atom_name]
        with open(outfile, "w") as out_pdb:
            for chain in self:
                if not chain_name is None and chain.get_id() not in chain_name:
                    continue
                for residue in chain:
                    for atom in residue:
                        if atom_name is None or atom.get_name() in atom_name:
                            if atom.occupancy is not None and atom.temp_factor is not None:
                                out_pdb.write("%-6s%5s %4s %3s %s%4s    %8.3f%8.3f%8.3f%6.2f%6.2f\n" %('ATOM', line_num, atom.name,
                                                            residue.name, chain.id, residue.num,round(atom.coords[0],3),
                                                            atom.coords[1], atom.coords[2],atom.occupancy, atom.temp_factor))
                            else:
                                out_pdb.write("%-6s%5s %4s %3s %s%4s    %8.3f%8.3f%8.3f\n" % ('ATOM', line_num, atom.name,
                                                            residue.name, chain.id, residue.num, round(atom.coords[0], 3),
                                                            atom.coords[1], atom.coords[2]))
                            #out_pdb.write("{:6s}{:5s} {:^4s}{:1s}{:3s} {:1s}{:4s}{:1s}   {:8s}{:8s}{:8s}{:6s}{:6s}\n".format('ATOM', atom.num, atom.name, residue.name, chain.id, residue.num,round(atom.coords[0],3), round(atom.coords[1],3), round(atom.coords[2],3), round(atom.occupancy, 2), round(atom.temp_factor, 2)))
                        line_num +=1
                out_pdb.write("TER\n")
class Chain(BASE):
    """Chain class in the typical hierarchical structure:
               Structure
                   ·Chain
                       ·Residue
                           ·Atom
        It inherits the attributes from BASE with some changes: its childs are a list of residues objects
        additionally its sequence is a ProteinSequence object"""
    def __init__(self, id, chain_dict):
        self.childs = self._init_residues(chain_dict)
        BASE.__init__(self, id)
        self.sequence = self._obtain_sequenceobj(chain_dict)
        self.mw = self.sequence.get_mw()
    def __str__(self):
        return self.sequence.get_sequence()
    def _init_residues(self, chain_dict):
        """Private method to generate and return the child residue objects for initialization"""
        r = []
        for res in chain_dict:
            r.append(Residue(res, chain_dict[res]))
        return r
    def _obtain_sequenceobj(self, chain_dict):
        """Private method to generate the ProteinSequence object for initialization"""
        seq = ''
        dna = {'DA' : 'A', 'DC': 'C', 'DG': 'G', 'DT': 'T'}
        rna = {'A': 'A', 'C': 'C', 'G': 'G', 'U': 'U'}
        flag = 'prot'
        for res_ky in chain_dict:
            try:
                seq += protein_letters_3to1[res_ky[1]]
            except:
                if res_ky[1] in dna:
                    seq += dna[res_ky[1]]
                    flag = 'dna'
                elif res_ky[1] in rna:
                    seq += rna[res_ky[1]]
                    flag = 'rna'
        if flag == 'prot':
            return ProteinSequence(self.id, seq)
        elif flag == 'dna':
            return DNASequence(self.id, seq)
        elif flag == 'rna':
            return RNASequence(self.id, seq)
    def get_mw(self):
        """Return molecular weight"""
        return self.mw
    def get_residues(self):
        """Returns a list of residues"""
        return self.childs
    def get_residue_by_num(self, number):
        for res in self:
            if res.num == number:
                return res
    def get_atoms_list(self):
        """Returns a list of atoms"""
        r = list()
        for residue in self.get_residues():
            r.extend(residue.childs)
        return r
    def iter_atoms(self):
        """Returns a generator of residues"""
        for residue in self.get_residues():
            for atom in residue:
                yield atom
    def get_sequence_str(self):
        """Returns string with the chain sequence."""
        return self.sequence.get_sequence()
    def get_sequence(self):
        """Returns sequence objct."""
        return self.sequence
    def renumber(self, ini):
        """
        Renumbers the residues of the chain starting on the @ini and the others consecutively
        WARNING!: It will mask the gaps.
        """
        i = ini
        for residue in self:
            residue.num = i
            residue.id = (str(i), residue.name)
            i += 1
        self.child_dict = self._get_childs_dict(self.childs)
    def interacting_residues(self, other_chain, dist = 3.5):
        """
        Iterates through all the possible pair of atoms (one from self and the other from other_chain)
        and returns a list of the residue numbers from self that interact with other_chain.

        :param other_chain: a chain object to be compared with
        :param dist: distance of interaction in Armstrong (default = 4)
        :return: List of intreacting residues
        """
        ns = NeighborSearch(other_chain.get_atoms_list())
        interacting = list()
        for res in self:
            for atom in res:
                anything_close = ns.search(center= atom.get_coord(), radius=dist)
                if len(anything_close) > 0:
                    interacting.append(res.num)
                    break
        if len(interacting)>0:
            return interacting
class Residue(BASE):
    """Residue class in the typical hierarchical structure:
                   Structure
                       ·Chain
                           ·Residue
                               ·Atom
            It inherits the attributes from BASE with some changes: its childs are a list of Atom object.
            additionally its sequence is a ProteinSequence object"""
    def __init__(self, id_tupple, res_list):
        self.childs = self._ini_atoms(res_list)
        BASE.__init__(self, id_tupple)
        self.num = id_tupple[0]
        self.name = id_tupple[1]
    def _ini_atoms(self, res_list):
        """Private method to generate and return the child atom objects for initialization"""
        a = []
        for aa in res_list:
            a.append(Atom(aa))
        return a
    def backbone(self):
        """Returns a list with the backbone atoms of the residue (Carbon, Nitrogen, Oxigen and Alfa Carbon)"""
        return [self['C'], self['N'], self['O'], self['CA']]
    def get_atoms(self):
        """Returns a list of atoms"""
        return self.childs
    def get_specific_atom(self, atom_name):
        """Returns the atom object with the specified name."""
        return self[atom_name]
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
        self.childs = []
        BASE.__init__(self, info[1])
        self.num = info[0]
        self.name = info[1]
        self.coords = info[2]
        self.occupancy = None
        self.temp_factor = None
        if len(info) > 4:
            self.temp_factor = info[4]
        if len(info) > 3:
            self.occupancy = info[3]
    def __sub__(self, other):
        """
        Calculate distance between two atoms.

        Example:
            · distance = atom1-atom2
        """
        if self.__class__ == other.__class__:
            diff = numpy.array(self.coords) - numpy.array(other.coords)
            return numpy.sqrt(numpy.dot(diff, diff))
        else:
            raise ArithmeticError('Impossible to calculate the distance from an Atom to a %s' %other.__class__.__name__)
    def get_name(self):
        """Return atom name."""
        return self.name
    def get_coords(self):
        """Returns a tupple of coords"""
        return tuple(self.coords)
    def get_coord(self):
        """Returns a numpy array of coords"""
        return numpy.array(self.coords)
    def get_occupancy(self):
        """Returns the occupancy of the atom"""
        return self.occupancy
    def get_temperature_factor(self):
        """Returns the temperature factor of the atom"""
        return self.temp_factor
    def transform(self, rot = numpy.array([[1,0,0],[0,1,0],[0,0,1]]), tran = numpy.array([0,0,0])):
        """Apply rotation and translation to the atomic coordinates.
        @param rot: A right multiplying rotation matrix
        @type rot: 3x3 Numeric array

        @param tran: the translation vector
        @type tran: size 3 Numeric array
        """
        self.coords = numpy.dot(self.coords, rot) + tran

if __name__ == '__main__':
    '''El sitio de las pruebas <3'''
