#!bin/env/ python3


class IncorrectSequenceLetter(Exception):
	def __init__(self, value, id, class_name):
		self.value = value
		self.raised_by = class_name
		self.id = id
	def __str__(self):
		return "The sequence item %s located in %s is not found in the alphabet of %s." %(self.value, self.id, self.raised_by)
class Sequence(object):
	'''
	A sequence abstract class, not ment to be used which constitudes the base for the other 
	sequence child classes line Protein and Nucleotide classes. It has the private attributes
	__identifier, and __sequence. As this attributes are private the get_identifier() and 
	get_sequence() are ment to be used to obtain the value of this attributes. Also the general
	method of get_length(), get_subsequence(another_sequence_object) and get_mw() are disponible.
	'''
	alphabet = set()
	def __init__(self, ident, seq):
		self.__identifier = ident
		self.__sequence = seq
		for letter in ''.join(set(seq)):
			if letter not in self.alphabet:
			 raise IncorrectSequenceLetter(letter, self.get_identifier(),self.__class__.__name__)
		self.mweight = None
	def get_identifier(self):
		return(self.__identifier)
	def get_sequence(self):
		return(self.__sequence)
	def get_mw(self):
		if self.mweight == None:
			self.mweight = 0
			for ky in self.mw_dict:
				time = self.__sequence.count(ky)
				self.mweight += time*self.mw_dict[ky]
		return(round(self.mweight,4))
	def get_length(self):
		return(len(self.__sequence))
	def has_subsequence(self, subseq_object):
		if subseq_object.get_sequence() in self.get_sequence():
			return True
		else:
			return False
	def __len__(self):
		return len(self.get_sequence())
	def __eq__(self, other_sequence):
		return self.get_sequence() == other_sequence.get_sequence()
	def __ne__(self, other_sequence):
		return self.get_sequence() != other_sequence.get_sequence()
	def __lt__(self, other_sequence):
		return self.get_mw() < other_sequence.get_mw()
	def __hash__(self):
		return hash(self.get_sequence() + self.get_identifier())
	def __add__(self, other_sequence):
		if self.__class__ == other_sequence.__class__:
			return self.__class__("%s+%s" %(self.get_identifier(), other_sequence.get_identifier()),
			 self.get_sequence()+other_sequence.get_sequence())
		else:
			raise ArithmeticError('Impossible to add sequences of different class')
	def __getitem__(self, key):
		return self.get_sequence()[key]
	def __contains__(self, insider):
		return str(insider) in self.get_sequence()
	def __str__(self):
		return self.get_sequence()
class NucleotideSequence(Sequence):
	'''
	A nucleotide class that has Sequence as parent and inherits the same atributes as the Sequence
	class. NucleotideSequence has also the translate method which is also generic and is ment to
	be used by the child classes like DNASequence and RNASequence and returns the aminoacid sequence
	based on the translation table.
	The get_complementary_sequence is also a generic method for the child clases that returns the
	complementary sequence of the sequence attribute.
	'''
	alphabet = ''
	def translate(self):
		start = list()
		prot = ''
		seq = self.get_sequence()
		i = 0
		while i < (len(seq)-3):
			if seq[i:i+3] in self.start_codons:
				break
			else:
				i +=1
		while i < (len(seq)-3):
			if seq[i:i+3] not in self.stop_codons:
				prot += self.table[seq[i:i+3]]
			else:
				break
			i += 3
		return prot
	def get_complementary_sequence(self):
		return "".join(self.complement.get(base, base) for base in self.get_sequence())
class ProteinSequence(Sequence):
	'''
	A protein class that has Sequence as parent and inherits the same atributes and methods as 
	Sequence.
	'''
	alphabet = 'ACEDGFIHKMLNQPSRTWVY'
	mw_dict = {'A': 89.09, 'C': 121.16, 'E': 147.13, 'D': 133.1, 'G': 75.07, 'F': 165.19, 
				   'I': 131.18, 'H': 155.16, 'K': 146.19, 'M': 149.21, 'L': 131.18, 'N': 132.12,
				   'Q': 146.15, 'P': 115.13, 'S': 105.09, 'R': 174.2, 'T': 119.12, 'W': 204.23,
				   'V': 117.15, 'Y': 181.19}
class DNASequence(NucleotideSequence):
	'''
	A DNA class that has NucleotideSequence and inherits the same and methods as its parent. Also, 
	it has the method transcribe, which returns the RNA version of the sequence atribute.
	'''
	alphabet = 'GATC'
	mw_dict = {'A': 347.0, 'C': 323.0, 'T': 322.0, 'G': 363.0}
	stop_codons = ['TAA', 'TAG', 'TGA']
	start_codons = ['TTG', 'CTG', 'ATG']
	table = {'CTT': 'L', 'ATG': 'M', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'ATA': 'I', 'AGG': 'R', 'CCT': 'P', 'AGC': 'S', 'AGA': 'R', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'ACT': 'T', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'GGG': 'G', 'GGA': 'G', 'GGC': 'G', 'TAC': 'Y', 'CGT': 'R', 'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GAC': 'D', 'GAA': 'E', 'AAG': 'K', 'AAA': 'K', 'AAC': 'N', 'CTC': 'L', 'CAT': 'H', 'AAT': 'N', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'TGT': 'C', 'TCT': 'S', 'GAT': 'D', 'TTT': 'F', 'TGC': 'C', 'TGG': 'W', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TCA': 'S', 'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A'}
	table_back = {'A': 'GCT', 'C': 'TGT', None: 'TAA', 'E': 'GAG', 'D': 'GAT', 'G': 'GGT', 'F': 'TTT', 'I': 'ATT', 'H': 'CAT', 'K': 'AAG', 'M': 'ATG', 'L': 'TTG', 'N': 'AAT', 'Q': 'CAG', 'P': 'CCT', 'S': 'TCT', 'R': 'CGT', 'T': 'ACT', 'W': 'TGG', 'V': 'GTT', 'Y': 'TAT'}
	complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
	def transcribe(self):
		return self.get_sequence().replace('T', 'U')
class RNASequence(NucleotideSequence):
	'''
	A RNA class that has NucleotideSequence and inherits the same and methods as its parent. Also, 
	it has the method transcribe, which returns the DNA version of the sequence atribute.
	'''
	alphabet = 'GAUC'
	mw_dict = {'A': 363.0, 'C': 339.0, 'U': 340.0, 'G': 379.0}
	stop_codons = ['UAA', 'UAG', 'UGA']
	start_codons = ['UUG', 'CUG', 'AUG']
	table = {'GUC': 'V', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GUU': 'V', 'AAC': 'N', 'AGG': 'R', 'UGG': 'W', 'AGC': 'S', 'AUC': 'I', 'AGA': 'R', 'AAU': 'N', 'ACU': 'T', 'CAC': 'H', 'GUG': 'V', 'CCG': 'P', 'CCA': 'P', 'AGU': 'S', 'CCC': 'P', 'GGU': 'G', 'UCU': 'S', 'GCG': 'A', 'CGA': 'R', 'CAG': 'Q', 'CGC': 'R', 'UAU': 'Y', 'CGG': 'R', 'UCG': 'S', 'CCU': 'P', 'GGG': 'G', 'GGA': 'G', 'GGC': 'G', 'GAG': 'E', 'UCC': 'S', 'UAC': 'Y', 'CGU': 'R', 'GAA': 'E', 'AUA': 'I', 'GCA': 'A', 'CUU': 'L', 'UCA': 'S', 'AUG': 'M', 'CUG': 'L', 'AUU': 'I', 'CAU': 'H', 'CUA': 'L', 'GCC': 'A', 'AAA': 'K', 'AAG': 'K', 'CAA': 'Q', 'UUU': 'F', 'GAC': 'D', 'GUA': 'V', 'UGC': 'C', 'GCU': 'A', 'UGU': 'C', 'CUC': 'L', 'UUG': 'L', 'UUA': 'L', 'GAU': 'D', 'UUC': 'F'}
	table_back = {'A': 'GCU', 'C': 'UGU', None: 'UAA', 'E': 'GAG', 'D': 'GAU', 'G': 'GGU', 'F': 'UUU', 'I': 'AUU', 'H': 'CAU', 'K': 'AAG', 'M': 'AUG', 'L': 'UUG', 'N': 'AAU', 'Q': 'CAG', 'P': 'CCU', 'S': 'UCU', 'R': 'CGU', 'T': 'ACU', 'W': 'UGG', 'V': 'GUU', 'Y': 'UAU'}
	complement = {'A': 'U', 'C': 'G', 'G': 'C', 'U': 'A'}
	def transcribe(self):
		return self.get_sequence().replace('U', 'T')
