# BioMacromplex

BioMacromplex is a python3 package to build Protein, DNA or RNA macrocomplex based on a set of interactions. It also includes a module to use and manipulate pdb (PDB) files (alternative to the Biopython one, although it still works with the other Biopython classes or functions)

Aditionally, it includes a stand-alone script *Reconstruct_macrocomplex.py* and a Graphical intereface version *InteractMakerGui.py* to reconstruct macrocomplex.
## How to install

Download the git folder on the directory you want

```bash
   $ git clone https://github.com/OriolGraCar/InteractMaker.git
```

Then check the BioMacromplex directory is here and also a script called setup.py. 

To be able to use the script *Reconstruct_macrocomplex.py* anywhere you'll need to install the package in your python site-packages.

```bash
   $ python3 setup.py install
```

Notice that the package has the following dependencies:
	
	- biopython
	- numpy
	- TKinter
	- pymol

To install the package of pymol as an API for python there are many ways; the easiest one is having the conda package manager and installing it through the followin repo.

```bash
   $ conda install -c samoturk pymol
```

## Usage

#### Stand-alone

```bash
    $ ./Reconstruct_macrocomplex.py -h

usage: Reconstruct_macrocomplex.py [-h] -i FOLDER [-o OUTPUT_PDB] [-v] [-s]
                                   [-q]

Reconstruct_macrocomplex.py is a script that uses the BioMacromplex module to
reconstruct Protein, DNA or RNA complexes through a set of interactions in pdb
format.

optional arguments:
  -h, --help     show this help message and exit
  -i FOLDER      Path to a folder where pdb are stored
  -o OUTPUT_PDB  Output file where the macrocomplex pdb will be stored
  -v, --verbose  Print the progress of the program and the log
  -s, --steps    Save a temporary pdb in tmp/ each time a chain is added to
                 track the process
  -q, --quiet    Stores standard output and error in /dev/null/

```

#### GUI

## Examples

Here we'll present a some examples macrocomplex that we have used as test subjects for the realization of the package and scripts (they are not the only examples we have used, tho).

#### Proteasoma (1pma.pdb) in command line

The first step is to obtain the pdb with the interactions between the pair of chains. You can do it either manually chosing the chains that will go in each pdb file or use our tools for it (either the GUI with PDB_split or the PDB_split.py as a script instead of a module).

In this case we will use the latest:
```bash
    $ python3 BioMacromplex/PDB_split.py test/1pma.pdb test/proteasoma
    
    $ ls test/proteasoma

1U.pdb  2O.pdb  2S.pdb  FE.pdb  SE.pdb  UB.pdb
```

As you can see we have created six pdb in the *test/proteasoma* directory, each one defining an interaction (in fact not all six interactions are needed for the program to work, in this case with 1U.pdb, 2O.pdb & FE.pdb would suffice ). 

Then we can proceed to *Reconstruct the macrocomplex*:
```bash
    $ python3 Reconstruct_macrocomplex.py -i test/proteasoma -o test/reconstructed_proteasoma.pdb

WARNING!: Chain O in pdb 2O doesn't start in residue 1
WARNING!: Chain E in pdb SE doesn't start in residue 1
WARNING!: Chain E in pdb FE doesn't start in residue 1
WARNING!: Chain F in pdb FE doesn't start in residue 1
The END

```

With this the process will be finished and the macrocomplex in our specified directory.

#### El que vulguis, oriol with the GUI



#### Other examples in test/

In the test folder you can find also some other test macrocomplexes that can be used as toy models with the interactions already extracted:

	- A Phosphate dehydratase (2f1d.pdb) in test/phosphate/
	- A Nucleosoma (3kuy.pdb) in test/nucl
	- An ATP-syntasa (5kw1.pdb) in test/ATP
	- An Hemoglobin (1a3n.pdb) in test/hemoglobin
	- A Mosaic Virus Capside (3j7l.pdb) in test/mosaic_virus

## References

Hidrogen bond distance taken as 3.5A so it will be the distance managed by the scripts in the pagackes to find and check interactions between different molecules.

Martz, Eric; Help, Index & Glossary for Protein Explorer, http://www.umass.edu/microbio/chime/pe_beta/pe/protexpl/igloss.htm
Jeffrey, George A.; An introduction to hydrogen bonding, Oxford University Press, 1997.
