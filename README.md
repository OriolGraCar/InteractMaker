# BioMacromplex

BioMacromplex is a python3 package to build Protein, DNA or RNA macrocomplex based on a set of interactions. It also includes a module to use and manipulate pdb (PDB) files (alternative to the Biopython one, although it still works with the other Biopython classes or functions)

Aditionally, it includes a stand-alone script *Reconstruct_macrocomplex.py* and a Graphical intereface version *InteractMakerGui.py* to reconstruct macrocomplex.
## How to install

Download the git folder on the directory you want

```bash
    git clone https://github.com/OriolGraCar/InteractMaker.git
```

Then check the BioMacromplex directory is here and also a script called setup.py. 

To be able to use the script *Reconstruct_macrocomplex.py* anywhere you'll need to install the package in your python site-packages.

```bash
    python3 setup.py install
```

Notice that the package has the following dependencies:
	
	- biopython
	- numpy
	- TKinter
	- pymol
	- argparse

To install the package of pymol as an API for python there are many ways; the easiest one is having the conda package manager and installing it through the followin repo.

```bash
    conda install -c samoturk pymol
```

## Usage

#### Stand-alone
```bash
    $ ./Reconstruct_macrocomplex.py -h


    usage: Reconstruct_macrocomplex.py [-h] -i FOLDER [-o OUTPUT_PDB] [-v] [-s]

./Reconstruct_macrocomplex.py is a script that uses the BioMacromplex module
to reconstruct Protein, DNA or RNA complexes through a set of interactions in
pdb format.

optional arguments:
  -h, --help     show this help message and exit
  -i FOLDER      Path to a folder where pdb are stored
  -o OUTPUT_PDB  Output file where the macrocomplex pdb will be stored
  -v, --verbose  Print the progress of the program and the log
  -s, --steps    Save a temporary pdb in tmp/ each time a chain is added to
                 track the process

```


## References

Hidrogen bond distance as 3.5 A

  Martz, Eric; Help, Index & Glossary for Protein Explorer, http://www.umass.edu/microbio/chime/pe_beta/pe/protexpl/igloss.htm
  Jeffrey, George A.; An introduction to hydrogen bonding, Oxford University Press, 1997.
