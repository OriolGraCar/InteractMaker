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

To install the package of pymol as an API for python there are many ways; the easiest one is having the conda package manager and installing it through the followin repo.

```bash
    conda install -c samoturk pymol
```




Based on // References

Hidrogen bond distance as 3.5 A

  Martz, Eric; Help, Index & Glossary for Protein Explorer, http://www.umass.edu/microbio/chime/pe_beta/pe/protexpl/igloss.htm
  Jeffrey, George A.; An introduction to hydrogen bonding, Oxford University Press, 1997.
