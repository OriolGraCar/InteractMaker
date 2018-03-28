# BioMacromplex


##How to install

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
