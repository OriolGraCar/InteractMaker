# Interactmaker

InteractMaker is a python3 program developed by Oriol Gracia Carmona and Alvaro Serrano Moras for the subject structural bioinformatics (SBI) of the Msc Bioinformatics for health sciences at the UPF.
This program is able to build different kinds of Macrocomplexes,(protein, DNA and RNA macrocomplexes) using its interactions pairs. 
It comes with a stand-alone command line script *Reconstruct_macrocomplex.py*, a Graphical User Interface (GUI) *InteractMakerGui.py* and a package called BioMacromplex, that contains all the functions and objects required to manipulate, load and build a pdb.

The theoretical approach used, and an in depth explanation of the usage of both the program and the package can be found in the *Documentation.md* that we also provide.

# Dependencies
The software versions stated here are the ones for which our program has been tested. However it may also work with older versions.

The following dependencies are the basic ones for the program to work.

	- python 3.5 or newer
	- biopython
	- numpy

For the GUI the following ones are also necessary:

	- TKinter 
	- pymol

# Installation and Download

Download the git folder on the directory you want

```bash
   $ git clone https://github.com/OriolGraCar/InteractMaker.git
```

Then check the BioMacromplex directory is there and also a script called setup.py. 

To be able to use the scripts provided anywhere you'll need to install the package BioMacromplex in your python site-packages.

```bash
   $ python3 setup.py install
```
Be sure to have the dependencies previously stated.

# Pymol installation

Since the installation of pymol may be tricky, here we provide an easy way to install it using the the conda package manager.
You just need to install it trough the following repository:

```bash
   $ conda install -c samoturk pymol
```

Despite that this installation is enough to run our program is not yet completed since it lacks the pmw package, necesary to run pymol comand window. to install it trough conda, use:

```bash
$ conda install -c fable pmw
```
The pmw package is also avaliable in pip


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
### GUI

To use the GUI just launch the script.

```bash

$ ./InteractMakerGui.py

```
For a detailed explanation of how to use the GUI check the *Documentation.md*

The program always returns one single structure that is the one most likely to be correct given the input data. A more detailed of the reason behing can be foun also in the *Documentation.md*

## Examples

The program also comes with some tested examples that we know they work properly.
All the examples have particularities that makes them different, this is so because we wanted to test and show how the program behaves in different situations.

The examples provided are:

	- A Phosphate dehydratase (2f1d.pdb) in test/phosphate/
	- A Nucleosoma (3kuy.pdb) in test/nucl/
	- An ATP-syntasa (5kw1.pdb) in test/ATP/
	- An Hemoglobin (1a3n.pdb) in test/hemoglobin/
	- A Mosaic Virus Capside (3j7l.pdb) in test/mosaic_virus/
	- An Aquaporin (2rc2.pdb) in test/aquaporin/
	- A RNA Polymerase I (5lmx.pdb) in test/RNA_polI/
	- A RNA exosome (2nn6.pdb) in test/RNA_exosome/
	- A Glutamine synthetase (1fpy.pdb) in test/glutamine_synth/

Particular characteritics of some of those examples:

- Hemoglobin: This one is a simple example of an small symetric protein but with the particulary of having hetatoms (not convencional atoms) in this case the iron ion.
- Proteosome: Example of a huge and higly symetric proteic macro-complex
- Nucleosome : Example of a protein - DNA macro-complex
- ATP-syntasa : Example of a big proteic complex with both symmetric and asymmetric parts.
- Phosphate dehydratase: Example of a big symetric macro-complex with lots of hetatoms.
- Virus capside 1: This is a enormous proteic macrocomplex made of multiple chains, that XXXXXX. Warning running this example may take a while do to its complexity.


# Command line example (Tutorial)

*Step 1:*
First we need the pdb files with the interactions. The program needs enough interactions to know how to place everything, for all the cases, the correct amount is all the different interactions found in the macro-complex.
However, it may still work with fewer interactions pairs as long as they are enough to know how the different proteins are placed together. If more interactions than the strictively necessary are provided the program will still work properly.
This files can be obtained in different ways, we also provide inside the package an script to split the pdbs called *pdb_split.py*

```bash
    $ python3 BioMacromplex/PDB_split.py path_of_the_input_file path_of_the_output_folder
 ```

This script will split the pdb in all the different interactions pairs found.

*Step 2:*
Once you have all the pairs of interactions just launch the Reconstruct_macrocomplex.py.

 ```bash
    $ python3 Reconstruct_macrocomplex.py -i path_of_the_folder_with_the_interactions_files -o outputfile_path

 ```
 You can add aditional arguments such as -v to have the vebrose or -s to have a pdb file of each one of the steps.
 Warning, if the pdb that is being constructed as more chains than possible letters, some chains will have the same name despite beign different chains. This happens for example when bulding a virus capside.

# Particular Example (*the proteosome*)

```bash
    $ python3 BioMacromplex/PDB_split.py test/1pma.pdb test/proteasoma
    
    $ ls test/proteasoma

1U.pdb  2O.pdb  2S.pdb  FE.pdb  SE.pdb  UB.pdb
```

```bash
    $ python3 Reconstruct_macrocomplex.py -i test/proteasoma -o test/reconstructed_proteasoma.pdb

WARNING!: Chain O in pdb 2O doesn't start in residue 1
WARNING!: Chain E in pdb SE doesn't start in residue 1
WARNING!: Chain E in pdb FE doesn't start in residue 1
WARNING!: Chain F in pdb FE doesn't start in residue 1
The END

```

# GUI example

*Step 1:*
To use the GUI first launch the InteractMakerGui.py.
Then a pop-up window will be opened.

To split a pdb into it's interactions, go to: menu bar > modules > PDB splitter.

This will load the PDB splitter interface.

the go to: menu bar > File > load and select the pdb that you want to split.

Once the pdb is loaded press, the Split PDB button. 
When the job finishes all the different interactions will be loaded in the unique pairs of interactions list.
Here you can examine the interactions obtained and rename them. 
Once you have finished checking the output you can save it pressing the save button or just send them to the PDB reconstruct by pressing the Send to PDB reconstruct button.
Once you have saved or/and send the output to the PDB reconstruct you will have to change the interface again back to the PDB reconstruct one.

To do so go to: menu bar > modules > PDB reconstruct

 *Step 2:*
if you have send the interactions in the previous step you will se that all the interactions are now loaded, if not go to:
menu bar > File > Open Folder and select the folder in which the interactions are saved. 

if you prefer to load the interactions files one by one, you can press the load button and then select the files that you want to load.
Once all the files are loaded press Run and the job will start.
A new protein called Macro-complex will appear in the list, this is the protein of the current step. There you can check how the protein is being constructed. 
In any moment you can pause the job pressing the stop button or kill it pressing the kill button. Once the job is finished go to: menu bar > File > Save to save the obtained macro-complex.

## FAQS

Question: When I press the Image or the Open in pymol button it gives an error but i have pymol installed and working. Why does this happens?

Answer: Sadly this happens because the pymol cmd is trash, and the python that should have pymol installed is the same that will be opened when you write python in the terminal. You can change the python interpreter that pymol will use by making a virtual environment or by changing the python path.

Question: How does the program knows how the pdb should be build?

Answer: Check the Documentation.md to have a detailed explanation of the approach used.

Question: What will happen if I give incorrect interactions?

Answer: The program assumes that all the input is valid input, if incorrect interactions are given, the resulting protein may be incorrect or incomplete. If the incorrect interaction are two protein that doesn't appear in any other interaction, the program will run fine as long as this interaction is not the first given. (because the program uses the first interaction as starting template).

Question: What will happen if the protein macro-complex that I'm trying to build has more than one possible solution?

Answer: If the macro-complex can have possible solutions the program will give the first one that it finds. Thats because trying all the possible combination will require a lot of recursivety and a lot of time, and we wanted to provide a fast an easy to understand approach.

Question: Is the program able to build a macro-complex that contains not common atoms such as metallic ions or small organic compounds?

Answer: Yes, the pdb loader that we have developed is able to recognise hetatoms and asign them to the chain they're interacting with. Thanks to that, our program is able to correctly place the hetatoms in the reconstructed structure.

For any other question contact with us:

XXXXXMAILSXXXXX

## References:

Hidrogen bond distance taken as 3.5A so it will be the distance managed by the scripts in the pagackes to find and check interactions between different molecules.

Martz, Eric; Help, Index & Glossary for Protein Explorer, http://www.umass.edu/microbio/chime/pe_beta/pe/protexpl/igloss.htm
Jeffrey, George A.; An introduction to hydrogen bonding, Oxford University Press, 1997.

For more detailed information consult the Documentation.



