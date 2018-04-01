# Intro to do

# Index

	1- Theoretical approach used

	2- Program Limitations

	3- BioMacromplex Package explanation
		3.1 - PDB.py
		3.2 - PDBaligner.py
		3.3 - PDB_split.py
		3.4 - pymolmanager.py
		3.5 - Sequences.py

	4- Gui explanaion



# 1. Theoretical approach:

The aim of this project is to build a protein Macro-complex using as input the pair of interactions that happen in the inner core of the protein. These Macro-complexes can also contain DNA, RNA or other compounds such as ions and small organic molecules.

The idea behind is the following: if we choose one interaction pair as the starting template, we can then place the rest of the chains by superposition. Due to the fact that, in a macro-complex, at least another interaction pair has the same chain, then we can align those chains and after that, move the different chain responsible of the interaction to the template that we are building, check image 1. Just like a puzzle, if we keep aligning the similar proteins and moving the chains from the interaction to the template we can eventually obtain the correct macro-complex.

To obtain the correct structure using this approach we need to know which proteins have to be aligned and in which order, to avoid making unnecessary clashes by adding more than one time the same chain, and also know when to stop adding chains to that protein because it already has all its interactions.

This can be achieved in different ways. The easiest one would be to use a brute-force approach, which means trying all the combinations and removing those combinations that generate clashes. However, this approach is too computationally dependent, especially when trying to build huge Macro-complexes since the amount of possible combinations grows exponentially. That's why we wanted to find another approach that could be faster and easier to understand, based on the biological properties of the proteins.

Our approach is based on the following: Proteins are usually on solution, and they obtain their folds and interact with each other in order to cover the hydrophobic regions, because by hiding the hydrophobic parts the entropy of the water decreases and the system is more stable. That means that, if two proteins are interacting by one side, that side of the protein will probably contain hydrophobic residues and so, it needs to be covered. That means that in the final macro-complex all the chains will probably have the sides that are interacting always covered.

Knowing that, the program loads into the memory all the residues that are interacting in each different protein and then starts building the macro-complex by forcing the condition that all the residues found to be interacting in a given protein must be interacting with something.
Thanks to that, the program can check in each step which proteins are lacking interactions and which of the interaction pairs provided in the input can be used to cover those residues. By doing that we decrease the number of combinations that we have to try.
Finally, when all the proteins have all the interactions, that means that the macro-complex has been build and the program knows that it have to finish.

This approach works also with non proteic chains such as DNA and RNA, Because this ensures that all the nucleotides maintain their Watson and Crick interactions and because it also ensures that the protein is interacting properly with the DNA/RNA strand (interacting in the correct place of the major groove for example).

Another advantage of our approach is that if the user gives redundant interactions pairs or unnecessary ones (not incorrect) our program will still go as fast as if those interactions pairs where not provided. That happens because, as stated before, our program doesn't try all the possible combinations, so if two residues are found to be interacting twice, the program will just ignore one of them, and if more than the necessary interactions are given, those interactions will probably be fulfilled when placing another chain, so the program will not need to make a new superposition.

# 2. Program limitations

Despite the fact that our program is quite versatile(it can deal with different kinds of chains, protein symmetry and even excessive input without affecting it's performance) in some special cases the program fails to give the correct solution or gives a partial correct one, those cases are:

*1. When there is more than one correct structure given the data.*
Our program ensures to give a good looking structure given the data, because it makes sure that all the interactions are on their place, which is very useful when there is only one correct structure. However, when the interactions provided can lead to more than one correct structure our program will give only one of them. For example with the case of XXXXX, the program will give any of the two depending in the order in which the interactions pairs are provided.

*2. When incorrect interactions are provided*

Another limitation happens when the user give incorrect data, everything is build assuming that the data is correct, so if the data have problems is probable that the output will be incorrect. This will depend in which kind of incorrect data is given and in which position in the order of the input is placed. 
The only situation in which the incorrect data will not affect the output is when the interaction pair is made of proteins that are not part of the final structure and it's not the first input interaction given (because our program uses the first one as template for the reconstruction).

*3. When the structure can grow unlimited such as a microtubul*

Our program will be able to correctly build this structure but as the termination condition will never be meted the program will never end. This can be solved by saving the steps and just killing the program when the desired length is achieved or by using the GUI and stopping the job when the structure looks fine to you.

