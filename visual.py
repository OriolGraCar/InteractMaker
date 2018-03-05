import tkinter
import pymol
from pymol import cmd

pymol.finish_launching(['pymol','-qc'])
cmd.load('pdb/junto.pdb')
cmd.png('pdb/junto.png')
