import tkinter
import pymol
from pymol import cmd

pymol.finish_launching(['pymol', '-qc'])

def photo_pdb(pdb_file, out_file):
    cmd.reinitialize()
    cmd.load(pdb_file)
    cmd.show_as('cartoon', 'all')
    cmd.util.cbc()
    cmd.png(filename=out_file, width=1024, height=926, ray=0, dpi=300)

photo_pdb('pdb/vira.pdb', 'pdb/vira.png')