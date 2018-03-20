import sys

if __name__ == "__main__":
    try:
        import pymol
        from pymol import cmd


        def photo_pdb(pdb_file, out_file):
            pymol.finish_launching(['pymol', '-qc'])
            cmd.reinitialize()
            cmd.load(pdb_file)
            cmd.show_as('cartoon', 'all')
            cmd.util.cbc()
            cmd.png(filename=out_file, width=350, height=350, ray=0, dpi=300)
            cmd.quit()


    except:
        exit(1)


    photo_pdb(sys.argv[1], "tmp/tmp_img")

