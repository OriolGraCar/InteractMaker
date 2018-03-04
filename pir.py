from PDB4 import ProteinStructure as PS


def superimpose_to_pir(finalpdb,pdb_used,outputfile):
    """
    :param finalpdb: pdb object of the final structure
    :param pdb_used: list of tuples with the pdb file name in the first and in the second a string with the names
                     of the chains in the final pdb.
    :param outputfile: name to write the pir alignment
    :return: pir format alignment file"""

    pir_file = open(outputfile, 'w')
    pir_file.write(">%s" % finalpdb.get_id())
    chains_order = []
    for chain in finalpdb.getchains():
        chain_id = chain.get_id()
        chain_seq = chain.get_sequence_str()
        chains_order.append(chain_id, len(chain))
        pir_file.write("\n%s/" % chain_seq)

    pir_file.write("*\n")

    for pdb in pdb_used:
        current_pdb = PS("PR", pdb[0])
        pir_file.write(">%s" % current_pdb.get_id())
        chains = current_pdb.getchains()
        placed_chains = 0
        for i in range(len(chains_order)):
            if chains_order[i][0] == pdb[1][placed_chains]:
                pir_file.write("\n%s/" % chains[placed_chains].get_sequence_str())
            else:
                string_to_write = '-' * chains_order[i][1]
                pir_file.write("\n%s" % string_to_write)

    pir_file.write("*\n")
