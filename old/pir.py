def superimpose_to_pir(finalpdb,pdb_used,outputfile, chain_id_dict):
    """
    :param finalpdb: pdb object of the final structure
    :param pdb_used: list of tuples with the pdb file name in the first and in the second a string with the names
                     of the chains in the final pdb. Or a list of pdb object.
    :param outputfile: name to write the pir alignment
    :param chain_id_dict: dict with the equivalencies of names
    :return: pir format alignment file"""

    pir_file = open(outputfile, 'w')
    pir_file.write(">%s" % finalpdb.get_id())
    chains_order = []
    for chain in finalpdb:
        chain_id = chain.get_id()
        chain_seq = chain.get_sequence_str()
        chains_order.append((chain_id, len(chain)))
        pir_file.write("\n%s/" % chain_seq)

    pir_file.write("*")

    for pdb in pdb_used:
        pir_file.write("\n>%s" % pdb.get_id())
        for i in range(len(chains_order)):
            if chain_id_dict[chains_order[i][0]] in pdb.child_dict:
                pir_file.write("\n%s/" % pdb.child_dict[chain_id_dict[chains_order[i][0]]].get_sequence_str())
            else:
                string_to_write = '-' * chains_order[i][1]
                pir_file.write("\n%s" % string_to_write)

    pir_file.write("*\n")
    pir_file.close()
