import sys

def load_pdb(file_name, chain1, chain2, n):
    pdb = dict()
    with open(file_name, "r") as file:
        for line in file:
            line = line.rstrip()
            if line.startswith('ATOM'):
                elemnts_lst = line.split()
                if elemnts_lst[4] == chain1 or elemnts_lst[4] == chain2:
                    elemnts_lst[6] = str(round(float(elemnts_lst[6]) + n,3))
                    elemnts_lst[7] = str(round(float(elemnts_lst[7]) + n,3))
                    elemnts_lst[8] = str(round(float(elemnts_lst[8]) + n,3))
                    print("%-6s%5s %4s %3s %s%4s    %8s%8s%8s"
                          %(elemnts_lst[0], elemnts_lst[1], elemnts_lst[2],
                            elemnts_lst[3], elemnts_lst[4], elemnts_lst[5],
                            elemnts_lst[6], elemnts_lst[7], elemnts_lst[8]))

load_pdb(sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]))