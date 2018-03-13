import subprocess
import sys
import re


def topology_build(script_template,pdb_name):
    script_name = pdb_name + script_template
    scriptfile = open(script_name,'w')
    with open ("/amberscripts/basictleap.in") as fl:
        for line in fl:
            line = line.rstrip()
            new_line = re.sub('complex_name',pdb_name,line)
            scriptfile.write(new_line + "\n")
    scriptfile.close()
    try:
        subprocess.run(["tleap", "-f %s" % outputfile_name], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except FileNotFoundError:
        sys.stderr.write("AmberTools is not properly installed")
        exit(1)


def runsander(script_template,pdb_name):

    try:
        subprocess.run(["sander", "-O", "-i %s " % script_template, "-o %s_energys.out" % pdb_name, "-p %s.prmtop" % pdb_name, "-c %s.inpcrd" % pdb_name,
                        "-r %s.restrt" % pdb_name], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except FileNotFoundError:
        sys.stderr.write("AmberTools is not properly installed")
        exit(1)


def createnewpdb(pdb_name,output_name):
    file = open(output_name,'w')
    try:
        subprocess.run(["ambpdb", "-p %s.prmtop" % pdb_name, "-c %s.rst" % pdb_name], stdout=file, stderr=subprocess.PIPE)
    except FileNotFoundError:
        sys.stderr.write("AmberTools is not properly installed")
        exit(1)
    file.close()