# -*- coding: utf-8 -*-
"""
Spyder Editor

This script will convert the DIMACS .graph files to .edg files for the zero-forcing codes.
"""

if __name__ == "__main__":
    num_nodes = ""
    num_edges = ""
    for basefilename in ["adjnoun", "celegansneural", "chesapeake", "dolphins", "football", "jazz", "karate", "lesmis", "polbooks"]:
        infilename = basefilename + ".graph"
        outfilename = basefilename + ".edg"
        counter = 0
        with open(infilename, 'r') as infile:
            line = infile.readline()
            line_split = line.split()
            num_nodes = line_split[0]
            num_edges = line_split[1]
            with open(outfilename, 'w') as outfile:
                outfile.write(num_nodes + " " + num_edges + "\n")
                line = infile.readline()
                while(line):
                    line_split = line.split()
                    for i in range(len(line_split)):
                        if int(line_split[i]) - 1 > counter:
                            outfile.write(str(counter) + " " + str(int(line_split[i]) - 1) + "\n")
                    line = infile.readline()
                    counter += 1
            