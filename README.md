# zero_forcing
This repository contains code for finding minimum zero forcing sets and minimum connected zero forcing sets of graphs.  The methods for the zero forcing problem are in the basic_forcing directory, and the methods for the connected zero forcing problem are in the connected_forcing directory.

The .edg files for various test instances are in the test_graphs directory.  The .edg format is currently the only format supported.  These files are basically an edge list where the first line of the file gives the number of nodes and then the number of edges in the graph.

The best way to run these methods is to download the whole repository and follow the compilation and usage instruction in the source file for each program.  However, the methods can also be downloaded and run individually with the caveat that the integer programming methods rely on the forcing_lib.hpp header file.
