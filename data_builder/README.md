## Data Builder
In this folder the user can find the code to create the data of the graphs and their geometric properties.  
Everything is organized in two folders:
  - `/graph_maker`: containing the code to build those graphs
  - `/data`: containing the output files storing the informations of the graphs

### Graph Maker
This sub-directories contains two tipes of files:
  - `*.m files`: Those files are used to build in the `../data/folder` the `*.txt` containing the infos about the graphs 
  we want to build. Those output will read on every line an edge containg source node and target node coordinates and then
  their boundary conditions.  
  The files are named after the shape of the graph. To change the settings and their dimensions you have to read the comment
  at the start of those files and then change the variables that they tell you to.  
  Most files create planar graphs but you can chose the option to embed them in 3 dimensions by simply putting the third
  coordinate to 0.  
  Each of those file builds even the adjacency matrix.  
  In the folder `/matlab_functions` the user can find all the auxiliary functions used to build the graphs of every shape.
  - `*.cpp/*.hpp files`: These files are `ptsbuilder_N.cpp`,`ptsbuilder_h.cpp` and the header `reader_femg.hpp`.  
  This is a C++ code that takes the `*.txt` files and it gives as an output the meshed `*.pts` code.  
  The first one divides every edge in N sub-edges, the second one takes as an input the length of the step h and divides every
  edge by it (taking the closest length that divides it uniformly).  
  This is based on the `BGLGeom Library` by Ilaria Speranza and Mattia Tantardini.  
  After compiling the makefile you should run in the command line:  
  `$ ./ptsbuilder_N (int number_of_interval_in_an_edge) (../data/graph.txt) (../data/graph_number.pts)`  
  EXAMPLE: `$ ./ptsbuilder_N 5 ../data/star.txt ../data/star_5.pts`  
  Or:  
  `$ ./ptsbuilder_h (double length_of_the_step) (../data/graph.txt) (../data/graph_number.pts)`  
  EXAMPLE: `$ ./ptsbuilder_h 0.1 ../data/star.txt ../data/star_01.pts` 
- `Makefile`: Makefile to compile the code. It should not be called by users from this folder,
as it is called from the Makefile contained in the upper level directory (also for partial compilations).
- `Makefile.inc`: auxiliary Makefile. Remember to change the directories in this file and in the Makefile.inc outside.  

### Data
This folder is used to store the data containing the geometrical informations about the graphs.
Files in this folder are:
  - `*.m`: these matlab files build the file storing the radii of the edges taking as input the `graph.txt` and having as
  output the `radii.txt`
  - `radii.txt`: list of dimensions of the radii of the edges
  - `graph.txt`: vector of the edges of the graph (source and target coordinates and then their BCs)
  - `graph.pts`: pts files describing the graphs
  
