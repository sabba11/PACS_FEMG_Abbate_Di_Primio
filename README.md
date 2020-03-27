# FEMG library (Finite Element Methods on Graphs)
Coding a library to define, analyze and solve differential problems on graphs.

This work is the final project for the course in Advanced Programming for Scientific Computing (a.y. 2018/19) held at Politecnico di Milano (Milan, Italy).

The repository contains a report with some detailed theoretical and numerical analysis made by the two authors on differential problem on graphs. Particular attention is given to eigenvalue problems.

Authors:  
Stefano Abbate (stefano.abbate@mail.polimi.it)  
Andrea Di Primio (andrea.diprimio@mail.polimi.it)

## Dependencies
The library depends on several external libraries.
##### **GetFEM++ (ver. 5.3 or higher)**
The GetFEM++ library can be downloaded [here.](http://getfem.org/download.html "GetFEM++ download page")  
For export routines, the GetFEM++\-MATLAB interface is necessary.

##### **Boost Graph Library (ver. 1.63 or higher)**
The Boost Graph Library (BGL) can be downloaded [here.](https://www.boost.org/doc/libs/1_63_0/libs/graph/doc/index.html "BGL download page")

##### **LAPACK**
The Linear Algebra PACKage (LAPACK) can be downloaded [here.]( http://www.netlib.org/lapack/#_software "LAPACK download page")  
On Linux, it can be installed using the Package Manager executing  
`$ sudo apt-get install liblapack-dev`  
on the terminal.  

##### **BGLgeom**
The BGLgeom library is not strictly necessary to build the project, although it is used to create data files (denoted by the format `.pts`). In all test cases, `.pts` files are pre-generated and read at runtime. To handle custom graphs data files refer to the README.md in the folder data_builder/graph_maker.
The BGLgeom repository can be found [here.](https://github.com/lformaggia/Pacs_BGLgeom_Ilaria_Mattia "BGLgeom repository")

## Installation
To build the library, firstly fill the files Makefile.inc and data_builder/graph_maker/Makefile.inc as indicated therein. Only empty fields have to be modified. Then, open the terminal and run  
`$ make all`  
to compile the whole project in the specified installation folder. Partial compilations can be made, type   
`$ make help`  
for the list of possible keywords.

To compile the code in verbose mode, use  
`$ make CPPFLAGS+=-DFEMG_VERBOSE_ all` 

Code documentation will be generated in both html and LaTeX format (LaTeX compiler eventually needed). To this end, alongside `doxygen`, the `graphviz` package is also needed.

#### Installation issues
A known installation issue in Linux involves BGL library linking. A workaround is to modify the Linux environment variable `LD_LIBRARY_PATH`, so to include the path to the Boost `lib` folder. This can be achieved through the command  
`$ export LD_LIBRARY_PATH=path/to/boost/lib:$LD_LIBRARY_PATH`  
every time the project is compiled from a new terminal session. To modify the variable once and for all, execute from the terminal  
`$ sudo nano ~/.bashrc`  
and add the `export` command above (without the initial `$` sign) as the last line of the file.

A known issue involves setting up the Makefile of the BGLgeom library. First, modify the following in Makefile.inc:
1. `libvtk*.so.1` in place of `libvtk*.so.5.10.1` at line 69 (if not using VTK 5.10);
2. `$(basename $(basename $(VTK_LIBS2)))` in place of `$(basename $(basename $(basename $(basename $(VTK_LIBS2)))))` at line 70 (if not using VTK 5.10);

Then, in libBGLgeom/Makefile:

3. `-soname` in place of `-$(SONAME)` at line 62.

A known issue involves the installation of the GetFEM++-MATLAB interface. In the file getfem-5.3/interface/src/matlab/gfm_commmon.c comment lines 45 and 107, both referring to `mxSPARSE_CLASS`. This only concerns the execution of the MATLAB routines provided in this repository. 

A known issue involves the usage of modules to load libraries. The code has been tested installing the external dependencies locally.
## Execution and postprocessing
For details on how to execute the code and postprocess the results, refer to the `README.md` in the `test_problem` folder.

## Dev settings
The code has been tested on Linux Ubuntu 18.04. 


