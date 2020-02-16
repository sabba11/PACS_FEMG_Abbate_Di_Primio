## Contents of the test_problem folder
The test_problem folder contains the implementation of several test cases.
The directory is organized as follows:
 - `/eigen` subfolder: contains all test cases concerning elliptic eigenvalue problems.
 - `/elliptic` subfolder: contains all test cases concerning elliptic boundary value problems.
 - `/matlab_functions` subfolder: auxiliary MATLAB routines to postprocess the results.
 
The lower level directories of the `/eigen` and `/elliptic` folders are:
 - `/laplacian` subfolder: contains code and setup files to treat the standard Laplacian case.
 - `/hamiltonian` subfolder: contains code and setup files to treat the simple Hamiltonian case (unitary potential).
 - `/potential` subfolder: contains code and setup files to treat the general Hamiltonian case.
These subfolders contain Makefiles for compilation and the source file `main.cpp` which, after compilation, generates the executable file, called `main`. Users should use the external Makefile routine to compile the code.

In order to run the executables from a terminal session, move in one of the lower level directories of the `/eigen` and `/elliptic` folders and execute the command

`$ ./main path/to/param/file`.

Some setup files are ready to be used in third-level subfolders (named after the graph topologies). To postprocess the results in MATLAB, open the `.m` file in the `/eigen` or `/elliptic` folders and follow the indications provided therein.
