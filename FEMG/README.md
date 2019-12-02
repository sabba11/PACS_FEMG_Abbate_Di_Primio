## Contents of the FEMG folder
All the C++ core code of the FEMG library is contained in the FEMG folder.  
The directory is organized as follows:
 - `/include` subfolder: contains all `.hpp` C++ headers. Class declarations are coded therein.
 - `/src` subfolder: contains all `.cpp` C++ source code files. Method definitions are coded therein.
 - `Makefile`: Makefile to compile the code. It should not be called by users from this folder, as it is called from the 
 Makefile contained in the upper level directory (also for partial compilations).
 - `Makefile.inc`: auxiliary Makefile. This file should not be modified.
