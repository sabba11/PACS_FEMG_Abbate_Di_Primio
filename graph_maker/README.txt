In this folder we find two tipes of file:


*.m files: Those files are used to build in the ../data/folder
           the *.txt containing the info about the graphs we want
           to build.
           To change the settings you have to read the comment at the
           start of those files and to change the variables that
           they tell you to.
           For now we only create 3d graphs even if they are on a plane
           by putting the third coordinate at 0.
           If it is useful in these files even the adjacency matrix is
           built.

Everything else
(the c++ coding part): Those files are writer.cpp, the header reader_femg and
                       the make files. (REMEMBER TO SET UP THE MAKEFILE.INC!)
                       This is a c++ code that takes the *.txt files and it gives
                       as an output the meshed *.pts code.
                       This is based on the BGLGeom Library by Ilaria Speranza and
                       Mattia Tantardini.

                       After running the main you should run in the command line:
                       ./writer (int number_of_interval_in_an_edge)
                       (../data/graph.txt)  (../data/graph_number.pts)
                       EXAMPLE:
                       ./writer 5 ../data/simple_star.txt ../data/simple_star_5.pts
