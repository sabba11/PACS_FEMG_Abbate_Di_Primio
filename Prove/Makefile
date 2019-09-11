#Including Makefile.inc
-include Makefile.inc
#==============================================================================
CXX = g++
CXXFLAGS += -std=c++11

INCLUDES = -I$(BGLGEOM_DIR)/include -I. -I$(mkBGLInc) -I$(EIGEN_INC_DIR) -I$(VTK_INC_DIR) -I$(GETFEM_INC_DIR)

CPPFLAGS += $(INCLUDES) -D FEMG_VERBOSE_

.PHONY = all clean distclean

.DEFAULT_GOAL = all

all: main

main: main.o quantum_graph_problem.o eigen_problem.o
	$(CXX) -o main main.o quantum_graph_problem.o eigen_problem.o $(LDLIBS)

main.o: main.cpp eigen_problem.hpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c main.cpp

quantum_graph_problem.o: quantum_graph_problem.cpp quantum_graph_problem.hpp node.hpp mesh1d_prova_new_order.hpp descr_qg.hpp type_aliases.hpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c quantum_graph_problem.cpp

eigen_problem.o: eigen_problem.cpp eigen_problem.hpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c eigen_problem.cpp

clean:
	rm *.o

distclean: clean
	rm main
