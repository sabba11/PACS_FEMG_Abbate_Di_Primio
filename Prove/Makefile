#Including Makefile.inc
-include Makefile.inc
#==============================================================================
CXX = g++
CXXFLAGS += -std=c++11

INCLUDES = -I$(BGLGEOM_DIR)/include -I. -I$(mkBGLInc) -I$(EIGEN_INC_DIR) -I$(VTK_INC_DIR) -I$(GETFEM_INC_DIR)

CPPFLAGS += $(INCLUDES)

.PHONY = all clean distclean

.DEFAULT_GOAL = all

all: main

main: main.o 
	$(CXX) -o main main.o $(LDLIBS)

main.o: main.cpp reader_femg.hpp 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c main.cpp

clean: 
	rm *.o

distclean: clean
	rm main






