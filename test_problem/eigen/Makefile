#Including Makefile.inc
-include ../../Makefile.inc
#==============================================================================
CXX = g++
CXXFLAGS += -std=c++11

INCLUDES = -I$(FEMG_DIR)/include -I. -I$(mkBGLInc) -I$(GETFEM_INC_DIR)

CPPFLAGS += $(INCLUDES) -D FEMG_VERBOSE_ -D GMM_USES_LAPACK

.PHONY = all clean distclean

.DEFAULT_GOAL = all

all: main

main: main.o
	$(CXX) -o main main.o $(LDLIBS)

main.o: main.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c main.cpp

clean:
	rm -f *.o

distclean: clean
	rm -f main
