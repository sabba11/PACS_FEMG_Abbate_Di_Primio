#Including Makefile.inc
-include Makefile.inc
#==============================================================================
CXX = g++
CXXFLAGS += -std=c++11

INCLUDES = -I$(BGLGEOM_DIR)/include -I. -I$(EIGEN_INC_DIR) -I$(GETFEM_INC_DIR)

CPPFLAGS += $(INCLUDES)

.PHONY = all clean distclean

.DEFAULT_GOAL = all

all: writer

writer: writer.o 
	$(CXX) -o writer writer.o $(LDLIBS)

writer.o: writer.cpp reader_femg.hpp 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c writer.cpp

clean: 
	rm -f *.o

distclean: clean
	rm -f writer






