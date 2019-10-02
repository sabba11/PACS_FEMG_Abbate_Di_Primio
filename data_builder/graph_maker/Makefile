#Including Makefile.inc
-include Makefile.inc
#==============================================================================
CXX = g++
CXXFLAGS += -std=c++11

INCLUDES += -I$(BGLGEOM_DIR)/include

override CPPFLAGS := $(CPPFLAGS) $(INCLUDES)

.PHONY = all clean distclean

.DEFAULT_GOAL = all

all: ptsbuilder_N ptsbuilder_h

ptsbuilder_N: ptsbuilder_N.o
	@echo " "
	@echo "--- Linking ---" 
	$(CXX) -o ptsbuilder_N ptsbuilder_N.o $(LDLIBS)

ptsbuilder_N.o: ptsbuilder_N.cpp reader_femg.hpp 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c ptsbuilder_N.cpp

ptsbuilder_h: ptsbuilder_h.o
	@echo " "
	@echo "--- Linking ---" 
	$(CXX) -o ptsbuilder_h ptsbuilder_h.o $(LDLIBS)

ptsbuilder_h.o: ptsbuilder_h.cpp reader_femg.hpp 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c ptsbuilder_h.cpp

clean:
	@echo " "
	@echo "--- Cleaning graph makers objects ---"
	rm -f *.o

distclean: clean
	@echo " "
	@echo "--- Cleaning graph makers executable ---"
	rm -f ptsbuilder_h
	rm -f ptsbuilder_N





