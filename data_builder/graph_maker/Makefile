#Including Makefile.inc
-include Makefile.inc
#==============================================================================
CXX = g++
CXXFLAGS += -std=c++11

INCLUDES += -I$(BGLGEOM_DIR)/include

override CPPFLAGS := $(CPPFLAGS) $(INCLUDES)

.PHONY = all clean distclean

.DEFAULT_GOAL = all

all: uniform_writer equisplit_writer

equisplit_writer: equisplit_writer.o
	@echo " "
	@echo "--- Linking ---" 
	$(CXX) -o equisplit_writer equisplit_writer.o $(LDLIBS)

equisplit_writer.o: equisplit_writer.cpp reader_femg.hpp 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c equisplit_writer.cpp

uniform_writer: uniform_writer.o
	@echo " "
	@echo "--- Linking ---" 
	$(CXX) -o uniform_writer uniform_writer.o $(LDLIBS)

uniform_writer.o: uniform_writer.cpp reader_femg.hpp 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c uniform_writer.cpp

clean:
	@echo " "
	@echo "--- Cleaning graph makers objects ---"
	rm -f *.o

distclean: clean
	@echo " "
	@echo "--- Cleaning graph makers executable ---"
	rm -f uniform_writer
	rm -f equisplit_writer





