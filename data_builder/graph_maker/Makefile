#Including Makefile.inc
-include Makefile.inc
#==============================================================================
CXX = g++
CXXFLAGS += -std=c++11

INCLUDES += -I$(BGLGEOM_DIR)/include

override CPPFLAGS := $(CPPFLAGS) $(INCLUDES)

.PHONY = all clean distclean

.DEFAULT_GOAL = all

all: writer

writer: writer.o
	@echo " "
	@echo "--- Linking ---" 
	$(CXX) -o writer writer.o $(LDLIBS)

writer.o: writer.cpp reader_femg.hpp 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c writer.cpp

clean:
	@echo " "
	@echo "--- Cleaning graph makers objects ---"
	rm -f *.o

distclean: clean
	@echo " "
	@echo "--- Cleaning graph makers executable ---"
	rm -f writer






