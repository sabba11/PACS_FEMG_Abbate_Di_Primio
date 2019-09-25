#Including Makefile.inc
-include ../../Makefile.inc
#==============================================================================
CXX = g++
CXXFLAGS += -std=c++11

INCLUDES += -I$(FEMG_DIR)/include

override CPPFLAGS := $(CPPFLAGS) $(INCLUDES) -DGMM_USES_LAPACK

.PHONY = all clean distclean

.DEFAULT_GOAL = all

all: main

main: main.o
	@echo " "
	@echo "--- Linking ---"
	$(CXX) -o main main.o $(LDLIBS)

main.o: main.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c main.cpp

clean:
	@echo " "
	@echo "--- Cleaning eigen objects ---"
	rm -f *.o

distclean: clean
	@echo " "
	@echo "--- Cleaning eigen executable ---"
	rm -f main
