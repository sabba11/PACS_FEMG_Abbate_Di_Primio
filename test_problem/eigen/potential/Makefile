#Including Makefile.inc
-include ../../../Makefile.inc
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
	@echo "--- Linking eigen simple potential ---"
	$(CXX) -o main main.o $(LDLIBS)
	@echo "--- COMPLETED ---"

main.o: main.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c main.cpp

clean:
	@echo " "
	@echo "--- Cleaning eigen simple potential objects ---"
	rm -f *.o
	@echo "--- COMPLETED ---"

distclean: clean
	@echo " "
	@echo "--- Cleaning eigen simple potential executable ---"
	rm -f main
	@echo "--- COMPLETED ---"
