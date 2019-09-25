#Including Makefile.inc
-include ../../Makefile.inc
#==============================================================================
CXX = g++
CXXFLAGS += -std=c++11

INCLUDES += -I$(FEMG_DIR)/include

override CPPFLAGS := $(CPPFLAGS) $(INCLUDES) -D GMM_USES_LAPACK

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
