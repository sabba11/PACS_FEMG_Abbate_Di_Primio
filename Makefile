#Including Makefile.inc
-include Makefile.inc
#==============================================================================
CXX = g++
CXXFLAGS += -std=c++11

override CPPFLAGS := $(CPPFLAGS) $(INCLUDES) -D GMM_USES_LAPACK

.PHONY = all clean distclean

.DEFAULT_GOAL = all

help:
	@echo "make help: ------- Prints this help"
	@echo "make all: -------- Makes libraries, test problems and the data builder"
	@echo "make library: ---- Makes both static and dynamic libraries"
	@echo "make eigen: ------ Makes the libraries and then builds eigen problems"
	@echo "make graph_maker:  Makes the code that can be used to build *.pts files"
	@echo "make clean: ------ Cleans all object files"
	@echo "make distclean: -- Cleans all"

all: graph_maker library
	@echo " "
	@echo "--- Building eigen ---"
	$(MAKE) -C test_problem/eigen all
	@echo "--- COMPLETED ---"

library:
	@echo " "
	@echo "--- Building library ---"
	$(MAKE) -C $(FEMG_DIR) all
	@echo "--- COMPLETED ---"

eigen: library
	@echo " "
	@echo "--- Building eigen ---"
	$(MAKE) -C test_problem/eigen all
	@echo "--- COMPLETED ---"

graph_maker:
	@echo " "
	@echo "--- Building graph maker ---"
	$(MAKE) -C data_builder/graph_maker all
	@echo "--- COMPLETED ---"

clean:
	$(MAKE) -C test_problem/eigen clean
	$(MAKE) -C data_builder/graph_maker clean

distclean:
	$(MAKE) -C test_problem/eigen distclean
	$(MAKE) -C $(FEMG_DIR) distclean
	$(MAKE) -C data_builder/graph_maker distclean
