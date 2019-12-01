#Including Makefile.inc
-include Makefile.inc
#==============================================================================

override CPPFLAGS := $(CPPFLAGS) $(INCLUDES) -DGMM_USES_LAPACK

.PHONY = all clean distclean

.DEFAULT_GOAL = all 

help:
	@echo "make help: ------- Prints this help"
	@echo "make all: -------- Makes libraries, test problems and the data builder"
	@echo "make library: ---- Makes both static and dynamic libraries"
	@echo "make eigen: ------ Builds eigen problems - library needed"
	@echo "make elliptic: --- Builds elliptic problems - library needed"
	@echo "make graph_maker:  Makes the code that can be used to build *.pts files"
	@echo "make resultclean:  Cleans all the exported files"
	@echo "make clean: ------ Cleans all object files"
	@echo "make distclean: -- Cleans all"

all: graph_maker library eigen elliptic

library:
	@echo " "
	@echo "--- Building library ---"
	$(MAKE) -C $(FEMG_DIR) all
	@echo "--- COMPLETED ---"

eigen:
	@echo " "
	@echo "--- Building eigen ---"
	$(MAKE) -C test_problem/eigen all
	@echo "--- COMPLETED ---"

elliptic:
	@echo " "
	@echo "--- Building elliptic ---"
	$(MAKE) -C test_problem/elliptic all
	@echo "--- COMPLETED ---"

graph_maker: 
	@echo " "
	@echo "--- Building graph maker ---"
	$(MAKE) -C data_builder/graph_maker all
	@echo "--- COMPLETED ---"

resultclean:
	$(MAKE) -C test_problem/elliptic resultclean
	$(MAKE) -C test_problem/eigen resultclean

clean:
	$(MAKE) -C test_problem/elliptic clean
	$(MAKE) -C test_problem/eigen clean
	$(MAKE) -C data_builder/graph_maker clean

distclean:
	$(MAKE) -C test_problem/elliptic distclean
	$(MAKE) -C test_problem/eigen distclean
	$(MAKE) -C $(FEMG_DIR) distclean
	$(MAKE) -C data_builder/graph_maker distclean
