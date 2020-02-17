#Including Makefile.inc
-include Makefile.inc
#==============================================================================

override CPPFLAGS := $(CPPFLAGS) $(INCLUDES) -DGMM_USES_LAPACK

.PHONY = all clean distclean

.DEFAULT_GOAL = all 

help:
	@echo "Usage: make [option] keyword"
	@echo " "
	@echo "List of possible keywords:"
	@echo "help: ------- Prints this help"
	@echo "doc: -------- Builds Doxygen documentation"
	@echo "all: -------- Makes libraries, test problems and the data builder"
	@echo "library: ---- Makes both static and dynamic libraries"
	@echo "eigen: ------ Builds eigen problems - library needed"
	@echo "elliptic: --- Builds elliptic problems - library needed"
	@echo "graph_maker:  Makes the code that can be used to build *.pts files"
	@echo "resultclean:  Cleans all the exported files"
	@echo "clean: ------ Cleans all object files"
	@echo "distclean: -- Cleans all"
	@echo " "
	@echo "List of possible options:"
	@echo "CPPFLAGS+=-DFEMG_VERBOSE_: Let user compile in verbose mode"

all: graph_maker library eigen elliptic doc clean

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

doc:
	@echo ""
	@echo "--- Making documentation ---"
	@mkdir -p $(DOC_PATH)
	doxygen Doxyfile
	$(MAKE) -C $(DOC_PATH)/latex
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
	-rm -rv -f $(DOC_PATH) $(INSTALL_PATH) 
