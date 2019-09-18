#Including Makefile.inc
-include Makefile.inc
#==============================================================================
CXX = g++
CXXFLAGS += -std=c++11

INCLUDES = -I$(FEMG_DIR)/include -I. -I$(mkBGLInc) -I$(GETFEM_INC_DIR)

CPPFLAGS += $(INCLUDES) -D FEMG_VERBOSE_ -D GMM_USES_LAPACK

.PHONY = all clean distclean

.DEFAULT_GOAL = all

all:
	$(MAKE) -C $(FEMG_DIR) all
	$(MAKE) -C test_problem/eigen all

library:
	$(MAKE) -C $(FEMG_DIR) all

eigen:
	$(MAKE) -C $(FEMG_DIR) all
	$(MAKE) -C test_problem/eigen all


clean:
	$(MAKE) -C test_problem/eigen clean

distclean:
	$(MAKE) -C test_problem/eigen distclean
	$(MAKE) -C $(FEMG_DIR) distclean
