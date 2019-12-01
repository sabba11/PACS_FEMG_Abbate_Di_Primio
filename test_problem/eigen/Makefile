#Including Makefile.inc
-include Makefile.inc
-include ../../Makefile.inc

#==============================================================================

override CPPFLAGS := $(CPPFLAGS) $(INCLUDES) -DGMM_USES_LAPACK

.PHONY: all $(EIGEN_TEST)

all: $(EIGEN_TEST)

$(EIGEN_TEST):
	$(MAKE) -C $@  $(MAKECMDGOALS)

resultclean:
	rm -r -f $(RESULTS)

clean: $(EIGEN_TEST)

distclean: $(EIGEN_TEST)
#resultclean

