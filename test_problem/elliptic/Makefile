#Including Makefile.inc
-include Makefile.inc
-include ../../Makefile.inc

#==============================================================================

override CPPFLAGS := $(CPPFLAGS) $(INCLUDES) -DGMM_USES_LAPACK

.PHONY: all $(ELLIPTIC_TEST)

all: $(ELLIPTIC_TEST)

$(ELLIPTIC_TEST):
	$(MAKE) -C $@  $(MAKECMDGOALS)

resultclean:
	rm -r -f $(RESULTS)

clean: $(ELLIPTIC_TEST)

distclean: $(ELLIPTIC_TEST)
#resultclean
