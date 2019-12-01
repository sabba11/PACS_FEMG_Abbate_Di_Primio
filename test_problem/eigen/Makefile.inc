# Subdirectorier and result libraries
EIGEN_TEST = $(wildcard */.)
RESULT_TEST = $(wildcard */*/.)
RESULTS := $(foreach dir,$(RESULT_TEST),$(dir)/export)
