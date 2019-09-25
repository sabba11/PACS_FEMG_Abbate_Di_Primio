# ====================================================================
#   					"BGLgeom library"
#      Course on Advanced Programming for Scientific Computing
#                     Politecnico di Milano
#                         A.Y. 2015-2016
#
#		Copyright (C) 2017 Ilaria Speranza & Mattia Tantardini
#=====================================================================

# Including general Makefile.inc of the project
MAKEFILEH_DIR=../
include $(MAKEFILEH_DIR)/Makefile.inc

# Including Makefile.inc of the local folder
-include Makefile.inc

#=========================================================================

override CPPFLAGS := $(CPPFLAGS) $(INCLUDES) -D GMM_USES_LAPACK

.PHONY: help all clean distclean doc library static dynamic install alltest

.DEFAULT_GOAL: help

help:
	@echo "make help: ----- Prints this help"
	@echo "make all: ------ Makes libraries, tests and documentation"
	@echo "make static: --- Makes static library and executable"
	@echo "make dynamic: -- Makes dynamic library and executable"
	@echo "make library: -- Makes both static and dynamic libraries"
	@echo "make install: -- Installs the libaries"
	@echo "make clean: ---- Cleans all object files"
	@echo "make distclean:  Cleans all"
	
all: library install clean

library: static dynamic

static: $(LIB_OBJS_STATIC)
	@echo " "
	@echo "--- Building static library ---"
	@mkdir -p $(LIB_DIR)
	ar -r -s $(STATIC_LIB) $(LIB_OBJS_STATIC)
	ranlib $(STATIC_LIB)
	@echo "--- COMPLETED ---"
	
$(BUILD_DIR_STATIC)/%.o: $(SOURCE_DIR)/%.cpp $(INCLUDE_DIR)/%.hpp
	@echo " "
	@echo "--- Building static object ---"
	@mkdir -p $(BUILD_DIR_STATIC)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

	
dynamic: CXXFLAGS += -fPIC
dynamic: $(LIB_OBJS_DYNAMIC)
	@echo " "
	@echo "--- Building dynamic library ---"
	@mkdir -p $(LIB_DIR)
	$(CXX) -shared -Wl,-soname,$(DYNAMIC_LIBFILE).1 $(LIB_OBJS_DYNAMIC) -o $(DYNAMIC_LIB).1.0
	@echo "--- COMPLETED ---"
	
$(BUILD_DIR_DYNAMIC)/%.o: $(SOURCE_DIR)/%.cpp $(INCLUDE_DIR)/%.hpp
	@echo " "
	@echo "--- Building dynamic object ---"
	@mkdir -p $(BUILD_DIR_DYNAMIC)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@


install:
	@echo ""
	@echo "--- Installing library ---"
	@mkdir -p $(INSTALL_PATH)/include
	@cp $(HEADERS) $(INSTALL_PATH)/include
	@mkdir -p $(INSTALL_PATH)/lib
	@cp $(DYNAMIC_LIB).1.0 $(STATIC_LIB) $(INSTALL_PATH)/lib
	@echo "--- INSTALLED ---"


clean:
	@echo ""
	@echo "--- Removing objects used to build the library ---"	
	-rm -rv -f $(BUILD_DIR)


distclean: clean
	@echo ""
	@echo "--- Removing library ---"
	-rm -rv -f $(LIB_DIR) $(DATA_DIR)

