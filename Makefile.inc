# Directory where the project resides.
PROJECT_DIR = 

# Installation path  
INSTALL_PATH = $(PROJECT_DIR)/install

# Documentation path  
DOC_PATH = $(PROJECT_DIR)/doc

# The place where we do find the files that let us install the library
FEMG_DIR = $(PROJECT_DIR)/FEMG

# Path to GetFEM library
GETFEM_LIB_DIR   = 
# Path to GetFEM includes
GETFEM_INC_DIR   = 
# Path to BGL includes
BGL_INC_DIR      = $(mkBGLInc)

# Location of Boost library
export mkBGLInc = 

# Collect all includes
ALL_INC_DIR= $(GETFEM_INC_DIR) $(BGL_INC_DIR)
includes := $(foreach inc,$(ALL_INC_DIR),-I$(inc))
INCLUDES += -I. -I./include $(includes)

# Linking libraries
LDLIBS += -L$(mkBGLInc)/../lib -lboost_filesystem -L$(mkBGLInc)/../lib -lboost_system -L$(FEMG_DIR)/lib -lFEMG -lqhull -llapack -lblas -L$(GETFEM_LIB_DIR) -lgetfem -lblas -llapack -lqhull -L$(FEMG_DIR)/lib -lFEMG  -L$(mkBGLInc)/../lib -lboost_system -L$(mkBGLInc)/../lib -lboost_filesystem -L$(mkBGLInc)/../lib






