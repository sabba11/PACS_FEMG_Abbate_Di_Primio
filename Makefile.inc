# Directory where the project resides.
PROJECT_DIR = /home/sguillo/Desktop/PACS_new

# Installation path  
INSTALL_PATH = /home/sguillo/Desktop/PACS_new/install

# The place where we do find the files that let us install the library
FEMG_DIR = $(PROJECT_DIR)/FEMG

# Path to GetFEM library
GETFEM_LIB_DIR   = /usr/local/lib
# Path to GetFEM includes
GETFEM_INC_DIR   = /usr/local/include/getfem
# Path to BGL includes
BGL_INC_DIR      = $(mkBGLInc)

# Location of Boost library
export mkBGLInc = /u/sw/pkgs/toolchains/gcc-glibc/7/pkgs/boost/1.63.0/include

# Collect all includes
ALL_INC_DIR= $(EIGEN_INC_DIR) $(GETFEM_INC_DIR) $(BGL_INC_DIR)
includes := $(foreach inc,$(ALL_INC_DIR),-I$(inc))
INCLUDES += -I. -I./include $(includes)

# Including headers from the FEMG library
CPPFLAGS += -I$(FEMG_DIR)/include $(INCLUDES)

# Linking libraries
LDLIBS += -L$(mkBGLInc)/../lib -lboost_filesystem -L$(mkBGLInc)/../lib -lboost_system -L$(FEMG_DIR)/lib -lFEMG -lqhull -llapack -llapacke -lblas -lgetfem -lblas -llapacke -llapack -lqhull -L$(FEMG_DIR)/lib -lFEMG  -L$(mkBGLInc)/../lib -lboost_system -L$(mkBGLInc)/../lib -lboost_filesystem -L$(mkBGLInc)/../lib






