# Location of BGLgeom library
BGLGEOM_DIR = /home/sguillo/Documents/Vecchi_progetti/Pacs_BGLgeom_Ilaria_Mattia-master/Install

# Location of Boost library
export mkBGLInc = /u/sw/pkgs/toolchains/gcc-glibc/7/pkgs/boost/1.63.0/include

# Path to Eigen includes
export mkEigenInc =   /home/sguillo/Documents/eigen/eigen

# Path to GetFEM library
GETFEM_LIB_DIR   = /usr/local/lib
# Path to GetFEM includes
GETFEM_INC_DIR   = /usr/local/include/getfem
# Path to Eigen includes
EIGEN_INC_DIR    = $(mkEigenInc)
# Path to BGL includes
BGL_INC_DIR      = $(mkBGLInc)

# Collect all includes
ALL_INC_DIR= $(EIGEN_INC_DIR) $(BGL_INC_DIR) $(GETFEM_INC_DIR)
includes := $(foreach inc,$(ALL_INC_DIR),-I$(inc))
INCLUDES += -I. -I./include $(includes)

# Linking libraries
LDLIBS += -lqhull -llapack -llapacke -lblas -lgetfem -lblas -llapacke -llapack -lqhull -L$(BGLGEOM_DIR)/lib -lBGLgeom

# Compiler flags
CXX = g++
CXXFLAGS += -std=c++11 -Wall








