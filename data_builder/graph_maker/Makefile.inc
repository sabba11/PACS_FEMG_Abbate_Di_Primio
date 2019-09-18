# Location of BGLgeom library
BGLGEOM_DIR = /home/sguillo/Documents/Vecchi_progetti/Pacs_BGLgeom_Ilaria_Mattia-master/Install

# Location of Boost library
export mkBGLInc = /u/sw/pkgs/toolchains/gcc-glibc/7/pkgs/boost/1.63.0/include

# Collect all includes
ALL_INC_DIR= $(EIGEN_INC_DIR) $(BGL_INC_DIR)
includes := $(foreach inc,$(ALL_INC_DIR),-I$(inc))
INCLUDES += -I. -I./include $(includes)

# Including headers from the BGLgeom library
CPPFLAGS += -I$(BGLGEOM_DIR)/include $(INCLUDES)

# Linking libraries
LDLIBS += -lqhull -llapack -llapacke -lblas -lgetfem -lblas -llapacke -llapack -lqhull -L$(BGLGEOM_DIR)/lib -lBGLgeom

# Path to Eigen includes
export mkEigenInc = /u/sw/pkgs/toolchains/gcc-glibc/7/pkgs/eigen/3.3.3/include/eigen3

# Path to GetFEM library
GETFEM_LIB_DIR   = /usr/local/lib
# Path to GetFEM includes
GETFEM_INC_DIR   = /usr/local/include/getfem
# Path to Eigen includes
EIGEN_INC_DIR    = $(mkEigenInc)
# Path to BGL includes
BGL_INC_DIR      = $(mkBGLInc)




