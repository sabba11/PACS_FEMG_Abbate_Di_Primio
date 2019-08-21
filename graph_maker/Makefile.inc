# Local directories
TARGET_DIR  = bin
SOURCE_DIR  = src
BUILD_DIR   = obj
INCLUDE_DIR = include
DOC_DIR     = doc
DATA_DIR    = data

# File names
SRCS          = $(wildcard $(SOURCE_DIR)/*.cpp)
SRCS_MAIN     = $(filter $(SOURCE_DIR)/main%.cpp, $(SRCS))
SRCS_NOT_MAIN = $(filter-out $(SOURCE_DIR)/main%.cpp, $(SRCS))
OBJS          = $(addprefix $(BUILD_DIR)/, $(notdir $(SRCS:.cpp=.o)))
OBJS_MAIN     = $(addprefix $(BUILD_DIR)/, $(notdir $(SRCS_MAIN:.cpp=.o)))
OBJS_NOT_MAIN = $(addprefix $(BUILD_DIR)/, $(notdir $(SRCS_NOT_MAIN:.cpp=.o)))
EXEC          = $(addprefix $(TARGET_DIR)/, $(notdir $(SRCS_MAIN:.cpp=)))

# Location of BGLgeom library
BGLGEOM_DIR = /home/sguillo/Documents/Progetto_PACS/PACS_Project_Abbate_DiPrimio/Pacs_BGLgeom_Ilaria_Mattia-master/Install

# Location of Boost library
export mkBGLInc = /u/sw/pkgs/toolchains/gcc-glibc/7/pkgs/boost/1.63.0/include

# For all code
# LDLIBS = -L$(VTK_LIB_DIR) $(foreach lib, $(VTK_LIBS), -l$(lib))

# Collect all includes
ALL_INC_DIR= $(EIGEN_INC_DIR) $(VTK_INC_DIR) $(BGL_INC_DIR)
includes := $(foreach inc,$(ALL_INC_DIR),-I$(inc))
INCLUDES += -I. -I./include $(includes)

# Including headers from the BGLgeom library
CPPFLAGS += -I$(BGLGEOM_DIR)/include $(INCLUDES)

# Linking libraries
LDLIBS += -lqhull -llapack -llapacke -lblas -lgetfem -lblas -llapacke -llapack -lqhull -L$(BGLGEOM_DIR)/lib -lBGLgeom

# Path to VTK includes
export mkVtkInc   = /usr/local/include/vtk-8.2
export mkVtkLib   = /usr/local/lib
# Path to Eigen includes
export mkEigenInc = /u/sw/pkgs/toolchains/gcc-glibc/7/pkgs/eigen/3.3.3/include/eigen3

# Path to GetFEM library
GETFEM_LIB_DIR   = /usr/local/lib
# Path to GetFEM includes
GETFEM_INC_DIR   = /usr/local/include/getfem
# Path to MANworks library
MANWORKS_LIB_DIR = $(mkMANworksHome)/lib
# Path to Manworks includes
MANWORKS_INC_DIR = $(mkMANworksHome)/include
# Path to Eigen includes
EIGEN_INC_DIR    = $(mkEigenInc)
# Path to VTK includes
VTK_INC_DIR      = $(mkVtkInc)
# Path to VTK library
VTK_LIB_DIR      = $(mkVtkLib)
# Path to BGL includes
BGL_INC_DIR      = $(mkBGLInc)
# Path to MANworks include (fluid)
MANWORKS_INC_DIR2 = $(mkMANworksFluid)/include
# Path to MANworks library
MANWORKS_LIB_DIR2 = $(mkMANworksFluid)/lib



