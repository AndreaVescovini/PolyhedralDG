#--------------------------------------------
# Directories that have to be set by the user
#--------------------------------------------

# Enter the directory where the Eigen library headers are:
EIGEN_INC  = $(mkEigenInc)

# Enter the directory where the GetPot header is (needed only by the executables):
GETPOT_INC = /vagrant/pacs/Material/Examples/src/Utilities

# Enter the directory where you want to install the library
POLYDG_PATH = /vagrant/pacs/polyDG

#-------------------------------------------------------------------------------
#------------------- DO NOT MODIFY ANYTHING BELOW THIS LINE --------------------
#-------------------------------------------------------------------------------

#-------------------------------
# Compiler and compilation flags
#-------------------------------
CXX = g++

WFLAGS  += -Wall
STDFLAG += -std=c++11

ifdef RELEASE
 OPTFLAGS += -O3 -DNDEBUG
else
 OPTFLAGS += -g -DVERBOSITY
 WFLAGS   += -Wextra -pedantic
endif

CXXFLAGS += $(WFLAGS) $(STDFLAG) $(OPTFLAGS)
CPPFLAGS += -I$(EIGEN_INC)

LD_LIBS  += -lPolyDG

#------------
# Directories
#------------
POLYDG_INC = $(POLYDG_PATH)/include
POLYDG_LIB = $(POLYDG_PATH)/lib
