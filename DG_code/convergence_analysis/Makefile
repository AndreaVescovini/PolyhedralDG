#------------
# Directories
#------------
EIGEN_INC  = $(mkEigenInc)
GETPOT_INC = /vagrant/pacs/Material/Examples/src/Utilities
POLYDG_INC = /vagrant/pacs/progetto_codici/DG_code/libPolyDG/include
POLYDG_LIB = /vagrant/pacs/progetto_codici/DG_code/libPolyDG/lib

#------------------------
# Sources and executables
#------------------------
SRCS = $(wildcard *.cpp)
EXE  = $(basename $(SRCS))

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
CPPFLAGS += -I$(EIGEN_INC) -I$(GETPOT_INC) -I$(POLYDG_INC)

LD_FLAGS += -L$(POLYDG_LIB)
LD_LIBS  += -lPolyDG

#--------
# Targets
#--------

all: $(EXE)

clean:
	@$(RM) -r -v *.o

distclean: clean
	@$(RM) -r -v $(EXE) *.vtu

.DEFAULT_GOAL = all

.PHONY: all clean distclean

#------
# Rules
#------
$(EXE): %: %.o
	$(CXX) $(CXXFLAGS) $< $(LD_FLAGS) $(LD_LIBS) $(OUTPUT_OPTION)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< $(OUTPUT_OPTION)
