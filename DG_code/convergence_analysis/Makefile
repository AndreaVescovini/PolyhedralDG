#------------------------------
# Including the global Makefile
#------------------------------
THAT_DIR = ..
include $(THAT_DIR)/Makefile.inc

#------------------------
# Sources and executables
#------------------------
SRCS = $(wildcard *.cpp)
EXE  = $(basename $(SRCS))

#------------------
# Compilation flags
#------------------
CPPFLAGS += -I$(GETPOT_INC) -I$(POLYDG_INC)
LD_FLAGS += -L$(POLYDG_LIB)

#--------
# Targets
#--------

help:
	@echo "---------------------- CONVERGENCE ANALISYS -------------------------"
	@echo "make help: ----- Prints this help"
	@echo "make all: ------ Makes the executables"
	@echo "make clean: ---- Cleans all the temporary files"
	@echo "make distclean:  Cleans all"
	@echo ""
	@echo "REMEMBER: Specify always RELEASE=yes in order to enable optimizations"
	@echo "---------------------------------------------------------------------"


all: $(EXE)

clean:
	@$(RM) -r -v *.o

distclean: clean
	@$(RM) -r -v $(EXE) *.vtu

.DEFAULT_GOAL = help

.PHONY: all clean distclean help

#------
# Rules
#------
$(EXE): %: %.o
	$(CXX) $(CXXFLAGS) $< $(LD_FLAGS) $(LD_LIBS) $(OUTPUT_OPTION)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< $(OUTPUT_OPTION)
