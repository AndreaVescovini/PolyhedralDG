#http://make.mad-scientist.net/papers/advanced-auto-dependency-generation/#gcc

#------------------------------
# Including the global Makefile
#------------------------------
THAT_DIR = ..
include $(THAT_DIR)/Makefile.inc

#-------------------------------------------------------------------------------
#------------------- DO NOT MODIFY ANYTHING BELOW THIS LINE --------------------
#-------------------------------------------------------------------------------

#------------
# Directories
#------------
SOURCE_DIR  	 = src
INCLUDE_DIR 	 = include
INCLUDE_SUBDIR = expr
BUILD_DIR 		 = obj
LIB_DIR 			 = lib
TEST_DIR  	   = test
EXE_DIR 		   = bin
DEP_DIR 		   = .d
DOC_DIR 		   = doc

BUILD_DIR_STATIC 	= $(BUILD_DIR)/static
BUILD_DIR_DYNAMIC = $(BUILD_DIR)/dynamic
BUILD_DIR_TEST 		= $(BUILD_DIR)/test

#---------------------------------------------
# Tests sources, objects files and executables
#---------------------------------------------
TEST_SRCS = $(wildcard $(TEST_DIR)/*.cpp)
TEST_OBJS = $(addprefix $(BUILD_DIR_TEST)/, $(notdir $(TEST_SRCS:.cpp=.o)))
EXE 			= $(addprefix $(EXE_DIR)/, $(basename $(notdir $(TEST_SRCS))))

#----------------------------------
# Library sources, headers and objects files
#----------------------------------
SRCS 				 = $(wildcard $(SOURCE_DIR)/*.cpp)
HEADERS      = $(wildcard $(INCLUDE_DIR)/*.hpp)
HEADERS     += $(wildcard $(INCLUDE_DIR)/$(INCLUDE_SUBDIR)/*.hpp)
STATIC_OBJS  = $(addprefix $(BUILD_DIR_STATIC)/, $(notdir $(SRCS:.cpp=.o)))
DYNAMIC_OBJS = $(addprefix $(BUILD_DIR_DYNAMIC)/, $(notdir $(SRCS:.cpp=.o)))

#----------
# Libraries
#----------
STATIC_LIB  = $(LIB_DIR)/libPolyDG.a
DYNAMIC_LIB = $(LIB_DIR)/libPolyDG.so

#-----------------------------------------------------------------------
# Files containing the debenpencies of objects files, they are generated
# compiling with DEPFLAGS
#-----------------------------------------------------------------------
DEPS  = $(addprefix $(DEP_DIR)/, $(notdir $(SRCS:.cpp=.d)))
DEPS += $(addprefix $(DEP_DIR)/, $(notdir $(TEST_SRCS:.cpp=.d)))

#------------------
# Compilation flags
#------------------
CPPFLAGS += -I./$(INCLUDE_DIR)
DEPFLAGS  = -MT $@ -MMD -MP -MF $(DEP_DIR)/$*.Td

LD_FLAGS += -L./$(LIB_DIR)

#---------------------------
# Assure that DEP_DIR exists
#---------------------------
$(shell mkdir -p $(DEP_DIR))

#--------
# Targets
#--------

help:
	@echo "--------------------------------- LIBRARY POLYDG ---------------------------------------"
	@echo "make help: ----- Prints this help"
	@echo "make all: ------ Makes the libraries, the executables of the tests and the documentation"
	@echo "make static: --- Makes the static library"
	@echo "make dynamic: -- Makes the dynamic library"
	@echo "make library: -- Makes both static and dynamic libraries"
	@echo "make test: ----- Makes the executables of the tests"
	@echo "make install: -- Installs the libraries"
	@echo "make uninstall - Uninstall the libraries"
	@echo "make doc: ------ Makes the documentation"
	@echo "make clean: ---- Cleans all the temporary files"
	@echo "make distclean:  Cleans all"
	@echo ""
	@echo "REMEMBER: Specify always RELEASE=yes in order to enable optimizations"
	@echo "-----------------------------------------------------------------------------------------"

all: static dynamic test doc

static: $(STATIC_LIB)

dynamic: $(DYNAMIC_LIB)

dynamic: CXXFLAGS += -fPIC

library: static dynamic

test: $(EXE)

test: CPPFLAGS += -I$(GETPOT_INC)

install: clean
	@mkdir -p $(POLYDG_INC)
	cp $(INCLUDE_DIR)/* -r $(POLYDG_INC)
	@mkdir -p $(POLYDG_LIB)
	cp $(STATIC_LIB) $(DYNAMIC_LIB) $(POLYDG_LIB)

uninstall:
	@$(RM) -r -v $(POLYDG_PATH)

doc: $(HEADERS) $(SRCS) $(TEST_SRCS) mainpage.dox Doxyfile
	@mkdir -p $(DOC_DIR)
	doxygen Doxyfile

clean:
	@$(RM) -r -v $(BUILD_DIR) $(DEP_DIR)

distclean: clean
	@$(RM) -r -v $(EXE_DIR) $(LIB_DIR) $(DOC_DIR)

.DEFAULT_GOAL = help

#----------------
# Special targets
#----------------
.PHONY:		  help all static dynamic library test install uninstall doc clean distclean
.PRECIOUS:  $(DEP_DIR)/%.d
.SECONDARY: $(STATIC_OBJS) $(DYNAMIC_OBJS) $(TEST_OBJS) $(STATIC_LIB) $(DYNAMIC_LIB)

#------
# Rules
#------
$(STATIC_LIB): $(STATIC_OBJS)
	@mkdir -p $(LIB_DIR)
	$(AR) -r -s $@ $^
	ranlib $@

$(DYNAMIC_LIB): $(DYNAMIC_OBJS)
	@mkdir -p $(LIB_DIR)
	$(CXX) -shared $^ $(OUTPUT_OPTION)

$(BUILD_DIR_STATIC)/%.o: $(SOURCE_DIR)/%.cpp $(DEP_DIR)/%.d
	@mkdir -p $(BUILD_DIR_STATIC)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(DEPFLAGS) -c $< $(OUTPUT_OPTION)
	@mv -f $(DEP_DIR)/$*.Td $(DEP_DIR)/$*.d && touch $(word 2,$^) && touch $@

$(BUILD_DIR_DYNAMIC)/%.o: $(SOURCE_DIR)/%.cpp $(DEP_DIR)/%.d
	@mkdir -p $(BUILD_DIR_DYNAMIC)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(DEPFLAGS) -c $< $(OUTPUT_OPTION)
	@mv -f $(DEP_DIR)/$*.Td $(DEP_DIR)/$*.d && touch $(word 2,$^) && touch $@

$(EXE_DIR)/%: $(BUILD_DIR_TEST)/%.o
	@mkdir -p $(EXE_DIR)
	$(CXX) $(CXXFLAGS) $< $(LD_FLAGS) $(LD_LIBS) $(OUTPUT_OPTION)

$(BUILD_DIR_TEST)/%.o: $(TEST_DIR)/%.cpp
	@mkdir -p $(BUILD_DIR_TEST)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< $(OUTPUT_OPTION)

$(DEP_DIR)/%.d: ;

-include $(DEPS)
