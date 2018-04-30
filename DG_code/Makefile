#http://make.mad-scientist.net/papers/advanced-auto-dependency-generation/#gcc

#------------
# Directories
#------------
SOURCE_DIR = src
INCLUDE_DIR = include
BUILD_DIR = obj
LIB_DIR = lib
TEST_DIR = test
EXE_DIR = bin
DEP_DIR = .d

BUILD_DIR_STATIC = $(BUILD_DIR)/static
BUILD_DIR_DYNAMIC = $(BUILD_DIR)/dynamic
BUILD_DIR_TEST = $(BUILD_DIR)/test

#---------------------------------------------
# Tests sources, objects files and executables
#- -------------------------------------------
TEST_SRCS = $(wildcard $(TEST_DIR)/*.cpp)
TEST_OBJS = $(addprefix $(BUILD_DIR_TEST)/, $(notdir $(TEST_SRCS:.cpp=.o)))
EXE = $(addprefix $(EXE_DIR)/, $(basename $(notdir $(TEST_SRCS))))

#----------------------------------
# Library sources and objects files
#----------------------------------
SRCS = $(wildcard $(SOURCE_DIR)/*.cpp)
STATIC_OBJS = $(addprefix $(BUILD_DIR_STATIC)/, $(notdir $(SRCS:.cpp=.o)))
DYNAMIC_OBJS = $(addprefix $(BUILD_DIR_DYNAMIC)/, $(notdir $(SRCS:.cpp=.o)))

#----------
# Libraries
#----------
LIB_NAME = PolyDG
STATIC_LIB = $(LIB_DIR)/lib$(LIB_NAME).a
DYNAMIC_LIB = $(LIB_DIR)/lib$(LIB_NAME).so

#-----------------------------------------------------------------------
# Files containing the debendencies of objects files, they are generated
# compiling with DEPFLAGS
#-----------------------------------------------------------------------
DEPS = $(addprefix $(DEP_DIR)/, $(notdir $(SRCS:.cpp=.d)))
DEPS += $(addprefix $(DEP_DIR)/, $(notdir $(TEST_SRCS:.cpp=.d)))

#-------------------------------
# Compiler and compilation flags
#-------------------------------
CXX = g++

WFLAGS += -Wall
STDFLAG += -std=c++11

ifdef RELEASE
OPTFLAGS += -O3 -DNDEBUG -pg #-DINLINE
else
OPTFLAGS += -g #-DVERBOSITY
WFLAGS += -Wextra -pedantic
endif

CXXFLAGS += $(WFLAGS) $(STDFLAG) $(OPTFLAGS)
CPPFLAGS += -I./$(INCLUDE_DIR) -I$(mkEigenInc)
DEPFLAGS = -MT $@ -MMD -MP -MF $(DEP_DIR)/$*.Td

LD_FLAGS += -L./$(LIB_DIR)
LD_LIBS += -l$(LIB_NAME)

#---------------------------
# Assure that DEP_DIR exists
#---------------------------
$(shell mkdir -p $(DEP_DIR))

#--------
# Targets
#--------

all: static tests

static: $(STATIC_LIB)

dynamic: $(DYNAMIC_LIB)

tests: $(EXE)

# run: $(word 2, $(MAKECMDGOALS))
# 	@$<

clean:
	@$(RM) -r -v $(BUILD_DIR) $(DEP_DIR)

distclean: clean
	@$(RM) -r -v $(EXE_DIR) $(LIB_DIR)

.DEFAULT_GOAL = all

#----------------
# Special targets
#----------------
.PHONY: all clean distclean run
.PRECIOUS: $(DEP_DIR)/%.d
.SECONDARY: $(STATIC_OBJS) $(DYNAMIC_OBJS) $(TEST_OBJS) $(STATIC_LIB) $(DYNAMIC_LIB)

#------
# Rules
#------
$(EXE_DIR)/%: $(BUILD_DIR_TEST)/%.o $(STATIC_LIB)
	@mkdir -p $(EXE_DIR)
	$(CXX) $(CXXFLAGS) $< $(LD_FLAGS) $(LD_LIBS) $(OUTPUT_OPTION)

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
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -fPIC $(DEPFLAGS) -c $< $(OUTPUT_OPTION)
	@mv -f $(DEP_DIR)/$*.Td $(DEP_DIR)/$*.d && touch $(word 2,$^) && touch $@

$(BUILD_DIR_TEST)/%.o: $(TEST_DIR)/%.cpp $(DEP_DIR)/%.d
	@mkdir -p $(BUILD_DIR_TEST)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(DEPFLAGS) -c $< $(OUTPUT_OPTION)
	@mv -f $(DEP_DIR)/$*.Td $(DEP_DIR)/$*.d && touch $(word 2,$^) && touch $@

$(DEP_DIR)/%.d: ;

-include $(DEPS)
