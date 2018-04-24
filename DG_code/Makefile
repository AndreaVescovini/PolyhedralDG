#http://make.mad-scientist.net/papers/advanced-auto-dependency-generation/#gcc

#------------
# Directories
#------------
SOURCE_DIR = src
HEADER_DIR = include
BUILD_DIR = obj
TEST_DIR = test
EXE_DIR = bin
DEP_DIR = .d

#----------------------------------------------------
# Sources, objects files and executables of the tests
#- --------------------------------------------------
TEST_SRCS = $(wildcard $(TEST_DIR)/*.cpp)
TEST_OBJS = $(addprefix $(BUILD_DIR)/, $(notdir $(TEST_SRCS:.cpp=.o)))
EXE = $(addprefix $(EXE_DIR)/, $(basename $(notdir $(TEST_SRCS))))

#--------------------------------
# Other sources and objects files
#--------------------------------
SRCS = $(wildcard $(SOURCE_DIR)/*.cpp)
OBJS = $(addprefix $(BUILD_DIR)/, $(notdir $(SRCS:.cpp=.o)))

#-----------------------------------------------------------------------
# Files containing the debendencies of objects files, they are generated
# compiling with DEPFLAGS
#-----------------------------------------------------------------------
DEPS = $(addprefix $(DEP_DIR)/, $(notdir $(OBJS:.o=.d)))
DEPS += $(addprefix $(DEP_DIR)/, $(notdir $(TEST_OBJS:.o=.d)))

#---------
# Compiler
#---------
CXX = g++

#------------------
# Compilation flags
#------------------
WFLAGS += -Wall -Wextra -pedantic
STDFLAG += -std=c++11

ifdef RELEASE
OPTFLAGS += -O3 -DNDEBUG # -pg #-DINLINE
else
OPTFLAGS += -g #-DVERBOSITY
endif

CXXFLAGS += $(WFLAGS) $(STDFLAG) $(OPTFLAGS)
CPPFLAGS += -I./$(HEADER_DIR) -I$(mkEigenInc)
DEPFLAGS = -MT $@ -MMD -MP -MF $(DEP_DIR)/$*.Td

# LD_FLAGS =
# LD_LIBS =

#---------------------------
# Assure that DEP_DIR exists
#---------------------------
$(shell mkdir -p $(DEP_DIR))

all: $(EXE)

.DEFAULT_GOAL = all

#----------------
# Special targets
#----------------
.PHONY: all clean distclean run
.PRECIOUS: $(DEP_DIR)/%.d
.SECONDARY: $(OBJS) $(TEST_OBJS)

#------
# Rules
#------
$(EXE_DIR)/%: $(BUILD_DIR)/%.o $(OBJS)
	@mkdir -p $(EXE_DIR)
	$(CXX) $(CXXFLAGS) $^ $(LD_FLAGS) $(LD_LIBS) $(OUTPUT_OPTION)

$(BUILD_DIR)/%.o: $(SOURCE_DIR)/%.cpp $(DEP_DIR)/%.d
	@mkdir -p $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(DEPFLAGS) -c $< $(OUTPUT_OPTION)
	@mv -f $(DEP_DIR)/$*.Td $(DEP_DIR)/$*.d && touch $(word 2,$^) && touch $@

$(BUILD_DIR)/%.o: $(TEST_DIR)/%.cpp $(DEP_DIR)/%.d
	@mkdir -p $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(DEPFLAGS) -c $< $(OUTPUT_OPTION)
	@mv -f $(DEP_DIR)/$*.Td $(DEP_DIR)/$*.d && touch $(word 2,$^) && touch $@

$(DEP_DIR)/%.d: ;

run: $(word 2, $(MAKECMDGOALS))
	@$<

clean:
	$(RM) -r $(BUILD_DIR) $(DEP_DIR)

distclean: clean
	$(RM) -r $(EXE_DIR)

-include $(DEPS)
