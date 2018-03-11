#http://make.mad-scientist.net/papers/advanced-auto-dependency-generation/#gcc

SOURCE_DIR = src
HEADER_DIR = include
BUILD_DIR = obj
TEST_DIR = test
EXE_DIR = bin
MESH_DIR = ../meshes/
DEP_DIR = .d

TEST_SRCS = $(wildcard $(TEST_DIR)/*.cpp)
TEST_OBJS = $(addprefix $(BUILD_DIR)/, $(notdir $(TEST_SRCS:.cpp=.o)))

SRCS = $(wildcard $(SOURCE_DIR)/*.cpp)
OBJS = $(addprefix $(BUILD_DIR)/, $(notdir $(SRCS:.cpp=.o)))

DEPS = $(addprefix $(DEP_DIR)/, $(notdir $(OBJS:.o=.d)))
DEPS += $(addprefix $(DEP_DIR)/, $(notdir $(TEST_OBJS:.o=.d)))

CXX = g++
CXXFLAGS += -Wall -std=c++11
CPPFLAGS += -I./$(HEADER_DIR) -I$(mkEigenInc)
# LD_FLAGS =
# LD_LIBS =
DEPFLAGS = -MT $@ -MMD -MP -MF $(DEP_DIR)/$*.Td

EXE = $(addprefix $(EXE_DIR)/, $(basename $(notdir $(TEST_SRCS))))

$(shell mkdir -p $(DEP_DIR) >/dev/null)

all: $(EXE)

.PHONY: all clean distclean run
.DEFAULT_GOAL = all
.PRECIOUS: $(DEP_DIR)/%.d

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

run: $(EXE)
	@$<

clean:
	$(RM) -r $(BUILD_DIR) $(DEP_DIR)

distclean: clean
	$(RM) -r $(EXE_DIR)

-include $(DEPS)
