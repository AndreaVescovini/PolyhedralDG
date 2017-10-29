SOURCE_DIR = src
HEADER_DIR = include
BUILD_DIR = obj
TEST_DIR = test
MESH_DIR = ../meshes/

TEST_SRCS = $(wildcard $(TEST_DIR)/*.cpp)
TEST_OBJS = $(addprefix $(BUILD_DIR)/, $(notdir $(TEST_SRCS:.cpp=.o)))

SRCS = $(wildcard $(SOURCE_DIR)/*.cpp)
SRCS += $(TEST_SRCS)
OBJS = $(addprefix $(BUILD_DIR)/, $(notdir $(SRCS:.cpp=.o)))
OBJS += $(TEST_OBJS)

CXX = g++
CXXFLAGS += -Wall -std=c++11
CPPFLAGS += -I./$(HEADER_DIR)
# LD_FLAGS =
# LD_LIBS =
EXE = $(basename $(notdir $(TEST_SRCS)))
DEPEND = make.dep

-include $(DEPEND)

all: $(DEPEND) $(EXE)

.PHONY: all clean distclean depend
.DEFAULT_GOAL = all

$(EXE): $(OBJS)
	$(CXX) $(CXXFLAGS) $^ $(LD_FLAGS) $(LD_LIBS) $(OUTPUT_OPTION)

$(BUILD_DIR)/%.o: $(SOURCE_DIR)/%.cpp
	@mkdir -p $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

$(BUILD_DIR)/%.o: $(TEST_DIR)/%.cpp
	@mkdir -p $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

clean:
	$(RM) -r $(BUILD_DIR)
	$(RM) $(DEPEND)

distclean: clean
	$(RM) $(EXE)

depend: $(DEPEND)

$(DEPEND): $(SRCS)
	@$(RM) $(DEPEND)
	@for ff in $(SRCS); do \
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -MM $$ff >> $(DEPEND); \
	done;
