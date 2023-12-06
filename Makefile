# Makefile

# Absolute path of the project directory
PROJECT_DIR := $(shell pwd)

# Directories
INCLUDE_DIR := $(PROJECT_DIR)/include
SRC_DIR := $(PROJECT_DIR)/src
LIB_DIR := $(PROJECT_DIR)/lib

# Files
HEADERS := $(wildcard $(INCLUDE_DIR)/*.h)
SOURCES := $(wildcard $(SRC_DIR)/*.cpp)
OBJECTS := $(patsubst $(SRC_DIR)/%.cpp, $(LIB_DIR)/%.o, $(SOURCES))
LIBRARY := $(LIB_DIR)/libEfficiency.so

# Compilation flags
CXX := g++
CXXFLAGS := -Wall -O2 -I$(INCLUDE_DIR) `root-config --cflags --libs`

# Linking flags
LDFLAGS := -shared
LDLIBS := `root-config --libs`

# Build rules
all: $(LIBRARY)
	@rm -f $(OBJECTS)

$(LIBRARY): $(OBJECTS)
	@mkdir -p $(LIB_DIR)
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

$(LIB_DIR)/%.o: $(SRC_DIR)/%.cpp $(HEADERS)
	@mkdir -p $(LIB_DIR)
	$(CXX) $(CXXFLAGS) -c -o $@ $< -fPIC

clean:
	rm -rf $(LIB_DIR)

.PHONY: all clean
