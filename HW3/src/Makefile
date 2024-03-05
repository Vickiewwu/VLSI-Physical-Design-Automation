CXX = g++   
CXXFLAGS = -Wall -std=c++11 -Ofast 

SRC_DIR = ../src
BIN_DIR = ../bin

SRC_FILES = $(SRC_DIR)/main.cpp $(wildcard $(SRC_DIR)/*.cpp)
OBJ_FILES = $(patsubst $(SRC_DIR)/%.cpp,$(BIN_DIR)/%.o,$(SRC_FILES))
TARGET = $(BIN_DIR)/hw3

$(TARGET): $(OBJ_FILES)
	$(CXX) $(CXXFLAGS) -o $@ $^

$(BIN_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

.PHONY: all
all: $(TARGET)

.PHONY: clean
clean:
	rm -f $(TARGET) $(OBJ_FILES)

	
