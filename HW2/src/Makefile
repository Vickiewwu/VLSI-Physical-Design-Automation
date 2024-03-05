CXX = g++   
CXXFLAGS = -Wall -std=c++11 -Ofast 

SRC_DIR = ../src
BIN_DIR = ../bin

SRC = $(SRC_DIR)/hw2.cpp
TARGET = $(BIN_DIR)/hw2
HYPERFILES = $(BIN_DIR)/hyper.hgr $(BIN_DIR)/hyper.hgr.part.2

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $@ $<
.PHONY: all
all: $(TARGET)

.PHONY: clean
clean:
	rm -f $(TARGET) $(HYPERFILES)

	
