#Compiler
# we dont want other compiler to be used yet.
CXX = g++

#Flags
CXX_ASSEMBLER_FLAGS := -std=c++11 $(ADD_G++_FLAGS) 

# INCLUDE_FLAGS = -I src/
# SHARED_LINK_FLAGS = 

DEBUG ?= 0
ifeq ($(DEBUG), 1)
    CXX_ASSEMBLER_FLAGS +=-g
	GCC_ASSEMBLER_FLAGS +=-g
endif


# Directories
SRC_DIR = src
OBJ_DIR = build
BIN_DIR = bin

OUTPUT_DIR = $(OBJ_DIR) $(BIN_DIR)

.PHONY: all
all:
	make directories
	make $(BIN_DIR)/seq $(BIN_DIR)/pt $(BIN_DIR)/omp

.PHONY: directories
directories:
	mkdir -p $(OUTPUT_DIR)
	touch $(BIN_DIR)/clusters.txt
	touch $(BIN_DIR)/centorids.txt

# BIN/Executable Rules
# seq exec
$(BIN_DIR)/seq: $(OBJ_DIR)/main_sequential.o $(OBJ_DIR)/lab1_sequential.o $(OBJ_DIR)/lab1_io.o
	$(CXX) $^ -o $@ $(SHARED_LINK_FLAGS) -lomp -lpthread

# pt exec
$(BIN_DIR)/pt: $(OBJ_DIR)/main_pthread.o $(OBJ_DIR)/lab1_pthread.o $(OBJ_DIR)/lab1_io.o
	$(CXX) $^ -o $@ $(SHARED_LINK_FLAGS) -lomp -lpthread

# omp exec
$(BIN_DIR)/omp: $(OBJ_DIR)/main_omp.o $(OBJ_DIR)/lab1_omp.o $(OBJ_DIR)/lab1_io.o
	$(CXX) $^ -o $@ $(SHARED_LINK_FLAGS) -lomp -lpthread

# OBJ/object Rules
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXX_ASSEMBLER_FLAGS) -c $^ -o $@ $(INCLUDE_FLAGS)

# omp special compile
$(OBJ_DIR)/lab1_omp.o: $(SRC_DIR)/lab1_omp.cpp
	$(CXX) -fopenmp $(CXX_ASSEMBLER_FLAGS) -c $^ -o $@ $(INCLUDE_FLAGS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	gcc -fopenmp $(GCC_ASSEMBLER_FLAGS) -c $^ -o $@ $(INCLUDE_FLAGS)

.PHONY: clean
clean:
	rm -rf bin/ build/