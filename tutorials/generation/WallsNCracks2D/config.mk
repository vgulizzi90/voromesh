SRC_DIR := ../../../src
VORO_DIR := $(SRC_DIR)/voro_old
MESH_DIR := $(SRC_DIR)/mesh

I_VORO := -I$(VORO_DIR)/src
L_VORO := -L$(VORO_DIR)/src
l_VORO := -lvoro++

I_MESH := -I$(MESH_DIR)/src
L_MESH := -L$(MESH_DIR)/src
l_MESH := -lvoromesh++


BUILD_DIR := ./build

CXX := g++
CXX_FLAGS := -Wall -ansi -pedantic -Wshadow -O3 -std=c++11
CXX_FLAGS_DBG := -Wall -ansi -pedantic -Wshadow -O0 -std=c++11
CXX_DEFS := 