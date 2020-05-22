SRC_DIR := .
BUILD_DIR := .

VORO_DIR := ../../voro

CXX := g++
CXX_FLAGS := -Wall -ansi -pedantic -O3 -std=c++11
CXX_DEFS := 

I_VORO := -I$(VORO_DIR)/src
L_VORO := -L$(VORO_DIR)/src
l_VORO := -lvoro++