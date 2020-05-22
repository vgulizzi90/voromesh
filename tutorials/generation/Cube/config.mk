SRC_DIR := ../../../src
BUILD_DIR := ./build

CXX := g++
CXX_FLAGS := -Wall -ansi -pedantic -O3 -std=c++11
CXX_DEFS := 

I_VORO := -I../../../src/voro/src
L_VORO := -L../../../src/voro/src
l_VORO := -lvoro++

I_MESH := -I../../../src/mesh/src
L_MESH := -L../../../src/mesh/src
l_MESH := -lvoromesh++