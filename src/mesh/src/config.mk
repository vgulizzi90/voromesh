SRC_DIR := .
BUILD_DIR := .

VORO_DIR := ../../voro_old
TRIANGLE_DIR := ../triangle
TETGEN_DIR := ../tetgen

TRIANGLE_DEFS := -DNO_TIMER -DTRILIBRARY
TETGEN_DEFS := -DTETLIBRARY

CC := gcc
CC_FLAGS := -Wall -ansi -pedantic -O3
CC_DEFS :=

CXX := g++
CXX_FLAGS := -Wall -ansi -pedantic -O3 -std=c++11
CXX_FLAGS_DBG := -Wall -ansi -pedantic -O0 -std=c++11
CXX_DEFS := 

I_VORO := -I$(VORO_DIR)/src
L_VORO := -L$(VORO_DIR)/src
l_VORO := -lvoro++