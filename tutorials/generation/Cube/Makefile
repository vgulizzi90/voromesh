
include config.mk

exec := cube.o

all:
	$(MAKE) -C $(SRC_DIR)
	$(MAKE) example_cube

example_cube: $(exec)

$(exec): main.cpp
	$(CXX) $(CXX_FLAGS) $(CXX_DEFS) $(I_VORO) $(L_VORO) $(I_MESH) $(L_MESH) -o $@ $< $(l_VORO) $(l_MESH)

clean:
	rm -rf $(exec)
	$(MAKE) -C $(SRC_DIR) clean

.PHONY: all example_cube clean