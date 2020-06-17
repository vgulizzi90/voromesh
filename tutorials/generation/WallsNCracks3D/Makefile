
include config.mk

exec := walls_n_cracks.o

all:
	$(MAKE) -C $(SRC_DIR)
	$(MAKE) example_walls_n_cracks

example_walls_n_cracks: $(exec)

$(exec): main.cpp $(VORO_DIR)/src/libvoro++.a $(MESH_DIR)/src/libvoromesh++.a
	$(CXX) $(CXX_FLAGS) $(CXX_DEFS) $(I_VORO) $(L_VORO) $(I_MESH) $(L_MESH) -o $@ $< $(l_VORO) $(l_MESH)

clean:
	rm -rf $(exec)
	$(MAKE) -C $(SRC_DIR) clean

.PHONY: all example_walls_n_cracks clean