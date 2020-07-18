
include config.mk

exec := cracks.o

all:
	$(MAKE) -C $(SRC_DIR)
	$(MAKE) example_cracks

example_cracks: $(exec)

$(exec): main.cpp $(VORO_DIR)/src/libvoro++.a $(MESH_DIR)/src/libvoromesh++.a
	$(CXX) $(CXX_FLAGS) $(CXX_DEFS) $(I_VORO) $(L_VORO) $(I_MESH) $(L_MESH) -o $@ $< $(l_VORO) $(l_MESH)

clean:
	rm -rf $(exec)
	$(MAKE) -C $(SRC_DIR) clean

.PHONY: all example_cracks clean