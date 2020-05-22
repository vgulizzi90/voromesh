
include config.mk

all:
	$(MAKE) -C $(VORO_DIR)/src/ libvoro++.a
	$(MAKE) -C $(MESH_DIR)/src/ libvoromesh++.a

clean:
	$(MAKE) -C $(MESH_DIR)/src/ clean
	$(MAKE) -C $(VORO_DIR)/src/ clean

.PHONY: all clean