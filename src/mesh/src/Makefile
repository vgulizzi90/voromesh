
include config.mk

cpp_src := utils_io.cpp \
	utils_vtk.cpp \
	geometry.cpp
cpp_obj := $(patsubst %.cpp,%.o,$(cpp_src))

obj := $(addprefix $(BUILD_DIR)/,$(cpp_obj))

all:
	$(MAKE) -C $(VORO_DIR)/src/ libvoro++.a
	$(MAKE) libvoromesh++.a

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXX_FLAGS) $(CXX_DEFS) $(I_VORO) $(L_VORO) -c $< $(l_VORO) -o $@

libvoromesh++.a: $(obj)
	rm -f libvoromesh++.a
	ar rs libvoromesh++.a $^

clean:
	rm -f $(obj) libvoromesh++.a
	$(MAKE) -C $(VORO_DIR)/src/ clean

.PHONY: all clean