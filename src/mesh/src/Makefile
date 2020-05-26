
include config.mk

c_src := $(TRIANGLE_DIR)/triangle.c \
	$(TRIANGLE_DIR)/tricall.c
c_obj := $(patsubst %.c,%.o,$(c_src))

cxx_src := $(TETGEN_DIR)/predicates.cxx \
	$(TETGEN_DIR)/tetgen.cxx
cxx_obj := $(patsubst %.cxx,%.o,$(cxx_src))

cpp_src := utils_io.cpp \
	utils_vtk.cpp \
	utils_math.cpp \
	utils_mesh.cpp \
	mesh.cpp
cpp_obj := $(patsubst %.cpp,%.o,$(cpp_src))

obj := $(addprefix $(TRIANGLE_DIR)/,$(c_obj)) $(addprefix $(TETGEN_DIR)/,$(cxx_obj)) $(addprefix $(BUILD_DIR)/,$(cpp_obj))

all:
	$(MAKE) -C $(VORO_DIR)/src/ libvoro++.a
	$(MAKE) libvoromesh++.a

$(TRIANGLE_DIR)/%.o: $(TRIANGLE_DIR)/%.c
	$(CC) $(CC_FLAGS) $(CC_DEFS) $(TRIANGLE_DEFS) -c $< -o $@

$(TETGEN_DIR)/%.o: $(TETGEN_DIR)/%.cxx
	$(CXX) $(CXX_FLAGS) $(CXX_DEFS) $(TETGEN_DEFS) -c $< -o $@

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp $(VORO_DIR)/src/libvoro++.a
	$(CXX) $(CXX_FLAGS) $(CXX_DEFS) $(I_VORO) $(L_VORO) -c $< $(l_VORO) -o $@

libvoromesh++.a: $(obj)
	rm -f libvoromesh++.a
	ar rs libvoromesh++.a $^

clean:
	rm -f $(obj) libvoromesh++.a
	$(MAKE) -C $(VORO_DIR)/src/ clean

.PHONY: all clean