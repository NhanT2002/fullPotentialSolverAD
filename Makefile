# Compiler
CHPL ?= chpl
CHPL_LLVM_BIN ?= /home/user/chapel-2.5.0/third-party/llvm/install/linux64-x86_64/bin
OPT ?= $(CHPL_LLVM_BIN)/opt
LLVM_AS ?= $(CHPL_LLVM_BIN)/llvm-as

# Directories
SRC_DIR := src
BIN_DIR := bin
BUILD_DIR := build

# Program name (without .chpl)
PROGRAM := main

# Source and target
SRC := $(wildcard $(SRC_DIR)/*.chpl)
TARGET := $(BIN_DIR)/$(PROGRAM)
CHPL_BC := $(BUILD_DIR)/chpl__module-opt1.bc
CHPL_LL := $(BUILD_DIR)/chpl__module-opt1.ll

# Default rule
all: $(TARGET)

# Module
CGNS_MOD_DIR := /home/user/test/fullPotentialSolver/champs/EXT_LIBS
CGNS_MOD_DIR_SRC := /home/user/test/fullPotentialSolver/champs/EXT_LIBS/src
COMMON_MOD_DIR := /home/user/test/fullPotentialSolver/champs/common/src

# PETSc was built against its own MPICH install. If PETSc exposes an install
# prefix, use that for MPI instead of a conflicting system MPIROOT.
PETSC_PREFIX := $(strip $(shell sed -n 's/^prefix=//p' $(PETSCROOT)/lib/pkgconfig/PETSc.pc 2>/dev/null))
ifneq ($(PETSC_PREFIX),)
MPIROOT := $(PETSC_PREFIX)
endif

MPI_LIB_DIR := $(if $(wildcard $(MPIROOT)/lib/release),$(MPIROOT)/lib/release,$(MPIROOT)/lib)
INCLUDE_DIRS := $(HDF5ROOT)/include $(MKLROOT)/include $(CGNSROOT)/include $(METISROOT)/include $(PETSCROOT)/include $(PETSC_DIR)/include $(SLEPCROOT_INCLUDE) $(SLEPCDIR_INCLUDE) $(MPIROOT)/include


HDF5_LIB := -L$(HDF5ROOT)/lib -lhdf5
CGNS_LIB := -L$(CGNSROOT)/lib -lcgns -lhdf5
MKL_LIB := -L$(MKLROOT)/lib/intel64 -lmkl_rt -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5
METIS_LIB := -L$(METISROOT)/lib -lmetis
MPI_LIB := -L$(MPI_LIB_DIR) -lmpi
PETSC_LIB := -L$(PETSCROOT)/lib -lpetsc

HDF5_INCLUDE := $(HDF5ROOT)/include hdf5.h
MKL_INCLUDE := $(MKLROOT)/include mkl_types.h mkl.h
CGNS_INCLUDE := $(CGNSROOT)/include cgnslib.h
METIS_INCLUDE := $(METISROOT)/include metis.h
MPI_INCLUDE := $(MPIROOT)/include mpi.h
PETSCROOT_INCLUDE := $(PETSCROOT)/include petscksp.h petsc.h petscmat.h
PETSCDIR_INCLUDE := $(PETSC_DIR)/include petscksp.h petsc.h petscmat.h

ALL_INCLUDES := $(addprefix -I,$(filter-out ,$(INCLUDE_DIRS)))
ALL_LIBS := $(MKL_LIB) $(HDF5_LIB) $(CGNS_LIB) $(METIS_LIB) $(PETSC_LIB) $(SLEPC_LIB) $(MPI_LIB)
CHPL_FLAGS := -M$(CGNS_MOD_DIR) -M$(CGNS_MOD_DIR_SRC) -M$(COMMON_MOD_DIR) $(ALL_INCLUDES) $(ALL_LIBS) --set blasImpl=off --local --fast

ENZYME_DIR ?= /home/user/Enzyme/enzyme/build-llvm19/Enzyme
ENZYME_PLUGIN ?= $(firstword $(sort $(wildcard $(ENZYME_DIR)/LLVMEnzyme*.so)))

# Compile Chapel program
$(TARGET): $(SRC) | $(BIN_DIR) $(BUILD_DIR)
	test -n "$(ENZYME_PLUGIN)"
	test -f "$(ENZYME_PLUGIN)"
	$(CHPL) $(CHPL_FLAGS) --driver-compilation-phase --driver-tmp-dir $(BUILD_DIR) $(SRC)
	$(OPT) $(CHPL_BC) -load-pass-plugin=$(ENZYME_PLUGIN) -passes=enzyme -o $(CHPL_LL) -S
	$(OPT) $(CHPL_LL) -O2 -o $(CHPL_LL) -S
	$(LLVM_AS) $(CHPL_LL) -o $(CHPL_BC)
	$(CHPL) $(CHPL_FLAGS) --driver-makebinary-phase --driver-tmp-dir $(BUILD_DIR) $(SRC) -o $(TARGET)

# Create bin directory if it does not exist
$(BIN_DIR):
	mkdir -p $(BIN_DIR)

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# Clean target
clean:
	rm -rf $(TARGET) $(BUILD_DIR)

# Phony targets
.PHONY: all clean
