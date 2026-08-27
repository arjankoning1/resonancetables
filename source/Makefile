# RESONANCETABLES Makefile
#
# Examples:
#
#   make
#   make FC=gfortran FFLAGS="-w -O3 -ffp-contract=off"
#   make FC=ifx FFLAGS="-O3"
#   make clean
#   make check

# Compiler-specific flags can be supplied by the user.
ifeq ($(origin FC), default)
  FC = gfortran
endif

ifeq ($(origin FFLAGS), undefined)
  ifneq ($(findstring gfortran,$(FC)),)
    FFLAGS = -w -O3 -ffp-contract=off
  else
    FFLAGS =
  endif
endif

# Installation directories.

RESONANCETABLES_DIR := $(abspath ..)
BIN_DIR             := $(RESONANCETABLES_DIR)/bin
TARGET              := $(BIN_DIR)/resonancetables

# Sibling data directories required for a complete reconstruction.

EXFOR_DIR_CHECK   := $(abspath $(RESONANCETABLES_DIR)/../exfortables/special)
RESBASE_DIR_CHECK := $(abspath $(RESONANCETABLES_DIR)/../libraries/resbase)
TENDL_DIR_CHECK   := $(abspath $(RESONANCETABLES_DIR)/../tendl)

# RESONANCETABLES source files.

SOURCES := $(sort $(wildcard *.f) $(wildcard *.f90))

.PHONY: all install clean check

all install: $(TARGET)

$(TARGET): $(SOURCES)
	@echo "Compiler: $(FC)"
	@echo "Flags:    $(FFLAGS)"
	@mkdir -p "$(BIN_DIR)"
	$(FC) $(FFLAGS) -o "$(TARGET)" $(SOURCES)
	@rm -f *.o *.mod
	@echo "Installed: $(TARGET)"

clean:
	rm -f *.o *.mod resonancetables "$(TARGET)"

check: all
	@test -d "$(EXFOR_DIR_CHECK)" || { \
	  echo "RESONANCETABLES check error: directory not found:" >&2; \
	  echo "  $(EXFOR_DIR_CHECK)" >&2; \
	  exit 1; \
	}
	@test -d "$(RESBASE_DIR_CHECK)" || { \
	  echo "RESONANCETABLES check error: directory not found:" >&2; \
	  echo "  $(RESBASE_DIR_CHECK)" >&2; \
	  exit 1; \
	}
	@test -d "$(TENDL_DIR_CHECK)" || { \
	  echo "RESONANCETABLES check error: directory not found:" >&2; \
	  echo "  $(TENDL_DIR_CHECK)" >&2; \
	  exit 1; \
	}
	@set -e; \
	  tmpdir=$$(mktemp -d); \
	  trap 'rm -rf "$$tmpdir"' EXIT; \
	  echo "Running RESONANCETABLES check in $$tmpdir"; \
	  cd "$$tmpdir"; \
	  RESONANCETABLES_DIR="$(RESONANCETABLES_DIR)" "$(TARGET)" > resonancetables.out; \
	  echo "RESONANCETABLES check completed successfully"
