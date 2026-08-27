#!/usr/bin/env bash

set -euo pipefail

# Determine the RESONANCETABLES installation directory independently of
# where the script is called from.

resonancetables_dir=$(cd "$(dirname "$0")" && pwd)
source_dir="$resonancetables_dir/source"

# Verify that the expected directories and build files exist.

if [[ ! -d "$source_dir" ]]; then
  echo "RESONANCETABLES installation error: source directory not found:" >&2
  echo "  $source_dir" >&2
  exit 1
fi

if [[ ! -f "$source_dir/Makefile" ]]; then
  echo "RESONANCETABLES installation error: Makefile not found:" >&2
  echo "  $source_dir/Makefile" >&2
  exit 1
fi

data_file="$resonancetables_dir/files/atlas_2006.txt"

if [[ ! -f "$data_file" ]]; then
  echo "RESONANCETABLES installation error: data files missing or incomplete:" >&2
  echo "  $data_file" >&2
  exit 1
fi

echo
echo "Installing RESONANCETABLES"
echo "Installation directory: $resonancetables_dir"
echo

# Pass all command-line arguments directly to make. This permits, e.g.:
#
# ./install_resonancetables.bash FC=ifx FFLAGS="-O3"
# ./install_resonancetables.bash FC=gfortran FFLAGS="-w -O3 -ffp-contract=off"

make -C "$source_dir" clean
make -C "$source_dir" all "$@"

resonancetables_exe="$resonancetables_dir/bin/resonancetables"

if [[ ! -x "$resonancetables_exe" ]]; then
  echo "RESONANCETABLES installation error: executable not created:" >&2
  echo "  $resonancetables_exe" >&2
  exit 1
fi

echo
echo "RESONANCETABLES executable:"
echo "  $resonancetables_exe"
echo
echo "If not already done, add the following lines to your shell configuration:"
echo
echo "  export RESONANCETABLES_DIR=\"$resonancetables_dir\""
echo "  export PATH=\"\$RESONANCETABLES_DIR/bin:\$PATH\""
echo
echo "For a full database reconstruction, sibling exfortables, libraries"
echo "and tendl directories are also required."
echo
echo "Alternatively, edit code_dir in source/machine.f90 and rebuild RESONANCETABLES."
echo
