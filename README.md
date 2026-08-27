# RESONANCETABLES

RESONANCETABLES is a software package to unify existing databases for thermal cross sections, resonance integrals, Maxwellian-averaged cross sections (MACS), resonance parameters and average resonance parameters. The final database is obtained through a priority selection among the original databases.

Most users can use the database supplied with this repository directly. Installing and running the code is mainly relevant when the database needs to be reconstructed.

## Documentation and reference

A description of the code and its options can be found in the [RESONANCETABLES tutorial (pdf)](https://github.com/arjankoning1/resonancetables/blob/main/doc/resonancetables.pdf).

The reference to be used for RESONANCETABLES is:

D. Rochman, A.J. Koning, J.-Ch. Sublet, *A statistical analysis of evaluated neutron resonances with TARES for JEFF-3.3, JENDL-4.0, ENDF/B-VIII.0 and TENDL-2019*, Nuclear Data Sheets 163, 163 (2020).

## Installation

### Prerequisites

The following are the prerequisites for compiling RESONANCETABLES:

- git (only if the package is downloaded via GitHub)
- GNU make
- a recent Fortran compiler, such as GNU Fortran (gfortran)

Reconstructing the complete database also requires external data installed as sibling directories of `resonancetables/`:

```text
.../resonancetables/
.../exfortables/
.../libraries/
.../tendl/
```

In particular, RESONANCETABLES uses:

```text
.../exfortables/special/
.../libraries/resbase/
.../tendl/
```

The required external databases can be obtained from the TALYS/TENDL distribution area at:

```text
https://nds.iaea.org/talys/
```

### Downloads

#### 1. Download the tar file (frozen version)

```bash
curl -LO https://nds.iaea.org/talys/resonancetables.tar
tar zxf resonancetables.tar
```

#### 2. Using git (latest beta version)

```bash
git clone https://github.com/arjankoning1/resonancetables.git
```

### Installation instructions

#### 1. For the frozen tar version

```bash
cd resonancetables
./install_resonancetables.bash
```

An alternative is:

```bash
cd resonancetables/source
make
```

The frozen distribution retains its own installation scripts and settings.

#### 2. For the git version (latest beta version)

```bash
cd resonancetables
./install_resonancetables.bash
```

which automatically executes the `Makefile` in `resonancetables/source`.

An alternative is:

```bash
cd resonancetables/source
make
```

For the git version, the default compiler is `gfortran`. When `gfortran` is used and no `FFLAGS` are supplied, the Makefile uses:

```text
-w -O3 -ffp-contract=off
```

For other compilers, no default compiler flags are imposed.

Compiler and compilation options can be passed through `install_resonancetables.bash`, for example:

```bash
./install_resonancetables.bash FC=gfortran FFLAGS="-O3 -ffp-contract=off"
./install_resonancetables.bash FC=ifx FFLAGS="-O3"
```

The executable is installed as:

```text
resonancetables/bin/resonancetables
```

Set `RESONANCETABLES_DIR` to the RESONANCETABLES installation directory. For example:

```bash
export RESONANCETABLES_DIR="/Users/koning/resonancetables"
```

If you want to run `resonancetables` from anywhere, add its `bin` directory to `PATH`:

```bash
export PATH="$RESONANCETABLES_DIR/bin:$PATH"
```

These lines can be added to `~/.zshrc` or `~/.profile`.

If setting `RESONANCETABLES_DIR` is not possible, edit `code_dir` in `source/machine.f90` and rebuild RESONANCETABLES.

No user-name environment variable is needed by RESONANCETABLES.

For the modern git version, `code_build.bash` and `path_change.bash` are no longer required and can be removed after adopting the new installer, Makefile and `machine.f90`.

## The RESONANCETABLES package

The `resonancetables/` directory contains:

- `README.md` this README file
- `LICENSE` the license file
- `install_resonancetables.bash` installation script
- `source/` the Fortran source code and Makefile
- `bin/` the executable after successful installation
- `doc/` the tutorial
- `files/` original input databases
- `libs/` results extracted from nuclear data libraries
- `thermal/` produced thermal cross-section tables
- `macs/` produced Maxwellian-averaged cross-section tables
- `resonance/` produced resonance and average-resonance-parameter tables

## Reconstructing and checking the database

A successful installation can be tested by reconstructing the database. The full run takes roughly 10 minutes.

**Important:** RESONANCETABLES removes and recreates directories named `thermal/`, `macs/` and `resonance/` in the directory from which it is run. Therefore, run it in an empty or dedicated working directory unless you deliberately want those directories replaced.

For example:

```bash
mkdir resonancetables_test
cd resonancetables_test
resonancetables > resonancetables.out
```

assuming that `resonancetables/bin` has been added to `PATH`.

From the top-level RESONANCETABLES directory, the same check can be started safely with:

```bash
make -C source check
```

`make check` verifies that the required sibling `exfortables/special/`, `libraries/resbase/` and `tendl/` directories exist, then runs RESONANCETABLES in a temporary empty directory. The temporary output is removed after the check, so the database directories shipped in the repository are not overwritten.

## License and Copyright

This software is distributed and copyrighted according to the [LICENSE](LICENSE) file.
