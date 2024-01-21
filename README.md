
# RESONANCETABLES
RESONANCETABLES is a software package to unify databases for thermal cross sections, resonance integrals, MACS and other average resonance parameters. The final database is the result of a priority setting of various original databases.
The provided database can be used directly, or the code can be installed for user-specified reproduction of the database.

## Documentation and reference
A description of the code and its options can be found in the [RESONANCETABLES tutorial (pdf)](https://github.com/arjankoning1/resonancetables/blob/main/doc/resonancetables.pdf).
The reference to be used for RESONANCETABLES is

D. Rochman, A.J. Koning, J.-Ch. Sublet, A statistical analysis of evaluated neutron resonances with TARES for JEFF-3.3, JENDL-4.0, ENDF/B-VIII.0 and TENDL-2019, Nuclear Data Sheets 163, 163 (2020).

## Installation

### Prerequisites:

The following are the prerequisites for compiling RESONANCETABLES:
  - git (only if the package is downloaded via Github)
  - a recent Fortran compiler such as gcc (gfortran)

### Downloads:

To download RESONANCEtables, you can use one of the following options:
#### 1. Download the tar file:
```
https://nds.iaea.org/talys/resonancetables.tar
tar zxf resonancetables.tar
```

#### 2. Using git:
```
git clone https://github.com/arjankoning1/resonancetables.git
```
### Installation instructions:

To install RESONANCETABLES, you can use one of the following options:
#### 1. Using make:
```
git clone https://github.com/arjankoning1/resonancetables.git
cd resonancetables/source
make
```
#### 2. Using the code_build script:
```
git clone https://github.com/arjankoning1/resonancetables.git
cd resonancetables
code_build resonancetables
```

The above instructions will produce a *resonancetables* executable in the *resonancetables/bin* directory. 
The compiler and its flags can be set in either the ource/Makefile* or in *code_build*.

## The RESONANCETABLES package

The *resonancetables/* directory contains the following directories and files:

+ `README.md` is this README file
+ `LICENSE` is the License file
+ `code_build` and `path_change` installation scripts
+ `source/` contains the Fortran source code of RESONANCETABLES and the Makefile
+ `bin/` contains the executable after successful installation
+ `doc/` contains the tutorial in pdf format
+ `files/` contains the input files for RESONANCETABLES: the original databases
+ `libs/` contains the results from the available nuclear data libraries
+ `macs/` contains the produced Maxwellian-averaged cross section (MACS) tables
+ `thermal/` contains the produced thermal cross section tables
+ `resonance/` contains the produced average resonance parameter tables

In total, you will need about 20 Mb of free disk space to install RESONANCETABLES.

## Sample cases

A successful installation can be verified by simply running *resonancetables* in an empty directory. This will take about 30 seconds.

## License and Copyright
This software is distributed and copyrighted according to the [LICENSE](LICENSE) file.
