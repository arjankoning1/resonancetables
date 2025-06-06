# This file is part of RESONANCETABLES
# 
#    Author: Arjan Koning Copyright: (C) 2024
#
# RESONANCETABLES is free software; you can redistribute it and/or modify it
# under the terms of the MIT License; see LICENSE file for more details.
#
# Basic installation (verified with the sample cases)
#
# FC="gfortran " FFLAGS=" "
#
# Distribution FC & FFLAGS (options provided by J-C Sublet)
#
# FC="gfortran " FFLAGS=" -Ofast "
# FC="ifort    " FFLAGS=" -Ofast "
# FC="nagfor   " FFLAGS=" -w     "
#
# Development FC & FFLAGS  (options provided by J-C Sublet)
#
# FC="gfortran " FFLAGS=" -Wall -fcheck=all -Og -g -fbacktrace   "
# FC="ifort    " FFLAGS=" -O0 -g -traceback -check all -debug all"
# FC="nagfor   " FFLAGS=" -C=all -O0 -g -gline                   "
#
CWD = $(shell pwd)
PREFIX ?= $(CWD)/../

FC = gfortran
FFLAGS = -w

fsub = $(shell echo [A-z]*.f90)
obj = $(fsub:.f=.o)
src = $(fsub)

all: change resonancetables install clean

change : 
	../path_change.bash

resonancetables : $(obj)
	$(FC) $(FFLAGS) -o  resonancetables  $(obj)

install:
	mkdir -p $(PREFIX)bin
	mv resonancetables $(PREFIX)bin  

clean:
	rm -f *.o *.mod

