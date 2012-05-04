#!/bin/bash

# Local Alignment Kernel
echo =================================
echo installing Local Alignment Kernel
echo =================================
cd lib/LAKernel/
make clean
make
cd ../..

# Mismatch Kernel
echo ==========================
echo installing Mismatch Kernel
echo ==========================
cd lib/MismatchKernel/src
make clean
make
cd ../../..

# PSIPRED
echo ==================
echo installing PSIPRED
echo ==================
cd lib/PSIPRED/src
make clean
make
make install
cd ../../..

# MoSta
echo ================
echo installing MoSta
echo ================
cd lib/MoSta/code
make clean
make all
cd ../../..

# GNU Scientific Library
echo ==============
echo installing GSL
echo ==============
cd lib/GSL/
tar -xzvf gsl-1.10.tar.gz
cd gsl-1.10
./configure --disable-shared --prefix=`pwd`
make
make install
cd ../../..

# STAMP
echo ================
echo installing STAMP
echo ================
export INSTALLDIR=`pwd`
cd lib/STAMP/code
./compileSTAMP $INSTALLDIR
cd ../../..