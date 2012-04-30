#!/bin/bash

# Local Alignment Kernel
echo =================================
echo installing Local Alignment Kernel
echo =================================
cd LAKernel/
make clean
make
cd ..

# Mismatch Kernel
echo ==========================
echo installing Mismatch Kernel
echo ==========================
cd MismatchKernel/src
make clean
make
cd ../..

# PSIPRED
echo ==================
echo installing PSIPRED
echo ==================
cd PSIPRED/src
make clean
make
make install
cd ../..

# MoSta
echo ================
echo installing MoSta
echo ================
cd MoSta/code
make clean
make all
cd ../..

# GNU Scientific Library
echo ==============
echo installing GSL
echo ==============
cd GSL/
tar -xzvf gsl-1.10.tar.gz
cd gsl-1.10
./configure --disable-shared --prefix=`pwd`
make
make install
cd ..
cd ..

# STAMP
echo ================
echo installing STAMP
echo ================
export INSTALLDIR=`pwd`
cd STAMP/code
./compileSTAMP $INSTALLDIR
cd ../..