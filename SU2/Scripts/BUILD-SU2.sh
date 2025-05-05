#! /bin/bash -e

export CXXFLAGS="-march=native -mtune=native"

./meson.py build \
  --prefix=/opt/SU2/8.1.0 \
  -Denable-autodiff=true \
  -Denable-directdiff=true \
  -Denable-pywrapper=false \
  -Dwith-mpi=disabled \
  -Dwith-omp=true \
  -Denable-cgns=false \
  -Denable-tecio=true \
  -Denable-mkl=true \
  -Dmkl_root=/opt/intel/oneapi/mkl/latest \
  -Denable-openblas=false \
  -Denable-pastix=false  \
  -Denable-mpp=false \
  -Denable-mixedprec=false \
  -Denable-mlpcpp=true \
  -Denable-gprof=false \
  --buildtype=release  \
  --optimization=3

sudo ./ninja -C build install
