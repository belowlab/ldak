# ldak
A Fork of Doug Speed's LDAK, Linkage-Disequilibrium Adjusted Kinships. The website of the original project with _very_ complete documetation can be found at http://dougspeed.com/.

## Included in this repository:
* LDAK source .c files
* Static library builds of libqsopt, https://www.math.uwaterloo.ca/~bico/qsopt/. Included are x86_64 Linux, x86_64 Mac, and Apple Silicon
* This documentation

## Building
Right now, there's no build system. LDAK's original documentation suggests that it can be built with 
```
gcc -O3 -o ldak5.1 ldak.c libqsopt_ex.a -lblas -llapack -lm -lz
```

### Alex's Getting It Working on arm64 notes:
* Went to go build a static libqsopt library from a modern fork: https://github.com/jonls/qsopt-ex
* Had to install GMP, which I did from homebrew
* Had to add homebrew to the include path via CPATH and LIBRARY_PATH
```shell
CPATH=/opt/homebrew/include/
LIBRARY_PATH=/opt/homebrew/lib/
```
* Had to install that qsopt-ex fork with `sudo make install`
* Updated ldak.c to include it with `#include <qsopt_ex/QSopt_ex.h>`

# Getting it building in CMake
In order to work around what I think is a CMAKE bug finding generic lapack, I had to do this: `cmake -DCMAKE_BUILD_TYPE=Release -DLAPACK_LIBRARIES=/usr/lib/x86_64-linux-gnu/liblapack.a ..`. On Ubuntu 20.04, I got it working
with the CMakeLists.txt that I'm checking in now, both dynamically linked and statically. Not only that, but it's using the openBLAS, which should be faster than generic BLAS or Intel BLAS on AMD CPUs.

# Dependencies to build on Ubuntu:
`sudo apt install cmake libopenblas-dev zlib1g-dev gfortran`

# Original documentation:
  To run LDAK using Linux:

1 - Download and unzip the Linux executable file (link at the bottom).

2 - Open a terminal window (on my Linux, I go to Applications / System Tools / Terminal)

3 - Type the name of the file; for example
/home/doug/Downloads/ldak5.1.linux
or if you are in the same folder as the executable file, you can simply type
./ldak5.1.linux

If compatible with your system, this should produce the LDAK welcome screen. Note that if your computer tells you that you do not have permission to run the file, then first run
chmod a+x ldak5.1.linux
_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _

The pre-compiled Linux version uses the Intel MKL Libraries. The command I used was

gcc -O3 -o ldak5.1.linux source/ldak.c source/libqsopt.linux.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lgomp -lpthread -lm -ldl -m64 -I${MKLROOT}/include -lz -fopenmp -static

Should you wish to compile a Linux version yourself, please download and unzip the source code, then from inside that folder run a command similar to

gcc -O3 -o ldak5.1 ldak.c libqsopt_ex.a -lblas -llapack -lm -lz
chmod a+x ldak5.1

The exact command will depend on which libraries you have installed. This should take less than a minute to complete. Note that if you do not have Intel MKL Libraries installed, you will be required to turn off the MKL libraries, by editing Line 63 of ldak.c (replace #define MKL 1 with #define MKL 0). For this reason, a self-compiled version will likely be slower than the pre-compiled version.
