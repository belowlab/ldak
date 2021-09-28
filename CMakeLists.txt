cmake_minimum_required(VERSION 3.16)
project(ldak C)

set(CMAKE_C_STANDARD 99)

SET(BLA_STATIC ON) # Find static LAPACK and BLAS.
find_package(LAPACK REQUIRED)
find_package(BLAS REQUIRED)

SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a") # Find static zlib
find_package(ZLIB REQUIRED)

add_executable(ldak5.1 ldak.c)

target_link_libraries(ldak5.1 ${ZLIB_LIBRARIES})
target_link_libraries(ldak5.1 ${LAPACK_LIBRARIES})
target_link_libraries(ldak5.1 ${BLAS_LIBRARIES})
target_link_libraries(ldak5.1 ${CMAKE_SOURCE_DIR}/libs/libqsopt.linux.a)
target_link_libraries(ldak5.1 gfortran) # gfortran for LAPACK
target_link_libraries(ldak5.1 quadmath) # quadmath for gfortran
target_link_libraries(ldak5.1 m)

target_link_libraries(ldak5.1 -static) # Static build. Not possible on MacOS because static libc is not available.

