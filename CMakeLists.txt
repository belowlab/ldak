cmake_minimum_required(VERSION 3.20)
project(ldak C)

set(CMAKE_C_STANDARD 99)

find_package(ZLIB REQUIRED)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

add_executable(ldak5.1 ldak.c)

target_link_libraries(ldak5.1 m)
target_link_libraries(ldak5.1 ${ZLIB_LIBRARIES})
target_link_libraries(ldak5.1 ${BLAS_LIBRARIES})
target_link_libraries(ldak5.1 ${LAPACK_LIBRARIES})
target_link_libraries(ldak5.1 ${CMAKE_SOURCE_DIR}/libs/qsopt.macarm64.a)
