cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
project(zstr)

# zlib is required
FIND_PACKAGE(ZLIB 1.2.3 REQUIRED)
include_directories(BEFORE ${ZLIB_INCLUDE_DIRECTORIES})
GET_FILENAME_COMPONENT(ZLIB_BASEPATH ${ZLIB_INCLUDE_DIRS} PATH)
LINK_DIRECTORIES(${ZLIB_BASEPATH}/lib64 ${ZLIB_BASEPATH}/lib)
MESSAGE("Adding ZLIB for ${PROJECT_NAME}: -I${ZLIB_INCLUDE_DIRS} and -L${ZLIB_BASEPATH}/lib64 -L${ZLIB_BASEPATH}/lib)")

include_directories(src)
