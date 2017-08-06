cmake_minimum_required(VERSION 3.0.2)
enable_language(Fortran C)

add_library(cgns STATIC IMPORTED GLOBAL)
message("add ${CMAKE_CURRENT_SOURCE_DIR}/cgnslib_2.5/LINUX/libcgns.a as CGNS usage")
set_property(TARGET cgns PROPERTY IMPORTED_LOCATION ${CMAKE_CURRENT_SOURCE_DIR}/cgnslib_2.5/LINUX/libcgns.a)
set_property(TARGET cgns PROPERTY IMPORTED_LINK_INTERFACE_LANGUAGES Fortran)
set_property(TARGET cgns PROPERTY INTERFACE_INCLUDE_DIRECTORIES $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/cgnslib_2.5>)
