cmake_minimum_required(VERSION 2.6)
enable_language(Fortran)

set(
	MAG_FSRCS
	cputim.F
	maggie.F
	wkstn.F
)

add_executable(maggie ${MAG_FSRCS})
target_include_directories(maggie PUBLIC "${PROJECT_SOURCE_DIR}/header")
target_link_libraries(maggie PUBLIC common)
target_link_libraries(maggie PUBLIC dist)


