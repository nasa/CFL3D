cmake_minimum_required(VERSION 2.6)
enable_language(Fortran)

set(
PRERON_FSRCS
main_pre.F
sizer.F
termn8.F
usrint.F
)

set(
RON_FSRCS
main.F
ronnie.F
setup.F
sizer.F
termn8.F
usrint.F
../cfl3d/dist/reass.F
)
add_executable(ronnie ${RON_FSRCS})
target_link_libraries(ronnie PUBLIC common)
target_link_libraries(ronnie PUBLIC dist)

add_executable(preronnie ${PRERON_FSRCS})
target_link_libraries(preronnie PUBLIC common)
target_link_libraries(preronnie PUBLIC dist)

