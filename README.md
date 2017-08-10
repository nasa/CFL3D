## CFL3D

CFL3D 
is a structured-grid, cell-centered, upwind-biased, Reynolds-averaged Navier-Stokes (RANS) code. It can be run
in parallel on multiple grid zones with point-matched, patched, overset, or embedded connectivities. Both
multigrid and mesh sequencing are available in time-accurate or steady-state modes.

The most up-to-date information can be found on the web at:

https://cfl3d.larc.nasa.gov

-------------

Copyright 2001 United States Government as represented by the Administrator
of the National Aeronautics and Space Administration. All Rights Reserved.

The CFL3D platform is licensed under the Apache License, Version 2.0 
(the "License"); you may not use this file except in compliance with the 
License. You may obtain a copy of the License at 
http://www.apache.org/licenses/LICENSE-2.0. 

Unless required by applicable law or agreed to in writing, software 
distributed under the License is distributed on an "AS IS" BASIS, WITHOUT 
WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the 
License for the specific language governing permissions and limitations 
under the License.

------------
## About the CMake branch
As a practice of cmake, I had ported the make system of cfl3d to CMake, a modern industry-level build system. Which is much more easy to use, understand, modify and maintain.

Now the CMake branch has been merged into master

### Compilation
```shell
## prerequirement (tested in Ubuntu 16.04)
sudo apt-get install libopenmpi-dev gfortran build-essential git cmake -y
# in RHEL/CentOS, you may try this:
# sudo yum install openmpi-devel gcc gcc-gfortran git cmake

## download  
git clone https://git.coding.net/chengdi123000/CFL3D.git

## manually build cgns 2.5 without hdf5 and mpi
# you can also use your own cgns library if you put your cgnslib_f.h and libcgns.a in /usr/share/include and /usr/share/lib respectively
cd CFL3D/external/cgns
./build_cgns

## build cfl3d executables out of source
cd ../../build
cmake ..
make
# or you may benefit from parallel run, try:
# make -j`nproc`
# There will be a lot of warnings because I do not turn off all warnings as the original CFL3D build system do.
```

### Test
First set the PATH
```shell
## When completed, the executables will be in CFL3D/build/bin directory.
# suggest to use this line of script to add executables to PATH
# export PATH=$PWD/bin:$PATH
# You can use following line of script to add a command `use_cfl3d` to load executables ready to run after re-login
# echo "alias use_cfl3d='export PATH=$PWD/bin:\$PATH'">>$HOME/.bashrc
```
Then you need to download testcases
```shell
cd CFL3D/testcase
./download_testcases #wget is needed.
cd cfl3d.larc.nasa.gov/Cfl3dv6/2DTestcases/Flatplate
tar xvf Flatplate.tar.Z
# use_cfl3d
splitter < split.inp_1blk 
cfl3d_seq < grdflat5.inp &
tail -f cfl3d.out #monitor the steps and residuals
```
