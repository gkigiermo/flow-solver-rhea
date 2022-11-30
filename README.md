# FlowSolverRHEA

RHEA: an open-source Reproducible Hybrid-architecture flow solver Engineered for Academia

Rhea was the Titaness great Mother of the ancient Greek Gods, and goddess of female fertility, motherhood, and generation. Her name means "flow" and "ease", representing the eternal flow of time and generations with ease.

The flow solver RHEA solves the conservation equations of fluid motion for compressible flows by means of second-order numerical schemes in combination with kinetic energy preserving or different all-speed Harten-Lax-van-Leer-type (HLL) Riemann solvers, and utilizes explicit Runge-Kutta methods for time integration. RHEA is written in C++, using object-oriented programming, utilizes YAML and HDF5 for input/output operations, and targets hybrid supercomputing architectures.

INSTALLATION:
- Requisites: C++ compiler and MPI library, YAML (yaml-cpp version 1.2) & HDF5 libraries/modules (compatible with MPI)
- OpenACC requisites (if CPU-GPU available): NVIDIA HPC SDK compiler (version 21.5)
- Clone/Download/Copy repository into working directory
- Set environment variable with RHEA's path: export RHEA_PATH='path_to_rhea_root_directory' 
- Adapt Makefile (flags & paths: CXXFLAGS, INC_LIB_YAML, INC_DIR_YAML, INC_LIB_HDF5, INC_DIR_HDF5) to the computing system. Examples:
   - CPU flags: CXXFLAGS = -O3 -Wall -std=c++0x -Wno-unknown-pragmas -Wno-unused-variable -I$(PROJECT_PATH)
   - CPU-GPU flags with managed OpenACC: CXXFLAGS = -fast -acc -ta=tesla:managed -Minfo=accel -O3 -Wall -std=c++0x -I$(PROJECT_PATH)
   - CPU-GPU flags with non-managed OpenACC: CXXFLAGS= -fast -acc -ta=tesla,pinned -Minfo=accel -O3 -Wall -std=c++0x -I$(PROJECT_PATH)
   - Ubuntu - Linux paths:
      INC_LIB_YAML =
      INC_DIR_YAML =
      INC_LIB_HDF5 = -L/usr/lib/x86_64-linux-gnu/hdf5/openmpi 
      INC_DIR_HDF5 = -I/usr/include/hdf5/openmpi
   - MAC - OS X paths:
      INC_LIB_YAML = -L/usr/local/lib
      INC_DIR_YAML = -I/usr/local/include
      INC_LIB_HDF5 =
      INC_DIR_HDF5 =

COMPILATION:
- In myRHEA.cpp, overwrite setInitialConditions and calculateSourceTerms
- Compile flow solver by executing: $ make

EXECUTION:
- Set simulation parameters in configuration file (YAML)
- Execute simulation by running: $ RHEA.exe configuration_file.yaml
- Post-process output data using HDF5 and xdmf reader (optional)

EXAMPLE TESTS:
- Change directory to any of the examples provided in the tests folder
- Adapt Makefile to the computing system (or copy the Makefile previously modified)
- Compile main file (myRHEA.cpp) by executing: $ make
- Start test by running (4 CPUs by default): $ ./execute.sh

SUPPORT, ISSUES & CONTRIBUTION
- Documentation of the classes, functions and variables is available in doxygen/html
- Additional support is provided via the "Service Desk" of the project
- Issues can be reported utilizing the "Issues" functionality menu
- Contributions can be proposed, discussed and incorporated through "Merge requests"

--------------------------------------------------
STEP-BY-STEP GUIDE: Ubuntu 22.04.1 LTS (Jammy Jellyfish)

The lines below provide detailed instructions for installing, compiling, running and post-processing the 3d_turbulent_channel_flow test on CPU-based computing architectures.

INSTALLATION:
- Install (if not already available) C++ compiler & toolchain:
$ sudo apt install build-essential make
- Install (required) openmpi, yaml-cpp and hdf5 libraries:
$ sudo apt install libopenmpi-dev libyaml-cpp-dev libhdf5-openmpi-dev
- Install (auxiliary) post-processing packages & libraries:
$ sudo apt install texlive-full python3.8 python3-pip paraview
- Install (auxiliary) python packages & libraries:
$ sudo pip3 install numpy scipy python3-matplotlib h5py
- Clone repository into working directory (e.g., $HOME):
$ git clone https://gitlab.com/ProjectRHEA/flowsolverrhea.git
- Write RHEA's path in ~/.bashrc (e.g., $HOME directory):
export RHEA_PATH='$HOME'
- Reload bashrc settings:
$ source ~/.bashrc
- Makefile files are already prepared for this installation

COMPILATION:
- Change directory to the 3d_turbulent_channel_flow test:
$ cd $RHEA_PATH/tests/3d_turbulent_channel_flow
- Compile test case by executing:
$ make

EXECUTION:
- Modify configuration_file.yaml file to output data every 1000 iterations:
output_frequency_iter: 1000
- Execute simulation on 4 CPUs by running:
$ mpirun -np 4 ./RHEA.exe configuration_file.yaml
- Visualize output files with paraview:
$ paraview
- Post-process output data using python script:
$ python3 post_process_plot_script.py 
--------------------------------------------------
