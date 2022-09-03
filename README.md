# FlowSolverRHEA

RHEA: an open-source Reproducible Hybrid-architecture flow solver Engineered for Academia

Rhea was the Titaness great Mother of the ancient Greek Gods, and goddess of female fertility, motherhood, and generation. Her name means "flow" and "ease", representing the eternal flow of time and generations with ease.

The flow solver RHEA solves the conservation equations of fluid motion for compressible flows by means of second-order numerical schemes in combination with kinetic energy preserving or different all-speed Harten-Lax-van-Leer-type (HLL) Riemann solvers, and utilizes explicit Runge-Kutta methods for time integration. RHEA is written in C++, using object-oriented programming, utilizes YAML and HDF5 for input/output operations, and targets hybrid supercomputing architectures.

INSTALLATION:
- Requisites: C++ compiler and MPI library, YAML (yaml-cpp version 1.2) & HDF5 libraries/modules (compatible with MPI)
- OpenACC requisites (if CPU-GPU available): NVIDIA HPC SDK compiler (version 21.5)
- Clone/Download/Copy repository into working directory
- Set environment variable with RHEA's path: export RHEA_PATH='path_to_rhea_root_directory' 
- Adapt Makefile (compiler, flags & paths: CXX, CXXFLAGS, LDFLAGS, INC_LIB_YAML, INC_DIR_YAML, INC_LIB_HDF5, INC_DIR_HDF5) to the computing system

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
