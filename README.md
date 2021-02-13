# FlowSolverRHEA

RHEA - an open-source Reproducible Hybrid-architecture flow solver Engineered for Academia

Rhea was the Titaness great Mother of the Gods, and goddess of female fertility, motherhood, and generation. Her name means "flow" and "ease", representing the eternal flow of time and generations with ease.

The flow solver RHEA solves the conservation equations of fluid motion for single-phase compressible flows by means of second-order central finite differences in combination with different all-speed Harten-Lax-van-Leer-type (HLL) Riemann solvers, and utilizes an explicit third-order strong-stability-preserving Runge-Kutta (SSP-RK3) method for time integration. RHEA is written in C++, using object-oriented programming, utilizes YAML and HDF5 for input/output operations, and targets hybrid supercomputing architectures.

INSTALLATION:
- Requisites: C++ compiler and MPI, YAML & HDF5 libraries/modules
- Clone/Download/Copy repository into working directory
- Modify Makefile (paths & flags) according to the computing system

COMPILATION:
- In myRHEA.cpp, overwrite setInitialConditions and calculateSourceTerms
- Compile flow solver by executing: $ make

EXECUTION:
- Set simulation parameters in configuration file (YAML)
- Execute simulation by running: $ RHEA.exe configuration_file.yaml
- Post-process output data using HDF5 and xdmf reader (optional)
