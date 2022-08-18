#mpirun -np 4 --hostfile my-hostfile ./RHEA.exe configuration_file.yaml
mpirun --mca btl_base_warn_component_unused 0 -np 4 --hostfile my-hostfile ./RHEA.exe configuration_file.yaml
