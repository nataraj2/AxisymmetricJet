# Compressible flow solver for axisymmetric, nozzle jets

# Installation and compilation 
The following are the instructions for compiling on Stampede. 
The gcc compilers in `/opt/apps/gcc9_1/mvapich2/2.3.7/bin` have to be used.
The intel compilers have some issue in the plot3d file reading. The compile script 
`run_compile.sh` uses the gcc compilers. 
```
git clone https://github.com/nataraj2/AxisymmetricJet.git
cd AxisymmetricJet
sh compile.sh 
ibrun -n <nprocs> ./run_AxiJet
```
Make sure `nprocs` is the product of the integers in `dims` in `ModuleVariables.f90`
`noutput` in ModuleVariables.f90 is the frequency of writing the output solution files.

# Running

# Visualization



