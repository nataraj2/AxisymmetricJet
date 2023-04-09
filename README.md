# Compressible flow solver for axisymmetric, nozzle jets
This repository contains a MPI parallelized direct numerical simualtion flow solver for 
aeroacoustic analysis of axisymmetric nozzles jets. The compressible Navier-Stokes 
equations are solved using fourth-order finite difference discretization with summation-by-parts (SBP) 
operators and the simulataneous approximation term (SAT) approach to implement the boundary conditions.

## Mach 1.5 axisymmetric nozzle jet
<img src="Images/TimeAvgNoControl_Final.gif?raw=true&v=100" alt="your_alternative_text" width="100%" height="100%" loop="true" autoplay="true">

# Installation and compilation 
The following are the instructions for compiling on Stampede. 
The gcc compilers in `/opt/apps/gcc9_1/mvapich2/2.3.7/bin` have to be used.
The intel compilers have some issue in the plot3d file reading. The compile script 
`run_compile.sh` uses the gcc compilers. 
```
git clone https://github.com/nataraj2/AxisymmetricJet.git
cd AxisymmetricJet
sh run_compile.sh 
idev -p development -N 2 -n 128 -m 150
ibrun -n <nprocs> ./run_AxiJet
```
Make sure `nprocs` is the product of the integers in `dims` in `ModuleVariables.f90`. 
`noutput` in `ModuleVariables.f90` is the frequency of writing the output solution files.

# Running

# Visualization



