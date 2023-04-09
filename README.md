# Compressible flow solver for axisymmetric, nozzle jets
This repository contains a MPI parallelized direct numerical simualtion flow solver for 
aeroacoustic analysis of axisymmetric nozzles jets. The compressible Navier-Stokes 
equations are solved using fourth-order finite difference discretization with summation-by-parts (SBP) 
operators and the simulataneous approximation term (SAT) approach to implement the boundary conditions.

## Mach 1.5 axisymmetric nozzle jet
<img src="Images/TimeAvgNoControl_Final.gif?raw=true&v=100" alt="your_alternative_text" width="50%" height="50%" loop="true" autoplay="true">

# Installation and compilation 
```
git clone https://github.com/nataraj2/AxisymmetricJet.git
cd AxisymmetricJet
sh run_AxiJetSolver_Nonlinear.sh 
```

# Running

# Visualization



