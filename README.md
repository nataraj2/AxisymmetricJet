# Compressible flow solver for axisymmetric, nozzle jets
This repository contains a MPI parallelized direct numerical simualtion flow solver for 
aeroacoustic analysis of axisymmetric nozzles jets. The compressible Navier-Stokes 
equations are solved using fourth-order finite difference discretization with summation-by-parts (SBP) 
operators and the simulataneous approximation term (SAT) approach to implement the boundary conditions.
The governing equations in cylindrical polar coordinates for the compressible Navier-Stokes equations are  
 
$\frac{\partial\rho}{\partial t} + \Bigg(V_r\frac{\partial\rho}{\partial r} + \frac{V_\theta}{r}\frac{\partial\rho}{\partial\theta} + V_z\frac{\partial\rho}{\partial z}\Bigg) + \rho\Bigg(\frac{\partial V_r}{\partial r} + \frac{V_r}{r} + \frac{1}{r}\frac{\partial V_\theta}{\partial\theta} + \frac{\partial V_z}{\partial z}\Bigg) = 0 $  
$\rho\frac{\partial V_r}{\partial t} + \rho\Bigg(V_r\frac{\partial V_r}{\partial r} + \frac{V_\theta}{r}\frac{\partial V_r}{\partial\theta} - \frac{V_\theta^2}{r} + V_z\frac{\partial V_r}{\partial z}\Bigg) + \frac{1}{\gamma }\Bigg(\rho\frac{\partial T}{\partial r}+T\frac{\partial \rho}{\partial r}\Bigg) = S_r$    
$\rho\frac{\partial V_\theta}{\partial t} + \rho\Bigg(V_r\frac{\partial V_\theta}{\partial r} + \frac{V_\theta}{r}\frac{\partial V_\theta}{\partial\theta} + \frac{V_rV_\theta}{r} + V_z\frac{\partial V_\theta}{\partial z}\Bigg) + \frac{1}{\gamma }\Bigg(\frac{\rho}{r}\frac{\partial T}{\partial \theta}+\frac{T}{r}\frac{\partial\rho}{\partial\theta}\Bigg)= S_\theta$  
$\rho\frac{\partial V_z}{\partial t} + \rho\Bigg(V_r\frac{\partial V_z}{\partial r} + \frac{V_\theta}{r}\frac{\partial V_z}{\partial\theta} + V_z\frac{\partial V_z}{\partial z}\Bigg) + \frac{1}{\gamma }\Bigg(\rho\frac{\partial T}{\partial z} + T\frac{\partial\rho}{\partial z}\Bigg)= S_z$  
$\rho\frac{\partial T}{\partial t} + \rho\Bigg(V_r\frac{\partial T}{\partial r} + \frac{V_\theta}{r}\frac{\partial T}{\partial\theta} + V_z\frac{\partial T}{\partial z}\Bigg) + (\gamma-1)\rho T\Bigg(\frac{\partial V_r}{\partial r} + \frac{V_r}{r} + \frac{1}{r}\frac{\partial V_\theta}{\partial\theta} + \frac{\partial V_z}{\partial z}\Bigg)= S_T$ 

where  
 
$S_r =\frac{1}{Re}\Bigg[\mu\Bigg(\Bigg(2+\cfrac{\lambda}{\mu}\Bigg)\frac{\partial^2V_r}{\partial r^2} + \Bigg(2+\cfrac{\lambda}{\mu}\Bigg)\frac{1}{r}\frac{\partial V_r}{\partial r} + \Bigg(1+\cfrac{\lambda}{\mu}\Bigg)\frac{1}{r}\frac{\partial^2V_r}{\partial r\partial\theta} + \Bigg(1+\cfrac{\lambda}{\mu}\Bigg)\frac{\partial^2V_z}{\partial z\partial r}$$ + \frac{1}{r^2}\frac{\partial V_r^2}{\partial\theta^2}
-\Bigg(3+\cfrac{\lambda}{\mu}\Bigg)\frac{1}{r^2}\frac{\partial V_\theta}{\partial\theta}+\frac{\partial^2V_r}{\partial z^2}-\Bigg(2+\cfrac{\lambda}{\mu}\Bigg)\frac{V_r}{r^2}\Bigg)
+\frac{\partial\mu}{\partial r}\Bigg(\Bigg(2+\cfrac{\lambda}{\mu}\Bigg)\frac{\partial V_r}{\partial r}+\cfrac{\lambda}{\mu}\frac{V_r}{r}+\cfrac{\lambda}{\mu}\frac{1}{r}\frac{\partial V_\theta}{\partial\theta} +\cfrac{\lambda}{\mu}\frac{\partial V_z}{\partial z}\Bigg)+\frac{\partial\mu}{\partial\theta}\Bigg(\frac{1}{r}\frac{\partial V_\theta}{\partial r} + \frac{1}{r^2}\frac{\partial V_r}{\partial \theta}-\frac{V_\theta}{r^2}\Bigg) + \frac{\partial\mu}{\partial z}\Bigg(\frac{\partial V_r}{\partial z}+\frac{\partial V_z}{\partial r}\Bigg)\Bigg]$
 


## Mach 1.5 axisymmetric nozzle jet - vorticity and noise
<img src="Images/TimeAvgNoControl_Final.gif?raw=true&v=100" alt="your_alternative_text" width="100%" height="100%" loop="true" autoplay="true">

# Installation, compilation and running
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

# Visualization





