# Compressible flow solver for axisymmetric, nozzle jets

## Mach 1.5 axisymmetric nozzle jet - vorticity and noise
<img src="Images/TimeAvgNoControl_Final.gif?raw=true&v=100" alt="your_alternative_text" width="100%" height="100%" loop="true" autoplay="true">

## Introduction
This repository contains a MPI parallelized direct numerical simualtion flow solver for 
aeroacoustic analysis of axisymmetric nozzles jets. The compressible Navier-Stokes 
equations are solved using fourth-order Runge-Kutta for the temporal discretization, 
fourth-order finite difference discretization for the spatial derivatives with summation-by-parts (SBP) 
operators, and the simulataneous approximation term (SAT) approach to implement the boundary conditions.

## Governing equations

The governing equations in cylindrical polar coordinates for the compressible Navier-Stokes equations are  
 
$\frac{\partial\rho}{\partial t} + \Bigg(V_r\frac{\partial\rho}{\partial r} + \frac{V_\theta}{r}\frac{\partial\rho}{\partial\theta} + V_z\frac{\partial\rho}{\partial z}\Bigg) + \rho\Bigg(\frac{\partial V_r}{\partial r} + \frac{V_r}{r} + \frac{1}{r}\frac{\partial V_\theta}{\partial\theta} + \frac{\partial V_z}{\partial z}\Bigg) = 0 $  
$\rho\frac{\partial V_r}{\partial t} + \rho\Bigg(V_r\frac{\partial V_r}{\partial r} + \frac{V_\theta}{r}\frac{\partial V_r}{\partial\theta} - \frac{V_\theta^2}{r} + V_z\frac{\partial V_r}{\partial z}\Bigg) + \frac{1}{\gamma }\Bigg(\rho\frac{\partial T}{\partial r}+T\frac{\partial \rho}{\partial r}\Bigg) = S_r$    
$\rho\frac{\partial V_\theta}{\partial t} + \rho\Bigg(V_r\frac{\partial V_\theta}{\partial r} + \frac{V_\theta}{r}\frac{\partial V_\theta}{\partial\theta} + \frac{V_rV_\theta}{r} + V_z\frac{\partial V_\theta}{\partial z}\Bigg) + \frac{1}{\gamma }\Bigg(\frac{\rho}{r}\frac{\partial T}{\partial \theta}+\frac{T}{r}\frac{\partial\rho}{\partial\theta}\Bigg)= S_\theta$  
$\rho\frac{\partial V_z}{\partial t} + \rho\Bigg(V_r\frac{\partial V_z}{\partial r} + \frac{V_\theta}{r}\frac{\partial V_z}{\partial\theta} + V_z\frac{\partial V_z}{\partial z}\Bigg) + \frac{1}{\gamma }\Bigg(\rho\frac{\partial T}{\partial z} + T\frac{\partial\rho}{\partial z}\Bigg)= S_z$  
$\rho\frac{\partial T}{\partial t} + \rho\Bigg(V_r\frac{\partial T}{\partial r} + \frac{V_\theta}{r}\frac{\partial T}{\partial\theta} + V_z\frac{\partial T}{\partial z}\Bigg) + (\gamma-1)\rho T\Bigg(\frac{\partial V_r}{\partial r} + \frac{V_r}{r} + \frac{1}{r}\frac{\partial V_\theta}{\partial\theta} + \frac{\partial V_z}{\partial z}\Bigg)= S_T$ 

where $S_r$, $S_\theta$, $S_z$ and $S_T$ include all viscous and heat transfer effects, and are given by 
 
$S_r =\frac{1}{Re}\Bigg[\mu\Bigg(\Bigg(2+\cfrac{\lambda}{\mu}\Bigg)\frac{\partial^2V_r}{\partial r^2} + \Bigg(2+\cfrac{\lambda}{\mu}\Bigg)\frac{1}{r}\frac{\partial V_r}{\partial r} + \Bigg(1+\cfrac{\lambda}{\mu}\Bigg)\frac{1}{r}\frac{\partial^2V_r}{\partial r\partial\theta} + \Bigg(1+\cfrac{\lambda}{\mu}\Bigg)\frac{\partial^2V_z}{\partial z\partial r} +\frac{1}{r^2}\frac{\partial V_r^2}{\partial\theta^2}-\Bigg(3+\cfrac{\lambda}{\mu}\Bigg)\frac{1}{r^2}\frac{\partial V_\theta}{\partial\theta}+\frac{\partial^2V_r}{\partial z^2}$  
$-\Bigg(2+\cfrac{\lambda}{\mu}\Bigg)\frac{V_r}{r^2}\Bigg)+\frac{\partial\mu}{\partial r}\Bigg(\Bigg(2+\cfrac{\lambda}{\mu}\Bigg)\frac{\partial V_r}{\partial r}+\cfrac{\lambda}{\mu}\frac{V_r}{r}+\cfrac{\lambda}{\mu}\frac{1}{r}\frac{\partial V_\theta}{\partial\theta} +\cfrac{\lambda}{\mu}\frac{\partial V_z}{\partial z}\Bigg)+\frac{\partial\mu}{\partial\theta}\Bigg(\frac{1}{r}\frac{\partial V_\theta}{\partial r} + \frac{1}{r^2}\frac{\partial V_r}{\partial \theta}-\frac{V_\theta}{r^2}\Bigg) + \frac{\partial\mu}{\partial z}\Bigg(\frac{\partial V_r}{\partial z}+\frac{\partial V_z}{\partial r}\Bigg)\Bigg]$  

$S_\theta=\frac{1}{Re}\Bigg[\mu\Bigg(\frac{\partial^2V_\theta}{\partial r^2}+\Bigg(1+\cfrac{\lambda}{\mu}\Bigg)\frac{1}{r}\frac{\partial^2V_\theta}{\partial\theta\partial r}+\Bigg(2+\cfrac{\lambda}{\mu}\Bigg)\frac{1}{r^2}\frac{\partial^2V_\theta}{\partial\theta^2}+ \Bigg(3+\cfrac{\lambda}{\mu}\Bigg)\frac{1}{r^2}\frac{\partial V_r}{\partial\theta}+\Bigg(1+\cfrac{\lambda}{\mu}\Bigg)\frac{1}{r}\frac{\partial^2V_z}{\partial\theta\partial z} + \frac{\partial^2V_\theta}{\partial z^2} + \frac{1}{r}\frac{\partial V_\theta}{\partial r} - \frac{V_\theta}{r^2}\Bigg)+$
$\frac{\partial\mu}{\partial r}\Bigg(\frac{\partial V_\theta}{\partial r} + \frac{1}{r}\frac{\partial V_r}{\partial\theta} - \frac{V_\theta}{r}\Bigg)$ 
$+\frac{\partial\mu}{\partial\theta}\Bigg(\Bigg(2+\cfrac{\lambda}{\mu}\Bigg)\frac{1}{r^2}\frac{\partial V_\theta}{\partial\theta} + \cfrac{(2\mu+\lambda)}{\mu}\frac{V_r}{r^2}+\cfrac{\lambda}{\mu}\frac{1}{r}\frac{\partial V_r}{\partial r} +\cfrac{\lambda}{\mu}\frac{1}{r}\frac{\partial V_z}{\partial z}\Bigg)+\frac{\partial\mu}{\partial z}\Bigg(\frac{\partial V_\theta}{\partial z} + \frac{1}{r}\frac{\partial V_z}{\partial\theta}\Bigg)\Bigg]$    
  

$S_z=\frac{1}{Re}\Bigg[\mu\Bigg(\frac{\partial^2V_z}{\partial r^2} + \frac{1}{r^2}\frac{\partial^2}{\partial\theta^2} + \Bigg(2+\cfrac{\lambda}{\mu}\Bigg)\frac{\partial^2V_z}{\partial z^2}+\Bigg(1+\cfrac{\lambda}{\mu}\Bigg)\frac{\partial^2 V_r}{\partial r\partial z} + \frac{1}{r}\frac{\partial V_z}{\partial r}+\Bigg(1+\cfrac{\lambda}{\mu}\Bigg)\frac{1}{r}\frac{\partial V_r}{\partial z} + \Bigg(1+\cfrac{\lambda}{\mu}\Bigg)\frac{1}{r}\frac{\partial^2V_\theta}{\partial\theta\partial z}\Bigg)$  
$+\frac{\partial\mu}{\partial r}\Bigg(\frac{\partial V_z}{\partial r} + \frac{\partial V_r}{\partial z}\Bigg) + \frac{\partial\mu}{\partial\theta}\Bigg(\frac{1}{r^2}\frac{\partial V_z}{\partial\theta} + \frac{1}{r}\frac{\partial V_\theta}{\partial z}\Bigg)+\frac{\partial\mu}{\partial z}\Bigg(\Bigg(2+\cfrac{\lambda}{\mu}\Bigg)\frac{\partial V_z}{\partial z}+\cfrac{\lambda}{\mu}\frac{\partial V_r}{\partial r} +\cfrac{\lambda}{\mu}\frac{V_r}{r} +\cfrac{\lambda}{\mu}\frac{1}{r}\frac{\partial V_\theta}{\partial\theta}\Bigg)\Bigg]$   

$S_T=\frac{\gamma}{RePr}\Bigg[\kappa\Bigg(\frac{\partial^2T}{\partial r^2}+\frac{1}{r}\frac{\partial T}{\partial r}+\frac{1}{r^2}\frac{\partial^2T}{\partial\theta^2}+\frac{\partial^2T}{\partial z^2}\Bigg) + \frac{\partial\kappa}{\partial r}\frac{\partial T}{\partial r} + \frac{1}{r^2}\frac{\partial\kappa}{\partial\theta}\frac{\partial T}{\partial\theta}+\frac{\partial\kappa}{\partial z}\frac{\partial T}{\partial z}\Bigg]$
$+\frac{\gamma(\gamma-1)}{Re}\mu\Bigg[\frac{\partial V_r}{\partial r}\Bigg((2\mu+\lambda)\frac{\partial V_r}{\partial r} +\lambda\frac{V_r}{r} +\lambda\frac{1}{r}\frac{\partial V_\theta}{\partial\theta}+\lambda\frac{\partial V_z}{\partial z}\Bigg)$
$+\Bigg(\frac{\partial V_r}{\partial\theta}+\frac{\partial V_\theta}{\partial r} - \frac{V_\theta}{r}\Bigg)\Bigg(\frac{\partial V_\theta}{\partial r} + \frac{1}{r}\frac{\partial V_r}{\partial\theta} - \frac{V_r}{r}\Bigg)$
$+\Bigg(\frac{1}{r}\frac{\partial V_\theta}{\partial\theta}+\frac{V_r}{r}\Bigg)\Bigg((2\mu+\lambda)\frac{1}{r}\frac{\partial V_\theta}{\partial\theta}+(2\mu+\lambda)\frac{V_r}{r}+\lambda\frac{\partial V_r}{\partial r} +\lambda\frac{\partial V_z}{\partial z}\Bigg)$
$+\Bigg(\frac{\partial V_\theta}{\partial z}+\frac{1}{r}\frac{\partial V_z}{\partial\theta}\Bigg)\Bigg(\frac{\partial V_\theta}{\partial z}+\frac{1}{r}\frac{\partial V_z}{\partial\theta}\Bigg)$
$+\frac{\partial V_z}{\partial z}\Bigg((2\mu+\lambda)\frac{\partial V_z}{\partial z}+\lambda\frac{\partial V_r}{\partial r} +\lambda\frac{V_r}{r}+\lambda\frac{1}{r}\frac{\partial V_\theta}{\partial\theta}\Bigg) + \Bigg(\frac{\partial V_r}{\partial z}+\frac{\partial V_z}{\partial r}\Bigg)\Bigg(\frac{\partial V_z}{\partial r}+\frac{\partial V_r}{\partial z}\Bigg)\Bigg]$  
 
A standard power law describes the temperature dependence of the fluid viscosity, $\mu$, and thermal conductivity, $\kappa$ as $\mu =\kappa=T^{2/3}$. The bulk viscosity is $\mu_B=\lambda+2/3\mu=0.6\mu$, where $\lambda$ is the second coefficient of viscosity and the Prandtl number $Pr=\mu C_p/\kappa=0.72$, where $C_p$ is the specific heat at constant pressure. For this nondimensionalization, the equation of state is $p = \rho T$.  

# Installation, compilation and running
The following are the instructions for compiling on the Stampede2 supercomputer at Texas Advanced Supercomputing Center (TACC). 
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
The I/O format used for the files is PLOT3D, which is a standard format for curvilinear, structured meshes and can be read into standard visualization 
packages such as VisIt, ParaView and Tecplot.




