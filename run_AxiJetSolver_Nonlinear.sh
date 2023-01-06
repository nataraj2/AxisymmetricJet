rm -rf out
rm -rf patch_data.*
rm -rf ModuleVariables.mod ModInputOutput.mod ModCommunicate.mod ModAllocate.mod ModAssign.mod ModRungeKutta.mod ModSolver.mod ModDecomp.mod ModOperator.mod ModDeriv.mod ModTGfilter.mod ModShockCapturing.mod ModControlForcing.mod ModRHS.mod ModAxis.mod ModBC.mod ModPatch.mod

mpif90 -c -ffree-line-length-512 ModuleVariables.f90
mpicc -c plot3d_format.c
mpif90 -ffree-line-length-512 ModuleVariables.f90  ModDecomp.f90 ModPLOT3D_IO.f90 ModInputOutput.f90 ModCommunicate.f90 ModAllocate.f90 ModAssign.f90  ModAxis.f90  ModBC.f90 ModOperator.f90 ModTGfilter.f90 ModDeriv.f90 ModControlForcing.f90 ModShockCapturing.f90 ModPatch.f90 ModRHS.f90 ModRungeKutta.f90 ModSolver.f90 plot3d_format.o -o run_AxiJet
#ibrun -np 30 ./out
