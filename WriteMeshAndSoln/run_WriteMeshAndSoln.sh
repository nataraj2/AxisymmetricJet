rm -rf out plot3d_format.o
mpicc -c plot3d_format.c
mpif90 ModPLOT3D_IO.f90 write_soln_file.f90 plot3d_format.o -o out
./out 
