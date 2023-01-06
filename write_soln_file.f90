PROGRAM write_soln_file


	
	USE module_variables
	USE ModGlobal
	USE ModPlot3D_IO 

	IMPLICIT NONE


	CALL Read_Grid(NDIM, ngrid, ND, G, IBLANK, prec, gf, vf, ib, GRID_FILENAME, OF)        


	Nr = ND(1,2)
	Nz = ND(1,1)

	ALLOCATE(rho_rhoj(Nr))
	ALLOCATE(rho_rhoinf(Nr))
	ALLOCATE(u_Uj(Nr))
	ALLOCATE(u_cinf(Nr))
	ALLOCATE(T_Tinf(Nr))
	ALLOCATE(rad(Nr))


       ALLOCATE(X3(1,ND(1,1),ND(1,2),ND(1,3),5)) 	

       X3 = 0.0D0

       gam     = 1.4D0	
       Mj      = 1.5D0
       Tj_Tinf = 0.6D0
       Uj_cinf = Mj*(Tj_Tinf)**0.5D0	
       theta = 0.1


	DO i = 1, Nr

	rad(i) = G(1,1,i,1,2)
	!PRINT*, r(i)
	
	ENDDO


      
	DO i = 1, Nr
           
             u_Uj(i) = (1.0D0/2.0D0)*( 1 - tanh(1/(4*theta)*(rad(i) - 1.0D0/(rad(i)))) )		
	     u_cinf(i) = u_Uj(i)*Uj_cinf

	ENDDO


	   DO i = 1, Nr

	    term1 = (1.0D0/2.0D0)*(gam -1)*u_Uj(i)*(1 - u_Uj(i))*Mj**2 
    
	    term2 = u_Uj(i) 
    
	    term3 = (1.0D0/Tj_Tinf) * (1 - u_Uj(i)) 
    
	    rho_rhoj(i)   = 1.0D0/(term1 + term2 + term3)

	    rho_rhoinf(i)  = rho_rhoj(i) / Tj_Tinf
    
	  ENDDO   
  

	DO j = 1, Nz
	   DO i = 1, Nr
	      DO k = 1, 1	

		 X3(1,j,i,k,1) = rho_rhoinf(i)		
		 X3(1,j,i,k,2) = 0.0D0!0.2D0*u_cinf(j)
		 X3(1,j,i,k,3) = 0.0D0
		 X3(1,j,i,k,4) = u_cinf(i)
		 X3(1,j,i,k,5) = 1.0D0/rho_rhoinf(i)

	
	      ENDDO
	   ENDDO
	ENDDO





	tau = 0.0D0

	CALL Write_Soln(NDIM, ngrid, ND, X3, tau, prec, gf, vf, EULER_SOLUTION_FILE)




END PROGRAM write_soln_file
