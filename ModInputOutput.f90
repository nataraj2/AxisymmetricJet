MODULE ModInputOutput

CONTAINS

SUBROUTINE GridRead

USE ModuleVariables
USE ModPlot3D_IO
USE ModDecomp
IMPLICIT NONE

INTEGER(KIND=8) :: datasize
		

	IF(rank .eq. 0)THEN
	CALL Read_Grid(NDIM, ngrid, ND, G, IBLANK, prec, gf, vf, ib, GRID_FILENAME, OF)        	
	Nz_grid = ND(1,1)
	Nr_grid = ND(1,2)	
	PRINT*, "Read grid done"		
	ENDIF

	neqns = 5

	
       DO procrank = 1, numtasks-1

        IF(rank .eq. 0)THEN 
        CALL MPI_Send(ND, 3, MPI_Integer8, procrank, 1, MPI_COMM_WORLD, ierr)
        ELSEIF(rank .eq. procrank)THEN
	ALLOCATE(ND(1,3))
        CALL MPI_Recv(ND, 3, MPI_Integer8, 0, 1, MPI_COMM_WORLD,stat, ierr)             
	Nz_grid = ND(1,1)
	Nr_grid = ND(1,2)
        ENDIF

       ENDDO

	PRINT*, "I have", rank, ND(1,1), ND(1,2), ND(1,3)

       DO procrank = 1, numtasks-1

        IF(rank .eq. 0)THEN     
  
	datasize = Nz_grid*Nr_grid*3
        CALL MPI_Send(G,datasize,MPI_Real8, procrank, 1, MPI_COMM_WORLD, ierr)
	CALL MPI_Send(IBLANK,Nz_grid*Nr_grid,MPI_Integer, procrank, 1, MPI_COMM_WORLD, ierr)

        ELSEIF(rank .eq. procrank)THEN	

        ALLOCATE(G(1,ND(1,1),ND(1,2),ND(1,3),3))
	ALLOCATE(IBLANK(1,ND(1,1),ND(1,2),ND(1,3)))

	datasize = Nz_grid*Nr_grid*3
        CALL MPI_Recv(G,datasize,MPI_Real8, 0, 1, MPI_COMM_WORLD,stat, ierr)     
        CALL MPI_Recv(IBLANK,Nz_grid*Nr_grid,MPI_Integer, 0, 1, MPI_COMM_WORLD,stat, ierr)

        ENDIF
        
       ENDDO

	PRINT*, "DONE!"
	
	
	!CALL StopTheCode
	

END SUBROUTINE GridRead


SUBROUTINE AssignGridCoordinates

USE ModuleVariables
IMPLICIT NONE

	DO i = 1, Nz
	   DO j = 1, Nr
	      DO k = 1, 1

		i_indx = istart + (i-1)
		j_indx = jstart + (j-1)
		
		zvec(i,j,k) = G(1,i_indx,j_indx,k,1)
		rvec(i,j,k) = G(1,i_indx,j_indx,k,2)
		ibval(i,j,k) = IBLANK(1,i_indx,j_indx,k)		

	      ENDDO
	   ENDDO
	ENDDO

END SUBROUTINE AssignGridCoordinates

SUBROUTINE SolnRead

USE ModuleVariables
USE ModPLOT3D_IO
IMPLICIT NONE

INTEGER(KIND=8) :: datasize

	
	IF(rank .eq. 0)THEN	
	CALL Read_Soln(NDIM, ngrid, ND, X_target, tau, prec, gf, vf, TARGET_FILE, OF)      	
	CALL Read_Soln(NDIM, ngrid, ND, X, tau, prec, gf, vf, INITIAL_CONDITION, OF)
	PRINT*, "DONE READING SOLUTION"
	ENDIF


	DO procrank = 1, numtasks-1

        IF(rank .eq. 0)THEN

	datasize = ND(1,1)*ND(1,2)*5
	CALL MPI_Send(tau(1), 1, MPI_Real8, procrank, 1, MPI_COMM_WORLD, ierr)
	CALL MPI_Send(tau(4), 1, MPI_Real8, procrank, 1, MPI_COMM_WORLD, ierr)
        CALL MPI_Send(X, datasize, MPI_Real8, procrank, 1, MPI_COMM_WORLD, ierr)
	CALL MPI_Send(X_target, datasize, MPI_Real8, procrank, 1, MPI_COMM_WORLD, ierr)

        ELSEIF(rank .eq. procrank)THEN

        datasize = ND(1,1)*ND(1,2)*5
        ALLOCATE(X(1,ND(1,1),ND(1,2),ND(1,3),5))
	ALLOCATE(X_target(1,ND(1,1),ND(1,2),ND(1,3),5))	
	CALL MPI_Recv(tau(1), 1, MPI_Real8, 0, 1, MPI_COMM_WORLD,stat, ierr)
	CALL MPI_Recv(tau(4), 1, MPI_Real8, 0, 1, MPI_COMM_WORLD,stat, ierr)
        CALL MPI_Recv(X, datasize, MPI_Real8, 0, 1, MPI_COMM_WORLD,stat, ierr)
        CALL MPI_Recv(X_target, datasize, MPI_Real8, 0, 1, MPI_COMM_WORLD,stat, ierr)

        ENDIF

       ENDDO



	DO i = 1, Nz
	   DO j = 1, Nr
	      DO k = 1, 1	
		i_indx = istart + (i-1)
		j_indx = jstart + (j-1)	        
	        DO eqn =  1, neqns
		   pv(i,j,k,eqn) = X(1,i_indx,j_indx,k,eqn)	
		   pv_target(i,j,k,eqn) = X_target(1,i_indx,j_indx,k,eqn)
		  		
		ENDDO
	      ENDDO
	   ENDDO
	ENDDO

	
END SUBROUTINE SolnRead

SUBROUTINE GatherStartEndIndices

USE ModuleVariables

IMPLICIT NONE
	
	ALLOCATE(startenddata(0:numtasks-1,5))

	IF(rank .eq. 0)THEN
	startenddata(0,1) = rank
	startenddata(0,2) = istart
	startenddata(0,3) = iend
	startenddata(0,4) = jstart
	startenddata(0,5) = jend
	ENDIF
	
	
       DO procrank = 1, numtasks-1
	
	IF(rank .eq. procrank)THEN	
	dummyvec(1) = rank
	dummyvec(2) = istart
	dummyvec(3) = iend
	dummyvec(4) = jstart
	dummyvec(5) = jend
	CALL MPI_Send(dummyvec, 5, MPI_Integer8, 0, 1, MPI_COMM_WORLD,ierr)
	ELSEIF(rank .eq. 0)THEN
	CALL MPI_Recv(dummyvec, 5, MPI_Integer8, procrank, 1, MPI_COMM_WORLD, stat, ierr)	
	startenddata(procrank,:) = dummyvec
	ENDIF

	!IF(rank .eq. 0)THEN
	!write(*,*), "Received data", dummyvec(1), dummyvec(2), dummyvec(3), dummyvec(4), dummyvec(5)
	!ENDIF

        ENDDO

	!IF(rank .eq. 0)THEN
	!DO procrank = 1, numtasks-1
	!PRINT*, startenddata(procrank,1:5)
	!ENDDO
	!ENDIF
	
	

30 format('rank= ',I0,'startend = ',I0,I0,I0,I0)


END SUBROUTINE GatherStartEndIndices


SUBROUTINE WriteOutput
USE ModuleVariables
USE ModPLOT3D_IO
IMPLICIT NONE

INTEGER(KIND=8) :: i_indx_start, i_indx_end, j_indx_start, j_indx_end, datasize
	    	
	IF(rank .eq. 0)THEN		
	ALLOCATE(Xout(1,Nz_grid,Nr_grid,1,5))  
	Xout = 0.0D0
	tau(1) = iteration
	IF(iteration .ne. iteration_start)THEN		
	tau(4) = time!time_start + (iteration-iteration_start)*delt
	ENDIF
	PRINT*, "TOTAL TIME=", tau(1), iteration, iteration_start, tau(4)

	ENDIF
	
	DO procrank = 1, numtasks-1

	  IF(rank .eq. procrank)THEN		   	
    	    CALL MPI_Send(pv(1:Nz,1:Nr,1,1:5), Nz*Nr*5, MPI_Real8, 0, 1, MPI_COMM_WORLD,ierr)	   
 	  ELSEIF(rank .eq. 0)THEN
	    i_indx_start = startenddata(procrank,2)
	    i_indx_end   = startenddata(procrank,3)
	    j_indx_start = startenddata(procrank,4)
	    j_indx_end   = startenddata(procrank,5)	
	    datasize = (i_indx_end-i_indx_start+1)*(j_indx_end-j_indx_start+1)*5		
	    CALL MPI_Recv(Xout(1,i_indx_start:i_indx_end,j_indx_start:j_indx_end,1,1:5), datasize, MPI_Real8, procrank, 1, MPI_COMM_WORLD, stat, ierr)	
 	  ENDIF	
	
	ENDDO	
		
	IF(rank .eq. 0)THEN

	PRINT*, "GOING TO WRITE SOLUTION", iteration, tau(1)
	
	DO i = 1, Nz
	   DO j = 1, Nr
	     DO k = 1, 1
		DO eqn  = 1, neqns
		 Xout(1,i,j,k,eqn) = pv(i,j,k,eqn)
		ENDDO		

	     ENDDO
	   ENDDO
	ENDDO 

	write (Solution_File_Name,'(A10,I8.8,A2)') 'RocFlo-CM.', iteration, SOLN_FILETYPE
        Call Write_Soln(NDIM, ngrid, ND, Xout, tau, prec, gf, vf, Solution_File_Name)

	PRINT*, "Density min/max",minval(Xout(1,:,:,:,1)), maxval(X(1,:,:,:,1))
		
	DEALLOCATE(Xout)
	
	
	ENDIF

	CALL MPI_Barrier(MPI_COMM_WORLD,ierr)



END SUBROUTINE WriteOutput


SUBROUTINE WriteRHS
USE ModuleVariables
USE ModPLOT3D_IO
IMPLICIT NONE

INTEGER(KIND=8) :: i_indx_start, i_indx_end, j_indx_start, j_indx_end, datasize
	    	
	IF(rank .eq. 0)THEN		
	ALLOCATE(Xout(1,Nz_grid,Nr_grid,1,5))  
	Xout = 0.0D0
	tau(1) = iteration
	tau(4) = iteration*delt
	PRINT*, "TOTAL TIME=", tau(4)

	ENDIF
	
	DO procrank = 1, numtasks-1

	  IF(rank .eq. procrank)THEN		   	
    	    CALL MPI_Send(rhs(1:Nz,1:Nr,1,1:5), Nz*Nr*5, MPI_Real8, 0, 1, MPI_COMM_WORLD,ierr)	   
 	  ELSEIF(rank .eq. 0)THEN
	    i_indx_start = startenddata(procrank,2)
	    i_indx_end   = startenddata(procrank,3)
	    j_indx_start = startenddata(procrank,4)
	    j_indx_end   = startenddata(procrank,5)	
	    datasize = (i_indx_end-i_indx_start+1)*(j_indx_end-j_indx_start+1)*5		
	    CALL MPI_Recv(Xout(1,i_indx_start:i_indx_end,j_indx_start:j_indx_end,1,1:5), datasize, MPI_Real8, procrank, 1, MPI_COMM_WORLD, stat, ierr)	
 	  ENDIF	
	
	ENDDO
		
	IF(rank .eq. 0)THEN
	
	DO i = 1, Nz
	   DO j = 1, Nr
	     DO k = 1, 1
		DO eqn  = 1, neqns
		Xout(1,i,j,k,eqn) = rhs(i,j,k,eqn)
		ENDDO		
	     ENDDO
	   ENDDO
	ENDDO 

	write (Solution_File_Name,'(A10,I8.8,A2)') 'RocFlo-CM.', iteration, RHS_SOLN_FILETYPE 	
        Call Write_Soln(NDIM, ngrid, ND, Xout, tau, prec, gf, vf, Solution_File_Name)
		
	DEALLOCATE(Xout)
	
	ENDIF



END SUBROUTINE WriteRHS


SUBROUTINE WriteVariableOutput(var1,var2,var3,var4,var5)
USE ModuleVariables
USE ModPLOT3D_IO
IMPLICIT NONE

INTEGER(KIND=8) :: i_indx_start, i_indx_end, j_indx_start, j_indx_end, datasize
REAL(KIND=8) :: var1(-3:Nz+4,-3:Nr+4,1), var2(-3:Nz+4,-3:Nr+4,1), var3(-3:Nz+4,-3:Nr+4,1), var4(-3:Nz+4,-3:Nr+4,1), var5(-3:Nz+4,-3:Nr+4,1)
	    	
	IF(rank .eq. 0)THEN		
	ALLOCATE(Xout(1,Nz_grid,Nr_grid,1,5))  
	Xout = 0.0D0
	tau(1) = iteration
	IF(iteration .ne. iteration_start)THEN		
	tau(4) = time !time_start + (iteration-iteration_start)*delt
	ENDIF	
	ENDIF
	
	
	DO procrank = 1, numtasks-1

	  IF(rank .eq. procrank)THEN		   	
    	    CALL MPI_Send(var1(1:Nz,1:Nr,1), Nz*Nr, MPI_Real8, 0, 1, MPI_COMM_WORLD,ierr)
	    CALL MPI_Send(var2(1:Nz,1:Nr,1), Nz*Nr, MPI_Real8, 0, 1, MPI_COMM_WORLD,ierr)
	    CALL MPI_Send(var3(1:Nz,1:Nr,1), Nz*Nr, MPI_Real8, 0, 1, MPI_COMM_WORLD,ierr)
	    CALL MPI_Send(var4(1:Nz,1:Nr,1), Nz*Nr, MPI_Real8, 0, 1, MPI_COMM_WORLD,ierr)
	    CALL MPI_Send(var5(1:Nz,1:Nr,1), Nz*Nr, MPI_Real8, 0, 1, MPI_COMM_WORLD,ierr)			   
 	  ELSEIF(rank .eq. 0)THEN
	    i_indx_start = startenddata(procrank,2)
	    i_indx_end   = startenddata(procrank,3)
	    j_indx_start = startenddata(procrank,4)
	    j_indx_end   = startenddata(procrank,5)	
	    datasize = (i_indx_end-i_indx_start+1)*(j_indx_end-j_indx_start+1)		
	    CALL MPI_Recv(Xout(1,i_indx_start:i_indx_end,j_indx_start:j_indx_end,1,1), datasize, MPI_Real8, procrank, 1, MPI_COMM_WORLD, stat, ierr)	
    	    CALL MPI_Recv(Xout(1,i_indx_start:i_indx_end,j_indx_start:j_indx_end,1,2), datasize, MPI_Real8, procrank, 1, MPI_COMM_WORLD, stat, ierr)	
       	    CALL MPI_Recv(Xout(1,i_indx_start:i_indx_end,j_indx_start:j_indx_end,1,3), datasize, MPI_Real8, procrank, 1, MPI_COMM_WORLD, stat, ierr)	
    	    CALL MPI_Recv(Xout(1,i_indx_start:i_indx_end,j_indx_start:j_indx_end,1,4), datasize, MPI_Real8, procrank, 1, MPI_COMM_WORLD, stat, ierr)	
    	    CALL MPI_Recv(Xout(1,i_indx_start:i_indx_end,j_indx_start:j_indx_end,1,5), datasize, MPI_Real8, procrank, 1, MPI_COMM_WORLD, stat, ierr)		

 	  ENDIF	
	
	ENDDO
		
	IF(rank .eq. 0)THEN
	
	   DO i = 1, Nz
	      DO j = 1, Nr
	        DO k = 1, 1
		 Xout(1,i,j,k,1) = var1(i,j,k)
   	         Xout(1,i,j,k,2) = var2(i,j,k)
		 Xout(1,i,j,k,3) = var3(i,j,k)
		 Xout(1,i,j,k,4) = var4(i,j,k)
		 Xout(1,i,j,k,5) = var5(i,j,k) 
		ENDDO		
	     ENDDO
	   ENDDO

	 
	  write (Solution_File_Name,'(A10,I8.8,A2)') 'RocFlo-CM.', iteration, ART_COEFF_SOLN_FILETYPE
          Call Write_Soln(NDIM, ngrid, ND, Xout, tau, prec, gf, vf, Solution_File_Name)

	  DEALLOCATE(Xout)

	ENDIF

	CALL MPI_Barrier(MPI_COMM_WORLD,ierr)

END SUBROUTINE WriteVariableOutput

SUBROUTINE ComputeTimeStep
USE ModuleVariables
USE ModDecomp
IMPLICIT NONE


	DO i = 1, Nz	
	  DO j = 1, Nr
	    DO k = 1, 1
	
		IF(ibval(i,j,k) .eq. 0)THEN
		deltvec(i,j,k) = 100.0D0
		CYCLE
		ENDIF
	
		delt_inv_term = ABS( Vzvec(i,j,k)*dxidz(i,j,k)+Vrvec(i,j,k)*detadr(i,j,k) ) + &
				DSQRT(T)*DSQRT( dxidz(i,j,k)**2+detadr(i,j,k)**2 )

		delt_visc_term = 2.0D0*MAX( muvec(i,j,k)/Re,lambdavec(i,j,k)/Re,kappavec(i,j,k)/(Re*Pr) )*&
				 ( ABS(dxidz(i,j,k)) + ABS(detadr(i,j,k)) )**2

		delt_denom = delt_inv_term + delt_visc_term

		deltvec(i,j,k) = CFL/delt_denom

	    ENDDO
	  ENDDO
	ENDDO

	min_delt = MINVAL(deltvec)

	!PRINT*, rank, "Mindelt=", MINVAL(deltvec)

	Call MPI_AllReduce(min_delt, delt, 1, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, ierr)

	IF(rank .eq. 0)THEN
	IF(delt .eq. 100.0D0)THEN
	PRINT*,"Time step blow up"
	CALL StopTheCode
	ENDIF
	ENDIF
	

	!PRINT*, rank, delt

	



END SUBROUTINE ComputeTimeStep

END MODULE ModInputOutput
