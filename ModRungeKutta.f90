MODULE ModRungeKutta

CONTAINS

SUBROUTINE RK4
USE ModuleVariables
USE ModInputOutput
USE ModAssign
USE ModAxis
USE ModRHS
USE ModPlot3D_IO

IMPLICIT NONE

      DO eqn = 1, neqns
	DO i = 1, Nz
	   DO j = 1, Nr
	      DO k = 1, 1
		pv_old(i,j,k,eqn) = pv(i,j,k,eqn)				
	      ENDDO
	   ENDDO	
	 ENDDO
	ENDDO
	
	! RK Step 1

	rkStep = 1
	
	 CALL AssignValues
         CALL ComputeRHS
	
	DO eqn = 1, neqns
	 DO i = 1, Nz
	   DO j = 1, Nr
	      DO k = 1, 1
		k1(i,j,k,eqn) = delt*rhs(i,j,k,eqn)				
	      ENDDO
	   ENDDO	
	 ENDDO
	ENDDO

	
	DO eqn = 1, neqns
	 DO i = 1, Nz
	   DO j = 1, Nr
	      DO k = 1, 1
		pv(i,j,k,eqn) = pv_old(i,j,k,eqn) + 0.5D0*k1(i,j,k,eqn)				
	      ENDDO
	   ENDDO	
	 ENDDO
	ENDDO

	CALL axis_conditions

	! RK Step 2

	rkStep = 2
	
	 CALL AssignValues
         CALL ComputeRHS

       DO eqn = 1, neqns
	 DO i = 1, Nz
	   DO j = 1, Nr
	      DO k = 1, 1
		k2(i,j,k,eqn) = delt*rhs(i,j,k,eqn)				
	      ENDDO
	    ENDDO	
	  ENDDO
	ENDDO

	DO eqn = 1, neqns	
	 DO i = 1, Nz
	   DO j = 1, Nr
	      DO k = 1, 1
		pv(i,j,k,eqn) = pv_old(i,j,k,eqn) + 0.5D0*k2(i,j,k,eqn)			
	      ENDDO
	   ENDDO	
	 ENDDO
	ENDDO

	CALL axis_conditions

	! RK Step 3

	rkStep = 3

	 CALL AssignValues
         CALL ComputeRHS

       DO eqn = 1, neqns	
	DO i = 1, Nz
	   DO j = 1, Nr
	      DO k = 1, 1
		k3(i,j,k,eqn) = delt*rhs(i,j,k,eqn)				
	      ENDDO
	   ENDDO	
	 ENDDO
	ENDDO
	

	DO eqn = 1, neqns
	 DO i = 1, Nz
	   DO j = 1, Nr
	      DO k = 1, 1
		pv(i,j,k,eqn) = pv_old(i,j,k,eqn) + k3(i,j,k,eqn)				
	      ENDDO
	   ENDDO	
	 ENDDO
	ENDDO

	CALL axis_conditions

	! RK Step 4

	rkStep = 4

	 CALL AssignValues
         CALL ComputeRHS

       DO eqn = 1, neqns
	DO i = 1, Nz
	   DO j = 1, Nr
	      DO k = 1, 1
		k4(i,j,k,eqn) = delt*rhs(i,j,k,eqn)				
	      ENDDO
	   ENDDO	
	 ENDDO
	ENDDO
	

	! Final update

       DO eqn = 1, neqns
	DO i = 1, Nz
	   DO j = 1, Nr
	      DO k = 1, 1
     
                IF(ibval(i,j,k) .eq. 0)CYCLE

		pv(i,j,k,eqn) =  pv_old(i,j,k,eqn) + 1.0D0/6.0D0*k1(i,j,k,eqn) + 1.0D0/3.0D0*(k2(i,j,k,eqn)+k3(i,j,k,eqn)) + 1.0D0/6.0D0*k4(i,j,k,eqn)
				
	      ENDDO
	   ENDDO	
	 ENDDO
	ENDDO

	CALL axis_conditions	

	CALL AssignValues	

	IF(rank .eq. 0)THEN
	time = time + delt
	WRITE(*,201) iteration, delt, time
	ENDIF	

201 FORMAT("Iteration=",I0,2X, "dt=",E17.9, 2X, "time=",E17.9)



	IF(mod(iteration,noutput) .eq. 0)THEN	
	IF(rank .eq. 0)PRINT*,"Calling output"
	CALL WriteOutput
	CALL MPI_Barrier(MPI_COMM_WORLD,ierr)
	ENDIF

	

	

END SUBROUTINE RK4

END MODULE ModRungeKutta
