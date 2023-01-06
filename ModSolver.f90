PROGRAM ModSolver

USE ModuleVariables
USE ModInputOutput
USE ModAllocate
USE ModAssign
USE ModDecomp
USE ModCommunicate
USE ModDeriv
USE ModTGfilter
USE ModOperator
USE ModRungeKutta
USE ModPatch
USE ModPlot3D_IO 
USE ModShockCapturing

IMPLICIT NONE
	
	CALL MPI_INIT(ierr)   
        CALL MPI_COMM_SIZE(MPI_COMM_WORLD, numtasks, ierr)        
        CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
	
	CALL GridRead	

	CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
	CALL CreateDecomposition		
	CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

	CALL DistributeGridDims	

	CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)	

	CALL DistributeGridPoints	
	CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)	

	CALL GatherStartEndIndices		
	CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)	


	CALL AllocateVariables

	CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
		if(rank==0)then
	print*, "Reached here "
	endif

	
	CALL AssignGridCoordinates		
	CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

	!CALL SolnRead	        	

	CALL CreateAllPatches	

	CALL SBP_operator

	CALL ProcessorCommunication(zvec)
	CALL ComputeMetric(zvec,1,dzdxi)

	CALL ProcessorCommunication(rvec)
	CALL ComputeMetric(rvec,2,drdeta)

	DO i = 1, Nz
	   DO j = 1, Nr
	      DO k = 1, 1
		dxidz(i,j,k) = 1.0D0/dzdxi(i,j,k)
		detadr(i,j,k) = 1.0D0/drdeta(i,j,k)
	      ENDDO	
	   ENDDO
	ENDDO

	CALL ProcessorCommunication(dxidz)
	CALL ComputeMetric(dxidz,1,ddxidxidz)
	
	CALL ProcessorCommunication(detadr)
	CALL ComputeMetric(detadr,2,ddetadetadr)
	
	CALL SolnRead

	CALL AssignValues		

	!CALL ProcessorCommunication(rhovec)

	!!!!!!!!! Check the filtering !!!!!!!!!!!!!!!
	
	!rhovecTGfilt = rhovec

	!CALL ProcessorCommunication(rhovecTGfilt)
	
	!DO dummy_int = 1, 10		
        !CALL PerformTGfilter(rhovecTGfilt,1,rhovecTGfilt)	
	!CALL ProcessorCommunication(rhovecTGfilt)
	!ENDDO
	
	!CALL WriteVariableOutput(rhovecTGfilt,rhovec,rhovec,rhovec,rhovec)

	!CALL StopTheCode

	!!!!!!!! Filtering checking ends here !!!!!!!!!!!	

	! Write the initial condition file		
	
	iteration = tau(1)	
	iteration_start = tau(1)	
	time_start = tau(4)
	time = time_start

	CALL WriteOutput

	! Determine the Time Step

	CALL ComputeTimeStep

	!CALL StopTheCode
				
	no_iterations = 10000

	!PRINT*, "iteration start=", rank, iteration_start, time_start
	!CALL StopTheCode
	

	DO iteration = iteration_start+1, iteration_start+no_iterations
	CALL ComputeTimeStep
	CALL RK4
	ENDDO


	CALL StopTheCode


END PROGRAM ModSolver


