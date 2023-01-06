MODULE ModCommunicate
CONTAINS

SUBROUTINE ProcessorCommunication(var)
USE ModuleVariables
USE ModDecomp
USE ModPLOT3D_IO

IMPLICIT NONE

INTEGER(KIND=8) :: i_indx_start, i_indx_end, j_indx_start, j_indx_end	
REAL(KIND=8) :: var(-3:Nz+4,-3:Nr+4,1)

	DO i = 1, 4

	IF(nbrs(i) .ge. 0)THEN
	
	IF(i .eq. 1)THEN
	CALL MPI_Send(var(1:Nz,Nr-3:Nr,1), 4*Nz, MPI_Real8, nbrs(i), 1, MPI_COMM_WORLD,ierr)
	ELSEIF(i .eq. 2)THEN
	CALL MPI_Send(var(1:Nz,1:4,1), 4*Nz, MPI_Real8, nbrs(i), 1, MPI_COMM_WORLD,ierr)
	ELSEIF(i .eq. 3)THEN
	CALL MPI_Send(var(1:4,1:Nr,1), 4*Nr, MPI_Real8, nbrs(i), 1, MPI_COMM_WORLD,ierr)	
	ELSEIF(i .eq. 4)THEN	
	CALL MPI_Send(var(Nz-3:Nz,1:Nr,1), 4*Nr, MPI_Real8, nbrs(i), 1, MPI_COMM_WORLD,ierr)			
	ENDIF

	ENDIF

	ENDDO

	DO i = 1, 4

	IF(nbrs(i) .ge. 0)THEN

	IF(i .eq. 1)THEN
	CALL MPI_Recv(var(1:Nz,Nr+1:Nr+4,1), 4*Nz, MPI_Real8, nbrs(i), 1, MPI_COMM_WORLD, stat, ierr)		
	ELSEIF(i .eq. 2)THEN
	CALL MPI_Recv(var(1:Nz,-3:0,1), 4*Nz, MPI_Real8, nbrs(i), 1, MPI_COMM_WORLD, stat, ierr)		
	ELSEIF(i .eq. 3)THEN
	CALL MPI_Recv(var(-3:0,1:Nr,1), 4*Nr, MPI_Real8, nbrs(i), 1, MPI_COMM_WORLD, stat, ierr)			
	ELSEIF(i .eq. 4)THEN
	CALL MPI_Recv(var(Nz+1:Nz+4,1:Nr,1), 4*Nr, MPI_Real8, nbrs(i), 1, MPI_COMM_WORLD, stat, ierr)		
	ENDIF

	ENDIF

	ENDDO

	

END SUBROUTINE ProcessorCommunication



END MODULE ModCommunicate


