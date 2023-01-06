MODULE ModCommunicate
CONTAINS

SUBROUTINE ProcessorCommunication
USE module_variables
USE ModDecomp
USE ModPLOT3D_IO

IMPLICIT NONE

INTEGER(KIND=8) :: i_indx_start, i_indx_end, j_indx_start, j_indx_end	


	DO i = 1, 4

	IF(nbrs(i) .ge. 0)THEN
	
	IF(i .eq. 1)THEN
	CALL MPI_Send(pv(1:Nz,Nr-1:Nr,1,1:5), 5*Nz, MPI_Real8, nbrs(i), 1, MPI_COMM_WORLD,ierr)
	ELSEIF(i .eq. 2)THEN
	CALL MPI_Send(pv(1:Nz,1:2,1,1:5), 5*Nz, MPI_Real8, nbrs(i), 1, MPI_COMM_WORLD,ierr)
	ELSEIF(i .eq. 3)THEN
	CALL MPI_Send(pv(1:2,1:Nr,1,1:5), 5*Nr, MPI_Real8, nbrs(i), 1, MPI_COMM_WORLD,ierr)
	ELSEIF(i .eq. 4)THEN
	CALL MPI_Send(pv(Nz-1:Nz,1:Nr,1,1:5), 5*Nr, MPI_Real8, nbrs(i), 1, MPI_COMM_WORLD,ierr)	
	ENDIF

	ENDIF

	ENDDO

	DO i = 1, 4

	IF(nbrs(i) .ge. 0)THEN

	IF(i .eq. 1)THEN
	CALL MPI_Recv(pv(1:Nz,Nr+1:Nr+2,1,1:5), 5*Nz, MPI_Real8, nbrs(i), 1, MPI_COMM_WORLD, stat, ierr)		
	ELSEIF(i .eq. 2)THEN
	CALL MPI_Recv(pv(1:Nz,-1:0,1,1:5), 5*Nz, MPI_Real8, nbrs(i), 1, MPI_COMM_WORLD, stat, ierr)		
	ELSEIF(i .eq. 3)THEN
	CALL MPI_Recv(pv(-1:0,1:Nr,1,1:5), 5*Nr, MPI_Real8, nbrs(i), 1, MPI_COMM_WORLD, stat, ierr)		
	ELSEIF(i .eq. 4)THEN
	CALL MPI_Recv(pv(Nz+1:Nz+2,1:Nr,1,1:5), 5*Nr, MPI_Real8, nbrs(i), 1, MPI_COMM_WORLD, stat, ierr)		
	ENDIF

	ENDIF

	ENDDO


	
	IF(rank .eq. 5)THEN
	ALLOCATE(Xout(1,Nz_grid,Nr_grid,1,5))

	i_indx_start = -1
	i_indx_end = Nz+2
	j_indx_start = -1
	j_indx_end = Nr+2

	IF(nbrs(1) .lt. 0)THEN	
	 j_indx_end = Nr
	ENDIF
	IF(nbrs(2) .lt. 0)THEN	
	 j_indx_start = 1
	ENDIF
	IF(nbrs(3) .lt. 0)THEN	
	 i_indx_start = 1
	ENDIF
	IF(nbrs(4) .lt. 0)THEN	
	 i_indx_end = Nz
	ENDIF
	

	DO i = i_indx_start, i_indx_end
	   DO j = j_indx_start, j_indx_end
	     DO k = 1, 1
		DO eqn  = 1, neqns		 
		 Xout(1,istart+i-1,jstart+j-1,k,eqn) = pv(i,j,k,eqn)
		ENDDO		
	     ENDDO
	   ENDDO
	ENDDO 

	tau(4) = 1.0D0
	iteration = 3
	write (Solution_File_Name,'(A10,I8.8,A2)') 'RocFlo-CM.', iteration, SOLN_FILETYPE
        Call Write_Soln(NDIM, ngrid, ND, Xout, tau, prec, gf, vf, Solution_File_Name)

	ENDIF
	



END SUBROUTINE ProcessorCommunication



END MODULE ModCommunicate


