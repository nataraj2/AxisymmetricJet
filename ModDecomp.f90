MODULE ModDecomp

CONTAINS

SUBROUTINE CreateDecomposition

USE ModuleVariables
IMPLICIT NONE
  		
      	 		
      !CALL MPI_INIT(ierr)   
      !CALL MPI_COMM_SIZE(MPI_COMM_WORLD, numtasks, ierr)	
      !call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
      dimsforcartcomm(1) = dims(2)
      dimsforcartcomm(2) = dims(1)  
      call MPI_CART_CREATE(MPI_COMM_WORLD, 2, dimsforcartcomm, periods, reorder, cartcomm, ierr) 	
      call MPI_COMM_RANK(cartcomm, rank, ierr)   
      call MPI_CART_COORDS(cartcomm, rank, 2, coordsforcartcomm, ierr)
      coords(1) = coordsforcartcomm(2)
      coords(2) = coordsforcartcomm(1)	
      call MPI_CART_SHIFT(cartcomm, 0, 1, nbrs(DOWN), nbrs(UP), ierr)
      call MPI_CART_SHIFT(cartcomm, 1, 1, nbrs(LEFT), nbrs(RIGHT), ierr)	      

      !write(*,20) rank,coords(1),coords(2),nbrs(UP),nbrs(DOWN), nbrs(LEFT),nbrs(RIGHT)      	         

!20 format('rank= ',I3,' coords= ',I2,I2,' neighbors(u,d,l,r)= ',I3,I3,I3,I3 )

  !CALL StopTheCode

END SUBROUTINE CreateDecomposition

SUBROUTINE DistributeGridDims
USE ModuleVariables
IMPLICIT NONE

INTEGER(KIND=8), ALLOCATABLE :: NDinproc(:)
INTEGER(KIND=8) :: dir
REAL(KIND=8) :: dummy_real

	ALLOCATE(NDinproc(2))

       DO dir = 1, 2

         IF(mod(ND(1,dir),dims(dir)) .eq. 0)THEN
            NDinproc(dir) = ND(1,dir)/dims(dir)
         ELSE
               
	    dummy_real = ND(1,dir)/dims(dir)        
            IF(coords(dir) .lt. dims(dir)-1)THEN              
            NDinproc(dir) = NINT(dummy_real)	
            ELSEIF(coords(dir) .eq. dims(dir)-1)THEN	    	
            NDinproc(dir) = ND(1,dir) - (dims(dir)-1)*NINT(dummy_real)
            ENDIF

         ENDIF
       ENDDO 
            
	Nz = NDinproc(1)
	Nr = NDinproc(2)

	!PRINT*, rank, Nz, Nr

	!CALL StopTheCode

END SUBROUTINE DistributeGridDims


SUBROUTINE DistributeGridPoints

USE ModuleVariables
IMPLICIT NONE
REAL(KIND=8) :: dummy_real	
	
	    
	    dummy_real = ND(1,1)/dims(1)

	    IF(nbrs(LEFT) .ge. 0)THEN	    
	    istart = coords(1)*NINT(dummy_real) + 1	 
	    iend = istart + Nz - 1
   	    ELSEIF(nbrs(LEFT) .lt. 0)THEN
	    istart = 1
	    iend = istart + Nz - 1
	    ENDIF	    
	    
	    dummy_real = ND(1,2)/dims(2)

	    IF(nbrs(DOWN) .ge. 0)THEN		       
	    jstart = coords(2)*NINT(dummy_real) + 1	 
	    jend = jstart + Nr - 1	        	
   	    ELSEIF(nbrs(DOWN) .lt. 0)THEN
	    jstart = 1
	    jend = jstart + Nr - 1
	    ENDIF
		
	    CALL MPI_Barrier(MPI_COMM_WORLD,ierr)

	    write(*,40) rank, istart, iend, jstart, jend

40 format('rank= ',I0,'startend = ',I0,I0,I0,I0 )

	       CALL MPI_Barrier(MPI_COMM_WORLD,ierr)

	    !CALL StopTheCode	     

END SUBROUTINE DistributeGridPoints



SUBROUTINE StopTheCode

USE ModuleVariables
IMPLICIT NONE

	CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)	
	!CALL MPI_ABORT(MPI_COMM_WORLD,ierr)
	CALL MPI_FINALIZE(ierr)
	STOP

END SUBROUTINE StopTheCode

END MODULE ModDecomp
