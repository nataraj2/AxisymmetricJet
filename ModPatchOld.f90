MODULE ModPatch

CONTAINS

SUBROUTINE AllocatePatches

USE module_variables
USE ModDecomp
IMPLICIT NONE	

	DO i = 1, Nz
	   DO j = 1, Nr
	      DO k = 1, 1

		is_it_i_patch(i,j,k) = 0
		is_it_j_patch(i,j,k) = 0 
		which_i_patch(i,j,k) = 0 
		which_j_patch(i,j,k) = 0 

	    ENDDO
	  ENDDO
	ENDDO
	
	total_no_bcs = 0
	
	reason = 0

	bcno = 0

	OPEN(UNIT=10,FILE=bc_file)	

	DO WHILE (abs(reason) .eq. 0)	
	  READ(10,*,IOSTAT=reason)dummy_int, dummy_int, dummy_int, dummy_int, dummy_int, dummy_int
	  total_no_bcs = total_no_bcs + 1		 	 
	ENDDO

	CLOSE(10)

	total_no_bcs = total_no_bcs - 1
	
	no_bcs = 0

	OPEN(UNIT=10,FILE=bc_file)	

	DO bcno = 1, total_no_bcs
	  
	  READ(10,*)dummy_int, bctype, bcdir, bcistart, bciend, bcjstart, bcjend, dummy_int, dummy_int	 
	  	  !PRINT*, bcno, bcdir	
		IF(bcistart .eq. -1)THEN
			bcistart = Nz_grid
		ENDIF
		IF(bcjstart .eq. -1)THEN
			bcjstart = Nr_grid
		ENDIF

		IF(bciend .eq. -1)THEN
			bciend = Nz_grid
		ENDIF
		IF(bcjend .eq. -1)THEN
			bcjend = Nr_grid
		ENDIF
	


       IF(abs(bcdir) .eq. 1)THEN	  	  	
	  IF(bcistart .ge. istart .and. bcistart .le. iend)THEN
              IF(bcjstart .ge. jstart .and. bcjstart .le. jend)THEN 		
	        IF(bcjend .ge. jstart .and. bcjend .le. jend)THEN 
		  PRINT*, rank, bcno			  		 
		  no_bcs = no_bcs+1
	        ELSE 	
		   PRINT*, rank, bcno			
		  no_bcs = no_bcs+1
		ENDIF
	      ELSEIF(bcjend .ge. jstart .and. bcjend .le. jend)THEN
		   PRINT*, rank, bcno			
		  no_bcs = no_bcs + 1
	      ENDIF
	   ENDIF	 
	   IF(istart .eq. bcistart .or. iend .eq. bcistart)THEN
	    IF(jstart .gt. bcjstart .and. jend .lt. bcjend)THEN
	  	PRINT*, rank, bcno			  		 
		no_bcs = no_bcs+1
	    ENDIF
	   ENDIF
	    
        ENDIF
		
	 IF(abs(bcdir) .eq. 2)THEN	  
	  IF(bcjstart .ge. jstart .and. bcjstart .le. jend)THEN	      
              IF(bcistart .ge. istart .and. bcistart .le. iend)THEN 		
	        IF(bciend .ge. istart .and. bciend .le. iend)THEN 	
		   PRINT*, rank, bcno				
		  no_bcs = no_bcs+1
	        ELSE 	
		   PRINT*, rank, bcno		
		  no_bcs = no_bcs+1
		ENDIF
	      ELSEIF(bciend .ge. istart .and. bciend .le. iend)THEN
		  PRINT*, rank, bcno		
		  no_bcs = no_bcs + 1
	      ENDIF
	   ENDIF
	   
	   IF(jstart .eq. bcjstart .or. jend .eq. bcjstart)THEN
	    IF(istart .gt. bcistart .and. iend .lt. bciend)THEN
		PRINT*, rank, bcno			  		 
		no_bcs = no_bcs+1
	    ENDIF
	   ENDIF	
		

	 ENDIF	 	 

	 

	ENDDO

	CLOSE(10)

	!PRINT*, rank, no_bcs

	CALL MPI_Barrier(MPI_COMM_WORLD,ierr)

	!CALL StopTheCode

	ALLOCATE(patch(no_bcs))

END SUBROUTINE AllocatePatches

SUBROUTINE CreateAllPatches

USE module_variables
USE ModDecomp

IMPLICIT NONE

	CALL AllocatePatches

	CALL MPI_Barrier(MPI_COMM_WORLD,ierr)

	!PRINT*, "Rank and bcs", rank, no_bcs	

	OPEN(UNIT=10,FILE=bc_file)	

	bcno = 1

	DO i = 1, total_no_bcs

		READ(10,*)dummy_int, bctype, bcdir, bcistart, bciend, bcjstart, bcjend, dummy_int, dummy_int

	
		IF(bcistart .eq. -1)THEN
			bcistart = Nz_grid
		ENDIF
		IF(bcjstart .eq. -1)THEN
			bcjstart = Nr_grid
		ENDIF

		IF(bciend .eq. -1)THEN
			bciend = Nz_grid
		ENDIF
		IF(bcjend .eq. -1)THEN
			bcjend = Nr_grid
		ENDIF

       IF(abs(bcdir) .eq. 1)THEN	  	  	
	  IF(bcistart .ge. istart .and. bcistart .le. iend)THEN
              IF(bcjstart .ge. jstart .and. bcjstart .le. jend)THEN 		
	        IF(bcjend .ge. jstart .and. bcjend .le. jend)THEN 
		  procistart = bcistart
		  prociend = bcistart
		  procjstart = bcjstart
		  procjend = bcjend		  
		  CALL PatchCreation
		  bcno = bcno + 1
	        ELSE 	
		  procistart = bcistart
		  prociend = bcistart
		  procjstart = bcjstart
		  procjend = jend	 
		  CALL PatchCreation
		  bcno = bcno + 1
		ENDIF
	      ELSEIF(bcjend .ge. jstart .and. bcjend .le. jend)THEN
		  procistart = bcistart
		  prociend = bcistart
		  procjstart = jstart
		  procjend = bcjend
		  CALL PatchCreation
		   bcno = bcno + 1
	      ENDIF
	   ENDIF	 
	   IF(istart .eq. bcistart .or. iend .eq. bcistart)THEN
	    IF(jstart .gt. bcjstart .and. jend .lt. bcjend)THEN
	  	  procistart = bcistart
		  prociend = bcistart
		  procjstart = jstart
		  procjend = jend
		  CALL PatchCreation
		  bcno = bcno + 1
	    ENDIF
	   ENDIF
	    
        ENDIF
		
	 IF(abs(bcdir) .eq. 2)THEN	  
	  IF(bcjstart .ge. jstart .and. bcjstart .le. jend)THEN	      
              IF(bcistart .ge. istart .and. bcistart .le. iend)THEN 		
	        IF(bciend .ge. istart .and. bciend .le. iend)THEN 	
		  procistart = bcistart
		  prociend = bciend
		  procjstart = bcjstart
		  procjend = bcjstart	
		  CALL PatchCreation
		  bcno = bcno + 1
	        ELSE 	
		  procistart = bcistart
		  prociend = iend
		  procjstart = bcjstart
		  procjend = bcjstart
		  CALL PatchCreation
		  bcno = bcno + 1
		ENDIF
	      ELSEIF(bciend .ge. istart .and. bciend .le. iend)THEN
		  procistart = istart
		  prociend = bciend
		  procjstart = bcjstart
		  procjend = bcjstart
		  CALL PatchCreation
		  bcno = bcno + 1
	      ENDIF
	   ENDIF
	   
	   IF(jstart .eq. bcjstart .or. jend .eq. bcjstart)THEN
	    IF(istart .gt. bcistart .and. iend .lt. bciend)THEN
		  procistart = istart
		  prociend = iend
		  procjstart = bcjstart
		  procjend = bcjstart
		  CALL PatchCreation
		  bcno = bcno + 1
	    ENDIF
	   ENDIF	
		

	 ENDIF	 	 

	ENDDO

	CLOSE(10)	



	WRITE(sol_file,101)rank
        101 FORMAT("patch_data.",I3.3)

	OPEN(UNIT=10,FILE=sol_file)
	DO bcno = 1, no_bcs

		WRITE(10,1)bcno, patch(bcno)%bc_dir, patch(bcno)%i_start, patch(bcno)%i_end, patch(bcno)%j_start, patch(bcno)%j_end	

	ENDDO
	CLOSE(10)


1 FORMAT(6(I))


END SUBROUTINE CreateAllPatches

SUBROUTINE PatchCreation

USE module_variables
IMPLICIT NONE

! bc, bc_dir, i_start, i_end, j_start, j_end		

		patch(bcno)%bc = bctype
		patch(bcno)%bc_dir = bcdir

		istartval = procistart - istart + 1		
		iendval = prociend - istart + 1	
		jstartval = procjstart - jstart + 1
		jendval = procjend - jstart + 1

		
		IF(bcdir .eq. 1)THEN
		patch(bcno)%i_start = istartval
		patch(bcno)%i_end = istartval+3
		patch(bcno)%j_start = jstartval
		patch(bcno)%j_end = jendval
		ENDIF
		IF(bcdir .eq. -1)THEN
		patch(bcno)%i_start = iendval-3
		patch(bcno)%i_end = iendval
		patch(bcno)%j_start = jstartval
		patch(bcno)%j_end = jendval
		ENDIF
		IF(bcdir .eq. 2)THEN
		patch(bcno)%i_start = istartval
		patch(bcno)%i_end = iendval
		patch(bcno)%j_start = jstartval
		patch(bcno)%j_end = jstartval+3
		ENDIF
		IF(bcdir .eq. -2)THEN
		patch(bcno)%i_start = istartval
		patch(bcno)%i_end = iendval
		patch(bcno)%j_start = jendval-3
		patch(bcno)%j_end = jendval
		ENDIF

	



END SUBROUTINE PatchCreation



END MODULE ModPatch
