MODULE ModPatch

CONTAINS

SUBROUTINE AllocatePatches

USE ModuleVariables
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
	  
	  exitvalue = 0	
	
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
	
		DO i = 1, Nz
		   DO j = 1, Nr
			i_indx = istart + i - 1 
			j_indx = jstart + j - 1	
	
			IF(bcdir .eq. 1)THEN	
			  patchistart = bcistart
			  patchiend = bcistart + 3
			  patchjstart = bcjstart
			  patchjend = bcjend
			ELSEIF(bcdir .eq. -1)THEN
			  patchistart = bciend - 3
			  patchiend = bciend
			  patchjstart = bcjstart
			  patchjend = bcjend
			ELSEIF(bcdir .eq. 2)THEN
			  patchistart = bcistart
			  patchiend = bciend
			  patchjstart = bcjstart
			  patchjend = bcjstart + 3
			ELSEIF(bcdir .eq. -2)THEN
			  patchistart = bcistart
			  patchiend = bciend
			  patchjstart = bcjend - 3
			  patchjend = bcjend
			ENDIF
			value1 = (i_indx-patchistart)*(i_indx-patchiend)
			value2 = (j_indx-patchjstart)*(j_indx-patchjend)
			
			IF(value1 .le. 0 .and. value2 .le. 0)THEN
			no_bcs = no_bcs + 1
			exitvalue = 1
			GOTO 100
			ENDIF			
		   ENDDO
		100 IF(exitvalue .eq. 1)EXIT
		ENDDO	

	ENDDO

	CLOSE(10)

	!PRINT*, rank, no_bcs

	CALL MPI_Barrier(MPI_COMM_WORLD,ierr)

	!CALL StopTheCode

	ALLOCATE(patch(no_bcs))

END SUBROUTINE AllocatePatches

SUBROUTINE CreateAllPatches

USE ModuleVariables
USE ModDecomp

IMPLICIT NONE

	CALL AllocatePatches
	
	OPEN(UNIT=20,FILE=bc_file)		

	bcno = 1

	DO bccount = 1, total_no_bcs

		READ(20,*)dummy_int, bctype, bcdir, bcistart, bciend, bcjstart, bcjend, dummy_int, dummy_int

		exitvalue = 0	
		procbcistart = 0 
	        procbciend = 0
	        procbcjstart = 0
		procbcjend = 0

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

		DO i = 1, Nz
		   DO j = 1, Nr
			i_indx = istart + i - 1 
			j_indx = jstart + j - 1	
	
			IF(bcdir .eq. 1)THEN	
			  patchistart = bcistart
			  patchiend = bcistart + 3
			  patchjstart = bcjstart
			  patchjend = bcjend
			ELSEIF(bcdir .eq. -1)THEN
			  patchistart = bciend - 3
			  patchiend = bciend
			  patchjstart = bcjstart
			  patchjend = bcjend
			ELSEIF(bcdir .eq. 2)THEN
			  patchistart = bcistart
			  patchiend = bciend
			  patchjstart = bcjstart
			  patchjend = bcjstart + 3
			ELSEIF(bcdir .eq. -2)THEN
			  patchistart = bcistart
			  patchiend = bciend
			  patchjstart = bcjend - 3
			  patchjend = bcjend
			ENDIF

			value1 = (i_indx-patchistart)*(i_indx-patchiend)
			value2 = (j_indx-patchjstart)*(j_indx-patchjend)

			IF(exitvalue .eq. 0)THEN
			  IF(value1 .le. 0 .and. value2 .le. 0)THEN			
			   procbcistart = i_indx
			   procbcjstart = j_indx			
			   exitvalue = 1
			  ENDIF
			ELSEIF(exitvalue .eq. 1)THEN
			  IF(value1 .le. 0 .and. value2 .le. 0)THEN			
			   procbciend = i_indx
			   procbcjend = j_indx								   
			  ENDIF
			ENDIF						
		   ENDDO	
		ENDDO			
		
		!IF(exitvalue .ne. 0)THEN
		!write(*,201), rank, bcno, procbcistart, procbciend, procbcjstart, procbcjend
		!ENDIF
			
		IF(exitvalue .ne. 0)THEN		
		CALL PatchCreation
		bcno = bcno + 1
		ENDIF

	ENDDO

201 FORMAT(6(I0)2X)
	CLOSE(20)	


      DO patch_no = 1, no_bcs	

	i_start = patch(patch_no)%i_start
	i_end = patch(patch_no)%i_end
	j_start = patch(patch_no)%j_start
	j_end = patch(patch_no)%j_end	
	
         DO i = i_start, i_end
	   DO j = j_start, j_end
	      DO k = 1, 1	    	

	    IF(patch(patch_no)%bc_dir .eq. -1 .or. patch(patch_no)%bc_dir .eq. 1)THEN
	    is_it_i_patch(i,j,k) = 1
	    which_i_patch(i,j,k) = patch_no			    
	    ENDIF
		
	    IF(patch(patch_no)%bc_dir .eq. -2 .or. patch(patch_no)%bc_dir .eq. 2)THEN
	    is_it_j_patch(i,j,k) = 1
	    which_j_patch(i,j,k) = patch_no			    	
	    ENDIF
	    	

	   ENDDO	 
         ENDDO
	ENDDO

	

      ENDDO	




	WRITE(sol_file,101)rank
        101 FORMAT("patch_data.",I3.3)

	OPEN(UNIT=10,FILE=sol_file)
	DO bcno = 1, no_bcs

		WRITE(10,1)bcno, patch(bcno)%bc, patch(bcno)%bc_dir, patch(bcno)%i_start, patch(bcno)%i_end, patch(bcno)%j_start, patch(bcno)%j_end	

	ENDDO
	CLOSE(10)


1 FORMAT(7(I0))


END SUBROUTINE CreateAllPatches

SUBROUTINE PatchCreation

USE ModuleVariables
IMPLICIT NONE

! bc, bc_dir, i_start, i_end, j_start, j_end		

		patch(bcno)%bc = bctype
		patch(bcno)%bc_dir = bcdir

		istartval = procbcistart - istart + 1		
		iendval = procbciend - istart + 1	
		jstartval = procbcjstart - jstart + 1
		jendval = procbcjend - jstart + 1

		
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
