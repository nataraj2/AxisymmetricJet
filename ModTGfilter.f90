MODULE ModTGfilter

CONTAINS

SUBROUTINE PerformTGfilter(var,dirn,out,op_flag)

USE ModuleVariables
IMPLICIT NONE

	REAL(KIND=8) :: var(-3:Nz+4,-3:Nr+4,1), out(-3:Nz+4,-3:Nr+4,1)
	INTEGER	:: dirn
	CHARACTER(LEN=8), OPTIONAL :: op_flag

	!out = 0.0D0

      IF(dirn .eq. 1)THEN

      	
	  DO i = 1, Nz	     
	     DO j = 1, Nr
		DO k = 1, 1

		IF(ibval(i,j,k) .eq. 0)CYCLE

		 IF(is_it_i_patch(i,j,k) .eq. 1)THEN
	          patch_no = which_i_patch(i,j,k)		  		  
		  
		  IF(patch(patch_no)%bc_dir .eq. 1)THEN

		  i_start = patch(patch_no)%i_start
		  SBP_row = i - i_start
  		  out(i,j,k) = DOT_PRODUCT( TGfilt(SBP_row,0:10),var(i_start:i_start+10,j,k) )

		  ELSEIF(patch(patch_no)%bc_dir .eq. -1)THEN

		  i_end = patch(patch_no)%i_end
		  SBP_row = 8 + i - i_end
  		  out(i,j,k) = DOT_PRODUCT( TGfilt(SBP_row,0:10),var(i_end-10:i_end,j,k) )
		 
	          ENDIF

		 ELSE

		  out(i,j,k) = DOT_PRODUCT( TGfilt(4,0:8),var(i-4:i+4,j,k) )		
		   
		 ENDIF
				
	
		ENDDO
	     ENDDO		
	  ENDDO

      ENDIF	


        IF(dirn .eq. 2)THEN
     	
	  DO i = 1, Nz	     
	     DO j = 1, Nr
		DO k = 1, 1

		IF(ibval(i,j,k) .eq. 0)CYCLE

		 IF(is_it_j_patch(i,j,k) .eq. 1)THEN
	          patch_no = which_j_patch(i,j,k)		  		  
		  
		  IF(patch(patch_no)%bc_dir .eq. 2)THEN

		  j_start = patch(patch_no)%j_start
		  SBP_row = j - j_start
  		  out(i,j,k) = DOT_PRODUCT( TGfilt(SBP_row,0:10),var(i,j_start:j_start+10,k) )

		  ELSEIF(patch(patch_no)%bc_dir .eq. -2)THEN

		  j_end = patch(patch_no)%j_end
		  SBP_row = 8 + j - j_end
  		  out(i,j,k) = DOT_PRODUCT( TGfilt(SBP_row,0:10),var(i,j_end-10:j_end,k) )
		 
	          ENDIF

		 ELSE

		  out(i,j,k) = DOT_PRODUCT( TGfilt(4,0:8),var(i,j-4:j+4,k) )		
		   
		 ENDIF
				
	
		ENDDO
	     ENDDO		
	  ENDDO

      ENDIF	


END SUBROUTINE PerformTGfilter


END MODULE ModTGfilter
