MODULE ModDeriv

CONTAINS

SUBROUTINE Compute1stDerivative(var,dirn,out,op_flag)

USE ModuleVariables
IMPLICIT NONE

	REAL(KIND=8) :: var(-3:Nz+4,-3:Nr+4,1), out(-3:Nz+4,-3:Nr+4,1)
	INTEGER	:: dirn
	CHARACTER(LEN=8), OPTIONAL :: op_flag

	out = 0.0D0

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
  		  out(i,j,k) = DOT_PRODUCT( SBP1(SBP_row,0:5),var(i_start:i_start+5,j,k) )

		  ELSEIF(patch(patch_no)%bc_dir .eq. -1)THEN

		  i_end = patch(patch_no)%i_end
		  SBP_row = 8 + i - i_end
  		  out(i,j,k) = DOT_PRODUCT( SBP1(SBP_row,0:5),var(i_end-5:i_end,j,k) )
		 
	          ENDIF

		 ELSE

		  out(i,j,k) = DOT_PRODUCT( SBP1(4,0:4),var(i-2:i+2,j,k) )		
		   
		 ENDIF
				
	
		ENDDO
	     ENDDO		
	  ENDDO

	   DO i = 1, Nz
	     DO j = 1, Nr
		DO k = 1, 1

		out(i,j,k) = out(i,j,k)*dxidz(i,j,k)	

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
  		  out(i,j,k) = DOT_PRODUCT( SBP1(SBP_row,0:5),var(i,j_start:j_start+5,k) )

		  ELSEIF(patch(patch_no)%bc_dir .eq. -2)THEN

		  j_end = patch(patch_no)%j_end
		  SBP_row = 8 + j - j_end
  		  out(i,j,k) = DOT_PRODUCT( SBP1(SBP_row,0:5),var(i,j_end-5:j_end,k) )
		 
	          ENDIF

		 ELSE

		  out(i,j,k) = DOT_PRODUCT( SBP1(4,0:4),var(i,j-2:j+2,k) )		
		   
		 ENDIF
				
	
		ENDDO
	     ENDDO		
	  ENDDO


	   DO i = 1, Nz
	     DO j = 1, Nr
		DO k = 1, 1

		out(i,j,k) = out(i,j,k)*detadr(i,j,k)	

		ENDDO
	     ENDDO
	  ENDDO		


       ENDIF	




END SUBROUTINE Compute1stDerivative


SUBROUTINE Compute2ndDerivative(var,dirn,first_deriv,out,op_flag)

USE ModuleVariables
IMPLICIT NONE

	REAL(KIND=8) :: var(-3:Nz+4,-3:Nr+4,1), out(-3:Nz+4,-3:Nr+4,1), first_deriv(-3:Nz+4,-3:Nr+4,1)
	INTEGER	:: dirn
	CHARACTER(LEN=8), OPTIONAL :: op_flag

	out = 0.0D0

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
  		  out(i,j,k) = DOT_PRODUCT( SBP2(SBP_row,0:5),var(i_start:i_start+5,j,k) )

		  ELSEIF(patch(patch_no)%bc_dir .eq. -1)THEN

		  i_end = patch(patch_no)%i_end
		  SBP_row = 8 + i - i_end
  		  out(i,j,k) = DOT_PRODUCT( SBP2(SBP_row,0:5),var(i_end-5:i_end,j,k) )
		 
	          ENDIF

		 ELSE

		  out(i,j,k) = DOT_PRODUCT( SBP2(4,0:4),var(i-2:i+2,j,k) )		
		   
		 ENDIF
				
	
		ENDDO
	     ENDDO		
	  ENDDO

	   DO i = 1, Nz
	     DO j = 1, Nr
		DO k = 1, 1

		out(i,j,k) = out(i,j,k)*dxidz(i,j,k)**2 + first_deriv(i,j,k)*ddxidxidz(i,j,k)		

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
  		  out(i,j,k) = DOT_PRODUCT( SBP2(SBP_row,0:5),var(i,j_start:j_start+5,k) )

		  ELSEIF(patch(patch_no)%bc_dir .eq. -2)THEN

		  j_end = patch(patch_no)%j_end
		  SBP_row = 8 + j - j_end
  		  out(i,j,k) = DOT_PRODUCT( SBP2(SBP_row,0:5),var(i,j_end-5:j_end,k) )
		 
	          ENDIF

		 ELSE

		  out(i,j,k) = DOT_PRODUCT( SBP2(4,0:4),var(i,j-2:j+2,k) )		
		   
		 ENDIF
				
	
		ENDDO
	     ENDDO		
	  ENDDO


	   DO i = 1, Nz
	     DO j = 1, Nr
		DO k = 1, 1

		out(i,j,k) = out(i,j,k)*detadr(i,j,k)**2 + first_deriv(i,j,k)*ddetadetadr(i,j,k)	

		ENDDO
	     ENDDO
	  ENDDO		


       ENDIF	



END SUBROUTINE Compute2ndDerivative


SUBROUTINE ComputeMetric(var,dirn,out)

USE ModuleVariables
IMPLICIT NONE

	REAL(KIND=8) :: var(-3:Nz+4,-3:Nr+4,1), out(-3:Nz+4,-3:Nr+4,1)
	INTEGER	:: dirn
	

	out = 0.0D0

      IF(dirn .eq. 1)THEN
	
	 DO i = 5, Nz-4	     
	     DO j = 1, Nr
		DO k = 1, 1			
		
		out(i,j,k) = DOT_PRODUCT( SBP1(4,0:4),var(i-2:i+2,j,k) )		
	
		ENDDO
	     ENDDO		
	  ENDDO

	 IF(nbrs(3) .lt. 0)THEN
	    DO i = 1, 4
	       DO j = 1, Nr
		  DO k = 1, 1

		   out(i,j,k) = DOT_PRODUCT( SBP1(i-1,0:5),var(1:6,j,k) )	

		  ENDDO
	       ENDDO
	    ENDDO
	 ELSE
	    DO i = 1, 4
	       DO j = 1, Nr
		  DO k = 1, 1

		   out(i,j,k) = DOT_PRODUCT( SBP1(4,0:4),var(i-2:i+2,j,k) )		

		  ENDDO
	       ENDDO
	    ENDDO
	  ENDIF


	IF(nbrs(4) .lt. 0)THEN	   				
	   DO i = Nz-3, Nz
	     DO j = 1, Nr
		DO k = 1, 1

		out(i,j,k) = DOT_PRODUCT( SBP1(8-Nz+i,0:5),var(Nz-5:Nz,j,k) )	

		ENDDO
	     ENDDO
	   ENDDO		 
	ELSE
	    DO i = Nz-3, Nz
	     DO j = 1, Nr
		DO k = 1, 1

		out(i,j,k) = DOT_PRODUCT( SBP1(4,0:4),var(i-2:i+2,j,k) )		

		ENDDO
	     ENDDO
	   ENDDO
	 ENDIF		 	

       ENDIF	

	 IF(dirn .eq. 2)THEN

	 DO i = 1, Nz	     
	     DO j = 5, Nr-4
		DO k = 1, 1		

		out(i,j,k) = DOT_PRODUCT( SBP1(4,0:4),var(i,j-2:j+2,k) )		
	
		ENDDO
	     ENDDO		
	  ENDDO

	IF(nbrs(2) .lt. 0)THEN

	  DO i = 1, Nz
	     DO j = 1, 4
		DO k = 1, 1

		out(i,j,k) = DOT_PRODUCT( SBP1(j-1,0:5),var(i,1:6,k) )	

		ENDDO
	     ENDDO
	  ENDDO		

	ELSE
	
	   DO i = 1, Nz
	     DO j = 1, 4
		DO k = 1, 1

		out(i,j,k) = DOT_PRODUCT( SBP1(4,0:4),var(i,j-2:j+2,k) )

		ENDDO
	     ENDDO
	  ENDDO	

	ENDIF

	IF(nbrs(1) .lt. 0)THEN
	 
		
	   DO i = 1, Nz
	     DO j = Nr-3, Nr
		DO k = 1, 1

		out(i,j,k) = DOT_PRODUCT( SBP1(8-Nr+j,0:5),var(i,Nr-5:Nr,k) )	

		ENDDO
	     ENDDO
	  ENDDO	

	ELSE
	  
	   DO i = 1, Nz
	     DO j = Nr-3, Nr
		DO k = 1, 1

		out(i,j,k) = DOT_PRODUCT( SBP1(4,0:4),var(i,j-2:j+2,k) )

		ENDDO
	     ENDDO
	  ENDDO	

	ENDIF




       ENDIF	


END SUBROUTINE ComputeMetric

 !!! Routines for computing the derivative in computational coordinates (dfdxi, not dfdx)


SUBROUTINE Compute2ndDerivativeComputational(var,dirn,out,op_flag)

USE ModuleVariables
IMPLICIT NONE

	REAL(KIND=8) :: var(-3:Nz+4,-3:Nr+4,1), out(-3:Nz+4,-3:Nr+4,1)
	INTEGER	:: dirn
	CHARACTER(LEN=8), OPTIONAL :: op_flag

	out = 0.0D0

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
  		  out(i,j,k) = DOT_PRODUCT( SBP2(SBP_row,0:5),var(i_start:i_start+5,j,k) )

		  ELSEIF(patch(patch_no)%bc_dir .eq. -1)THEN

		  i_end = patch(patch_no)%i_end
		  SBP_row = 8 + i - i_end
  		  out(i,j,k) = DOT_PRODUCT( SBP2(SBP_row,0:5),var(i_end-5:i_end,j,k) )
		 
	          ENDIF

		 ELSE

		  out(i,j,k) = DOT_PRODUCT( SBP2(4,0:4),var(i-2:i+2,j,k) )		
		   
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
  		  out(i,j,k) = DOT_PRODUCT( SBP2(SBP_row,0:5),var(i,j_start:j_start+5,k) )

		  ELSEIF(patch(patch_no)%bc_dir .eq. -2)THEN

		  j_end = patch(patch_no)%j_end
		  SBP_row = 8 + j - j_end
  		  out(i,j,k) = DOT_PRODUCT( SBP2(SBP_row,0:5),var(i,j_end-5:j_end,k) )
		 
	          ENDIF

		 ELSE

		  out(i,j,k) = DOT_PRODUCT( SBP2(4,0:4),var(i,j-2:j+2,k) )		
		   
		 ENDIF
				
	
		ENDDO
	     ENDDO		
	  ENDDO


       ENDIF	



END SUBROUTINE Compute2ndDerivativeComputational
















END MODULE ModDeriv
