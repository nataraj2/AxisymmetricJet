MODULE ModBC

CONTAINS

SUBROUTINE CallBoundaryConditions
USE module_variables
IMPLICIT NONE
		
	CALL BoundaryConditions	
	CALL sponge_bcs


END SUBROUTINE CallBoundaryConditions


SUBROUTINE BoundaryConditions
USE module_variables
IMPLICIT NONE

	! Inflow and outflow SAT boundary condition

     
	lam = 0.0D0
	  
    DO patch_no = 1, no_bcs                		  	
	
	  IF(patch(patch_no)%bc .eq. AXIS)CYCLE
	
	     IF(patch(patch_no)%bc_dir .eq. 1)THEN
                  i_start = patch(patch_no)%i_start
		  i_end = patch(patch_no)%i_start
		  j_start = patch(patch_no)%j_start
		  j_end = patch(patch_no)%j_end
	     ELSEIF(patch(patch_no)%bc_dir .eq. -1)THEN	
	     	  i_start = patch(patch_no)%i_end
		  i_end = patch(patch_no)%i_end
		  j_start = patch(patch_no)%j_start
		  j_end = patch(patch_no)%j_end
	    	
	     ELSEIF(patch(patch_no)%bc_dir .eq. 2)THEN	
		  i_start = patch(patch_no)%i_start
		  i_end = patch(patch_no)%i_end		
	     	  j_start = patch(patch_no)%j_start
		  j_end = patch(patch_no)%j_start
		
	     ELSEIF(patch(patch_no)%bc_dir .eq. -2)THEN	
	     	  i_start = patch(patch_no)%i_start
		  i_end = patch(patch_no)%i_end		
	     	  j_start = patch(patch_no)%j_end
		  j_end = patch(patch_no)%j_end
		
	     ENDIF
	 
       IF(patch(patch_no)%bc .eq. SAT_FARFIELD)THEN               
	                    		 	
          DO i = i_start, i_end
	     DO j = j_start, j_end
               DO k = 1, 1

		  IF(j .eq. 1)CYCLE
              	
		  IF(abs(patch(patch_no)%bc_dir) .eq. 1)THEN		   
		 
		    CALL eigen_matrices_in_z_direction
		    CALL compute_SAT_RHS	

		    DO eqn = 1, neqns		 		      		
		     rhs(i,j,k,eqn) = rhs(i,j,k,eqn) + sigma_1*dxidz(i,j,k)*one_by_h_44*rhs_SAT(eqn)
		    ENDDO		

		  ENDIF  		 

	       ENDDO !k
             ENDDO !j
 	   ENDDO !i	

	ENDIF ! SAT_FAR_FIELD

	
	IF(patch(patch_no)%bc .eq. SAT_NOSLIP_WALL)THEN               
	                    		 	
          DO i = i_start, i_end
	     DO j = j_start, j_end
               DO k = 1, 1
              	
		  IF(abs(patch(patch_no)%bc_dir) .eq. 1)THEN		   
		 
		    pv_targetI(i,j,k,1) = pv(i,j,k,1)	
		    pv_targetI(i,j,k,2) = pv(i,j,k,2)
		    pv_targetI(i,j,k,3) = pv(i,j,k,3)
		    pv_targetI(i,j,k,4) = 0.0D0
		    pv_targetI(i,j,k,5) = pv(i,j,k,5)		    
		
		    CALL eigen_matrices_in_z_direction
		    CALL compute_SAT_RHS	

		    DO eqn = 1, neqns		 		      		
		     rhs(i,j,k,eqn) = rhs(i,j,k,eqn) + sigma_1*dxidz(i,j,k)*one_by_h_44*rhs_SAT(eqn) 
		    ENDDO		

		    pv_targetV(i,j,k,1) = pv(i,j,k,1)	
		    pv_targetV(i,j,k,2:4) = 0.0D0	
		    pv_targetV(i,j,k,5) = pv_target(i,j,k,5) 		    

		    DO eqn = 1, neqns		 		      		
		     rhs(i,j,k,eqn) = rhs(i,j,k,eqn) + sigma_2*dxidz(i,j,k)*one_by_h_44*(pv(i,j,k,eqn)-pv_targetV(i,j,k,eqn))
		    ENDDO

		  ELSEIF(abs(patch(patch_no)%bc_dir) .eq. 2)THEN	

		    pv_targetI(i,j,k,1) = pv(i,j,k,1)	
		    pv_targetI(i,j,k,2) = 0.0D0
		    pv_targetI(i,j,k,3) = pv(i,j,k,3)
		    pv_targetI(i,j,k,4) = pv(i,j,k,4)
		    pv_targetI(i,j,k,5) = pv(i,j,k,5)

		    CALL eigen_matrices_in_r_direction
		    CALL compute_SAT_RHS	

		    DO eqn = 1, neqns		 		      		
		     rhs(i,j,k,eqn) = rhs(i,j,k,eqn) + sigma_1*detadr(i,j,k)*one_by_h_44*rhs_SAT(eqn) 
		    ENDDO	
		    
		    pv_targetV(i,j,k,1) = pv(i,j,k,1)
		    pv_targetV(i,j,k,2:4) = 0.0D0			
		    pv_targetV(i,j,k,5) = pv_target(i,j,k,5)			

		    DO eqn = 1, neqns		 		      		
		     rhs(i,j,k,eqn) = rhs(i,j,k,eqn) + sigma_2*detadr(i,j,k)*one_by_h_44*(pv(i,j,k,eqn)-pv_targetV(i,j,k,eqn))
		    ENDDO

		  ENDIF			

	       ENDDO !k
             ENDDO !j
 	   ENDDO !i	

	ENDIF ! SAT_NOSLIP_WALL



	IF(patch(patch_no)%bc .eq. NSCBC_NONREFLECTING)THEN               
	                    		 	
          DO i = i_start, i_end
	     DO j = j_start, j_end
               DO k = 1, 1
              	
		  IF(abs(patch(patch_no)%bc_dir) .eq. 1)THEN		   
		 
		    CALL eigen_matrices_in_z_direction	
		    CALL ComputeNSCBCRHS	

		  ENDIF 

		    DO eqn = 1, neqns		 		      		
		      rhs(i,j,k,eqn) = rhs(i,j,k,eqn) - rhs_NSCBC(eqn)    		   
		    ENDDO		

		 
	       ENDDO !k
             ENDDO !j
 	   ENDDO !i	

	ENDIF ! SAT_FAR_FIELD



     ENDDO ! patch	


END SUBROUTINE BoundaryConditions


SUBROUTINE eigen_matrices_in_z_direction

USE module_variables
IMPLICIT NONE
			 	   
	        rho = rhovec(i,j,k)
	        Vz = Vzvec(i,j,k)
	        T = Tvec(i,j,k)	    		

	        lam(1,1) = Vz
	        lam(2,2) = Vz
   	        lam(3,3) = Vz
	        lam(4,4) = (Vz-DSQRT(T))
	        lam(5,5) = (Vz+DSQRT(T))	
	 
	    	T_mat(1,1) = -1.0D0*rho/T
		T_mat(1,2) = 0.0D0
		T_mat(1,3) = 0.0D0
		T_mat(1,4) = rho/((gam-1.0D0)*T)
		T_mat(1,5) = rho/((gam-1.0D0)*T)
	 	  
		T_mat(2,1) = 0.0D0
		T_mat(2,2) = 0.0D0
		T_mat(2,3) = 1.0D0
		T_mat(2,4) = 0.0D0
		T_mat(2,5) = 0.0D0

		T_mat(3,1) = 0.0D0
		T_mat(3,2) = 1.0D0
		T_mat(3,3) = 0.0D0
		T_mat(3,4) = 0.0D0
		T_mat(3,5) = 0.0D0

		T_mat(4,1) = 0.0D0	
		T_mat(4,2) = 0.0D0
		T_mat(4,3) = 0.0D0
		T_mat(4,4) = -1.0D0/((gam-1.0D0)*DSQRT(T))
		T_mat(4,5) = 1.0D0/((gam-1.0D0)*DSQRT(T))	

		T_mat(5,1) = 1.0D0	
		T_mat(5,2) = 0.0D0
		T_mat(5,3) = 0.0D0
		T_mat(5,4) = 1.0D0
		T_mat(5,5) = 1.0D0	


		T_mat_inv(1,1) = -1.0D0*(gam-1.0D0)*T/(gam*rho)
		T_mat_inv(1,2) = 0.0D0
		T_mat_inv(1,3) = 0.0D0
		T_mat_inv(1,4) = 0.0D0
		T_mat_inv(1,5) = 1.0D0/gam

		T_mat_inv(2,1) = 0.0D0
		T_mat_inv(2,2) = 0.0D0
		T_mat_inv(2,3) = 1.0D0
		T_mat_inv(2,4) = 0.0D0
		T_mat_inv(2,5) = 0.0D0		

		T_mat_inv(3,1) = 0.0D0
		T_mat_inv(3,2) = 1.0D0
		T_mat_inv(3,3) = 0.0D0
		T_mat_inv(3,4) = 0.0D0
		T_mat_inv(3,5) = 0.0D0

		T_mat_inv(4,1) = (gam-1.0D0)*T/(2.0D0*gam*rho)
		T_mat_inv(4,2) = 0.0D0
		T_mat_inv(4,3) = 0.0D0
		T_mat_inv(4,4) = -1.0D0*(gam-1.0D0)*DSQRT(T)/2.0D0
		T_mat_inv(4,5) = (gam-1.0D0)/(2.0D0*gam)

		T_mat_inv(5,1) = (gam-1.0D0)*T/(2.0D0*gam*rho)
		T_mat_inv(5,2) = 0.0D0
		T_mat_inv(5,3) = 0.0D0
		T_mat_inv(5,4) = (gam-1.0D0)*DSQRT(T)/2.0D0
		T_mat_inv(5,5) = (gam-1.0D0)/(2.0D0*gam)


		PRINT*, MATMUL(MATMUL(T_mat,lam),T_mat_inv)
		PRINT*, "rho ......", rho, T
		STOP


			

END SUBROUTINE eigen_matrices_in_z_direction

SUBROUTINE eigen_matrices_in_r_direction

USE module_variables
IMPLICIT NONE
	 	   
	        rho = rhovec(i,j,k)
	        Vr = Vrvec(i,j,k)
	        T = Tvec(i,j,k)	    		

	        lam(1,1) = Vr
	        lam(2,2) = Vr
   	        lam(3,3) = Vr
	        lam(4,4) = (Vr-DSQRT(T))
	        lam(5,5) = (Vr+DSQRT(T))	
	 
	    	T_mat(1,1) = -1.0D0*rho/T
		T_mat(1,2) = 0.0D0
		T_mat(1,3) = 0.0D0
		T_mat(1,4) = rho/((gam-1.0D0)*T)
		T_mat(1,5) = rho/((gam-1.0D0)*T)
	 	  
		T_mat(2,1) = 0.0D0
		T_mat(2,2) = 0.0D0
		T_mat(2,3) = 0.0D0
		T_mat(2,4) = -1.0D0/((gam-1.0D0)*DSQRT(T))
		T_mat(2,5) = 1.0D0/((gam-1.0D0)*DSQRT(T))

		T_mat(3,1) = 0.0D0
		T_mat(3,2) = 0.0D0
		T_mat(3,3) = 1.0D0
		T_mat(3,4) = 0.0D0
		T_mat(3,5) = 0.0D0

		T_mat(4,1) = 0.0D0	
		T_mat(4,2) = 1.0D0
		T_mat(4,3) = 0.0D0
		T_mat(4,4) = 0.0D0
		T_mat(4,5) = 0.0D0

		T_mat(5,1) = 1.0D0	
		T_mat(5,2) = 0.0D0
		T_mat(5,3) = 0.0D0
		T_mat(5,4) = 1.0D0
		T_mat(5,5) = 1.0D0	


		T_mat_inv(1,1) = -1.0D0*(gam-1.0D0)*T/(gam*rho)
		T_mat_inv(1,2) = 0.0D0
		T_mat_inv(1,3) = 0.0D0
		T_mat_inv(1,4) = 0.0D0
		T_mat_inv(1,5) = 1.0D0/gam

		T_mat_inv(2,1) = 0.0D0
		T_mat_inv(2,2) = 0.0D0
		T_mat_inv(2,3) = 0.0D0
		T_mat_inv(2,4) = 1.0D0
		T_mat_inv(2,5) = 0.0D0		

		T_mat_inv(3,1) = 0.0D0
		T_mat_inv(3,2) = 0.0D0
		T_mat_inv(3,3) = 1.0D0
		T_mat_inv(3,4) = 0.0D0
		T_mat_inv(3,5) = 0.0D0

		T_mat_inv(4,1) = (gam-1.0D0)*T/(2.0D0*gam*rho)
		T_mat_inv(4,2) = -1.0D0*(gam-1.0D0)*DSQRT(T)/2.0D0
		T_mat_inv(4,3) = 0.0D0
		T_mat_inv(4,4) = 0.0D0
		T_mat_inv(4,5) = (gam-1.0D0)/(2.0D0*gam)

		T_mat_inv(5,1) = (gam-1.0D0)*T/(2.0D0*gam*rho)
		T_mat_inv(5,2) = (gam-1.0D0)*DSQRT(T)/2.0D0
		T_mat_inv(5,3) = 0.0D0
		T_mat_inv(5,4) = 0.0D0
		T_mat_inv(5,5) = (gam-1.0D0)/(2.0D0*gam)
			

END SUBROUTINE eigen_matrices_in_r_direction

SUBROUTINE compute_SAT_RHS
USE module_variables

		IF(patch(patch_no)%bc_dir .gt. 0)THEN			
		   lam = (abs(lam)+lam)/2.0D0
		ELSEIF(patch(patch_no)%bc_dir .lt. 0)THEN
		   lam = (abs(lam)-lam)/2.0D0
		ENDIF				

	        IF(patch(patch_no)%bc .eq. SAT_FARFIELD)THEN

		DO counter = 1, neqns		 
		  waveamp(counter) = pv(i,j,k,counter) - pv_target(i,j,k,counter)
		ENDDO

		ENDIF

		IF(patch(patch_no)%bc .eq. SAT_NOSLIP_WALL)THEN

		DO counter = 1, neqns		 
		  waveamp(counter) = pv(i,j,k,counter) - pv_targetI(i,j,k,counter)
		ENDDO

		ENDIF

		waveamp = MATMUL(T_mat_inv,waveamp)		 

		rhs_SAT = MATMUL(MATMUL(T_mat,lam),waveamp)


END SUBROUTINE compute_SAT_RHS

SUBROUTINE ComputeNSCBCRHS
USE module_variables
IMPLICIT NONE

		    drhodz = drhovecdz(i,j,k)
 		    dVrdz = dVrvecdz(i,j,k)
		    dVtdz = dVtvecdz(i,j,k)
		    dVzdz = dVzvecdz(i,j,k)
		    dTdz = dTvecdz(i,j,k)		

		    deriv_vector(1) = drhodz
		    deriv_vector(2) = dVrdz
		    deriv_vector(3) = dVtdz
		    deriv_vector(4) = dVzdz
		    deriv_vector(5) = dTdz 

		    waveamp = MATMUL(T_mat_inv,deriv_vector)		 
		   
		    IF (patch(patch_no)%bc_dir .gt. 0)THEN	
		      DO counter = 1, neqns
		      IF(lam(counter,counter) .gt. 0.0D0)THEN
		        waveamp(counter) = 0.0D0
	    	      ENDIF
		      ENDDO	 
                    ELSEIF(patch(patch_no)%bc_dir .lt. 0)THEN
	              DO counter = 1, neqns
		      IF(lam(counter,counter) .lt. 0.0D0)THEN
		        waveamp(counter) = 0.0D0
	    	      ENDIF
	              ENDDO 	   
		    ENDIF		    

		    rhs_NSCBC = MATMUL(MATMUL(T_mat,lam),waveamp)		    				    				  

END SUBROUTINE ComputeNSCBCRHS


SUBROUTINE sponge_bcs
USE module_variables

IMPLICIT NONE

	! Top sponge

       DO eqn = 1, neqns
	DO i = 1, Nz
	   DO j = Nr+sponge_start_top, Nr
	      DO k = 1, 1
 
                zeta = ( rvec(i,j,k) - rvec(i,Nr+sponge_start_top,k) )/( rvec(i,Nr,k) - rvec(i,Nr+sponge_start_top,k) )		
		rhs(i,j,k,eqn) = rhs(i,j,k,eqn) - A_sponge*zeta**n_sponge*(pv(i,j,k,eqn) - pv_target(i,j,k,eqn) ) 

	     ENDDO
           ENDDO
	 ENDDO
	ENDDO


	! Inflow sponge
	
      DO eqn = 1, neqns
	DO i = 1, sponge_end_inflow
	   DO j = 1, Nr
	      DO k = 1, 1
 
                zeta = ( zvec(sponge_end_inflow,j,k) - zvec(i,j,k) )/( zvec(sponge_end_inflow,j,k) - zvec(1,j,k) )		
		!rhs(i,j,k,eqn) = rhs(i,j,k,eqn) - A_sponge*zeta**n_sponge*(pv(i,j,k,eqn) - pv_target(i,j,k,eqn) ) 

	     ENDDO
           ENDDO
	 ENDDO
	ENDDO
	
	
	! Outflow sponge
	
      DO eqn = 1, neqns
	DO i = Nz+sponge_start_outflow, Nz
	   DO j = 1, Nr
	      DO k = 1, 1
 
                zeta =  ( zvec(i,j,k) - zvec(Nz+sponge_start_outflow,j,k) )/( zvec(Nz,j,k) - zvec(Nz+sponge_start_outflow,j,k) )				
		!rhs(i,j,k,eqn) = rhs(i,j,k,eqn) - A_sponge*zeta**n_sponge*(pv(i,j,k,eqn) - pv_target(i,j,k,eqn) ) 

	     ENDDO
           ENDDO
	 ENDDO
	ENDDO

END SUBROUTINE sponge_bcs


END MODULE ModBC
