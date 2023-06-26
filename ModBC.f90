MODULE ModBC

CONTAINS

SUBROUTINE CallBoundaryConditions
USE ModuleVariables
IMPLICIT NONE
		
	CALL BoundaryConditions	
	CALL sponge_bcs


END SUBROUTINE CallBoundaryConditions


SUBROUTINE BoundaryConditions
USE ModuleVariables
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

		  IF(nbrs(2) .lt. 0)THEN
		  IF(j .eq. 1)CYCLE
		  ENDIF
              	
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

		  IF(nbrs(2) .lt. 0)THEN
		  IF(j .eq. 1)CYCLE
		  ENDIF
              	
		 rho = rhovec(i,j,k)		 
		 Vr  = Vrvec(i,j,k)
		 Vt  = Vtvec(i,j,k)		
		 Vz  = Vzvec(i,j,k)
		 T   = Tvec(i,j,k)

		 drhodr = drhovecdr(i,j,k)
 		 dVrdr = dVrvecdr(i,j,k)
		 dVtdr = dVtvecdr(i,j,k)
		 dVzdr = dVzvecdr(i,j,k)
		 dTdr = dTvecdr(i,j,k)
		 dpdr = dpvecdr(i,j,k)

		 drhodz = drhovecdz(i,j,k)
 		 dVrdz = dVrvecdz(i,j,k)
		 dVtdz = dVtvecdz(i,j,k)
		 dVzdz = dVzvecdz(i,j,k)
		 dTdz = dTvecdz(i,j,k)
		 dpdz = dpvecdz(i,j,k)		

		  IF(abs(patch(patch_no)%bc_dir) .eq. 1)THEN		   		 		    

                    rhs(i,j,k,1) = rhs(i,j,k,1) + (Vz*drhodz + rho*dVzdz)
		    rhs(i,j,k,2) = rhs(i,j,k,2) + Vz*dVrdz
   	            rhs(i,j,k,3) = rhs(i,j,k,3) + Vz*dVtdz
		    rhs(i,j,k,4) = rhs(i,j,k,4) + (Vz*dVzdz + 1.0D0/(gam*rho)*dpdz)
		    rhs(i,j,k,5) = rhs(i,j,k,5) + (Vz*dTdz + (gam-1.0D0)*T*dVzdz)


		    CALL eigen_matrices_in_z_direction	
		    CALL ComputeNSCBCRHS	

		    DO eqn = 1, neqns		 		      		
		      rhs(i,j,k,eqn) = rhs(i,j,k,eqn) - rhs_NSCBC(eqn)    		   
		    ENDDO	
				    		 
		  ENDIF 
		
		  IF(abs(patch(patch_no)%bc_dir) .eq. 2)THEN		   				    	

	            rhs(i,j,k,1) = rhs(i,j,k,1) + (Vr*drhodr + rho*dVrdr)
		    rhs(i,j,k,2) = rhs(i,j,k,2) + (Vr*dVrdr + 1.0D0/(gam*rho)*dpdr)
   	            rhs(i,j,k,3) = rhs(i,j,k,3) + Vr*dVtdr
		    rhs(i,j,k,4) = rhs(i,j,k,4) + Vr*dVzdr
		    rhs(i,j,k,5) = rhs(i,j,k,5) + (Vr*dTdr + (gam-1.0D0)*T*dVrdr) 		   	 		      		

		    CALL eigen_matrices_in_r_direction	
		    CALL ComputeNSCBCRHS

		    DO eqn = 1, neqns		 		      		
		      rhs(i,j,k,eqn) = rhs(i,j,k,eqn) - rhs_NSCBC(eqn)    		   
		    ENDDO	


		  ENDIF 
	
		 
	       ENDDO !k
             ENDDO !j
 	   ENDDO !i	

	ENDIF ! NSCBC_NONREFLECTING

	


	IF(patch(patch_no)%bc .eq. NSCBC_NOSLIP_WALL)THEN               
	                    		 	
          DO i = i_start, i_end
	     DO j = j_start, j_end
               DO k = 1, 1

		 rho = rhovec(i,j,k)
		 Vr = Vrvec(i,j,k)			
		 Vz = Vzvec(i,j,k)
		 drhodr = drhovecdr(i,j,k)
 		 dVrdr = dVrvecdr(i,j,k)		 
		 drhodz = drhovecdz(i,j,k) 		
		 dVzdz = dVzvecdz(i,j,k)
              	
		  IF(abs(patch(patch_no)%bc_dir) .eq. 1)THEN		   
		 
		    CALL eigen_matrices_in_z_direction	
		    CALL ComputeNSCBCRHS	
		   	 		      		
		    rhs(i,j,k,1) = rhs(i,j,k,1) + (Vz*drhodz + rho*dVzdz)

		    DO eqn = 1, neqns		 		      		
		      rhs(i,j,k,eqn) = rhs(i,j,k,eqn) - rhs_NSCBC(eqn)    		   
		    ENDDO	
				    
		    rhs(i,j,k,2:5) = 0.0D0

		  ENDIF 
		
		  IF(abs(patch(patch_no)%bc_dir) .eq. 2)THEN		   
		 
		    CALL eigen_matrices_in_r_direction	
		    CALL ComputeNSCBCRHS	
		   	 		      		
		    rhs(i,j,k,1) = rhs(i,j,k,1) + (Vr*drhodr + rho*dVrdr)

		    DO eqn = 1, neqns		 		      		
		      rhs(i,j,k,eqn) = rhs(i,j,k,eqn) - rhs_NSCBC(eqn)    		   
		    ENDDO	

		    rhs(i,j,k,2:5) = 0.0D0

		  ENDIF 

		 
	       ENDDO !k
             ENDDO !j
 	   ENDDO !i	

	ENDIF ! NSCBC_NOSLIP_WALL






     ENDDO ! patch	


END SUBROUTINE BoundaryConditions


SUBROUTINE eigen_matrices_in_z_direction

USE ModuleVariables
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



END SUBROUTINE eigen_matrices_in_z_direction

SUBROUTINE eigen_matrices_in_r_direction

USE ModuleVariables
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
USE ModuleVariables

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
USE ModuleVariables
IMPLICIT NONE

		    drhodz = drhovecdz(i,j,k)
 		    dVrdz = dVrvecdz(i,j,k)
		    dVtdz = dVtvecdz(i,j,k)
		    dVzdz = dVzvecdz(i,j,k)
		    dTdz = dTvecdz(i,j,k)

		    drhodr = drhovecdr(i,j,k)
 		    dVrdr = dVrvecdr(i,j,k)
		    dVtdr = dVtvecdr(i,j,k)
		    dVzdr = dVzvecdr(i,j,k)
		    dTdr = dTvecdr(i,j,k)				
		   		 		   
	          IF(patch(patch_no)%bc .eq. NSCBC_NONREFLECTING)THEN

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
		
		  ENDIF

		 IF(patch(patch_no)%bc .eq. NSCBC_NOSLIP_WALL)THEN
		 
		   IF(abs(patch(patch_no)%bc_dir) .eq. 1)THEN	

		    deriv_vector(1) = drhodz
		    deriv_vector(2) = dVrdz
		    deriv_vector(3) = dVtdz
		    deriv_vector(4) = dVzdz
		    deriv_vector(5) = dTdz 

		   ELSEIF(abs(patch(patch_no)%bc_dir) .eq. 2)THEN			    		      		

		    deriv_vector(1) = drhodr
		    deriv_vector(2) = dVrdr
		    deriv_vector(3) = dVtdr
		    deriv_vector(4) = dVzdr
		    deriv_vector(5) = dTdr

		   ENDIF

		   waveamp = MATMUL(T_mat_inv,deriv_vector)		   

                   IF(patch(patch_no)%bc_dir .lt. 0)THEN
		     waveamp(4) = -1.0D0*waveamp(5)
		   ELSEIF(patch(patch_no)%bc_dir .gt. 0)THEN		             	   
		     waveamp(5) = -1.0D0*waveamp(4)
		   ENDIF	
		
		  ENDIF
	
		    rhs_NSCBC = MATMUL(MATMUL(T_mat,lam),waveamp)		    				    				  

END SUBROUTINE ComputeNSCBCRHS


SUBROUTINE sponge_bcs
USE ModuleVariables


IMPLICIT NONE

INTEGER(KIND=8) :: j_indx_start


        IF(nbrs(2) .lt. 0)THEN
        j_indx_start = 2
        ELSE
        j_indx_start = 1
        ENDIF

        ! Top sponge


        DO i = 1, Nz
           DO j = 1, Nr
              DO k = 1, 1
                IF(rvec(i,j,k) .ge. sponge_start_top)THEN
                DO eqn = 1, neqns
                zeta = ( rvec(i,j,k) - sponge_start_top )/(rvec(i,Nr,k) - sponge_start_top)
                rhs(i,j,k,eqn) = rhs(i,j,k,eqn) - 0.5D0*zeta**n_sponge*(pv(i,j,k,eqn) - pv_target(i,j,k,eqn) )
                ENDDO
                ENDIF
           ENDDO
         ENDDO
        ENDDO


	 ! Inflow sponge

        DO i = 1, Nz
           DO j = j_indx_start, Nr
              DO k = 1, 1
                IF(zvec(i,j,k) .le. 6.0D0)THEN
                DO eqn = 1, neqns
                zeta = ( 6.0D0 - zvec(i,j,k) )/( 6.0D0 )
                rhs(i,j,k,eqn) = rhs(i,j,k,eqn) - 5.0D0*zeta**n_sponge*(pv(i,j,k,eqn) - pv_target(i,j,k,eqn) )
                ENDDO
                ENDIF
             ENDDO
           ENDDO
         ENDDO


           ! Outflow sponge

        !DO i = 1, Nz
           !DO j = j_indx_start, Nr
              !DO k = 1, 1

                !IF(zvec(i,j,k) .ge. 55.0D0)THEN
                !DO eqn = 1, neqns
                !zeta = (zvec(i,j,k)-55.0D0 )/( 20.0D0 )
                !rhs(i,j,k,eqn) = rhs(i,j,k,eqn) - 0.5D0*zeta**n_sponge*(pv(i,j,k,eqn) - pv_target(i,j,k,eqn) )
                !ENDDO
                !ENDIF

             !ENDDO
           !ENDDO
         !ENDDO





END SUBROUTINE sponge_bcs


END MODULE ModBC
