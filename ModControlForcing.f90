MODULE ModControlForcing

CONTAINS


SUBROUTINE ControlForcingTerms
USE ModuleVariables
IMPLICIT NONE 



     IF(iteration .eq. iteration_start+1)THEN


       forcing_mat = 0.0D0          

       OPEN(UNIT=101,FILE='Vz_Control_Vr_Feedback.dat')        

        DO i = 1, 5
           DO j = 1, 5
             READ(101,*)forcing_mat(i,j)
           ENDDO
        ENDDO

        READ(101,*)z0
        READ(101,*)r0
        READ(101,*)lz
        READ(101,*)lr
       
        CLOSE(101)
        
        IF(rank .eq. 0)THEN
        DO i = 1, 5
            DO j = 1, 5
                PRINT*,i,j,forcing_mat(i,j)
            ENDDO
        ENDDO
        PRINT*, z0, r0, lz, lr
        ENDIF
 
     ENDIF



       DO i = 1, Nz
          DO j = 1, Nr           

           z_loc = G(1,istart+(i-1),jstart+(j-1),1,1)
           r_loc = G(1,istart+(i-1),jstart+(j-1),1,2)

           GaussianFactor = exp( -(z_loc-z0)**2/lz**2 - (r_loc-r0)**2/lr**2 )   

	   DO eqn = 1, neqns
             DO k = 1, neqns
  	       rhs(i,j,1,eqn) = rhs(i,j,1,eqn) + alpha*forcing_mat(eqn,k)*GaussianFactor*( pv(i,j,1,k) - pv_target(i,j,1,k) )
	     ENDDO
	   ENDDO

           
          ENDDO
       ENDDO 
       


END SUBROUTINE ControlForcingTerms

END MODULE ModControlForcing


