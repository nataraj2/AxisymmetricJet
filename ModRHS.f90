MODULE ModRHS

CONTAINS

SUBROUTINE ComputeRHS

USE ModuleVariables
USE ModInputOutput
USE ModDecomp
USE ModCommunicate
USE ModDeriv
USE ModControlForcing
USE ModBC
USE ModPLOT3D_IO
USE ModShockCapturing
IMPLICIT NONE
	
INTEGER(KIND=8) :: j_indx_start

	CALL ProcessorCommunication(rhovec)
	CALL ProcessorCommunication(Vrvec)
	CALL ProcessorCommunication(Vtvec)
	CALL ProcessorCommunication(Vzvec)
	CALL ProcessorCommunication(Tvec)
	CALL ProcessorCommunication(pvec)
	CALL ProcessorCommunication(muvec)
	CALL ProcessorCommunication(kappavec)
	
	CALL Compute1stDerivative(rhovec,1,drhovecdz)
	CALL Compute1stDerivative(Vrvec,1,dVrvecdz)
	CALL Compute1stDerivative(Vtvec,1,dVtvecdz)
	CALL Compute1stDerivative(Vzvec,1,dVzvecdz)
	CALL Compute1stDerivative(Tvec,1,dTvecdz)
	CALL Compute1stDerivative(pvec,1,dpvecdz)

	CALL Compute1stDerivative(rhovec,2,drhovecdr)
	CALL Compute1stDerivative(Vrvec,2,dVrvecdr)
	CALL Compute1stDerivative(Vtvec,2,dVtvecdr)
	CALL Compute1stDerivative(Vzvec,2,dVzvecdr)
	CALL Compute1stDerivative(Tvec,2,dTvecdr)	
	CALL Compute1stDerivative(pvec,2,dpvecdr)					
	
	CALL ProcessorCommunication(dVrvecdz)
	CALL ProcessorCommunication(dVzvecdz)

	IF(ShockCapturing .eq. 1)THEN
	CALL ComputeArtTransCoeffs
	ENDIF

	CALL ProcessorCommunication(muvec)
	CALL ProcessorCommunication(kappavec)

	!CALL WriteVariableOutput(d2Vrvecdzdr,d2Vrvecdzdr,d2Vtvecdr2,d2Vzvecdr2,d2Tvecdr2)
	!CALL WriteVariableOutput(d2Vrvecdzdr,d2Vrvecdzdr,d2Vtvecdr2,d2Vzvecdzdr,d2Tvecdr2)
	!CALL StopTheCode	

	rhs = 0.0D0
	
	IF(nbrs(2) .lt. 0)THEN
	j_indx_start = 2
	ELSE
	j_indx_start = 1
	ENDIF


	DO i = 1, Nz
	   DO j = j_indx_start, Nr
	      DO k = 1, 1		 
	
		 r = rvec(i,j,k)

		 rho = rhovec(i,j,k)
		 Vr = Vrvec(i,j,k)	
		 Vt = Vtvec(i,j,k)
		 Vz = Vzvec(i,j,k)
		 T = Tvec(i,j,k)

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
		
		 rhs(i,j,k,1) = Vr*drhodr + Vz*drhodz + rho*(dVrdr + Vr/r + dVzdz)
		 rhs(i,j,k,2) = Vr*dVrdr + Vz*dVrdz + 1.0D0/(gam*rho)*dpdr!(rho*dTdr+T*drhodr)!dpdr
		 rhs(i,j,k,3) = Vr*dVtdr + Vr*Vt/r + Vz*dVtdz
		 rhs(i,j,k,4) = Vr*dVzdr + Vz*dVzdz + 1.0D0/(gam*rho)*dpdz!(rho*dTdz+T*drhodz)!dpdz
		 rhs(i,j,k,5) = Vr*dTdr + Vz*dTdz + (gam-1.0D0)*T*(dVrdr + Vr/r + dVzdz)	
		
	      ENDDO
	   ENDDO
	ENDDO		


       DO eqn = 1, neqns
	DO i = 1, Nz
	   DO j = 1, Nr
	      DO k = 1, 1  

		rhs(i,j,k,eqn) = -1.0D0*rhs(i,j,k,eqn)

	      ENDDO
	   ENDDO
	 ENDDO	
	ENDDO


	IF(do_viscous_terms .eq. 1)THEN

	CALL Compute2ndDerivative(Vrvec,1,dVrvecdz,d2Vrvecdz2)
	CALL Compute2ndDerivative(Vtvec,1,dVtvecdz,d2Vtvecdz2)
	CALL Compute2ndDerivative(Vzvec,1,dVzvecdz,d2Vzvecdz2)
	CALL Compute2ndDerivative(Tvec,1,dTvecdz,d2Tvecdz2)
	CALL Compute1stDerivative(muvec,1,dmuvecdz)
	CALL Compute1stDerivative(kappavec,1,dkappavecdz)

	CALL Compute2ndDerivative(Vrvec,2,dVrvecdr,d2Vrvecdr2)
	CALL Compute2ndDerivative(Vtvec,2,dVtvecdr,d2Vtvecdr2)
	CALL Compute2ndDerivative(Vzvec,2,dVzvecdr,d2Vzvecdr2)
	CALL Compute2ndDerivative(Tvec,2,dTvecdr,d2Tvecdr2)	
	CALL Compute1stDerivative(muvec,2,dmuvecdr)	
	CALL Compute1stDerivative(kappavec,2,dkappavecdr)

	CALL Compute1stDerivative(dVrvecdz,2,d2Vrvecdzdr)
	CALL Compute1stDerivative(dVzvecdz,2,d2Vzvecdzdr)	

	DO i = 1, Nz
	   DO j = j_indx_start, Nr
	      DO k = 1, 1

		 r = rvec(i,j,k)
		 onebyr = 1.0D0/r 
		 rho = rhovec(i,j,k)
		 Vr = Vrvec(i,j,k)	
		 Vt = Vtvec(i,j,k)
		 Vz = Vzvec(i,j,k)
		 T = Tvec(i,j,k)
		 mu = muvec(i,j,k)
		 kappa = kappavec(i,j,k)

		 drhodr = drhovecdr(i,j,k)
 		 dVrdr = dVrvecdr(i,j,k)
		 dVtdr = dVtvecdr(i,j,k)
		 dVzdr = dVzvecdr(i,j,k)
		 dTdr = dTvecdr(i,j,k)
		 dpdr = dpvecdr(i,j,k)
		 dmudr = dmuvecdr(i,j,k)
		 dkappadr = dkappavecdr(i,j,k)

		 drhodz = drhovecdz(i,j,k)
 		 dVrdz = dVrvecdz(i,j,k)
		 dVtdz = dVtvecdz(i,j,k)
		 dVzdz = dVzvecdz(i,j,k)
		 dTdz = dTvecdz(i,j,k)
		 dpdz = dpvecdz(i,j,k)
		 dmudz = dmuvecdz(i,j,k)
		 dkappadz = dkappavecdz(i,j,k)
	
		 d2Vrdz2 = d2Vrvecdz2(i,j,k)
		 d2Vtdz2 = d2Vtvecdz2(i,j,k)
		 d2Vzdz2 = d2Vzvecdz2(i,j,k)
		 d2Tdz2 = d2Tvecdz2(i,j,k)

 		 d2Vrdr2 = d2Vrvecdr2(i,j,k)
		 d2Vtdr2 = d2Vtvecdr2(i,j,k)
		 d2Vzdr2 = d2Vzvecdr2(i,j,k)
		 d2Tdr2 = d2Tvecdr2(i,j,k)

		 d2Vrdzdr = d2Vrvecdzdr(i,j,k)
		 d2Vzdzdr = d2Vzvecdzdr(i,j,k)

		 IF(ShockCapturing .eq. 1)THEN

		 tmp = lambdavec(i,j,k)/muvec(i,j,k)
		
		 fourby3 = 2.0D0 + tmp
		 oneby3 = 1.0D0 + tmp
		 twoby3 = -1.0D0*tmp

		 ENDIF
 				
		 rhs(i,j,k,2) = rhs(i,j,k,2) + 1.0D0/(Re*rho)*( mu*(fourby3*d2Vrdr2 + fourby3*onebyr*dVrdr + oneby3*d2Vzdzdr + d2Vrdz2 - fourby3*Vr*onebyr**2)+&
					       dmudr*(fourby3*dVrdr - twoby3*Vr*onebyr - twoby3*dVzdz) +&
					       dmudz*(dVrdz+dVzdr) )
		 rhs(i,j,k,3) = rhs(i,j,k,3) + 1.0D0/(Re*rho)*(mu*(d2Vtdr2+d2Vtdz2+onebyr*dVtdr-Vt*onebyr**2) +&
					       dmudr*(dVtdr - Vt*onebyr) +&
					       dmudz*dVtdz )	
		 rhs(i,j,k,4) = rhs(i,j,k,4) + 1.0D0/(Re*rho)*( mu*(d2Vzdr2 + fourby3*d2Vzdz2 + oneby3*d2Vrdzdr + onebyr*dVzdr + oneby3*onebyr*dVrdz)+&
					       dmudr*(dVzdr+dVrdz) + &
					       dmudz*(fourby3*dVzdz-twoby3*dVrdr-twoby3*Vr*onebyr) )
		 
		 rhs(i,j,k,5) = rhs(i,j,k,5) + gam/(Re*Pr*rho)*( kappa*(d2Tdr2+onebyr*dTdr + d2Tdz2) + dkappadr*dTdr + dkappadz*dTdz ) + &
					       gam*(gam-1.0D0)/Re*mu/rho*( dVrdr*(fourby3*dVrdr-twoby3*Vr*onebyr-twoby3*dVzdz) + &
								       (dVtdr-Vt*onebyr)*(dVtdr-Vr*onebyr)+&
								       Vr*onebyr*(fourby3*Vr*onebyr-twoby3*dVrdr-twoby3*dVzdz)+&	
								       dVtdz*dVtdz + dVzdz*(fourby3*dVzdz-twoby3*dVrdr-twoby3*Vr*onebyr)+&
								       (dVrdz+dVzdr)*(dVzdr+dVrdz) )		
				
	      ENDDO
	   ENDDO
	ENDDO		


	ENDIF


	!IF(write_rhs_file .eq. 1)THEN	

	!CALL WriteRHS	
	!CALL StopTheCode
	
	!ENDIF


        IF(ControlForcing .eq. 1)THEN
        CALL ControlForcingTerms
        ENDIF


	CALL CallBoundaryConditions


	!IF(write_rhs_file .eq. 1)THEN	

	!CALL WriteRHS	
	!CALL StopTheCode
	
	!ENDIF


	IF(do_mms .eq. 1)THEN
	
	DO eqn = 1, neqns
	DO i = 1, Nz
	   DO j = 1, Nr
	      DO k = 1, 1

		rhs(i,j,k,eqn) = rhs(i,j,k,eqn) - rhs_initial(i,j,k,eqn)

	      ENDDO
	   ENDDO
	 ENDDO	
	ENDDO

	ENDIF

	


END SUBROUTINE ComputeRHS




END MODULE ModRHS
