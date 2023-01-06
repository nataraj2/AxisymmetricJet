MODULE ModShockCapturing

CONTAINS

SUBROUTINE ComputeArtTransCoeffs
USE ModuleVariables
USE ModInputOutput
USE ModCommunicate
USE ModDecomp
USE ModDeriv
USE ModTGFilter

IMPLICIT NONE


	! Compute magnitude of the Strain-Rate tensor

	CALL ComputeMgnStrRate

	! Compute second derivative of MgnStrRate d2Sdxi2
	
	CALL ProcessorCommunication(MgnStrRate)
        CALL ProcessorCommunication(Tvec)

	CALL Compute2ndDerivativeComputational(MgnStrRate,1,d2MgnStrRatedxi2)
	CALL Compute2ndDerivativeComputational(Tvec,1,d2Tvecdxi2)		

	CALL Compute2ndDerivativeComputational(MgnStrRate,2,d2MgnStrRatedeta2)
        CALL Compute2ndDerivativeComputational(Tvec,2,d2Tvecdeta2)
	
	! Compute fourth derivatives of MgnStrRate and Temperature

	CALL ProcessorCommunication(d2MgnStrRatedxi2)
	CALL ProcessorCommunication(d2Tvecdxi2)
	
	CALL ProcessorCommunication(d2MgnStrRatedeta2)
        CALL ProcessorCommunication(d2Tvecdeta2)

	CALL Compute2ndDerivativeComputational(d2MgnStrRatedxi2,1,d4MgnStrRatedxi4)
	CALL Compute2ndDerivativeComputational(d2Tvecdxi2,1,d4Tvecdxi4)

        CALL Compute2ndDerivativeComputational(d2MgnStrRatedeta2,2,d4MgnStrRatedeta4)
        CALL Compute2ndDerivativeComputational(d2Tvecdeta2,2,d4Tvecdeta4)
	
	CALL ComputeTobeFilteredQuantities

	! Perform the filtering and multiply by coefficients C_mu, C_beta or
	! C_kappa and (Re or RePr) to get the artificial transport coefficients

	CALL ProcessorCommunication(mustarQuantity)
	CALL ProcessorCommunication(kappastarQuantity)

	!!!!!!!!!!!!!!!!!!!!!!!!!
	!DO dummy_int = 1, 4
	!CALL PerformTGfilter(mustarQuantity,1,mustarQuantity)	
	!CALL ProcessorCommunication(mustarQuantity)
	!CALL PerformTGfilter(mustarQuantity,2,mustarQuantity)
	!CALL ProcessorCommunication(mustarQuantity)
	!ENDDO
	!!!!!!!!!!!!!!!!!!!!!!!!

	CALL PerformTGfilter(mustarQuantity,1,mustarQuantity)
        CALL ProcessorCommunication(mustarQuantity)
	CALL PerformTGfilter(mustarQuantity,2,mustar)
	mustar = C_mu*Re*mustar

	!!!!!!!!!!!!!!!!!!!!!!!!
	!DO dummy_int = 1, 4
        !CALL PerformTGfilter(kappastarQuantity,1,kappastarQuantity)
	!CALL ProcessorCommunication(kappastarQuantity)
	!CALL PerformTGfilter(kappastarQuantity,2,kappastarQuantity)
	!CALL ProcessorCommunication(kappastarQuantity)
	!ENDDO	
	!!!!!!!!!!!!!!!!!!!!!!!!!

	CALL PerformTGfilter(kappastarQuantity,1,kappastarQuantity)
        CALL ProcessorCommunication(kappastarQuantity)
	CALL PerformTGfilter(kappastarQuantity,2,kappastar)

	kappastar = C_kappa*Re*Pr*kappastar

	DO i = 1, Nz
           DO j = 1, Nr
              DO k = 1, 1

		tmp = ( 1.0D0+tanh( 5.0D0*(zvec(i,j,k)-10.0D0) ) )/2.0D0
                mustar(i,j,k) = mustar(i,j,k)*tmp
                kappastar(i,j,k) = kappastar(i,j,k)*tmp

              ENDDO
           ENDDO
        ENDDO
		
	CALL UpdateTransportCoeffs

	IF(mod(iteration,noutput) .eq. 0 .and. rkStep .eq. 4)THEN
	CALL WriteVariableOutput(mustar,betavec,kappastar,muvec,muvec)
        CALL MPI_Barrier(MPI_COMM_WORLD,ierr)
	ENDIF

	!CALL StopTheCode

	! RHS to be called after this
		

END SUBROUTINE ComputeArtTransCoeffs


SUBROUTINE ComputeMgnStrRate
USE ModuleVariables
IMPLICIT NONE

	DO i = 1, Nz
	   DO j = 1, Nr
	     DO k = 1, 1

		S11 = dVrvecdr(i,j,k)	
		S12 = 0.0D0
		S13 = 0.5D0*(dVrvecdz(i,j,k) + dVzvecdr(i,j,k) )
		S21 = 0.0D0
		S23 = 0.0D0
		S31 = S13
		S32 = 0.0D0
		S33 = dVzvecdz(i,j,k)

		IF(j .eq. 1)THEN	
		S22 = dVrvecdr(i,j,k)
		ELSE		
		S22 = Vrvec(i,j,k)/rvec(i,j,k)
		ENDIF

		MgnStrRate(i,j,k) = DSQRT(S11**2+S12**2+S13**2+S21**2+S22**2+S23**2+S31**2+S32**2+S33**2)

	      ENDDO
	   ENDDO
	ENDDO


END SUBROUTINE ComputeMgnStrRate


SUBROUTINE ComputeTobeFilteredQuantities

USE ModuleVariables
IMPLICIT NONE

	DO i = 1, Nz
	   DO j = 1, Nr
	     DO k = 1, 1
		
		IF(i .ge. 2 .and. i .le. Nz-1)THEN
		deltaz = (zvec(i+1,j,k)-zvec(i-1,j,k))/2.0D0
		ENDIF
		IF(i .eq. 1)THEN
		deltaz = zvec(2,j,k) - zvec(1,j,k)
		ENDIF
		IF(i .eq. Nz)THEN
		deltaz = zvec(Nz,j,k) - zvec(Nz-1,j,k)
		ENDIF

		IF(j .ge. 2 .and. j .le. Nr-1)THEN
                deltar = (rvec(i,j+1,k)-rvec(i,j-1,k))/2.0D0
                ENDIF
                IF(j .eq. 1)THEN
                deltar = rvec(i,2,k) - rvec(i,1,k)
                ENDIF
                IF(j .eq. Nr)THEN
                deltar = rvec(i,Nr,k) - rvec(i,Nr-1,k)
                ENDIF
		
		
		mustarQuantity(i,j,k) = rhovec(i,j,k)*ABS( ( d4MgnStrRatedxi4(i,j,k)*deltaz**2 + &
							     d4MgnStrRatedeta4(i,j,k)*deltar**2) )


		!betastarQuantity is the same as mustarQuantity

		kappastarQuantity(i,j,k) = rhovec(i,j,k)/(gam*DSQRT(Tvec(i,j,k)))*ABS( (d4Tvecdxi4(i,j,k)*deltaz + &
                                                                                        d4Tvecdeta4(i,j,k)*deltar) )

	      ENDDO
	   ENDDO
	ENDDO 


END SUBROUTINE ComputeTobeFilteredQuantities



SUBROUTINE UpdateTransportCoeffs

USE ModuleVariables
IMPLICIT NONE


	DO i = 1, Nz
	   DO j = 1, Nr
	      DO k = 1, 1

		muvec(i,j,k) = muvec(i,j,k) + mustar(i,j,k)
		betavec(i,j,k) = betavec(i,j,k) + mustar(i,j,k)/C_mu*C_beta
		lambdavec(i,j,k) = betavec(i,j,k) - (2.0D0/3.0D0)*muvec(i,j,k)
		kappavec(i,j,k) = kappavec(i,j,k) + kappastar(i,j,k)

	      ENDDO
	   ENDDO
	ENDDO

	
END SUBROUTINE UpdateTransportCoeffs

END MODULE ModShockCapturing
