MODULE ModAssign 

CONTAINS

SUBROUTINE AssignValues

USE ModuleVariables
IMPLICIT NONE

	DO i = 1, Nz
	   DO j = 1, Nr
	      DO k = 1, 1

		rhovec(i,j,k) = pv(i,j,k,1)
		Vrvec(i,j,k)  = pv(i,j,k,2)
		Vtvec(i,j,k)  = pv(i,j,k,3)
		Vzvec(i,j,k)  = pv(i,j,k,4)	
		Tvec(i,j,k)  = pv(i,j,k,5)
		pvec(i,j,k) = rhovec(i,j,k)*Tvec(i,j,k)
		muvec(i,j,k) = Tvec(i,j,k)**0.666D0
		kappavec(i,j,k) = muvec(i,j,k)
		betavec(i,j,k) = (3.0D0/5.0D0)*muvec(i,j,k)
		lambdavec(i,j,k) = betavec(i,j,k) - (2.0D0/3.0D0)*muvec(i,j,k)
				
	      ENDDO
	   ENDDO	
	 ENDDO

END SUBROUTINE AssignValues

END MODULE ModAssign
