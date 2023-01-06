MODULE ModAxis

CONTAINS

SUBROUTINE axis_conditions

USE ModuleVariables

	IF(nbrs(2) .lt. 0)THEN

	DO i = 1, Nz
	   DO k = 1, 1

		rhovec(i,1,k) = -1.0D0*DOT_PRODUCT(SBP1(0,1:5),rhovec(i,2:6,k)) /SBP1(0,0)
	        Vzvec(i,1,k) = -1.0D0*DOT_PRODUCT(SBP1(0,1:5),Vzvec(i,2:6,k)) /SBP1(0,0)
		Tvec(i,1,k) = -1.0D0*DOT_PRODUCT(SBP1(0,1:5),Tvec(i,2:6,k)) /SBP1(0,0)

	    ENDDO
	ENDDO

	DO i = 1, Nz
	   DO k = 1, 1

		pv(i,1,k,1) = rhovec(i,1,k)		
		pv(i,1,k,4) = Vzvec(i,1,k)
		pv(i,1,k,5) = Tvec(i,1,k)		

	    ENDDO
	ENDDO

	ENDIF



END SUBROUTINE axis_conditions

END MODULE ModAxis
