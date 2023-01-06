MODULE ModOperator

CONTAINS

SUBROUTINE SBP_operator

USE ModuleVariables

	ALLOCATE(SBP1(0:8,0:5))
	ALLOCATE(SBP2(0:8,0:5))	
	ALLOCATE(TGFilt(0:8,0:10))

	DO i = 1, 9
           DO j = 1, 6

		SBP1(i-1,j-1) = 0.0D0

	   ENDDO
	ENDDO	

	!Form the SBP1 ddx operator

	!Form the top block (4 x 6)

	SBP1(0,0:3) = (/ -24.0D0/17.0D0, 59.0D0/34.0D0, -4.0D0/17.0D0, -3.0D0/34.0D0 /)
	SBP1(1,0:2) = (/ -1.0/2.0D0, 0.0D0, 1.0D0/2.0D0 /)
	SBP1(2,0:4) = (/ 4.0D0/43.0D0, -59.0D0/86.0D0, 0.0D0, 59.0/86.0D0, -4.0D0/43.0D0 /)
	SBP1(3,0:5) = (/ 3.0D0/98.0D0, 0.0D0, -59.0D0/98.0D0, 0.0D0, 32.0D0/49.0D0, -4.0D0/49.0D0 /)

	!Form the bottom block (4 x 6)

	DO i = 9, 6, -1
	   DO j = 6, 1, -1		

	     SBP1(i-1,j-1) = -1.0D0*SBP1(10-i-1,7-j-1)

	   ENDDO
	ENDDO

	!Form the interior stencils

	 SBP1(4,0:4) = (/1.0D0/12.0D0, -2.0D0/3.0D0, 0.0D0, 2.0D0/3.0D0, -1.0D0/12.0D0 /)		


	! SBP second deriv operator

	

	DO i = 1, 9
           DO j = 1, 6

		SBP2(i-1,j-1) = 0.0D0

	   ENDDO
	ENDDO	

	!Form the SBP ddx operator

	!Form the top block (4 x 6)

	SBP2(0,0:3) = (/ 2.0D0, -5.0D0, 4.0D0, -1.0D0 /)
	SBP2(1,0:2) = (/ 1.0D0, -2.0D0, 1.0D0 /)
	SBP2(2,0:4) = (/ -4.0D0/43.0D0, 59.0D0/43.0D0, -110.0D0/43.0D0, 59.0/43.0D0, -4.0D0/43.0D0 /)
	SBP2(3,0:5) = (/ -1.0D0/49.0D0, 0.0D0, 59.0D0/49.0D0, -118.0D0/49.0D0, 64.0D0/49.0D0, -4.0D0/49.0D0 /)

	!Form the bottom block (4 x 6)

	DO i = 9, 6, -1
	   DO j = 6, 1, -1		

	     SBP2(i-1,j-1) = 1.0D0*SBP2(10-i-1,7-j-1)

	   ENDDO
	ENDDO

	!Form the interior stencils

	 SBP2(4,0:4) = (/-1.0D0/12.0D0, 4.0D0/3.0D0, -5.0D0/2.0D0, 4.0D0/3.0D0, -1.0D0/12.0D0 /)		

	


	!!! Stencils for Truncated Gaussian filter

	DO i = 1, 9
	   DO j = 1, 11
	      TGfilt(i-1,j-1) = 0.0D0
	   ENDDO
	ENDDO

	!Form the top block (4 x 11)

        TGfilt(0,0) = 1.0D0
        TGfilt(1,0:6)  = (/0.085777408970000,0.722371828476000,0.356848072173000,-0.223119093072000,0.057347064865000,0.000747264596000,0.000027453993000/)
        TGfilt(2,0:6)  = (/-0.032649010764000,0.143339502575000,0.726678822020000,0.294622121167000,-0.186711738069000,0.062038376258000,-0.007318073189000/)
        TGfilt(3,0:10) = (/0.000054596010000,-0.042124772446000,0.173103107841000,0.700384128648000,0.276543612935000,&
                          -0.131223506571000,0.023424966418000,-0.013937561779000,0.024565095706000,-0.013098287852000,0.002308621090000/)

        !Form the bottom block (4 x 11))

        DO i = 9, 6, -1
           DO j = 11, 1, -1

             TGfilt(i-1,j-1) = 1.0D0*TGfilt(10-i-1,12-j-1)

           ENDDO
        ENDDO

	! Form the interior stencils

        TGfilt(4,0:8)=(/ 107.0D0/103680.0D0, 149.0D0/12960.0D0, 1997.0D0/25920.0D0, 3091.0D0/12960.0D0,&
		       3565.0D0/10368.0D0, 3091.0D0/12960.0D0,1997.0D0/25920.0D0, 149.0D0/12960.0D0, 107.0D0/103680.0D0 /)	








END SUBROUTINE SBP_operator

END MODULE ModOperator

