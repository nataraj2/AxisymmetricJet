MODULE ModAllocate

CONTAINS

SUBROUTINE AllocateVariables

USE ModuleVariables

	ALLOCATE( rhs(Nz,Nr,1,neqns) )
	ALLOCATE( zvec(-3:Nz+4,-3:Nr+4,1) )
	ALLOCATE( rvec(-3:Nz+4,-3:Nr+4,1) )
	ALLOCATE( pv(Nz,Nr,1,neqns) )
	ALLOCATE( pv_old(Nz,Nr,1,neqns) )
	ALLOCATE( pv_target(Nz,Nr,1,neqns) )
	ALLOCATE( pv_targetI(Nz,Nr,1,neqns) )
	ALLOCATE( pv_targetV(Nz,Nr,1,neqns) )
	ALLOCATE( pvforSoln(Nz_grid,Nr_grid,1,neqns) )
	ALLOCATE( ibval(Nz,Nr,1) )

	ALLOCATE( dxidz(-3:Nz+4,-3:Nr+4,1) )
	ALLOCATE( detadr(-3:Nz+4,-3:Nr+4,1) )
	ALLOCATE( dzdxi(-3:Nz+4,-3:Nr+4,1) )
	ALLOCATE( drdeta(-3:Nz+4,-3:Nr+4,1) )
	ALLOCATE( ddxidxidz(-3:Nz+4,-3:Nr+4,1) )
	ALLOCATE( ddetadetadr(-3:Nz+4,-3:Nr+4,1) )

	ALLOCATE( rhovec(-3:Nz+4,-3:Nr+4,1) )
	ALLOCATE( Vrvec(-3:Nz+4,-3:Nr+4,1) )
	ALLOCATE( Vtvec(-3:Nz+4,-3:Nr+4,1) )
	ALLOCATE( Vzvec(-3:Nz+4,-3:Nr+4,1) )
	ALLOCATE( Tvec(-3:Nz+4,-3:Nr+4,1) )
	ALLOCATE( pvec(-3:Nz+4,-3:Nr+4,1) )
	ALLOCATE( muvec(-3:Nz+4,-3:Nr+4,1) )
	ALLOCATE( kappavec(-3:Nz+4,-3:Nr+4,1) )
	ALLOCATE( lambdavec(-3:Nz+4,-3:Nr+4,1) )
	ALLOCATE( betavec(-3:Nz+4,-3:Nr+4,1) )

	ALLOCATE( rhovecTGfilt(-3:Nz+4,-3:Nr+4,1) )
	
	ALLOCATE( drhovecdz(-3:Nz+4,-3:Nr+4,1))
	ALLOCATE( dVrvecdz(-3:Nz+4,-3:Nr+4,1))
	ALLOCATE( dVtvecdz(-3:Nz+4,-3:Nr+4,1))
	ALLOCATE( dVzvecdz(-3:Nz+4,-3:Nr+4,1))
	ALLOCATE( dTvecdz(-3:Nz+4,-3:Nr+4,1))
	ALLOCATE( dpvecdz(-3:Nz+4,-3:Nr+4,1))
	ALLOCATE( dmuvecdz(-3:Nz+4,-3:Nr+4,1))	
	ALLOCATE( dkappavecdz(-3:Nz+4,-3:Nr+4,1))		
		
	ALLOCATE( drhovecdr(-3:Nz+4,-3:Nr+4,1))
	ALLOCATE( dVrvecdr(-3:Nz+4,-3:Nr+4,1))
	ALLOCATE( dVtvecdr(-3:Nz+4,-3:Nr+4,1))
	ALLOCATE( dVzvecdr(-3:Nz+4,-3:Nr+4,1))
	ALLOCATE( dTvecdr(-3:Nz+4,-3:Nr+4,1))
	ALLOCATE( dpvecdr(-3:Nz+4,-3:Nr+4,1))
	ALLOCATE( dmuvecdr(-3:Nz+4,-3:Nr+4,1))
	ALLOCATE( dkappavecdr(-3:Nz+4,-3:Nr+4,1))
	
	ALLOCATE( d2Vrvecdz2(-3:Nz+4,-3:Nr+4,1))
	ALLOCATE( d2Vtvecdz2(-3:Nz+4,-3:Nr+4,1))
	ALLOCATE( d2Vzvecdz2(-3:Nz+4,-3:Nr+4,1))
	ALLOCATE( d2Tvecdz2(-3:Nz+4,-3:Nr+4,1))
		
	ALLOCATE( d2Vrvecdr2(-3:Nz+4,-3:Nr+4,1))
	ALLOCATE( d2Vtvecdr2(-3:Nz+4,-3:Nr+4,1))
	ALLOCATE( d2Vzvecdr2(-3:Nz+4,-3:Nr+4,1))
	ALLOCATE( d2Tvecdr2(-3:Nz+4,-3:Nr+4,1))
	
	ALLOCATE( d2Vrvecdzdr(-3:Nz+4,-3:Nr+4,1))	
	ALLOCATE( d2Vzvecdzdr(-3:Nz+4,-3:Nr+4,1))	

	ALLOCATE(k1(Nz,Nr,1,neqns) )
	ALLOCATE(k2(Nz,Nr,1,neqns) )
	ALLOCATE(k3(Nz,Nr,1,neqns) )
	ALLOCATE(k4(Nz,Nr,1,neqns) )

	ALLOCATE(is_it_i_patch(Nz,Nr,1))
	ALLOCATE(is_it_j_patch(Nz,Nr,1))
	ALLOCATE(which_i_patch(Nz,Nr,1))
	ALLOCATE(which_j_patch(Nz,Nr,1))


	! Allocate variables for shock capturing scheme

	ALLOCATE(MgnStrRate(-3:Nz+4,-3:Nr+4,1))
	ALLOCATE(d2MgnStrRatedxi2(-3:Nz+4,-3:Nr+4,1))
	ALLOCATE(d2MgnStrRatedeta2(-3:Nz+4,-3:Nr+4,1))
	ALLOCATE(d4MgnStrRatedxi4(-3:Nz+4,-3:Nr+4,1))
	ALLOCATE(d4MgnStrRatedeta4(-3:Nz+4,-3:Nr+4,1))

	ALLOCATE(d2Tvecdxi2(-3:Nz+4,-3:Nr+4,1))
	ALLOCATE(d2Tvecdeta2(-3:Nz+4,-3:Nr+4,1))
	ALLOCATE(d4Tvecdxi4(-3:Nz+4,-3:Nr+4,1))
	ALLOCATE(d4Tvecdeta4(-3:Nz+4,-3:Nr+4,1))


	ALLOCATE(dilvec(-3:Nz+4,-3:Nr+4,1))
	ALLOCATE(d2dilvecdxi2(-3:Nz+4,-3:Nr+4,1))
	ALLOCATE(d2dilvecdeta2(-3:Nz+4,-3:Nr+4,1))
	ALLOCATE(d4dilvecdxi4(-3:Nz+4,-3:Nr+4,1))
	ALLOCATE(d4dilvecdeta4(-3:Nz+4,-3:Nr+4,1))




	ALLOCATE(mustarQuantity(-3:Nz+4,-3:Nr+4,1))
	ALLOCATE(mustar(-3:Nz+4,-3:Nr+4,1))
	
	ALLOCATE(betastarQuantity(-3:Nz+4,-3:Nr+4,1))
	ALLOCATE(betastar(-3:Nz+4,-3:Nr+4,1))
		
	ALLOCATE(kappastarQuantity(-3:Nz+4,-3:Nr+4,1))
	ALLOCATE(kappastar(-3:Nz+4,-3:Nr+4,1))	

	ALLOCATE(deltvec(Nz,Nr,1))


END SUBROUTINE AllocateVariables

END MODULE ModAllocate
