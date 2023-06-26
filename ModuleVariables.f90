MODULE ModuleVariables
include 'mpif.h'

	!!!!! Inputs
	
    CHARACTER(LEN=200) :: GRID_FILENAME = 'axijet_mesh.xyz'
    CHARACTER(LEN=80)  :: INITIAL_CONDITION = 'Restart.q'
    CHARACTER(LEN=200) :: TARGET_FILE = 'axijet_soln.q'
    CHARACTER(LEN=200) :: bc_file = 'bc.dat'

	! Sponge	
	double precision :: sponge_start_top = 12.5d0, sponge_A_top = 0.5d0
	double precision :: sponge_end_left = 6.0d0, sponge_A_left = 5.0d0

	! Processor decomposition
	INTEGER dims(2)
	DATA dims /4,2/

	!!!!!!! Inputs end
	
	CHARACTER(LEN=200) :: RHS_FILENAME = 'jet_rhs.q'


	INTEGER                 :: NDIM, ngrid, ftype, grid_NO
        INTEGER,POINTER         :: IBLANK(:,:,:,:)
        INTEGER,POINTER         :: ND(:,:)    
        CHARACTER(2)            :: prec, gf, vf, ib
        INTEGER, PARAMETER      :: RFREAL = SELECTED_REAL_KIND(PRECISION(1.0D0))
        REAL(rfreal)            :: tau(4)
        REAL(rfreal), POINTER   :: G(:,:,:,:,:), X(:,:,:,:,:), rhs(:,:,:,:), rhs_initial(:,:,:,:), X_rhs(:,:,:,:,:), X_target(:,:,:,:,:), Xout(:,:,:,:,:), X2(:,:,:,:,:)
        REAL(rfreal), POINTER   :: Xoutfinal(:,:,:,:,:)
        INTEGER                 :: OF = 1    

	INTEGER(KIND=8)		:: i, j, k, Nz_grid, Nr_grid, eqn, neqns, SBP_row
	INTEGER(KIND=8) 	:: iteration, iteration_start
	CHARACTER(LEN=2)	:: SOLN_FILETYPE = '.q', ART_COEFF_SOLN_FILETYPE = '.v', RHS_SOLN_FILETYPE = '.r'
	CHARACTER(LEN=200) 	:: Solution_File_Name, flag = 'periodic'

	REAL(KIND=8), ALLOCATABLE :: zvec(:,:,:), rvec(:,:,:), pv(:,:,:,:), pv_old(:,:,:,:), pv_target(:,:,:,:), pv_targetI(:,:,:,:), pv_targetV(:,:,:,:), pvforSoln(:,:,:,:)
	REAL(KIND=8), ALLOCATABLE :: ibval(:,:,:)
	REAL(KIND=8), ALLOCATABLE :: k1(:,:,:,:), k2(:,:,:,:), k3(:,:,:,:), k4(:,:,:,:)
	REAL(KIND=8), ALLOCATABLE :: dxidz(:,:,:), detadr(:,:,:), dzdxi(:,:,:), drdeta(:,:,:), ddxidxidz(:,:,:), ddetadetadr(:,:,:)
	REAL(KIND=8), ALLOCATABLE :: rhovec(:,:,:), Vrvec(:,:,:), Vtvec(:,:,:), Vzvec(:,:,:), Tvec(:,:,:), pvec(:,:,:), muvec(:,:,:), kappavec(:,:,:)
	REAL(KIND=8), ALLOCATABLE :: lambdavec(:,:,:), betavec(:,:,:)
	REAL(KIND=8), ALLOCATABLE :: rhovecTGfilt(:,:,:)
	REAL(KIND=8), ALLOCATABLE :: drhovecdz(:,:,:), dVrvecdz(:,:,:), dVtvecdz(:,:,:), dVzvecdz(:,:,:), dTvecdz(:,:,:), dpvecdz(:,:,:), dmuvecdz(:,:,:)
	REAL(KIND=8), ALLOCATABLE :: drhovecdr(:,:,:), dVrvecdr(:,:,:), dVtvecdr(:,:,:), dVzvecdr(:,:,:), dTvecdr(:,:,:), dpvecdr(:,:,:), dmuvecdr(:,:,:)
	REAL(KIND=8), ALLOCATABLE :: dkappavecdz(:,:,:), dkappavecdr(:,:,:)
	REAL(KIND=8), ALLOCATABLE ::  d2Vrvecdz2(:,:,:), d2Vtvecdz2(:,:,:), d2Vzvecdz2(:,:,:), d2Tvecdz2(:,:,:)
	REAL(KIND=8), ALLOCATABLE ::  d2Vrvecdr2(:,:,:), d2Vtvecdr2(:,:,:), d2Vzvecdr2(:,:,:), d2Tvecdr2(:,:,:)
	REAL(KIND=8), ALLOCATABLE :: d2Vrvecdzdr(:,:,:), d2Vzvecdzdr(:,:,:)
	REAL(KIND=8), ALLOCATABLE :: drhovecdr_initial(:,:,:), dVzvecdr_initial(:,:,:), dTvecdr_initial(:,:,:)
	REAL(KIND=8) :: z, r, rho, Vr, Vt, Vz, T, mu, kappa
	REAL(KIND=8) :: drhodr, dVrdr, dVtdr, dVzdr, dTdr, dpdr, dmudr, dkappadr
	REAL(KIND=8) :: drhodz, dVrdz, dVtdz, dVzdz, dTdz, dpdz, dmudz, dkappadz
	REAL(KIND=8) :: d2Vrdr2, d2Vtdr2, d2Vzdr2, d2Tdr2
	REAL(KIND=8) ::  d2Vrdz2, d2Vtdz2, d2Vzdz2, d2Tdz2
	REAL(KIND=8) :: d2Vrdzdr, d2Vzdzdr



	REAL(KIND=8), ALLOCATABLE :: SBP1(:,:), SBP2(:,:), TGfilt(:,:)

	! Constants

	REAL(KIND=8) :: pi = 4.0D0*ATAN(1.0D0), gam = 1.4D0
	REAL(KIND=8) :: c = 1.0D0, total_time = 40.0D0, delt = 0.00025D0
	INTEGER(KIND=8) :: no_iterations
	INTEGER :: noutput = 200
	
	INTEGER(KIND=8) :: do_viscous_terms = 1, write_rhs_file = 1, do_mms = 0
	INTEGER(KIND=8) :: NSCBC=1
	REAL(KIND=8) :: Re = 1000.0D0, Pr = 0.72D0
	REAL(KIND=8) :: fourby3 = 29.0D0/15.0D0, oneby3 = 14.0D0/15.0D0, twoby3 = 1.0D0/15.0D0, sevenby3 = 44.0D0/15.0D0
	REAL(KIND=8) :: onebyr
 

	REAL(KIND=8) :: zeta
	REAL(KIND=8), PARAMETER :: A_sponge = 5.0D0, n_sponge = 2
		
	
	
	! Variables for NSCBC and SAT boundary conditions

	REAL(KIND=8) :: lam(5,5), T_mat(5,5), T_mat_inv(5,5), deriv_vector(5), waveamp(5), rhs_NSCBC(5), rhs_SAT(5)
	INTEGER(KIND=8) :: counter
        REAL(KIND=8), PARAMETER :: sigma_1 = -2.0D0, sigma_2 = -2.0D0
        REAL(KIND=8), PARAMETER :: one_by_h_11 = 48.0D0/17.0D0, one_by_h_44 = 48.0D0/49.0D0


	! Error analysis

	REAL(KIND=8) :: rms_error
	REAL(KIND=8), ALLOCATABLE :: error(:)

	! Patch variables

	INTEGER(KIND=8), ALLOCATABLE :: is_it_i_patch(:,:,:), is_it_j_patch(:,:,:), which_i_patch(:,:,:), which_j_patch(:,:,:)
	INTEGER(KIND=8) :: i_start, i_end, j_start, j_end, bc, bc_dir, no_bcs, patch_no
	INTEGER(KIND=8) :: dummy_int, reason

	TYPE patch_type

	INTEGER(KIND=8) :: i_start, i_end, j_start, j_end, bc, bc_dir

	END TYPE patch_type

	TYPE(patch_type), ALLOCATABLE :: patch(:)

	! Boundary conditions

	INTEGER(KIND=8) :: AXIS = 90
	INTEGER(KIND=8) :: SAT_FARFIELD = 91
	INTEGER(KIND=8) :: NSCBC_NONREFLECTING = 92
	INTEGER(KIND=8) :: SAT_NOSLIP_WALL = 93
	INTEGER(KIND=8) :: NSCBC_NOSLIP_WALL = 94
	INTEGER(KIND=8) :: SPONGE = 95
	

	
	! Variables for parallelization

	INTEGER SIZE, UP, DOWN, LEFT, RIGHT
   	PARAMETER(SIZE=16)
   	PARAMETER(UP=1)
    PARAMETER(DOWN=2)
    PARAMETER(LEFT=3)
    PARAMETER(RIGHT=4)
    INTEGER numtasks, rank, source, dest, outbuf, tag, ierr, inbuf(4), nbrs(4), dimsforcartcomm(2), &
            coords(2), coordsforcartcomm(2), stats(MPI_STATUS_SIZE, 8), reqs(8), cartcomm, periods(2), reorder, procrank
    DATA inbuf /MPI_PROC_NULL,MPI_PROC_NULL,MPI_PROC_NULL, MPI_PROC_NULL/, dimsforcartcomm /0,0/, tag /1/, periods /0,0/, reorder /1/ 
		

	INTEGER(KIND=8) :: Nz, Nr, istart, iend, jstart, jend, i_indx, j_indx, bcistart, bciend, bcjstart, bcjend, bctype, bcdir, bcno, total_no_bcs, bccount
	INTEGER(KIND=8) :: patchistart, patchiend, patchjstart, patchjend, value1, value2, exitvalue

	INTEGER(KIND=8) :: procbcistart, procbciend, procbcjstart, procbcjend
	INTEGER(KIND=8) :: istartval, iendval, jstartval, jendval
	CHARACTER(LEN=80) :: sol_file
	


	INTEGER(KIND=8), ALLOCATABLE :: startenddata(:,:)
	INTEGER(KIND=8) :: dummyvec(5)
	integer stat(MPI_STATUS_SIZE)
	

	REAL(KIND=8) :: time_start

	! Variables for Shock Capturing Scheme

	INTEGER(KIND=8) :: ShockCapturing = 0
	INTEGER(KIND=8) :: rkStep

	REAL(KIND=8), ALLOCATABLE :: MgnStrRate(:,:,:), d2MgnStrRatedxi2(:,:,:), d2MgnStrRatedeta2(:,:,:), d4MgnStrRatedxi4(:,:,:), d4MgnStrRatedeta4(:,:,:)
	REAL(KIND=8), ALLOCATABLE :: d2Tvecdxi2(:,:,:), d2Tvecdeta2(:,:,:),d4Tvecdxi4(:,:,:), d4Tvecdeta4(:,:,:)
	REAL(KIND=8), ALLOCATABLE :: dilvec(:,:,:), d2dilvecdxi2(:,:,:), d2dilvecdeta2(:,:,:),d4dilvecdxi4(:,:,:), d4dilvecdeta4(:,:,:)
	REAL(KIND=8), ALLOCATABLE :: mustarQuantity(:,:,:), mustar(:,:,:), kappastarQuantity(:,:,:), kappastar(:,:,:)
	REAL(KIND=8), ALLOCATABLE :: betastarQuantity(:,:,:), betastar(:,:,:)

	REAL(KIND=8) :: dilatation, minus_dilatation, vorticity, f_sw, H_of_minus_dilatation


	REAL(KIND=8), PARAMETER :: C_mu = 0.000D0, C_beta = 0.0D0, C_kappa = 0.0D0
	REAL(KIND=8) :: S11, S12, S13, S21, S22, S23, S31, S32, S33
	REAL(KIND=8) :: deltaz, deltar

	REAL(KIND=8) :: tmp

	REAL(KIND=8), ALLOCATABLE :: tmp_vector1(:,:,:), tmp_vector2(:,:,:), tmp_vector3(:,:,:)


	! Time Step Calculations

	REAL(KIND=8) ::	delt_inv_term, delt_visc_term, delt_denom, min_delt, time

	REAL(KIND=8) :: CFL = 0.5D0

	REAL(KIND=8), ALLOCATABLE :: deltvec(:,:,:)

	    ! Variables for control forcing

        REAL(KIND=8) :: forcing_mat(5,5), z_loc, r_loc, z0, r0, lz, lr, GaussianFactor
        REAL(KIND=8), PARAMETER :: alpha = 6.0D0
        INTEGER(KIND=8), PARAMETER :: ControlForcing = 0

	


END MODULE ModuleVariables

