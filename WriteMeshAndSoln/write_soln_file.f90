program write_soln_file

use ModPlot3D_IO 

implicit none

	integer		 :: NDIM, ngrid, ftype, grid_NO
	integer,pointer	 :: IBLANK(:,:,:,:)
	integer,pointer	 :: ND(:,:)
	character(2)		:: prec, gf, vf, ib
	integer, parameter	  :: RFREAL = SELECTED_REAL_KIND(PRECISION(1.0d0))
	REAL(rfreal)		:: tau(4)
	REAL(rfreal), pointer   :: G(:,:,:,:,:), X3(:,:,:,:,:) 
	integer		 :: OF = 1

	double precision, allocatable :: rho_rhoinf(:), u_Uj(:), u_cinf(:), T_Tinf(:), rad(:)
	double precision :: gam, Mj, Tj_Tinf, Uj_cinf, theta
	
	double precision :: term1, term2, term3, rho_rhoj

	integer(KIND=8)	 :: i, j, k, Nz, Nr

	character(LEN=200) :: GRID_FILENAME = "axijet_mesh.xyz"
	character(LEN=200) :: SOLUTION_FILENAME = "axijet_soln.q"
	
	! Create Grid
	
	call create_AxiJet_CartGrid(G, GRID_FILENAME)
	
	! Create initial condition

	call Read_Grid(NDIM, ngrid, ND, G, IBLANK, prec, gf, vf, ib, GRID_FILENAME, OF)	

	Nz = ND(1,1)
	Nr = ND(1,2)

	allocate(rho_rhoinf(Nr))
	allocate(u_Uj(Nr))
	allocate(u_cinf(Nr))
	allocate(T_Tinf(Nr))
	allocate(rad(Nr))


	allocate(X3(1,ND(1,1),ND(1,2),ND(1,3),5)) 	

	X3 = 0.0d0

	gam	 = 1.4d0	
	Mj	  = 0.8d0
	Tj_Tinf = 0.6d0
	Uj_cinf = Mj*(Tj_Tinf)**0.5d0	
	theta = 0.1


	do j = 1, Nr
		rad(j) = G(1,1,j,1,2)
	end do

	do i = 1, Nz
	   do j = 1, Nr
		 u_Uj(j) = (1.0d0/2.0d0)*( 1.0d0 - tanh(1.0d0/(4.0d0*theta)*(rad(j) - 1.0d0/(rad(j)))) )
		 u_Uj(j) = 2.0d0*u_Uj(j) - 1.0d0;
		 if(G(1,i,1,1,1) .gt. 6.0)then
		 	u_Uj(j) = u_Uj(j)/(G(1,i,1,1,1)-6.0+1.0)**1.5	
		 endif
		 if(G(1,1,j,1,2) .gt. 1.0d0)u_Uj(j) = 0.0 
		 u_cinf(j) = u_Uj(j)*Uj_cinf
		 term1 = (1.0d0/2.0d0)*(gam -1)*u_Uj(j)*(1 - u_Uj(j))*Mj**2 
		 term2 = u_Uj(j) 
		 term3 = (1.0d0/Tj_Tinf) * (1 - u_Uj(j)) 
		 rho_rhoj   = 1.0d0/(term1 + term2 + term3)
		 rho_rhoinf(j)  = rho_rhoj / Tj_Tinf

		  do k = 1, 1	

		 X3(1,i,j,k,1) = rho_rhoinf(j)		
		 X3(1,i,j,k,2) = 0.0d0
		 X3(1,i,j,k,3) = 0.0d0
		 X3(1,i,j,k,4) = u_cinf(j)
		 X3(1,i,j,k,5) = 1.0d0/rho_rhoinf(j)

	
		  end do
	   end do
	end do

	tau = 0.0d0

	call Write_Soln(NDIM, ngrid, ND, X3, tau, prec, gf, vf, SOLUTION_FILENAME)

	

contains 

subroutine compute_coord_val(i, nx, delta_x, pos_x, coord_min, coord_max, coord_val)

	integer, intent(in) :: i, nx
	double precision, intent(in) :: delta_x, coord_min, coord_max
	double precision :: pos_x, pos_x1

	double precision, intent(out) :: coord_val

	double precision :: xi, u1, u2, fac


    pos_x1 = pos_x/(coord_max-coord_min);

    xi = float(i)/float(nx);
    u1 = tanh(delta_x*(1.0-xi))*(1.0-pos_x1);
    u2 = (2.0-tanh(delta_x*xi))*pos_x1;
    fac = 1.0 - ((u1+u2)-pos_x1);

    coord_val = coord_min + fac*(coord_max-coord_min);

end subroutine compute_coord_val

subroutine create_AxiJet_CartGrid(G, GRID_FILENAME)
			
	integer :: i, j, k, ipos2;
    double precision :: xpos2 = 11.0;
    logical :: found_index = .false.;
	
	REAL(rfreal), pointer   :: G(:,:,:,:,:)

	integer :: nx = 150, ny = 50, nz = 1
	integer :: npoints_imp_wall_zone = 100	
 	double precision :: xmin = 0.0, xmax = 15.0;
	double precision :: ymin = 0.0, ymax = 15.0;

	double precision :: delta_y = 6.0;
    double precision :: pos_y = 1.0;

	double precision :: xval, yval

	character(len=200), intent(in) :: GRID_FILENAME


	integer		 :: NDIM, ngrid, ftype, grid_NO
	integer,pointer	 :: IBLANK(:,:,:,:)
	integer,pointer	 :: ND(:,:)
	character(2)		:: prec, gf, vf, ib
	integer, parameter	  :: RFREAL = SELECTED_REAL_KIND(PRECISION(1.0d0))
	REAL(rfreal)		:: tau(4)
	integer		 :: OF = 1

	integer :: jbottom, jtop, iright, iwall


	allocate(G(1,nx+npoints_imp_wall_zone,ny,nz,3))
	allocate(IBLANK(1,nx+npoints_imp_wall_zone,ny,nz))

    do i = 1, nx

        call compute_coord_val(i,nx,3.0d0,6.0d0,xmin,xmax,xval)

        if(xval .gt. xpos2)then
            ipos2 = i
            found_index = .true.
            exit
        end if

        do j= 1, ny
            call compute_coord_val(j,ny,delta_y,pos_y,ymin,ymax,yval);
            do k=1, nz
                G(1,i,j,k,1) = xval
                G(1,i,j,k,2) = yval
                G(1,i,j,k,3) = 0.0d0
            end do
        end do
    end do

    do i=ipos2, nx+npoints_imp_wall_zone

        call compute_coord_val(i-ipos2,nx+npoints_imp_wall_zone-ipos2,2.0d0,3.0d0,xpos2,xmax,xval);

        do j=1, ny
            call compute_coord_val(j,ny,delta_y,pos_y,ymin,ymax,yval);
            do k=1, nz
                G(1,i,j,k,1) = xval
                G(1,i,j,k,2) = yval
                G(1,i,j,k,3) = 0.0
			end do
		end do
	end do

	! write iblank data to mesh

	do i = 1, nx+npoints_imp_wall_zone
		do j=1, ny
			do k=1, nz
				xval = G(1,i,j,k,1)
				yval = G(1,i,j,k,2)
				if(yval .gt. 1.0d0 .and. yval .lt. 1.05d0 .and. xval .lt. 6.0d0)then
					IBLANK(1,i,j,k) = 0
				else	
                    IBLANK(1,i,j,k) = 1
                end if
			end do
		end do
	end do
			
	

	NDIM = 3;
	ngrid = 1;
	prec(1:1) = 'd'
	gf(1:1) = 'm'
	vf(1:1) = 'w'
	ib(1:1) = 'y'
	
	allocate(ND(1,3))
	ND(1,1) = nx+npoints_imp_wall_zone
	ND(1,2) = ny
	ND(1,3) = 1

 
	call Write_Grid(NDIM, ngrid, ND, G, IBLANK, prec, gf, vf, ib, GRID_FILENAME)	

	print*, "Done writing mesh" 
	
	! Write bc file

	do j=1, ND(1,2)
		yval = G(1,1,j,1,2)
		if(yval .gt. 1.0d0)then
			jbottom = j-1	
			exit
		endif
	end do

	do j= ND(1,2), 1, -1
		yval = G(1,1,j,1,2)
		if(yval .lt. 1.05d0)then
			jtop = j+1	
			exit
		endif
	end do

	do i = ND(1,1), 1, -1
		xval = G(1,i,1,1,1)
		if(xval .lt. 6.0d0)then
			iright = i+1
			exit
		endif
	end do


	open(unit=10,file="bc.dat")
	write(10,'(9(I0,2X))') 1, 91, 1, 1, 1, 1, jbottom, 1, 1
	write(10,'(9(I0,2X))') 1, 94, -2, 1, iright, jbottom, jbottom, 1, 1
    write(10,'(9(I0,2X))') 1, 94, 1, iright,  iright, jbottom, jtop, 1, 1
	write(10,'(9(I0,2X))') 1, 94, 2, 1, iright, jtop, jtop, 1, 1
	write(10,'(9(I0,2X))') 1, 91, 1, 1, 1, jtop,  -1, 1, 1
	write(10,'(9(I0,2X))') 1, 92, -2, 1, -1, -1, -1, 1, 1
	write(10,'(9(I0,2X))') 1, 94, -1, -1, -1, 1, -1, 1, 1
	write(10,'(9(I0,2X))') 1, 90, 2, 1 , -1, 1, 1, 1, 1 
	close(10)	


			


end subroutine create_AxiJet_CartGrid

end program write_soln_file
