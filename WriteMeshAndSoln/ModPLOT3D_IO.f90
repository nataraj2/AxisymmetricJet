! Copyright (c) 2015, University of Illinois at Urbana-Champaign, XPACC
! License: MIT, http://opensource.org/licenses/MIT
Module ModPLOT3D_IO

  Implicit None

Contains

  Subroutine Read_Grid_Size(NDIM, ngrid, ND, file, OUTPUT_FLAG, ib)

    Implicit None

    Character(LEN=2), Intent(out), Optional :: ib
    Real(KIND=8), Dimension(:,:,:,:,:), Pointer :: X
    Integer, Dimension(:,:,:,:), Pointer :: IBLANK
    Character(LEN=2) :: prec, gf, vf, iblocal
    Character(LEN=*) :: file
    Integer :: ngrid, I, J, K, L, Nx, Ny, Nz, NDIM, M, ftype
    Integer, Pointer :: ND(:,:)
    Integer :: OMap(3), NMap(3), ND_copy(3)
    Real(KIND=4), Dimension(:,:,:,:), Pointer :: fX
    Logical :: gf_exists
    Integer :: OUTPUT_FLAG

    Inquire(file=Trim(file),exist=gf_exists)
    If (.NOT.(gf_exists.eqv..true.)) Then
      Write (*,'(A,A,A)') 'File ', file(1:LEN_TRIM(file)), ' is not readable'
      Stop
    End If

    Call p3d_detect(LEN(Trim(file)),Trim(file),prec,NDIM,gf,vf,iblocal,ftype)
    If (OUTPUT_FLAG .gt. 0) Call output_format(file,prec,NDIM,gf,vf,iblocal,ftype)
    If (PRESENT(ib) .eqv. .true.) ib(1:1) = iblocal(1:1)

    Open (unit=10, file=trim(file), form='unformatted', status='old')

    !  Read number of grids
    If (gf(1:1) .eq. 's') Then
      ngrid = 1
    Else
      Read (10) ngrid
    End If
    Allocate(ND(ngrid,3))

    If (OUTPUT_FLAG > 0) Then
      Write (*,'(A)') ' '
      Write (*,'(A,I2)') 'RFLOCM: Number of grids: ', ngrid
    End If

    !  Read grid sizes
    ND(:,3) = 1
    Read (10) ((ND(I,J),J=1,NDIM),I=1,ngrid)
    If (OUTPUT_FLAG > 0) Then
      Do I = 1, ngrid
        Write (*,'(A,I2,A,3(I5,1X))') 'RFLOCM: Grid ', I, ' is ', (ND(I,J),J=1,NDIM)
      End Do
    End IF

    close (10)

  End Subroutine Read_Grid_Size

  Subroutine Read_Soln_Header(NDIM, ngrid, ND, fname, tau)

    Implicit None

    ! ... function call variables
    Integer :: NDIM, ngrid
    Integer, Dimension(:,:), Pointer :: ND
    Character(LEN=*) :: fname
    Real(KIND=8) :: tau(4)

    ! ... local variables
    Integer :: ibuf, I, J

    ! ... open the file
    Open (Unit=10, file=trim(fname), form='unformatted', status='old')

    ! ... read the number of grids, exit if error
    Read (10) ibuf
    ! If (ibuf /= ngrid) Then
    !   Write (*,*) 'RFLOCM: ERROR: assumed and read ngrid do not match for file '//trim(fname)
    !   Write (*,'(A,I4)') 'RFLOCM: ERROR: assumed ngrid = ', ngrid
    !   Write (*,'(A,I4)') 'RFLOCM: ERROR: read ngrid = ', ibuf
    !   STOP 'RFLOCM: ERROR: Stopping in read_soln_header()'
    ! End If

    ! ... read the size of each grid, throw value away
    Read (10) ((ibuf,J=1,NDIM),I=1,ngrid)

    ! ... read the header
    Read (10) tau

    ! ... close and exit
    Close (10)

  End Subroutine Read_Soln_Header

  Subroutine Read_Func_Header(NDIM, ngrid, fname, nvar)

    Implicit None

    Character(LEN=*) :: fname
    Integer :: ngrid, I, J, NDIM, nvar, ibuf
    Integer :: ND2(ngrid, NDIM+1)
    Integer, Parameter :: funit = 10

    Open (unit=funit, file=trim(fname), form='unformatted', status='old')

    !  Write number of grids
    Read (funit) ibuf

    !  Read grid sizes & nvars
    Read (funit) ((ND2(I,J),J=1,NDIM+1),I=1,ngrid)

    !  Save nvar
    nvar = ND2(ngrid,NDIM+1) 

    Close (funit)

  End Subroutine Read_Func_Header

  Subroutine Read_Grid(NDIM, ngrid, ND, X, IBLANK, prec, gf, vf, ib, file, OF)

    Implicit None

    Real(KIND=8), Dimension(:,:,:,:,:), Pointer :: X
    Integer, Dimension(:,:,:,:), Pointer :: IBLANK
    Character(LEN=2) :: prec, gf, vf, ib
    Character(LEN=80) :: file
    Integer :: ngrid, I, J, K, L, Nx, Ny, Nz, NDIM, M, ftype
    Integer, Dimension(:,:), Pointer :: ND
    Integer :: OMap(3), NMap(3), ND_copy(3)
    Real(KIND=4), Dimension(:,:,:,:), Pointer :: fX
    Logical :: gf_exists
    Integer :: OF

    Inquire(file=Trim(file),exist=gf_exists)
    If (.NOT.gf_exists) Then
      Write (*,'(A,A,A)') 'File ', file(1:LEN_TRIM(file)), ' is not readable'
      Stop
    End If

    Call p3d_detect(LEN(Trim(file)),Trim(file),prec,NDIM,gf,vf,ib,ftype)
    If (OF > 0) Call output_format(file,prec,NDIM,gf,vf,ib,ftype)

    Open (unit=10, file=trim(file), form='unformatted', status='old')

    !  Read number of grids
    If (gf(1:1) .eq. 's') Then
      ngrid = 1
    Else
      Read (10) ngrid
    End If
    Allocate(ND(ngrid,3))
    If (OF > 0) Write (*,'(A)') ' '
    If (OF > 0)Write (*,'(A,I2)') 'RFLOCM: Number of grids: ', ngrid

    !  Read grid sizes
    ND(:,3) = 1
    Read (10) ((ND(I,J),J=1,NDIM),I=1,ngrid)
    If (OF > 0) Then
      Do I = 1, ngrid
        Write (*,'(A,I2,A,3(I5,1X))') 'RFLOCM: Grid ', I, ' is ', (ND(I,J),J=1,NDIM)
      End Do
    End If

    !  Allocate grid memory
    Nx = MAXVAL(ND(:,1))
    Ny = MAXVAL(ND(:,2))
    Nz = MAXVAL(ND(:,3))
    Allocate(X(ngrid,Nx,Ny,Nz,NDIM))
    If (prec(1:1) .eq. 's') Allocate(fX(Nx,Ny,Nz,NDIM))
    Allocate(IBLANK(ngrid,Nx,Ny,Nz))

    !  Read grid values
    If (vf(1:1) .eq. 'w') Then
      If (prec(1:1) .eq. 'd') Then
        If (ib(1:1) .eq. 'n') Then
          Do L = 1, ngrid
            Read (10) ((((X(L,I,J,K,M),I=1,ND(L,1)),J=1,ND(L,2)),K=1,ND(L,3)),M=1,NDIM)
          End Do
        Else
          Do L = 1, ngrid
            Read (10) ((((X(L,I,J,K,M),I=1,ND(L,1)),J=1,ND(L,2)),K=1,ND(L,3)),M=1,NDIM), &
                 (((IBLANK(L,I,J,K),I=1,ND(L,1)),J=1,ND(L,2)),K=1,ND(L,3))
          End Do
        End If
      Else
        If (ib(1:1) .eq. 'n') Then
          Do L = 1, ngrid
            Read (10) ((((fX(I,J,K,M),I=1,ND(L,1)),J=1,ND(L,2)),K=1,ND(L,3)),M=1,NDIM)
            X(L,:,:,:,:) = Dble(fX(:,:,:,:))
          End Do
        Else
          Do L = 1, ngrid
            Read (10) ((((fX(I,J,K,M),I=1,ND(L,1)),J=1,ND(L,2)),K=1,ND(L,3)),M=1,NDIM), &
                 (((IBLANK(L,I,J,K),I=1,ND(L,1)),J=1,ND(L,2)),K=1,ND(L,3))
            X(L,:,:,:,:) = Dble(fX(:,:,:,:))
          End Do
        End If
      End If

    Else
      If (NDIM .EQ. 3) Then
        If (prec(1:1) .eq. 'd') Then
          If (ib(1:1) .eq. 'n') Then
            Do L = 1, ngrid
              Do K = 1, ND(L,3)
                Read (10) ((X(L,I,J,K,1),I=1,ND(L,1)),J=1,ND(L,2)), &
                     ((X(L,I,J,K,2),I=1,ND(L,1)),J=1,ND(L,2)), &
                     ((X(L,I,J,K,3),I=1,ND(L,1)),J=1,ND(L,2))
              End Do
            End Do
          Else
            Do L = 1, ngrid
              Do K = 1, ND(L,3)
                Read (10) ((X(L,I,J,K,1),I=1,ND(L,1)),J=1,ND(L,2)), &
                     ((X(L,I,J,K,2),I=1,ND(L,1)),J=1,ND(L,2)), &
                     ((X(L,I,J,K,3),I=1,ND(L,1)),J=1,ND(L,2)), &
                     ((IBLANK(L,I,J,K),I=1,ND(L,1)),J=1,ND(L,2))
              End Do
            End Do
          End If
        Else
          If (ib(1:1) .eq. 'n') Then
            Do L = 1, ngrid
              Do K = 1, ND(L,3)
                Read (10) ((fX(I,J,K,1),I=1,ND(L,1)),J=1,ND(L,2)), &
                     ((fX(I,J,K,2),I=1,ND(L,1)),J=1,ND(L,2)), &
                     ((fX(I,J,K,3),I=1,ND(L,1)),J=1,ND(L,2))
              End Do
              X(L,:,:,:,:) = Dble(fX(:,:,:,:))
            End Do
          Else
            Do L = 1, ngrid
              Do K = 1, ND(L,3)
                Read (10) ((fX(I,J,K,1),I=1,ND(L,1)),J=1,ND(L,2)), &
                     ((fX(I,J,K,2),I=1,ND(L,1)),J=1,ND(L,2)), &
                     ((fX(I,J,K,3),I=1,ND(L,1)),J=1,ND(L,2)), &
                     ((IBLANK(L,I,J,K),I=1,ND(L,1)),J=1,ND(L,2))
              End Do
              X(L,:,:,:,:) = Dble(fX(:,:,:,:))
            End Do
          End If
        End If
      Else
        If (prec(1:1) .eq. 'd') Then
          If (ib(1:1) .eq. 'n') Then
            Do L = 1, ngrid
              Do J = 1, ND(L,2)
                Read (10) (X(L,I,J,1,1),I=1,ND(L,1)), &
                     (X(L,I,J,1,2),I=1,ND(L,1))
              End Do
            End Do
          Else
            Do L = 1, ngrid
              Do J = 1, ND(L,2)
                Read (10) (X(L,I,J,1,1),I=1,ND(L,1)), &
                     (X(L,I,J,1,2),I=1,ND(L,1)), &
                     (IBLANK(L,I,J,1),I=1,ND(L,1))
              End Do
            End Do
          End If
        Else
          If (ib(1:1) .eq. 'n') Then
            Do L = 1, ngrid
              Do K = J, ND(L,2)
                Read (10) (fX(I,J,1,1),I=1,ND(L,1)), &
                     (fX(I,J,1,2),I=1,ND(L,1))
              End Do
              X(L,:,:,:,:) = Dble(fX(:,:,:,:))
            End Do
          Else
            Do L = 1, ngrid
              Do K = J, ND(L,2)
                Read (10) (fX(I,J,1,1),I=1,ND(L,1)), &
                     (fX(I,J,1,2),I=1,ND(L,1)), &
                     (IBLANK(L,I,J,1),I=1,ND(L,1))
              End Do
              X(L,:,:,:,:) = Dble(fX(:,:,:,:))
            End Do
          End If
        End If
      End If
    End If

    Close (10)

    If (prec(1:1) .eq. 's')  Deallocate(fX)
    If (ib(1:1) .eq. 'n') IBLANK = 1

    Return
  End Subroutine Read_Grid

  Subroutine Write_Grid(NDIM, ngrid, ND, X, IBLANK, prec, gf, vf, ib, file)

    Implicit None

    Real(KIND=8), Dimension(:,:,:,:,:), Pointer :: X
    Integer, Dimension(:,:,:,:), Pointer :: IBLANK
    Character(LEN=2) :: prec, gf, vf, ib
    Character(LEN=80) :: file
    Integer :: ngrid, I, J, K, L, Nx, Ny, Nz, NDIM, M, ftype
    Integer, Dimension(:,:), Pointer :: ND
    Integer :: OMap(3), NMap(3), ND_copy(3)
    Real(KIND=4), Dimension(:,:,:,:), Pointer :: fX

    ftype = 0
    Call output_format(file,prec,NDIM,gf,vf,ib,ftype)

    Open (unit=10, file=trim(file), form='unformatted', status='unknown')

    !  Write number of grids
    If (gf(1:1) .eq. 'm') Then
      Write (10) ngrid
    End If

    !  Write grid sizes
    Write (10) ((ND(I,J),J=1,NDIM),I=1,ngrid)

    !  Read grid values
    If (vf(1:1) .eq. 'w') Then
      If (prec(1:1) .eq. 'd') Then
        If (ib(1:1) .eq. 'n') Then
          Do L = 1, ngrid
            Write (10) ((((X(L,I,J,K,M),I=1,ND(L,1)),J=1,ND(L,2)),K=1,ND(L,3)),M=1,NDIM)
          End Do
        Else
          Do L = 1, ngrid
            Write (10) ((((X(L,I,J,K,M),I=1,ND(L,1)),J=1,ND(L,2)),K=1,ND(L,3)),M=1,NDIM), &
                     (((IBLANK(L,I,J,K),I=1,ND(L,1)),J=1,ND(L,2)),K=1,ND(L,3))
          End Do
        End If
      Else
        If (ib(1:1) .eq. 'n') Then
          Do L = 1, ngrid
            Write (10) ((((Real(X(L,I,J,K,M)),I=1,ND(L,1)),J=1,ND(L,2)),K=1,ND(L,3)),M=1,NDIM)
          End Do
        Else
          Do L = 1, ngrid
            Write (10) ((((Real(X(L,I,J,K,M)),I=1,ND(L,1)),J=1,ND(L,2)),K=1,ND(L,3)),M=1,NDIM), &
                      (((IBLANK(L,I,J,K),I=1,ND(L,1)),J=1,ND(L,2)),K=1,ND(L,3))
          End Do
        End If
      End If

    Else
      If (NDIM .EQ. 3) Then
        If (prec(1:1) .eq. 'd') Then
          If (ib(1:1) .eq. 'n') Then
            Do L = 1, ngrid
              Do K = 1, ND(L,3)
                Write (10) ((X(L,I,J,K,1),I=1,ND(L,1)),J=1,ND(L,2)), &
                     ((X(L,I,J,K,2),I=1,ND(L,1)),J=1,ND(L,2)), &
                     ((X(L,I,J,K,3),I=1,ND(L,1)),J=1,ND(L,2))
              End Do
            End Do
          Else
            Do L = 1, ngrid
              Do K = 1, ND(L,3)
                Write (10) ((X(L,I,J,K,1),I=1,ND(L,1)),J=1,ND(L,2)), &
                     ((X(L,I,J,K,2),I=1,ND(L,1)),J=1,ND(L,2)), &
                     ((X(L,I,J,K,3),I=1,ND(L,1)),J=1,ND(L,2)), &
                     ((IBLANK(L,I,J,K),I=1,ND(L,1)),J=1,ND(L,2))
              End Do
            End Do
          End If
        Else
          If (ib(1:1) .eq. 'n') Then
            Do L = 1, ngrid
              Do K = 1, ND(L,3)
                Write (10) ((Real(X(L,I,J,K,1)),I=1,ND(L,1)),J=1,ND(L,2)), &
                     ((Real(X(L,I,J,K,2)),I=1,ND(L,1)),J=1,ND(L,2)), &
                     ((Real(X(L,I,J,K,3)),I=1,ND(L,1)),J=1,ND(L,2))
              End Do
            End Do
          Else
            Do L = 1, ngrid
              Do K = 1, ND(L,3)
                Write (10) ((Real(X(L,I,J,K,1)),I=1,ND(L,1)),J=1,ND(L,2)), &
                     ((Real(X(L,I,J,K,2)),I=1,ND(L,1)),J=1,ND(L,2)), &
                     ((Real(X(L,I,J,K,3)),I=1,ND(L,1)),J=1,ND(L,2)), &
                     ((IBLANK(L,I,J,K),I=1,ND(L,1)),J=1,ND(L,2))
              End Do
            End Do
          End If
        End If
      Else
        If (prec(1:1) .eq. 'd') Then
          If (ib(1:1) .eq. 'n') Then
            Do L = 1, ngrid
              Do J = 1, ND(L,2)
                Write (10) (X(L,I,J,1,1),I=1,ND(L,1)), &
                     (X(L,I,J,1,2),I=1,ND(L,1))
              End Do
            End Do
          Else
            Do L = 1, ngrid
              Do J = 1, ND(L,2)
                Write (10) (X(L,I,J,1,1),I=1,ND(L,1)), &
                     (X(L,I,J,1,2),I=1,ND(L,1)), &
                     (IBLANK(L,I,J,1),I=1,ND(L,1))
              End Do
            End Do
          End If
        Else
          If (ib(1:1) .eq. 'n') Then
            Do L = 1, ngrid
              Do K = J, ND(L,2)
                Write (10) (Real(X(L,I,J,1,1)),I=1,ND(L,1)), &
                     (Real(X(L,I,J,1,2)),I=1,ND(L,1))
              End Do
            End Do
          Else
            Do L = 1, ngrid
              Do K = J, ND(L,2)
                Write (10) (Real(X(L,I,J,1,1)),I=1,ND(L,1)), &
                     (Real(X(L,I,J,1,2)),I=1,ND(L,1)), &
                     (IBLANK(L,I,J,1),I=1,ND(L,1))
              End Do
            End Do
          End If
        End If
      End If
    End If

    Close (10)

    Return
  End Subroutine Write_Grid


  Subroutine Read_Soln(NDIM, ngrid, ND, X, tau, prec, gf, vf, file, OF)

    Implicit None

    Real(KIND=8), Dimension(:,:,:,:,:), Pointer :: X
    Real(KIND=8) :: tau(4)
    Character(LEN=2) :: prec, gf, vf, ib
    Character(LEN=80) :: file
    Integer :: ngrid, I, J, K, L, M, Nx, Ny, Nz, NDIM, ftype
    Integer, Dimension(:,:), Pointer :: ND
    Real(KIND=4), Dimension(:,:,:,:), Pointer :: fX
    Real(KIND=4) :: ftau(4)
    Logical :: gf_exists
    Integer :: OF

    Inquire(file=Trim(file),exist=gf_exists)
    If (.NOT.gf_exists) Then
      Write (*,'(A,A,A)') 'File ', file(1:LEN_TRIM(file)), ' is not readable'
      Stop
    End If

    Call p3d_detect(LEN(Trim(file)),Trim(file),prec,NDIM,gf,vf,ib,ftype)
    If (OF > 0) Call output_format(file,prec,NDIM,gf,vf,ib,ftype)

    Open (unit=10, file=trim(file), form='unformatted', status='old')

    !  Read number of grids
    If (gf(1:1) .eq. 's') Then
      ngrid = 1
    Else
      Read (10) ngrid
    End If

    If (OF > 0) Then
      Write (*, *) ' '
      Write (*,'(A,I2)') 'RFLOCM: Number of grids: ', ngrid
    End If

    !  Read grid sizes
    Allocate(ND(ngrid,3))
    ND(:,3) = 1
    Read (10) ((ND(I,J),J=1,NDIM),I=1,ngrid)
    If (OF > 0) Then
      Do I = 1, ngrid
        Write (*,'(A,I2,A,3(I5,1X))') 'RFLOCM: Grid ', I, ' is ', (ND(I,J),J=1,NDIM)
      End Do
    End If

    !  Allocate soln memory
    Nx = MAXVAL(ND(:,1))
    Ny = MAXVAL(ND(:,2))
    If (NDIM .EQ. 3) Then
      Nz = MAXVAL(ND(:,3))
    Else
      Nz = 1
    End If
    If (associated(X)) deallocate(X)
    Allocate(X(ngrid,Nx,Ny,Nz,NDIM+2))
    If (prec(1:1) .eq. 's') Allocate(fX(Nx,Ny,Nz,NDIM+2))
    X = 0D0

    !  Read soln values
    If (vf(1:1) .eq. 'w') Then
      If (prec(1:1) .eq. 'd') Then
        Do L = 1, ngrid
          Read (10) tau(1), tau(2), tau(3), tau(4)
          Read (10) ((((X(L,I,J,K,M),I=1,ND(L,1)),J=1,ND(L,2)),K=1,ND(L,3)),M=1,NDIM+2)
        End Do
      Else
        Do L = 1, ngrid
          Read (10) ftau(1), ftau(2), ftau(3), ftau(4)
          Read (10) ((((fX(I,J,K,M),I=1,ND(L,1)),J=1,ND(L,2)),K=1,ND(L,3)),M=1,NDIM+2)
          X(L,:,:,:,:) = Dble(fX(:,:,:,:))
          tau = Dble(ftau)
        End Do
      End If

    Else
      If (prec(1:1) .eq. 'd') Then
        Do L = 1, ngrid
          Read (10) tau(1), tau(2), tau(3), tau(4)
          Do K = 1, ND(L,3)
			print*, K
            Read (10) (((X(L,I,J,K,M),I=1,ND(L,1)),J=1,ND(L,2)),M=1,NDIM+2)
          End Do
        End Do
      Else
        Do L = 1, ngrid
          Read (10) ftau(1), ftau(2), ftau(3), ftau(4)
          Do K = 1, ND(L,3)
            Read (10) (((fX(I,J,K,M),I=1,ND(L,1)),J=1,ND(L,2)),M=1,NDIM+2)
          End Do
          X(L,:,:,:,:) = Dble(fX(:,:,:,:))
          tau = Dble(ftau)
        End Do
      End If
    End If

    Close (10)

    If (prec(1:1) .eq. 's')  Deallocate(fX)

    Return
  End Subroutine Read_Soln

  Subroutine Write_Soln(NDIM, ngrid, ND, X, tau, prec, gf, vf, file)

    Implicit None

    Real(KIND=8), Dimension(:,:,:,:,:), Pointer :: X
    Real(KIND=8) :: tau(4)
    Character(LEN=2) :: prec, gf, vf
    Character(LEN=80) :: file, ib
    Integer :: ngrid, I, J, K, L, M, NDIM, ftype
    Integer, Dimension(:,:), Pointer :: ND
    Integer, Dimension(5) :: OMap(3), NMap(3), ND_copy(3)

    ftype = 1
    Call output_format(file,prec,NDIM,gf,vf,ib,ftype)

    Open (unit=10, file=trim(file), form='unformatted', status='unknown')

    !  Write number of grids
    If (gf(1:1) .eq. 'm') Then
      Write (10) ngrid
    End If

    !  Write grid sizes
    Write (10) ((ND(I,J),J=1,NDIM),I=1,ngrid)

    !  Write soln values
    If (vf(1:1) .eq. 'w') Then
      If (prec(1:1) .eq. 'd') Then
        Do L = 1, ngrid
          Write (10) ((tau(i)),i=1,4)
          Write (10) ((((X(L,I,J,K,M),I=1,ND(L,1)),J=1,ND(L,2)),K=1,ND(L,3)),M=1,NDIM+2)
        End Do
      Else
        Do L = 1, ngrid
          Write (10) (Real(tau(i)),i=1,4)
          Write (10) ((((Real(X(L,I,J,K,M)),I=1,ND(L,1)),J=1,ND(L,2)),K=1,ND(L,3)),M=1,NDIM+2)
        End Do
      End If

    Else
      If (prec(1:1) .eq. 'd') Then
        Do L = 1, ngrid
          Write (10) ((tau(i)),i=1,4)
          Do K = 1, ND(L,3)
            Write (10) (((X(L,I,J,K,M),I=1,ND(L,1)),J=1,ND(L,2)),M=1,NDIM+2)
          End Do
        End Do
      Else
        Do L = 1, ngrid
          Write (10) (Real(tau(i)),i=1,4)
          Do K = 1, ND(L,3)
            Write (10) (((Real(X(L,I,J,K,M)),I=1,ND(L,1)),J=1,ND(L,2)),M=1,NDIM+2)
          End Do
        End Do
      End If
    End If

    Close (10)

    Return
  End Subroutine Write_Soln

  Subroutine Read_Single_Grid(NDIM, ngrid, ND, X, IBLANK, prec, gf, vf, ib, file, gridID, OF)

    Implicit None

    Real(KIND=8), Dimension(:,:,:,:), Pointer :: X ! ... reduction in rank, JKim 04/2008
    Integer, Dimension(:,:,:), Pointer :: IBLANK ! ... reduction in rank, JKim 04/2008
    Character(LEN=2) :: prec, gf, vf, ib
    Character(LEN=80) :: file
    Integer :: ngrid, I, J, K, L, Nx, Ny, Nz, NDIM, M, ftype
    Integer, Dimension(:,:), Pointer :: ND
    Integer :: OMap(3), NMap(3), ND_copy(3)
    Real(KIND=4), Dimension(:,:,:,:), Pointer :: fX
    integer :: gridID ! ... specify which grid we want to read, JKim 04/2008
    Logical :: gf_exists
    Integer :: OF

    Inquire(file=Trim(file),exist=gf_exists)
    If (.NOT.gf_exists) Then
      Write (*,'(A,A,A)') 'RFLOCM: File ', file(1:LEN_TRIM(file)), ' is not readable'
      Stop
    End If

    Call p3d_detect(LEN(Trim(file)),Trim(file),prec,NDIM,gf,vf,ib,ftype)
    If (OF > 0) Call output_format(file,prec,NDIM,gf,vf,ib,ftype)

    Open (unit=10, file=trim(file), form='unformatted', status='old')

    !  Read number of grids
    If (gf(1:1) .eq. 's') Then
      ngrid = 1
    Else
      Read (10) ngrid
    End If
    Allocate(ND(ngrid,3))
    If (OF > 0) Write (*,'(A)') ' '
    If (OF > 0) Write (*,'(A,I2)') 'RFLOCM: Number of grids: ', ngrid

    !  Read grid sizes
    ND(:,3) = 1
    Read (10) ((ND(I,J),J=1,NDIM),I=1,ngrid)
    If (OF > 0) Then
      Do I = 1, ngrid
        Write (*,'(A,I2,A,3(I5,1X))') 'RFLOCM: Grid ', I, ' is ', (ND(I,J),J=1,NDIM)
      End Do
    End If

    !  Allocate grid memory
    Nx = MAXVAL(ND(:,1))
    Ny = MAXVAL(ND(:,2))
    Nz = MAXVAL(ND(:,3))
    Allocate(X(Nx,Ny,Nz,NDIM)) ! ... reduction in rank, JKim 04/2008
    If (prec(1:1) .eq. 's') Allocate(fX(Nx,Ny,Nz,NDIM))
    Allocate(IBLANK(Nx,Ny,Nz)) ! ... reduction in rank, JKim 04/2008

    !  Read grid values
    If (vf(1:1) .eq. 'w') Then
      If (prec(1:1) .eq. 'd') Then
        If (ib(1:1) .eq. 'n') Then
          Do L = 1, gridID ! ... stop after reading the desired grid, JKim 04/2008
            Read (10) ((((X(I,J,K,M),I=1,ND(L,1)),J=1,ND(L,2)),K=1,ND(L,3)),M=1,NDIM)
          End Do
        Else
          Do L = 1, gridID ! ... stop after reading the desired grid, JKim 04/2008
            Read (10) ((((X(I,J,K,M),I=1,ND(L,1)),J=1,ND(L,2)),K=1,ND(L,3)),M=1,NDIM), &
                 (((IBLANK(I,J,K),I=1,ND(L,1)),J=1,ND(L,2)),K=1,ND(L,3))
          End Do
        End If
      Else
        If (ib(1:1) .eq. 'n') Then
          Do L = 1, gridID ! ... stop after reading the desired grid, JKim 04/2008
            Read (10) ((((fX(I,J,K,M),I=1,ND(L,1)),J=1,ND(L,2)),K=1,ND(L,3)),M=1,NDIM)
            X(:,:,:,:) = Dble(fX(:,:,:,:))
          End Do
        Else
          Do L = 1, gridID ! ... stop after reading the desired grid, JKim 04/2008
            Read (10) ((((fX(I,J,K,M),I=1,ND(L,1)),J=1,ND(L,2)),K=1,ND(L,3)),M=1,NDIM), &
                 (((IBLANK(I,J,K),I=1,ND(L,1)),J=1,ND(L,2)),K=1,ND(L,3))
            X(:,:,:,:) = Dble(fX(:,:,:,:))
          End Do
        End If
      End If

    Else
      If (NDIM .EQ. 3) Then
        If (prec(1:1) .eq. 'd') Then
          If (ib(1:1) .eq. 'n') Then
            Do L = 1, gridID ! ... stop after reading the desired grid, JKim 04/2008
              Do K = 1, ND(L,3)
                Read (10) ((X(I,J,K,1),I=1,ND(L,1)),J=1,ND(L,2)), &
                     ((X(I,J,K,2),I=1,ND(L,1)),J=1,ND(L,2)), &
                     ((X(I,J,K,3),I=1,ND(L,1)),J=1,ND(L,2))
              End Do
            End Do
          Else
            Do L = 1, gridID ! ... stop after reading the desired grid, JKim 04/2008
              Do K = 1, ND(L,3)
                Read (10) ((X(I,J,K,1),I=1,ND(L,1)),J=1,ND(L,2)), &
                     ((X(I,J,K,2),I=1,ND(L,1)),J=1,ND(L,2)), &
                     ((X(I,J,K,3),I=1,ND(L,1)),J=1,ND(L,2)), &
                     ((IBLANK(I,J,K),I=1,ND(L,1)),J=1,ND(L,2))
              End Do
            End Do
          End If
        Else
          If (ib(1:1) .eq. 'n') Then
            Do L = 1, gridID ! ... stop after reading the desired grid, JKim 04/2008
              Do K = 1, ND(L,3)
                Read (10) ((fX(I,J,K,1),I=1,ND(L,1)),J=1,ND(L,2)), &
                     ((fX(I,J,K,2),I=1,ND(L,1)),J=1,ND(L,2)), &
                     ((fX(I,J,K,3),I=1,ND(L,1)),J=1,ND(L,2))
              End Do
              X(:,:,:,:) = Dble(fX(:,:,:,:))
            End Do
          Else
            Do L = 1, gridID ! ... stop after reading the desired grid, JKim 04/2008
              Do K = 1, ND(L,3)
                Read (10) ((fX(I,J,K,1),I=1,ND(L,1)),J=1,ND(L,2)), &
                     ((fX(I,J,K,2),I=1,ND(L,1)),J=1,ND(L,2)), &
                     ((fX(I,J,K,3),I=1,ND(L,1)),J=1,ND(L,2)), &
                     ((IBLANK(I,J,K),I=1,ND(L,1)),J=1,ND(L,2))
              End Do
              X(:,:,:,:) = Dble(fX(:,:,:,:))
            End Do
          End If
        End If
      Else
        If (prec(1:1) .eq. 'd') Then
          If (ib(1:1) .eq. 'n') Then
            Do L = 1, gridID ! ... stop after reading the desired grid, JKim 04/2008
              Do J = 1, ND(L,2)
                Read (10) (X(I,J,1,1),I=1,ND(L,1)), &
                     (X(I,J,1,2),I=1,ND(L,1))
              End Do
            End Do
          Else
            Do L = 1, gridID ! ... stop after reading the desired grid, JKim 04/2008
              Do J = 1, ND(L,2)
                Read (10) (X(I,J,1,1),I=1,ND(L,1)), &
                     (X(I,J,1,2),I=1,ND(L,1)), &
                     (IBLANK(I,J,1),I=1,ND(L,1))
              End Do
            End Do
          End If
        Else
          If (ib(1:1) .eq. 'n') Then
            Do L = 1, gridID ! ... stop after reading the desired grid, JKim 04/2008
              Do K = J, ND(L,2)
                Read (10) (fX(I,J,1,1),I=1,ND(L,1)), &
                     (fX(I,J,1,2),I=1,ND(L,1))
              End Do
              X(:,:,:,:) = Dble(fX(:,:,:,:))
            End Do
          Else
            Do L = 1, gridID ! ... stop after reading the desired grid, JKim 04/2008
              Do K = J, ND(L,2)
                Read (10) (fX(I,J,1,1),I=1,ND(L,1)), &
                     (fX(I,J,1,2),I=1,ND(L,1)), &
                     (IBLANK(I,J,1),I=1,ND(L,1))
              End Do
              X(:,:,:,:) = Dble(fX(:,:,:,:))
            End Do
          End If
        End If
      End If
    End If

    Close (10)

    If (prec(1:1) .eq. 's')  Deallocate(fX)
    If (ib(1:1) .eq. 'n') IBLANK = 1

    Return
  End Subroutine Read_Single_Grid

  Subroutine Read_Single_Soln(NDIM, ngrid, ND, X, tau, prec, gf, vf, file, gridID, OF)

    Implicit None

    Real(KIND=8), Dimension(:,:,:,:), Pointer :: X ! ... reduction in rank, JKim 04/2008
    Real(KIND=8) :: tau(4)
    Character(LEN=2) :: prec, gf, vf, ib
    Character(LEN=80) :: file
    Integer :: ngrid, I, J, K, L, M, Nx, Ny, Nz, NDIM, ftype
    Integer, Dimension(:,:), Pointer :: ND
    Real(KIND=4), Dimension(:,:,:,:), Pointer :: fX
    Real(KIND=4) :: ftau(4)
    Logical :: gf_exists
    Integer :: OF, gridID ! ... specify which grid we want to read, JKim 04/2008

    Inquire(file=Trim(file),exist=gf_exists)
    If (.NOT.gf_exists) Then
      Write (*,'(A,A,A)') 'File ', file(1:LEN_TRIM(file)), ' is not readable'
      Stop
    End If

    Call p3d_detect(LEN(Trim(file)),Trim(file),prec,NDIM,gf,vf,ib,ftype)
    If (OF > 0) Call output_format(file,prec,NDIM,gf,vf,ib,ftype)

    Open (unit=10, file=trim(file), form='unformatted', status='old')

    !  Read number of grids
    If (gf(1:1) .eq. 's') Then
      ngrid = 1
    Else
      Read (10) ngrid
    End If

    If (OF > 0) Then
      Write (*, *) ' '
      Write (*,'(A,I2)') 'RFLOCM: Number of grids: ', ngrid
    End If

    !  Read grid sizes
    Allocate(ND(ngrid,3))
    ND(:,3) = 1
    Read (10) ((ND(I,J),J=1,NDIM),I=1,ngrid)
    If (OF > 0) Then
      Do I = 1, ngrid
        Write (*,'(A,I2,A,3(I5,1X))') 'RFLOCM: Grid ', I, ' is ', (ND(I,J),J=1,NDIM)
      End Do
    End If

    !  Allocate soln memory
    Nx = MAXVAL(ND(:,1))
    Ny = MAXVAL(ND(:,2))
    If (NDIM .EQ. 3) Then
      Nz = MAXVAL(ND(:,3))
    Else
      Nz = 1
    End If
    Allocate(X(Nx,Ny,Nz,NDIM+2)) ! ... reduction in rank, JKim 04/2008
    If (prec(1:1) .eq. 's') Allocate(fX(Nx,Ny,Nz,NDIM+2))
    X = 0D0

    !  Read soln values
    If (vf(1:1) .eq. 'w') Then
      If (prec(1:1) .eq. 'd') Then
        Do L = 1, gridID ! ... stop after reading the desired grid, JKim 04/2008
          Read (10) tau(1), tau(2), tau(3), tau(4)
          Read (10) ((((X(I,J,K,M),I=1,ND(L,1)),J=1,ND(L,2)),K=1,ND(L,3)),M=1,NDIM+2)
        End Do
      Else
        Do L = 1, gridID ! ... stop after reading the desired grid, JKim 04/2008
          Read (10) ftau(1), ftau(2), ftau(3), ftau(4)
          Read (10) ((((fX(I,J,K,M),I=1,ND(L,1)),J=1,ND(L,2)),K=1,ND(L,3)),M=1,NDIM+2)
          X(:,:,:,:) = Dble(fX(:,:,:,:))
          tau = Dble(ftau)
        End Do
      End If

    Else
      If (prec(1:1) .eq. 'd') Then
        Do L = 1, gridID ! ... stop after reading the desired grid, JKim 04/2008
          Read (10) tau(1), tau(2), tau(3), tau(4)
          Do K = 1, ND(L,3)
            Read (10) (((X(I,J,K,M),I=1,ND(L,1)),J=1,ND(L,2)),M=1,NDIM+2)
          End Do
        End Do
      Else
        Do L = 1, gridID ! ... stop after reading the desired grid, JKim 04/2008
          Read (10) ftau(1), ftau(2), ftau(3), ftau(4)
          Do K = 1, ND(L,3)
            Read (10) (((fX(I,J,K,M),I=1,ND(L,1)),J=1,ND(L,2)),M=1,NDIM+2)
          End Do
          X(:,:,:,:) = Dble(fX(:,:,:,:))
          tau = Dble(ftau)
        End Do
      End If
    End If

    Close (10)

    If (prec(1:1) .eq. 's')  Deallocate(fX)

    Return
  End Subroutine Read_Single_Soln

  Subroutine output_format(file,prec,NDIM,gf,vf,ib,ftype)

    Character(LEN=80) :: file
    Character(LEN=2) :: prec, gf, vf, ib
    Integer :: NDIM, ftype

    Write (*,'(A)') ""
    Write (*,'(A)') 'File "'//file(1:LEN_TRIM(file))//'" is '

    If (ftype .eq. 0) Then
      Write (*,'(5X,A)') 'a grid file'
    Else
      Write (*,'(5X,A)') 'a solution file'
    End If

    Write (*,'(5X,I1,A)') NDIM, 'D'

    If (prec(1:1) .eq. 'd') Then
      Write (*,'(5X,A)') 'double precision'
    Else
      Write (*,'(5X,A)') 'single precision'
    End If

    If (gf(1:1) .eq. 's') Then
      Write (*,'(5X,A)') 'single grid'
    Else
      Write (*,'(5X,A)') 'multi-block'
    End If

    If (vf(1:1) .eq. 'w') Then
      Write (*,'(5X,A)') 'whole'
    Else
      Write (*,'(5X,A)') 'planes'
    End If

    If (ib(1:1) .eq. 'y') Then
      Write (*,'(5X,A)') 'with IBLANK'
    Else
      Write (*,'(5X,A)') 'without IBLANK'
    End If

    Return

  End Subroutine output_format

  Subroutine Write_Soln_Single_Grid(NDIM, ND, X, tau, fnamei, fnameleni, funit)

    Implicit None

    Real(KIND=8), Dimension(:,:,:,:), Pointer :: X
    Real(KIND=8) :: tau(4)
    Character(LEN=80) :: fnamei
    Integer :: I, J, K, L, M, NDIM, funit, fnameleni
    Integer, Dimension(:), Pointer :: ND

    ! ... open
    Open (unit=funit, file=fnamei(1:fnameleni), form='unformatted', status='unknown', position='append')

    ! ... Write soln values
    Write (funit) ((tau(i)),i=1,4)
    if (NDIM == 3) then
      Write (funit) ((((X(I,J,K,M),I=1,ND(1)),J=1,ND(2)),K=1,ND(3)),M=1,size(X,4))
    else if (NDIM == 2) Then
      Write (funit) ((((X(I,J,K,M),I=1,ND(1)),J=1,ND(2)),K=1,ND(3)),M=1,3), &
                    ((((0.0_8,I=1,ND(1)),J=1,ND(2)),K=1,ND(3)),M=1,1), &
                    ((((X(I,J,K,M),I=1,ND(1)),J=1,ND(2)),K=1,ND(3)),M=4,4)
    else if (NDIM == 1) Then
      Write (funit) ((((X(I,J,K,M),I=1,ND(1)),J=1,ND(2)),K=1,ND(3)),M=1,2), &
                    ((((0.0_8,I=1,ND(1)),J=1,ND(2)),K=1,ND(3)),M=1,2), &
                    ((((X(I,J,K,M),I=1,ND(1)),J=1,ND(2)),K=1,ND(3)),M=3,3)
    end if

    ! ... close
    Close (funit)

    Return

  End Subroutine Write_Soln_Single_Grid

  Subroutine Write_Func_Single_Grid(NDIM, ND, X, fnamei, fnameleni, funit)

    Implicit None

    Real(KIND=8), Dimension(:,:,:,:), Pointer :: X
    Character(LEN=80) :: fnamei
    Integer :: I, J, K, L, M, NDIM, funit, fnameleni
    Integer, Dimension(:), Pointer :: ND

    ! ... open
    Open (unit=funit, file=fnamei(1:fnameleni), form='unformatted', status='unknown', position='append')

    ! ... Write soln values
    Write (funit) ((((X(I,J,K,M),I=1,ND(1)),J=1,ND(2)),K=1,ND(3)),M=1,size(X,4))

    ! ... close
    Close (funit)

    Return

  End Subroutine Write_Func_Single_Grid

  Subroutine Write_Grid_Single_Grid(NDIM, ND, X, IB, fnamei, fnameleni, funit)

    Implicit None

    Real(KIND=8), Dimension(:,:,:,:), Pointer :: X
    Integer, Pointer :: IB(:,:,:)
    Real(KIND=8) :: tau(4)
    Character(LEN=80) :: fnamei
    Integer :: I, J, K, L, M, NDIM, funit, fnameleni
    Integer, Dimension(:), Pointer :: ND

    ! ... open
    Open (unit=funit, file=fnamei(1:fnameleni), form='unformatted', status='unknown', position='append')

    ! ... Write grid values
    if (NDIM == 3) then
      Write (funit) ((((X(I,J,K,M),I=1,ND(1)),J=1,ND(2)),K=1,ND(3)),M=1,size(X,4)), &
                    (((IB(I,J,K),I=1,ND(1)),J=1,ND(2)),K=1,ND(3))
    else if (NDIM == 2) then
      Write (funit) ((((X(I,J,K,M),I=1,ND(1)),J=1,ND(2)),K=1,ND(3)),M=1,2), &
                    ((((0.0_8,I=1,ND(1)),J=1,ND(2)),K=1,ND(3)),M=1,1), &
                    (((IB(I,J,K),I=1,ND(1)),J=1,ND(2)),K=1,ND(3))
    else if (NDIM == 1) then
      Write (funit) ((((X(I,J,K,M),I=1,ND(1)),J=1,ND(2)),K=1,ND(3)),M=1,1), &
                    ((((0.0_8,I=1,ND(1)),J=1,ND(2)),K=1,ND(3)),M=1,2), &
                    (((IB(I,J,K),I=1,ND(1)),J=1,ND(2)),K=1,ND(3))
    end if

    ! ... close
    Close (funit)

    Return

  End Subroutine Write_Grid_Single_Grid


  Subroutine Write_Grid_Soln_Header(NDIM, ngrid, ND, file, fnamelen, funit)

    Implicit None

    Character(LEN=80) :: file
    Integer :: ngrid, I, J, NDIM, funit, fnamelen
    Integer, Dimension(:,:), Pointer :: ND

    Open (unit=funit, file=file(1:fnamelen), form='unformatted', status='unknown')

    !  Write number of grids
    Write (funit) ngrid

    !  Write grid sizes
    Write (funit) ((ND(I,J),J=1,NDIM),I=1,ngrid)

    Close (funit)

  End Subroutine Write_Grid_Soln_Header

  Subroutine Write_Grid_Header_And_Pad(NDIM, ngrid, ND, file, fnamelen, funit)

    Implicit None

    Character(LEN=80) :: file
    Integer :: ngrid, I, J, NDIM, funit, fnamelen
    Integer, Dimension(:,:), Pointer :: ND

    Open (unit=funit, file=file(1:fnamelen), form='unformatted', status='unknown')

    !  Write number of grids
    Write (funit) ngrid

    !  Write grid sizes
    Write (funit) ((ND(I,J),J=1,NDIM),I=1,ngrid)

    ! ... write a bunch of zeros
    Do I = 1, ngrid
      Write (funit) (0.0_8, J=1,PRODUCT(ND(I,:))*3), (0, J=1,PRODUCT(ND(I,:)))
    End Do

    Close (funit)

  End Subroutine Write_Grid_Header_And_Pad


  Subroutine Write_Soln_Header_And_Pad(NDIM, ngrid, ND, file, fnamelen, funit)

    Implicit None

    Character(LEN=80) :: file
    Integer :: ngrid, I, J, NDIM, funit, fnamelen
    Integer, Dimension(:,:), Pointer :: ND

    Open (unit=funit, file=file(1:fnamelen), form='unformatted', status='unknown')

    !  Write number of grids
    Write (funit) ngrid

    !  Write grid sizes
    Write (funit) ((ND(I,J),J=1,NDIM),I=1,ngrid)

    ! ... write a bunch of zeros
    Do I = 1, ngrid
      Write (funit) (0.0_8, J=1,4)
      Write (funit) (0.0_8, J=1,PRODUCT(ND(I,:))*5)
    End Do

    Close (funit)

  End Subroutine Write_Soln_Header_And_Pad

  Subroutine Write_Func_Header(NDIM, ngrid, ND, nvar, file, fnamelen, funit)

    Implicit None

    Character(LEN=80) :: file
    Integer :: ngrid, I, J, NDIM, funit, nvar, fnamelen
    Integer, Dimension(:,:), Pointer :: ND
    Integer :: ND2(ngrid, NDIM+1)

    ND2(:,1:NDIM) = ND(:,:)
    ND2(:,NDIM+1) = nvar

    Open (unit=funit, file=file(1:fnamelen), form='unformatted', status='unknown')

    !  Write number of grids
    Write (funit) ngrid

    !  Write grid sizes
    Write (funit) ((ND2(I,J),J=1,NDIM+1),I=1,ngrid)

    Close (funit)

  End Subroutine Write_Func_Header

  Subroutine Write_Func_Header_And_Pad(NDIM, ngrid, ND, nvar, file, fnamelen, funit)

    Implicit None

    Character(LEN=80) :: file
    Integer :: ngrid, I, J, NDIM, funit, nvar, fnamelen
    Integer, Dimension(:,:), Pointer :: ND
    Integer :: ND2(ngrid, NDIM+1)

    ND2(:,1:NDIM) = ND(:,:)
    ND2(:,NDIM+1) = nvar

    Open (unit=funit, file=file(1:fnamelen), form='unformatted', status='unknown')

    !  Write number of grids
    Write (funit) ngrid

    !  Write grid sizes
    Write (funit) ((ND2(I,J),J=1,NDIM+1),I=1,ngrid)

    ! ... now write a bunch of zeros
    Do I = 1, ngrid
      Write (funit) (0.0_8,J=1,PRODUCT(ND2(I,:)))
    End Do

    Close (funit)

  End Subroutine Write_Func_Header_And_Pad

  Subroutine Read_Single_Grid_IBLANK(NDIM, ngrid, ND, IBLANK, prec, gf, vf, ib, file, gridID, OF)

    Implicit None

    Integer, Dimension(:,:,:), Pointer :: IBLANK ! ... reduction in rank, JKim 04/2008
    Character(LEN=2) :: prec, gf, vf, ib
    Character(LEN=80) :: file
    Integer :: ngrid, I, J, K, L, Nx, Ny, Nz, NDIM, M, ftype
    Integer, Dimension(:,:), Pointer :: ND
    Integer :: OMap(3), NMap(3), ND_copy(3)
    integer :: gridID ! ... specify which grid we want to read, JKim 04/2008
    Logical :: gf_exists
    Integer :: OF
    Real(KIND=8) :: dX
    Real(KIND=4) :: fX

    Inquire(file=Trim(file),exist=gf_exists)
    If (.NOT.gf_exists) Then
      Write (*,'(A,A,A)') 'File ', file(1:LEN_TRIM(file)), ' is not readable'
      Stop
    End If

    Call p3d_detect(LEN(Trim(file)),Trim(file),prec,NDIM,gf,vf,ib,ftype)
    If (OF > 0) Call output_format(file,prec,NDIM,gf,vf,ib,ftype)

    Open (unit=10, file=trim(file), form='unformatted', status='old')

    !  Read number of grids
    If (gf(1:1) .eq. 's') Then
      ngrid = 1
    Else
      Read (10) ngrid
    End If
    Allocate(ND(ngrid,3))
    If (OF > 0) Write (*,'(A)') ' '
    If (OF > 0) Write (*,'(A,I2)') 'Number of grids: ', ngrid

    !  Read grid sizes
    ND(:,3) = 1
    Read (10) ((ND(I,J),J=1,NDIM),I=1,ngrid)
    If (OF > 0) Then
      Do I = 1, ngrid
        Write (*,'(A,I2,A,3(I5,1X))') 'Grid ', I, ' is ', (ND(I,J),J=1,NDIM)
      End Do
    End If

    !  Allocate grid memory
    Nx = MAXVAL(ND(:,1))
    Ny = MAXVAL(ND(:,2))
    Nz = MAXVAL(ND(:,3))
    Allocate(IBLANK(Nx,Ny,Nz)) ! ... reduction in rank, JKim 04/2008

    If (ib(1:1) .eq. 'n') Then

      IBLANK(:,:,:) = 1
      return

    End If


    !  Read grid values
    If (vf(1:1) .eq. 'w') Then
      If (prec(1:1) .eq. 'd') Then
        If (ib(1:1) .eq. 'n') Then
          Do L = 1, gridID ! ... stop after reading the desired grid, JKim 04/2008
            Read (10) ((((dX,I=1,ND(L,1)),J=1,ND(L,2)),K=1,ND(L,3)),M=1,NDIM)
          End Do
        Else
          Do L = 1, gridID ! ... stop after reading the desired grid, JKim 04/2008
            Read (10) ((((dX,I=1,ND(L,1)),J=1,ND(L,2)),K=1,ND(L,3)),M=1,NDIM), &
                 (((IBLANK(I,J,K),I=1,ND(L,1)),J=1,ND(L,2)),K=1,ND(L,3))
          End Do
        End If
      Else
        If (ib(1:1) .eq. 'n') Then
          Do L = 1, gridID ! ... stop after reading the desired grid, JKim 04/2008
            Read (10) ((((fX,I=1,ND(L,1)),J=1,ND(L,2)),K=1,ND(L,3)),M=1,NDIM)
          End Do
        Else
          Do L = 1, gridID ! ... stop after reading the desired grid, JKim 04/2008
            Read (10) ((((fX,I=1,ND(L,1)),J=1,ND(L,2)),K=1,ND(L,3)),M=1,NDIM), &
                 (((IBLANK(I,J,K),I=1,ND(L,1)),J=1,ND(L,2)),K=1,ND(L,3))
          End Do
        End If
      End If

    Else
      If (NDIM .EQ. 3) Then
        If (prec(1:1) .eq. 'd') Then
          If (ib(1:1) .eq. 'n') Then
            Do L = 1, gridID ! ... stop after reading the desired grid, JKim 04/2008
              Do K = 1, ND(L,3)
                Read (10) ((dX,I=1,ND(L,1)),J=1,ND(L,2)), &
                     ((dX,I=1,ND(L,1)),J=1,ND(L,2)), &
                     ((dX,I=1,ND(L,1)),J=1,ND(L,2))
              End Do
            End Do
          Else
            Do L = 1, gridID ! ... stop after reading the desired grid, JKim 04/2008
              Do K = 1, ND(L,3)
                Read (10) ((dX,I=1,ND(L,1)),J=1,ND(L,2)), &
                     ((dX,I=1,ND(L,1)),J=1,ND(L,2)), &
                     ((dX,I=1,ND(L,1)),J=1,ND(L,2)), &
                     ((IBLANK(I,J,K),I=1,ND(L,1)),J=1,ND(L,2))
              End Do
            End Do
          End If
        Else
          If (ib(1:1) .eq. 'n') Then
            Do L = 1, gridID ! ... stop after reading the desired grid, JKim 04/2008
              Do K = 1, ND(L,3)
                Read (10) ((fX,I=1,ND(L,1)),J=1,ND(L,2)), &
                     ((fX,I=1,ND(L,1)),J=1,ND(L,2)), &
                     ((fX,I=1,ND(L,1)),J=1,ND(L,2))
              End Do
            End Do
          Else
            Do L = 1, gridID ! ... stop after reading the desired grid, JKim 04/2008
              Do K = 1, ND(L,3)
                Read (10) ((fX,I=1,ND(L,1)),J=1,ND(L,2)), &
                     ((fX,I=1,ND(L,1)),J=1,ND(L,2)), &
                     ((fX,I=1,ND(L,1)),J=1,ND(L,2)), &
                     ((IBLANK(I,J,K),I=1,ND(L,1)),J=1,ND(L,2))
              End Do
            End Do
          End If
        End If
      Else
        If (prec(1:1) .eq. 'd') Then
          If (ib(1:1) .eq. 'n') Then
            Do L = 1, gridID ! ... stop after reading the desired grid, JKim 04/2008
              Do J = 1, ND(L,2)
                Read (10) (dX,I=1,ND(L,1)), &
                     (dX,I=1,ND(L,1))
              End Do
            End Do
          Else
            Do L = 1, gridID ! ... stop after reading the desired grid, JKim 04/2008
              Do J = 1, ND(L,2)
                Read (10) (dX,I=1,ND(L,1)), &
                     (dX,I=1,ND(L,1)), &
                     (IBLANK(I,J,1),I=1,ND(L,1))
              End Do
            End Do
          End If
        Else
          If (ib(1:1) .eq. 'n') Then
            Do L = 1, gridID ! ... stop after reading the desired grid, JKim 04/2008
              Do K = J, ND(L,2)
                Read (10) (fX,I=1,ND(L,1)), &
                     (fX,I=1,ND(L,1))
              End Do
            End Do
          Else
            Do L = 1, gridID ! ... stop after reading the desired grid, JKim 04/2008
              Do K = J, ND(L,2)
                Read (10) (fX,I=1,ND(L,1)), &
                     (fX,I=1,ND(L,1)), &
                     (IBLANK(I,J,1),I=1,ND(L,1))
              End Do
            End Do
          End If
        End If
      End If
    End If

    Close (10)


    Return

  End Subroutine Read_Single_Grid_IBLANK

  Subroutine Read_Func(NDIM, ngrid, ND, X, file)

    Implicit None

    Character(LEN=80) :: file
    Integer :: ngrid, I, J, K, L, NDIM, NG, nvar, fnamelen
    Integer, Dimension(:,:), Pointer :: ND
    Integer :: ND2(ngrid, NDIM+1)
    REAL(KIND=8), Pointer :: X(:,:,:,:,:)
    Integer, Parameter :: funit = 10

    Open (unit=funit, file=trim(file), form='unformatted', status='unknown')

    !  read number of grids
    Read (funit) ngrid

    !  read grid sizes
    Allocate(ND(ngrid,4))
    Read (funit) ((ND(I,J),J=1,4),I=1,ngrid)

    ! ... read data
    Allocate(X(ngrid,MAXVAL(ND(:,1)),MAXVAL(ND(:,2)),MAXVAL(ND(:,3)),MAXVAL(ND(:,4))))
    Do ng = 1, ngrid
      Read (funit) ((((X(ng,I,J,K,L),I=1,ND(ng,1)),J=1,ND(ng,2)),K=1,ND(ng,3)),L=1,ND(ng,4))
    End Do

    Close (funit)

  End Subroutine Read_Func

  Subroutine Write_Func(ngrid, ND, X, file)

    Implicit None

    Real(KIND=8), Dimension(:,:,:,:,:), Pointer :: X
    Real(KIND=8) :: tau(4)
    Character(LEN=2) :: prec, gf, vf
    Character(LEN=80) :: file, ib
    Integer :: ngrid, I, J, K, L, M, NDIM, ftype, ng
    Integer, Dimension(:,:), Pointer :: ND, ND_aux


    Open (unit=10, file=trim(file), form='unformatted', status='unknown')

    !  Write number of grids
    Write (10) ngrid

    !  Write grid sizes and number of variables
    allocate(ND_aux(ngrid,4))
    ND_aux(:,1:3) = ND(:,:)
    ND_aux(:,4) = size(X,5) 
    Write (10) ((ND_aux(I,J),J=1,4),I=1,ngrid)
    deallocate(ND_aux)

    !  Write func values
    Do ng = 1, ngrid
      Write (10) ((((X(ng,I,J,K,M),I=1,ND(ng,1)),J=1,ND(ng,2)),K=1,ND(ng,3)),M=1,size(X,5))
    End Do

    Close (10)
    Return

  End Subroutine Write_Func

End Module ModPLOT3D_IO
