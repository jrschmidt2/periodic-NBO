MODULE rd_basis
  USE rd_wavefunction
  USE projection_shared
  IMPLICIT NONE
  PRIVATE
  PUBLIC read_basis,basis_set_info_out

  CONTAINS

! This subroutine reads in a basis input file
! The process is split into two parts
! First the basis set is read to see how many orbitals there are and make sure the same number of atom types is contained as the wavefunction input file
! Second the file is read again and specific quantities are assigned to basis functions
!  *** The information for each basis function is stored in the AO_function structures, there is one for each basis function
!  *** One for each basis function on each atomic center will be filled out in this subroutine


SUBROUTINE read_basis(basis_fn,AO_basis)
  IMPLICIT NONE

  CHARACTER(128),INTENT(IN)  ::  basis_fn

  TYPE(AO_function), ALLOCATABLE :: AO_basis(:)

  INTEGER  ::  nu
  INTEGER  ::  atom_shift,orb_sum,iatom
  INTEGER  ::  igauss,ityp,icart
  INTEGER  ::  level_shift

  INTEGER          ::  IO
  CHARACTER(2)     ::  dummy
  CHARACTER(2)     ::  orb_type
  CHARACTER(2)     ::  element_sym
  INTEGER          ::  orb_gauss
  INTEGER          ::  basis_per_atom


  ALLOCATE(atom_basis(n_atom+1),atom_symbol(n_atom),atomic_number(n_atom),atom_valence(n_atom))

  OPEN(20, file=basis_fn)
  CALL skip_header(20,1)

  ityp = 0
  s_dim = 0

  !First the basis file is pre-read to find the number of: 1)Types of atoms & 2)Basis functions per atom
  !The loops to read the file are open, since nothing is known about them before hand
  DO
     READ(20,*,IOSTAT=IO)element_sym
     IF( IO < 0 )EXIT
     !WRITE(6,*)element_sym

     ityp = ityp + 1
     basis_per_atom = 0

     DO
        READ(20,'(A2)',ADVANCE='NO')orb_type
        IF( orb_type == '**' )THEN
           READ(20,'(A2)')dummy
          EXIT
        ENDIF
        READ(20,*) orb_gauss
        !WRITE(6,*)orb_type,orb_gauss

        DO igauss=1,orb_gauss
           READ(20,*)
        ENDDO

        IF( orb_type=='S' .OR. orb_type=='s' )THEN
            basis_per_atom = basis_per_atom + 1
        ELSEIF( orb_type=='SP' .OR. orb_type=='sp' )THEN
            basis_per_atom = basis_per_atom + 4
        ELSEIF( orb_type=='P' .OR. orb_type=='p' )THEN
            basis_per_atom = basis_per_atom + 3
        ELSEIF( orb_type=='D' .OR. orb_type=='d' )THEN
            basis_per_atom = basis_per_atom + 5
        ELSEIF( orb_type=='F' .OR. orb_type=='f' )THEN
            basis_per_atom = basis_per_atom + 7
        ELSE
            WRITE(6,'(A)')"Unreadable basis function type in basis input file, pre-screen ", orb_type
            STOP
        ENDIF
        !WRITE(6,*)basis_per_atom
     ENDDO

     !WRITE(6,*)

     s_dim = s_dim + basis_per_atom*atom_types(ityp)
     !WRITE(6,*)'Current s_dim',s_dim

  ENDDO

  CLOSE(20)

  !WRITE(6,*)ityp,ntyp
  !The number of types found in the basis set input is compared to that from the wavefunction.dat file.  If they are different there is an error.
  IF( ityp /= ntyp )THEN
     WRITE(6,*)"Number of atom types in basis set does not match wavefunction input"
     WRITE(6,*)ityp,ntyp
     STOP
  ENDIF

  WRITE(6,*) '*** Reading Atomic Basis Set ***'
  WRITE(6,*)
  WRITE(6,*)'Total number of basis functions per unit cell',s_dim

  ALLOCATE(AO_basis(s_dim))

  !Now the input file must be read again for actual values of exponents and coefficients of basis functions
  !The assignment is non trivial as their can be multiple iterations of the same atom in the system, even though there is only one listing in the basis set file
  !To deal with that the file will be read multiple times (based on how many centers there are of a particular atom).
  OPEN(20, file=basis_fn)

  nu=1 !We can not loop over the basis functions directly, but this variable is used to keep track of what basis function is currently being interpreted

  level_shift=0  !Keeps track of how many unique sets of exponents are looked at.

  atom_shift=0 !Keeps track of what atom the basis function is being assigned to.  Independent of type

  DO ityp=1,ntyp !This index runs over the total number of different types of atoms; number of unique atomic basis sets.
     DO iatom=1,atom_types(ityp) !This runs over the number of unique centers of each particular atom type; input from rd_wavefunction
        atom_shift = atom_shift + 1

        !Skip to the appropriate atomic section of the basis input file
        CALL skip_header(20,ityp)

        !Keep track of information about each atom.
        !While not used in the projection code this information is necessary for NBO analysis
        READ(20,*)atom_symbol(atom_shift)
        CALL symbol2number(atom_symbol(atom_shift),atomic_number(atom_shift))
        atom_valence(atom_shift) = P(ityp)%val_corr
        atom_basis(atom_shift) = nu

        DO
           !It is not known how many orbitals there will be for a particular type, so we again have to read them in an open loop
           READ(20,'(A2)',ADVANCE='NO')orb_type
           IF( orb_type=='**' )THEN
              READ(20,*)
              EXIT
           ENDIF
           READ(20,*)orb_gauss

           !Actual basis functions will now be constructed from the read in orbital type and number of gaussians
           level_shift = level_shift + 1
           !WRITE(6,*)'nu on calling read_orbital',nu
           CALL read_orbital(20,orb_type,orb_gauss,AO_basis,nu,atom_shift,level_shift)
        ENDDO

        REWIND(20)

     ENDDO

  ENDDO

  atom_basis(atom_shift+1) = nu  !Index of 1 past end of basis set.  Used in NBO code

  IF( nu-1 .NE. s_dim )THEN
     WRITE(6,*)'An incorrect number of basis functions have been interpreted',nu-1,s_dim
     STOP
  ENDIF

  IF( n_atom .NE. atom_shift )THEN
     WRITE(6,*)'An incorrect number of atoms have been interpreted in rd_basis',n_atom,atom_shift
     STOP
  ENDIF

  WRITE(6,'(A,30A4)')' Atoms of system  ',atom_symbol


  !DO iatom=1,n_atom
  !   WRITE(6,*)'Atom number',iatom
  !   WRITE(6,*)'Symbol ',atom_symbol(iatom)
  !   WRITE(6,*)'Number ',atomic_number(iatom)
  !   WRITE(6,*)'Valence',SNGL(atom_valence(iatom))
  !   WRITE(6,*)'# Basis',atom_basis(iatom)
  !   WRITE(6,*)
  !ENDDO

  DO nu=1,s_dim
     !WRITE(6,*)'For basis function',nu
     !WRITE(6,*)'Atomic center',AO_basis(nu)%atom
     !WRITE(6,'(A,3F12.8)')' Atom Position ',AO_basis(nu)%pos
     !WRITE(6,*)'Exponent Level',AO_basis(nu)%level
     !WRITE(6,*)'number gauss',AO_basis(nu)%num_gauss
     !WRITE(6,*)'l- and m-  ',AO_basis(nu)%l,AO_basis(nu)%m
     !WRITE(6,*)' Exponent          Coefficient       Norm'
     !DO igauss=1,AO_basis(nu)%num_gauss
     !   WRITE(6,'(3E18.10)')AO_basis(nu)%alpha(igauss),AO_basis(nu)%coeff(igauss),AO_basis(nu)%norm(igauss)
     !ENDDO

     CALL interp_cart(AO_basis(nu))

     !WRITE(6,*)'Num. Cartesian',AO_basis(nu)%ncart
     !WRITE(6,*)'Cart. Coeff and Exponents'
     !DO icart=1,AO_basis(nu)%ncart
     !   WRITE(6,*)AO_basis(nu)%cart_coeff(icart),AO_basis(nu)%cart_mat(icart,:)
     !ENDDO
     !WRITE(6,*)

  ENDDO

  WRITE(6,*)


END SUBROUTINE read_basis



!The basis set input file is actually read and values are assigned to particular basis functions
!First the type is determined, then the appropriate number of basis function are allocated (1 for s, 4 for sp, 3 for p, 5 for d)
!Then the exponent and coefficient values are actually read in from the file.
!The next step is to assign those values to all the equivalent basis functions.  
SUBROUTINE read_orbital(unit,orb_type,orb_gauss,AO_basis,nu,iatom,ilevel)
  IMPLICIT NONE

  TYPE(AO_function), DIMENSION(:)  ::  AO_basis

  INTEGER,INTENT(IN)       ::  unit
  CHARACTER(2),INTENT(IN)  ::  orb_type
  INTEGER,INTENT(IN)       ::  orb_gauss
  INTEGER,INTENT(INOUT)    ::  nu
  INTEGER,INTENT(IN)       ::  iatom
  INTEGER,INTENT(IN)       ::  ilevel

  INTEGER  ::  nu_shift
  INTEGER  ::  igauss,j

  AO_basis(nu)%atom  = iatom
  AO_basis(nu)%pos   = atoms(:,iatom)
  AO_basis(nu)%level = ilevel

  !Based on the type the appropriate number of basis functions are allocated and the l-quantum number is set
  !ONLY FOR SP-ORBITALS is more required
  IF( orb_type == 'S' .OR. orb_type == 's' )THEN
     AO_basis(nu)%l = 0

  ELSEIF( orb_type == 'P' .OR. orb_type == 'p' )THEN
     AO_basis(nu)%l = 1

  !SP are more complicated to interpret, so everything is done for the s-type here.
  !Then the first p-orbital is set as the seed for use below
  ELSEIF( orb_type == 'SP' .OR. orb_type == 'sp' )THEN

     CALL array_allocation(AO_basis,nu,3,orb_gauss)
     !The SP must be specially read in, since there are two coefficients for each exponent
     DO igauss=1,AO_basis(nu)%num_gauss
        READ(unit,*)AO_basis(nu)%alpha(igauss),AO_basis(nu)%coeff(igauss),AO_basis(nu+1)%coeff(igauss)
     ENDDO
     AO_basis(nu)%l = 0
     AO_basis(nu)%m = 1
     CALL norm_calc(AO_basis(nu)%alpha,AO_basis(nu)%norm,AO_basis(nu)%l)

     !Set up the p-orbital lmn and norm vectors
     nu = nu + 1
     AO_basis(nu)%l = 1
     !Since the first p-orbital is the seed to filling, some other values need to be shifted from the s- to the first p-orbital
     AO_basis(nu)%alpha = AO_basis(nu-1)%alpha
     AO_basis(nu)%atom  = AO_basis(nu-1)%atom
     AO_basis(nu)%pos   = AO_basis(nu-1)%pos
     AO_basis(nu)%level = AO_basis(nu-1)%level

     GOTO 400  !This skips the general read in used below for single types

  ELSEIF( orb_type == 'D' .OR. orb_type == 'd' )THEN
     AO_basis(nu)%l = 2

  ELSEIF( orb_Type == 'F' .OR. orb_type == 'f' )THEN
     AO_basis(nu)%l = 3

  ELSE
     STOP 'unreadable orbital type in read_orbital subroutine' 
  ENDIF

  !For single type orbitals the coeff, norm, and alpha arrays for all equivalent orbitals are allocated
  CALL array_allocation(AO_basis,nu,2*AO_basis(nu)%l,orb_gauss)

  !All single orbital types can be read in using the same format
  DO igauss=1,AO_basis(nu)%num_gauss
     READ(unit,*)AO_basis(nu)%alpha(igauss),AO_basis(nu)%coeff(igauss)
  ENDDO

400 CONTINUE

  !All norms are calculated in the same way.
  !Note higher order cartesian effects on the norms are not include here.  
  !Those will instead be included in cart_coeff
  CALL norm_calc(AO_basis(nu)%alpha,AO_basis(nu)%norm,AO_basis(nu)%l)

  !The m-index just runs from 1 to 2*l+1 and therfore all first functions have 1 as their m-quantum number
  AO_basis(nu)%m = 1

  !For s-type orbitals, there is only one and no other equivalent orbitals to fill with alpha's or coeff's etc.
  !So nu is shifted and the subroutine exited
  IF( AO_basis(nu)%l == 0 )THEN
     nu = nu + 1
     GOTO 500
  ENDIF

  !Now all equivalent orbitals are filled with the same information already assigned to the nu-th basis function
  !This way all p- or d-type orbitals will contain the same information
  DO j=1,2*AO_basis(nu)%l

     AO_basis(nu+j)%alpha = AO_basis(nu)%alpha
     AO_basis(nu+j)%coeff = AO_basis(nu)%coeff
     AO_basis(nu+j)%norm  = AO_basis(nu)%norm
     AO_basis(nu+j)%atom  = AO_basis(nu)%atom
     AO_basis(nu+j)%pos   = AO_basis(nu)%pos
     AO_basis(nu+j)%level = AO_basis(nu)%level
     AO_basis(nu+j)%l     = AO_basis(nu)%l

     AO_basis(nu+j)%m = j + 1  !The m-index is also assigned

  ENDDO

  !nu_shift is necessary since it is probably not a good idea to change the value of nu using a variable indexed by nu
  nu_shift = 2*AO_basis(nu)%l+1
  nu = nu + nu_shift

500 CONTINUE


END SUBROUTINE  read_orbital


!For a GTO with alphas=expon and cartesian powers=lmn, the norm is--> norm_vec
!The following subroutine comes from: Clementi and Davis, Journal of Comp. Physics; 2,223-244(1967)
SUBROUTINE norm_calc(expon, norm_vec, l) 
  IMPLICIT NONE

  REAL*8, DIMENSION(:), INTENT(IN)  :: expon
  REAL*8, DIMENSION(:), INTENT(OUT) :: norm_vec
  INTEGER, INTENT(IN)               :: l

  REAL*8, PARAMETER  ::  pi=3.14159265358979323846264338d0

  !Note the actually double factorials we are looking to calculate are 2*l-1. There is a shift up by 2 due to the array indexing. dble_factor(1) = -1!!
  !norm_vec = 2.d0**(2*l+1.5) / ( pi**1.5 )

  !norm_vec = ( norm_veci * expon**(SUM(lmn)+1.5) )**0.5


  norm_vec = ( expon**(2*l+3) * 2.d0**(4*l+3) / pi**3 )**(0.25d0)


ENDSUBROUTINE norm_calc



!For a particular basis function type, all the vectors(norm,coeff,alpha) must be allocated
!Since multiple basis sets are read at once (i.e. px,py,pz), the arrays for multiple basis functions are allocated in a loop
!The 'shift' variable is input to determine how many additional basis functions are treated (2 for p-type, 4 for d-type)
SUBROUTINE array_allocation(AO_basis,nu,shift,orb_gauss)
  IMPLICIT NONE

  TYPE(AO_function), DIMENSION(:) :: AO_basis

  INTEGER  :: nu
  INTEGER  :: shift
  INTEGER  :: orb_gauss
  INTEGER  :: l

  INTEGER  :: j

  DO j=0,shift
     AO_basis(nu+j)%num_gauss=orb_gauss
     !WRITE(6,*)'Gaussian number for orbital',nu+j,AO_basis(nu+j)%num_gauss
     ALLOCATE(AO_basis(nu+j)%norm(orb_gauss))
     ALLOCATE(AO_basis(nu+j)%alpha(orb_gauss))
     ALLOCATE(AO_basis(nu+j)%coeff(orb_gauss))
  ENDDO

END SUBROUTINE array_allocation



!Based on the l- and m- indices of the basis function 'AO' the approproate cartesian component is assigned 
!The number of cartesian components as well as each's coefficient and exponents are all assigned
!There is not much fancy about this, just checking the indices and assigning value from there
!Wherever orbitals could be grouped together I have (same ncart's, same coeff's, etc.)
!
!The linear combinations used come from:
!Schlegel and Frisch, Intl. J. of Quant. Chem.; 54, 83-87 (1995)
SUBROUTINE interp_cart(AO)
  IMPLICIT NONE

  TYPE(AO_Function)    ::  AO

  IF( AO%l .LT. 2 )THEN   !s, px, py, pz

     AO%ncart = 1
     ALLOCATE(AO%cart_mat(AO%ncart,3),AO%cart_coeff(AO%ncart))
     AO%cart_coeff = 1.d0
     AO%cart_mat = 0
     AO%cart_mat(1,AO%m) = AO%l

  ELSEIF( AO%l .EQ. 2 )THEN   !d-orbitals

     IF( AO%m .LE. 3 )THEN    !yz, xz, xy
        AO%ncart = 1
        ALLOCATE(AO%cart_mat(AO%ncart,3),AO%cart_coeff(AO%ncart))
        AO%cart_coeff = 1.d0
        AO%cart_mat = 1
        AO%cart_mat(1,AO%m) = 0
     ELSEIF( AO%m .LE. 5 )THEN

        IF( AO%m .EQ. 4 )THEN    !x^2-y^2
           AO%ncart = 2
           ALLOCATE(AO%cart_mat(AO%ncart,3),AO%cart_coeff(AO%ncart))
           AO%cart_coeff(1) = 0.5d0
           AO%cart_coeff(2) = -0.5d0
        ELSE       ! z^2
           AO%ncart = 3
           ALLOCATE(AO%cart_mat(AO%ncart,3),AO%cart_coeff(AO%ncart))
           AO%cart_coeff = -1.d0 / (2.d0 * SQRT(3.d0))
           AO%cart_coeff(3) = 1.d0 / SQRT(3.d0)
        ENDIF

        AO%cart_mat = 0
        CALL fill_cart_mat(AO%cart_mat)

     ELSE
        STOP 'Improperly indexed d-orbital'
     ENDIF

  ELSEIF( AO%l .EQ. 3 )THEN   !f-orbitals

     IF( AO%m .EQ. 1 )THEN  !xyz orbital

        AO%ncart = 1
        ALLOCATE(AO%cart_mat(AO%ncart,3),AO%cart_coeff(AO%ncart))
        AO%cart_coeff = 1.d0
        AO%cart_mat = 1

     ELSEIF( AO%m .LE. 4 )THEN

        AO%ncart = 2
        ALLOCATE(AO%cart_mat(AO%ncart,3),AO%cart_coeff(AO%ncart))

        IF( AO%m .EQ. 2 )THEN  !z(x^2-y^2)
           AO%cart_coeff(1) = 0.5d0
           AO%cart_coeff(2) = -0.5d0
        ELSEIF( AO%m .EQ. 3 )THEN  !y(3x^2-y^2)
           AO%cart_coeff(1) = 0.5d0 * SQRT(1.5d0)
           AO%cart_coeff(2) = -0.5d0 / SQRT(6.d0)
        ELSE                     !x(x^2-3y^2)
           AO%cart_coeff(1) = 0.5d0 / SQRT(6.d0)
           AO%cart_coeff(2) = -0.5d0 * SQRT(1.5d0) 
        ENDIF

        AO%cart_mat = 0
        AO%cart_mat(:,5-AO%m) = 1
        CALL fill_cart_mat(AO%cart_mat)

     ELSEIF( AO%m .LE. 7 )THEN

        AO%ncart = 3
        ALLOCATE(AO%cart_mat(AO%ncart,3),AO%cart_coeff(AO%ncart))

        IF( AO%m .LE. 6 )THEN  !xz^2 and yz^2
           AO%cart_coeff(3) = 2.d0 / SQRT(10.d0)
           AO%cart_coeff(1:2) = -0.5d0 / SQRT(10.d0)
        ELSE       !z^3
           AO%cart_coeff(3) = 1.d0 / SQRT(15.d0)
           AO%cart_coeff(1:2) = -0.5d0 * SQRT(0.6d0) 
        ENDIF

        AO%cart_mat = 0
        AO%cart_mat(:,AO%m-4) = 1
        CALL fill_cart_mat(AO%cart_mat)

     ELSE
        WRITE(6,*)'Improper m-index given to an f-orbital',AO%m
        STOP
     ENDIF

  ELSE
     STOP 'can only interpret up to f-type orbitals in interp_cart'
  ENDIF


CONTAINS

  !For many orbitals the cartesian exponents have a squared factor across all components
  !This is an overly simple subroutine, but it greatly cleans up the assignments above.
  SUBROUTINE  fill_cart_mat(cart_mat) 
    IMPLICIT NONE
    INTEGER,INTENT(INOUT)   ::  cart_mat(:,:)
    INTEGER   ::  j
    DO j=1,SIZE(cart_mat,1)
       cart_mat(j,j) = cart_mat(j,j) + 2
    ENDDO 
  END SUBROUTINE



END SUBROUTINE




!This subroutine shifts the reading point in the input file to the desired atom type.
!The 'loops' variable determines which atom type to shift to. For loops=1 it is the first, loops=2 is the second etc.
SUBROUTINE skip_header(unit,loops)
  IMPLICIT NONE
  CHARACTER(2)  ::  dummy
  INTEGER,INTENT(IN)  ::  loops
  INTEGER       ::  unit,counter
  counter = loops
  DO
     READ(unit,*)dummy
     IF( dummy == '**' )THEN
        counter = counter - 1
     ENDIF
     IF( counter == 0 )EXIT
  ENDDO
END SUBROUTINE skip_header



!This subrotuine generates a file with information on each basis function
!This file can then easily be read in for visualization purposes.
SUBROUTINE basis_set_info_out(AO_basis)
  IMPLICIT NONE

  TYPE(AO_function), DIMENSION(:)  ::  AO_basis

  INTEGER :: j,nu

  OPEN(35, file='basis_set.info')

  WRITE(35,*)'info on basis set that was projected into, from BDD projection code'
  WRITE(35,*)
  WRITE(35,*)n_atom, '!number of atoms'
  WRITE(35,*)s_dim, '!number of basis functions'
  WRITE(35,*)
  DO j=1,3
     WRITE(35,*)a(:,j)
  ENDDO
  WRITE(35,*)
  DO j=1,n_atom
     WRITE(35,*)atomic_number(j)
     WRITE(35,*)atoms(:,j)
  ENDDO
  WRITE(35,*)
  DO nu=1,s_dim
     WRITE(35,'(I)')AO_basis(nu)%num_gauss
     WRITE(35,'(10D17.9)')AO_basis(nu)%alpha
     WRITE(35,'(10D17.9)')AO_basis(nu)%coeff
     WRITE(35,'(10D17.9)')AO_basis(nu)%norm
     WRITE(35,'(3D16.8)')AO_basis(nu)%pos
     !WRITE(35,'(3I)')AO_basis(nu)%lmn
     WRITE(35,'(I)')AO_basis(nu)%ncart
     DO j=1,AO_basis(nu)%ncart
        WRITE(35,'(D17.9,3I)')AO_basis(nu)%cart_coeff(j),AO_basis(nu)%cart_mat(j,:)
     ENDDO

     WRITE(35,*)
  ENDDO

  CLOSE(35)





END SUBROUTINE basis_set_info_out



!This subroutine takes an atomic symbol as input and outputs the appropriate atomic number
!It is nothing fancy just simple brute force scanning of all symbols.
SUBROUTINE symbol2number(symbol,atomic_number)
  IMPLICIT NONE
  CHARACTER(2), INTENT(IN)  ::  symbol
  INTEGER, INTENT(OUT)      ::  atomic_number

  IF( symbol == 'H' )THEN
    atomic_number = 1
  ELSEIF( symbol == 'He' )THEN
    atomic_number = 2
  ELSEIF( symbol == 'Li' )THEN
    atomic_number = 3
  ELSEIF( symbol == 'Be' )THEN
    atomic_number = 4
  ELSEIF( symbol == 'B' )THEN
    atomic_number = 5
  ELSEIF( symbol == 'C' )THEN
    atomic_number = 6
  ELSEIF( symbol == 'N' )THEN
    atomic_number = 7
  ELSEIF( symbol == 'O' )THEN
    atomic_number = 8
  ELSEIF( symbol == 'F' )THEN
    atomic_number = 9
  ELSEIF( symbol == 'Ne' )THEN
    atomic_number = 10
  ELSEIF( symbol == 'Na' )THEN
    atomic_number = 11
  ELSEIF( symbol == 'Mg' )THEN
    atomic_number = 12
  ELSEIF( symbol == 'Al' )THEN
    atomic_number = 13
  ELSEIF( symbol == 'Si' )THEN
    atomic_number = 14
  ELSEIF( symbol == 'P' )THEN
    atomic_number = 15
  ELSEIF( symbol == 'S' )THEN
    atomic_number = 16
  ELSEIF( symbol == 'Cl' )THEN
    atomic_number = 17
  ELSEIF( symbol == 'Ar' )THEN
    atomic_number = 18
  ELSEIF( symbol == 'K' )THEN
    atomic_number = 19
  ELSEIF( symbol == 'Ca' )THEN
    atomic_number = 20
  ELSEIF( symbol == 'Sc' )THEN
    atomic_number = 21
  ELSEIF( symbol == 'Ti' )THEN
    atomic_number = 22
  ELSEIF( symbol == 'V' )THEN
    atomic_number = 23
  ELSEIF( symbol == 'Cr' )THEN
    atomic_number = 24
  ELSEIF( symbol == 'Mn' )THEN
    atomic_number = 25
  ELSEIF( symbol == 'Fe' )THEN
    atomic_number = 26
  ELSEIF( symbol == 'Co' )THEN
    atomic_number = 27
  ELSEIF( symbol == 'Ni' )THEN
    atomic_number = 28
  ELSEIF( symbol == 'Cu' )THEN
    atomic_number = 29
  ELSEIF( symbol == 'Zn' )THEN
    atomic_number = 30
  ELSEIF( symbol == 'Ga' )THEN
    atomic_number = 31
  ELSEIF( symbol == 'Ge' )THEN
    atomic_number = 32
  ELSEIF( symbol == 'As' )THEN
    atomic_number = 33
  ELSEIF( symbol == 'Se' )THEN
    atomic_number = 34
  ELSEIF( symbol == 'Br' )THEN
    atomic_number = 35
  ELSEIF( symbol == 'Kr' )THEN
    atomic_number = 36
  ELSEIF( symbol == 'Rb' )THEN
    atomic_number = 37
  ELSEIF( symbol == 'Sr' )THEN
    atomic_number = 38
  ELSEIF( symbol == 'Y' )THEN
    atomic_number = 39
  ELSEIF( symbol == 'Zr' )THEN
    atomic_number = 40
  ELSEIF( symbol == 'Nb' )THEN
    atomic_number = 41
  ELSEIF( symbol == 'Mo' )THEN
    atomic_number = 42
  ELSEIF( symbol == 'Tc' )THEN
    atomic_number = 43
  ELSEIF( symbol == 'Ru' )THEN
    atomic_number = 44
  ELSEIF( symbol == 'Rh' )THEN
    atomic_number = 45
  ELSEIF( symbol == 'Pd' )THEN
    atomic_number = 46
  ELSEIF( symbol == 'Ag' )THEN
    atomic_number = 47
  ELSEIF( symbol == 'Cd' )THEN
    atomic_number = 48
  ELSEIF( symbol == 'In' )THEN
    atomic_number = 49
  ELSEIF( symbol == 'Sn' )THEN
    atomic_number = 50
  ELSEIF( symbol == 'Sb' )THEN
    atomic_number = 51
  ELSEIF( symbol == 'Te' )THEN
    atomic_number = 52
  ELSEIF( symbol == 'I' )THEN
    atomic_number = 53
  ELSEIF( symbol == 'Xe' )THEN
    atomic_number = 54
  ELSEIF( symbol == 'Cs' )THEN
    atomic_number = 55
  ELSEIF( symbol == 'Ba' )THEN
    atomic_number = 56
  ELSEIF( symbol == 'La' )THEN
    atomic_number = 57
  ELSEIF( symbol == 'Ce' )THEN
    atomic_number = 58
  ELSEIF( symbol == 'Pr' )THEN
    atomic_number = 59
  ELSEIF( symbol == 'Nd' )THEN
    atomic_number = 60
  ELSEIF( symbol == 'Pm' )THEN
    atomic_number = 61
  ELSEIF( symbol == 'Sm' )THEN
    atomic_number = 62
  ELSEIF( symbol == 'Eu' )THEN
    atomic_number = 63
  ELSEIF( symbol == 'Gd' )THEN
    atomic_number = 64
  ELSEIF( symbol == 'Tb' )THEN
    atomic_number = 65
  ELSEIF( symbol == 'Dy' )THEN
    atomic_number = 66
  ELSEIF( symbol == 'Ho' )THEN
    atomic_number = 67
  ELSEIF( symbol == 'Er' )THEN
    atomic_number = 68
  ELSEIF( symbol == 'Tm' )THEN
    atomic_number = 69
  ELSEIF( symbol == 'Yb' )THEN
    atomic_number = 70
  ELSEIF( symbol == 'Lu' )THEN
    atomic_number = 71
  ELSEIF( symbol == 'Hf' )THEN
    atomic_number = 72
  ELSEIF( symbol == 'Ta' )THEN
    atomic_number = 73
  ELSEIF( symbol == 'W' )THEN
    atomic_number = 74
  ELSEIF( symbol == 'Re' )THEN
    atomic_number = 75
  ELSEIF( symbol == 'Os' )THEN
    atomic_number = 76
  ELSEIF( symbol == 'Ir' )THEN
    atomic_number = 77
  ELSEIF( symbol == 'Pt' )THEN
    atomic_number = 78
  ELSEIF( symbol == 'Au' )THEN
    atomic_number = 79
  ELSEIF( symbol == 'Hg' )THEN
    atomic_number = 80
  ELSEIF( symbol == 'Tl' )THEN
    atomic_number = 81
  ELSEIF( symbol == 'Pb' )THEN
    atomic_number = 82
  ELSEIF( symbol == 'Bi' )THEN
    atomic_number = 83
  ELSEIF( symbol == 'Po' )THEN
    atomic_number = 84
  ELSEIF( symbol == 'At' )THEN
    atomic_number = 85
  ELSEIF( symbol == 'Rn' )THEN
    atomic_number = 86
  ELSE
    WRITE(6,*)'non recognized atomic symbol passed to symbol2number subroutine'
    STOP
  ENDIF


ENDSUBROUTINE symbol2number



END MODULE
