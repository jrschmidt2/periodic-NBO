MODULE projection_shared
  IMPLICIT NONE

  !For use in describing basis set
  INTEGER  :: s_dim  !Total number of atomic orbital basis functions per unit cell
  TYPE AO_function
     REAL*8, ALLOCATABLE, DIMENSION(:) :: norm, alpha, coeff
     INTEGER                :: num_gauss
     !INTEGER, DIMENSION(3)  :: lmn
     INTEGER                :: atom
     REAL*8, DIMENSION(3)   :: pos
     INTEGER                :: level  !Keeps track of what basis functions share the same sets of exponents, on the same atom.
     INTEGER                :: l,m     !l- and m-quantum numbers of basis function.  m actually runs from 0,2*l+1 and is onyl and index.

     INTEGER                :: ncart
     INTEGER,ALLOCATABLE    :: cart_mat(:,:)
     REAL*8,ALLOCATABLE     :: cart_coeff(:)

  END TYPE AO_function

  !Information about each atom in the system.  Filled in rd_basis MODULE
  INTEGER,ALLOCATABLE      ::  atom_basis(:)     !Index of the first basis function belonging to that atom. Includes additional term for one after last basis function.
  CHARACTER(2),ALLOCATABLE ::  atom_symbol(:)    !Atomic symbol of atom
  INTEGER,ALLOCATABLE      ::  atomic_number(:)  !Atomic number of atom
  REAL*8,ALLOCATABLE       ::  atom_valence(:)   !Effective valence of atom.  How many electrons are explicitly treated for this atom in PW calculation.

  LOGICAL   ::  surface

  !Used in storing various matrices used in bloch_overlap MODULE
  REAL*8, ALLOCATABLE        ::  per_s_mat(:,:,:)   !Contains the overlap matrix between all cell pairs in S{u,v}^{0,l}
  COMPLEX*16, ALLOCATABLE    ::  bloch_s_mat(:,:,:) !Fourier transform of the per_s_mat, gives the overlap between all bloch orbitals at different k-points
  COMPLEX*16, ALLOCATABLE    ::  bloch_s_inv(:,:,:) !Inverse of bloch_s_mat, unique inverse for each k-point
  REAL*8, ALLOCATABLE        ::  s_mat_tilde(:,:)   !Used specially for Gamma Point Calculations.

  COMPLEX*16,ALLOCATABLE   ::  bloch_density(:,:,:,:) !Density matrix in the projected AO basis
  COMPLEX*16,ALLOCATABLE   ::  bloch_fock(:,:,:,:)    !Fock matrix in the projected AO basis

  REAL*8, PARAMETER        ::  nu_g_screen=45.d0
  REAL*8, PARAMETER        ::  pi=3.141592653589793238462d0
  COMPLEX*16, PARAMETER    ::  sqrt_neg_one=(0.d0,1.d0)
  REAL*8, PARAMETER        ::  proj_2pi = (2.d0 * pi ) ** 0.75


CONTAINS



!This subroutine calculates the overlap for an arbitrary GTO basis functions (centered at the origin) and a planewave basis function
!This is equivalent to the Fourier transform of the basis function with a phase factor of the exponenet of the planewave
SUBROUTINE nu_g_overlap(AO,gk,pw_overlap)
  IMPLICIT NONE

  TYPE(AO_function),INTENT(IN)            ::   AO   !Basis function overlap is calculated for
  REAL*8,DIMENSION(3),INTENT(IN)          ::   gk         !Vector of planewave, includes k-point pahse factor
  COMPLEX*16,INTENT(OUT)                  ::   pw_overlap   !resulting overlap of the pw and gaussian basis functions

  REAL*8             ::   gsqr, screen
  INTEGER            ::   igauss,icart
  COMPLEX*16         ::   cartesian

  !REAL*8, PARAMETER        ::  nu_g_screen=45.d0
  !REAL*8, PARAMETER        ::  pi=3.141592653589793238462d0
  !COMPLEX*16, PARAMETER    ::  sqrt_neg_one=(0.d0,1.d0)

  !Start by initializing some values
  gsqr = DOT_PRODUCT(gk,gk)/4.d0
  pw_overlap = (0.d0,0.d0)

  !Calculate the component of the overlap that comes from the cartesian components of the AO
  !The linear combos of cartesian-GTOs make it so that all other constant terms cancel.
  !THIS MAY NOT BE TRUE FOR ALL ORBITAL TYPES
  cartesian = (0.d0,0.d0)
  DO icart=1,AO%ncart
     cartesian = cartesian + AO%cart_coeff(icart) * PRODUCT(gk**AO%cart_mat(icart,:))
  ENDDO

  !Then the contribution of all gaussian orbitals are calculated
  DO igauss=1,AO%num_gauss
     screen = gsqr / AO%alpha(igauss)
     IF( screen .GT. nu_g_screen )EXIT !Since the values of alpha only shrink, if one alpha is too small all following it will be as well
     pw_overlap = pw_overlap + AO%coeff(igauss) * proj_2pi / (AO%alpha(igauss))**(0.75) * (sqrt_neg_one/SQRT(AO%alpha(igauss)))**AO%l * EXP(-screen)
     
  ENDDO

  !Gaussian and cartesian are combined
  pw_overlap = pw_overlap * cartesian

ENDSUBROUTINE nu_g_overlap


!This subrotuine performs an Occupancy Weighted Symmetric Orthogonalization
!This orthogonalizes the vectors represented in 'band_coeff' so that 'overlap_mat' becomes an identity matrix
!This process DOES NOT orhtogonalize the atomic orbitals that make up the bands, simply the bands relative to each other
!The bands are in the same original AO basis
!The orhtogonalizing matrix is constructed as W(WSW)^-1/2, where S is the 'overlap_mat' and W is a diagonal matrix containing the occupancy of each band
SUBROUTINE do_OWSO(overlap_mat,band_occ,band_coeff)
  USE mkl95_BLAS
  USE mkl95_LAPACK
  IMPLICIT NONE

  COMPLEX*16,DIMENSION(:,:),INTENT(IN)    ::  overlap_mat  !Input overlap matrix, used to construct orthogonalizing matrix
  REAL*8,DIMENSION(:),INTENT(IN)          ::  band_occ     !Occupancies of bands.  OWSO preferentially preserves the character of occupied bands.
  COMPLEX*16,DIMENSION(:,:),INTENT(INOUT) ::  band_coeff   !On output represent the orhtogonalized bands in the AO basis

  COMPLEX*16,DIMENSION(SIZE(overlap_mat,1),SIZE(overlap_mat,1)) ::  sym_weight !Used in OWSO equals (WSW), where W is weight_mat. Also used as dummy for BLAS
  REAL*8,DIMENSION(SIZE(overlap_mat,1))   ::  weight_mat   !Diagonal matrix of occupancies to weight orthogonalization.  

  COMPLEX*16,DIMENSION(SIZE(overlap_mat,1),SIZE(overlap_mat,1)) ::  eigvec  !Eigenvectors used in matrix square root
  REAL*8,DIMENSION(SIZE(overlap_mat,1))   ::  eigval                        !Eigenvalues used in matrix square root
  INTEGER   ::  INFO  !Used to identify error on exit of LAPACK eigenvector subroutine

  COMPLEX*16,DIMENSION(SIZE(overlap_mat,1),SIZE(overlap_mat,1)) ::  sqrt_inv    !Inverse of the square root of the matrix (WSW) used in orthogonalization
  COMPLEX*16,DIMENSION(SIZE(band_coeff,1),SIZE(band_coeff,2))  ::  coeff_dummy  !Dummy matrix used in BLAS subroutine

  INTEGER   :: nbands,iband   !Number of bands and counter

  !Start by verifying dimensions of all matrices match
  IF( SIZE(overlap_mat,1) .NE. SIZE(overlap_mat,2) )STOP 'Input overlap matrix is not square in do_OWSO subroutine'
  nbands = SIZE(overlap_mat,1)
  IF( nbands .NE. SIZE(band_occ,1) )STOP 'Input overlap matrix and band occupancies are not the same dimension in do_OWSO subroutine'
  IF( nbands .NE. SIZE(band_coeff,2) )STOP 'Input overlap matrix and band coefficients are not the same dimension in do_OWSO subroutine'

  !Construct weight matrix.
  !Below some cutoff occupancy, 0.005 is used as the occupancy.  
  !If real occupancies are used, the WSW matrix will become linear dependent.  
  !This floor used is still far below full (2) so occupied and unoccupied will be differentiated
  DO iband=1,nbands
     IF( ABS(band_occ(iband)) > 5.d-3 )THEN
        weight_mat(iband) = ABS(band_occ(iband))
     ELSE
        weight_mat(iband) = 5.d-3
     ENDIF
  ENDDO

  !sym_weight = MATMUL(weight_mat, MATMUL(r_mat(:,:,ik),weight_mat))
  !Since weight_mat is a diagonal matrix, the matrix multiplications can be represented as a scaling of either the columns or rows.
  !Scaling columns is multiplication by a diagonal matrix from the right
  DO iband=1,nbands
     sym_weight(:,iband) = overlap_mat(:,iband) * weight_mat(iband)
  ENDDO
  !Scaling rows is multiplication by a diagonal matrix from the left
  DO iband=1,nbands
     sym_weight(iband,:) = weight_mat(iband) * sym_weight(iband,:)
  ENDDO

  !Now obtain the eigenvectors to obtain the square root matrix
  eigvec = sym_weight
  CALL ZHEEV_MKL95(eigvec,eigval,'V','U',INFO)
  IF( INFO /= 0 )THEN
     WRITE(6,*)'The eigenvalues of the overlap matrix could not be computed'
     WRITE(6,*)'INFO on exit of ZHEEV',INFO
     STOP
  ENDIF

  !sqrt_inv = MATMUL(eigvec, MATMUL(diag,CONJG(TRANSPOSE(eigvec))))
  DO iband=1,nbands
     sym_weight(iband,:) = CONJG(eigvec(:,iband)) / SQRT(eigval(iband))  !Sym_weight is now just a dummy matrix for use in the BLAS routine
  ENDDO
  CALL ZGEMM_MKL95(eigvec,sym_weight,sqrt_inv,'N','N',(1.d0,0.d0),(0.d0,0.d0))


  !To obtain the orthogonlaizing matrix an additional multiplication by the weight matrix must be applied from the left
  !sqrt_inv = MATMUL(weight_mat,sqr_rt_inv)
  DO iband=1,nbands
     sqrt_inv(iband,:) = weight_mat(iband) * sqrt_inv(iband,:)
  ENDDO

  !Now the bands can finally be orhtogonalized
  !Note that the original coefficients are lost. They are currently unused but this will have to be rewritten if they become necessary.
  coeff_dummy = band_coeff
  CALL ZGEMM_MKL95(coeff_dummy,sqrt_inv,band_coeff,'N','N',(1.d0,0.d0),(0.d0,0.d0))





END SUBROUTINE do_OWSO




!Now write out information about the system and k-space matrices for input into periodic NBO code of JRS and BDD
!Two files will be written to.
! 1) A header file whose name can be sent manually by the user. This contains information on basis functions and atoms, as well as l-vectors and k-points.
! 2) Unformatted file for the overlap, density and Fock matrices in the projected AO basis.  Unformatted file is used for memory considerations and speed of read/write
SUBROUTINE write_NBO_output(NBO_fn,mat_fn,AO_basis,index_l)
  USE rd_wavefunction
  IMPLICIT NONE

  CHARACTER(128),INTENT(IN)                  ::  NBO_fn
  CHARACTER(64),INTENT(IN)                   ::  mat_fn
  TYPE(AO_function),DIMENSION(:),INTENT(IN)  ::  AO_basis
  INTEGER,DIMENSION(:,:),INTENT(IN)          ::  index_l

  INTEGER  ::  l_half
  INTEGER  ::  i,j,ik,ispin,nu

  INTEGER  ::  nneigh  !Number of neighbor unit cells to print

  OPEN (55,file=NBO_fn)
  !Begin by writing out some general info about the system
  WRITE(55,*)'comment line, output from BDD projection code'
  WRITE(55,'(I8,A7)')n_atom,'#natom'
  WRITE(55,'(I8,A8)')s_dim,'#nbasis'
  WRITE(55,'(I8,A8)')nspins,'#nspins'

  !Not all l-vectors are used, only those corresponding to neighboring unit  cells 
  nneigh = 1
  DO j=1,3 !Don't want any neighbors if there is only one k-point along that direction though.
     IF( kdim(j).GT.1 )nneigh=nneigh*3
  ENDDO
  WRITE(55,'(I8,A4)')nneigh,'#ng'
  WRITE(55,'(I8,A4)')nkpts,'#nk'

  WRITE(55,*)
  
  !Now write out some info about basis function ordering and info about each basis function (ie. set of exponents and angular momentum)
  !The index of the first basis function on each atom.
  !There is an additional number written out for the next index AFTER the end of the basis set.
  DO j=1,n_atom+1
     WRITE(55,'(I5)', ADVANCE='no')atom_basis(j)
  ENDDO
  WRITE(55,'(A12)',ADVANCE='NO')'#ibasismap'
  WRITE(55,*)

  !Indices are the same for atoms with the same exponents.
  !This counts up over different atoms.
  DO j=1,s_dim
     WRITE(55,'(I5)',ADVANCE='NO')AO_basis(j)%level
  ENDDO
  WRITE(55,'(A11)',ADVANCE='NO')'#ishellmap'
  WRITE(55,*)

  !l-quantum number of each basis set
  DO j=1,s_dim
     WRITE(55,'(I2)',ADVANCE='NO')AO_basis(j)%l
  ENDDO
  WRITE(55,'(A7)',ADVANCE='NO')'#ilmap'
  WRITE(55,*)

  !Then some info about the atoms present in the system (symbol and atomic number)
  !Atomic symbol of each atom.
  DO j=1,n_atom
     WRITE(55,'(A3)',ADVANCE='NO')atom_symbol(j)
  ENDDO
  WRITE(55,'(A8)',ADVANCE='NO')'#symbols'
  WRITE(55,*)

  !Atomic number of each atom.
  DO j=1,n_atom
     WRITE(55,'(I3)',ADVANCE='NO')atomic_number(j)
  ENDDO
  WRITE(55,'(A8)',ADVANCE='NO')'#iatnum'
  WRITE(55,*)

  !Number of electrons fully represented by the pseudopotential for each atom.
  DO j=1,n_atom
     WRITE(55,'(F6.1)',ADVANCE='NO')atom_valence(j)
  ENDDO
  WRITE(55,'(A8)',ADVANCE='NO')'#iatval'
  WRITE(55,*)

  WRITE(55,*)

  !Lattice parameters
  DO j=1,3
     WRITE(55,'(3F18.10)',ADVANCE='NO')a(:,j)
     IF( j == 1 )WRITE(55,'(A)',ADVANCE='NO')' #lattice vectors'
     WRITE(55,*)
  ENDDO
  WRITE(55,*)

  !Atom positions
  DO j=1,n_atom
     WRITE(55,'(3F18.10)',ADVANCE='NO')atoms(:,j)
     IF( j == 1 )WRITE(55,'(A)',ADVANCE='NO')' #atom positions'
     WRITE(55,*)
  ENDDO
  WRITE(55,*)

  !Write out the l-vectors
  !The ordering of l_vectors expected by the NBO code is different than that used here, so they will be written out non-linearly
  !The first must be (/0,0,0/)
  !Then they can proceed in any order but must be printed in pairs of negatives.
  !Only those of neighboring cells will be included as those are all that will be used in the NBO analysis in searching for spanning bonds
  l_half = (SIZE(index_l,2)+1)/2
  WRITE(55,'(3I3,A)')index_l(:,l_half),' #real space unit cells to search'
  IF( .NOT. gamma_point )THEN
     DO j=1,l_half-1
        !IF( surface .AND. index_l(3,l_half-j) .NE. 0 )GOTO 55
        DO i=1,3 !Generalize old surface test to all dimensions
           IF( kdim(i).EQ.1 .AND. index_l(i,l_half-j).NE.0 )GOTO 55
        ENDDO
        IF( MAXVAL(ABS(index_l(:,l_half-j))) < 2 )THEN
           WRITE(55,'(3I3)')index_l(:,l_half - j)
           WRITE(55,'(3I3)')index_l(:,l_half + j)
        ENDIF
55      CONTINUE
     ENDDO
  ENDIF

  WRITE(55,*)

  !Write out the k-vectors
  !The ordering will be that output from VASP.  Gamma point is always first.
  !The density and overlap matrices will be written out in this same order.
  !The weight of each k-point is also included. Unless hihger symmetry is included, this is 1/nkpts for the Gamma point and 2/nkpts for all other.
  DO ik=1,nkpts
     WRITE(55,'(4F18.10)',ADVANCE='NO')kpt(:,ik),kpt_wt(ik)
     IF( ik == 1 )WRITE(55,'(A)',ADVANCE='NO')' #k-points'
     WRITE(55,*)
  ENDDO
  WRITE(55,*)

  WRITE(55,'(A64)')mat_fn  !'NBO_mat.out'
  WRITE(55,*)

  !Basis set information
  WRITE(55,*)'Basis set information. Used to produce output file for visualtion.'
  WRITE(55,*)
  DO nu=1,s_dim
     WRITE(55,'(I)')AO_basis(nu)%num_gauss
     WRITE(55,'(10D17.9)')AO_basis(nu)%alpha
     WRITE(55,'(10D17.9)')AO_basis(nu)%coeff
     WRITE(55,'(10D17.9)')AO_basis(nu)%norm
     WRITE(55,'(3D16.8)')AO_basis(nu)%pos
     !WRITE(35,'(3I)')AO_basis(nu)%lmn
     WRITE(55,'(I)')AO_basis(nu)%ncart
     DO j=1,AO_basis(nu)%ncart
        WRITE(55,'(D17.9,3I)')AO_basis(nu)%cart_coeff(j),AO_basis(nu)%cart_mat(j,:)
     ENDDO

     WRITE(55,*)
  ENDDO


  CLOSE(55)

  OPEN(66,file=mat_fn,FORM='UNFORMATTED')

  !Varios Bloch-space matrices are now written out for the bands in the projected AO basis.
  !All of these are written out in a lower triangular form
  !The ordering 1)overlap, 2)density, 3)fock is expected by the NBO code
  CALL write_sym_matrices(66,bloch_s_mat)

  DO ispin=1,nspins  
     CALL write_sym_matrices(66,bloch_density(:,:,:,ispin))
  ENDDO

  DO ispin=1,nspins
     CALL write_sym_matrices(66,bloch_fock(:,:,:,ispin))
  ENDDO


  CLOSE(66)


END SUBROUTINE write_NBO_output



!This subroutine writes out the contents of a Bloch space matrix
!At each k-point, the matrices are Hermitian so only the lower triangular part is written out.
!An unformatted output file is used for memory considerations.
SUBROUTINE write_sym_matrices(fn,matrix)
  IMPLICIT NONE

  INTEGER, INTENT(IN)     ::   fn   !File number index of file that needs to be written to
  COMPLEX*16, INTENT(IN)  ::   matrix(:,:,:)  !periodic matrix to write out
  INTEGER                 ::   nk   !Total number of k vectors to write out a matrix for
  INTEGER                 ::  i,j,ik

  IF( SIZE(matrix,1) .NE. SIZE(matrix,2) )STOP 'A non-square matrix has been passed to write_sym_matrices'

  nk = SIZE(matrix,3)

  DO ik=1,nk
     DO i=1,s_dim
       WRITE(fn)matrix(i,1:i,ik)
     ENDDO
  ENDDO

END SUBROUTINE write_sym_matrices







!This subroutine performs a Bloch transform of 'real_mat' to 'bloch_mat'
!The arrays 'l_index' and 'k-index' should coincide with 'real_mat' and 'bloch_mat' respectively,
!  where the n-th vector in 'l_index(:,n)' describes the n-th matrix of 'real_mat(:,:,n)'.
!This procedure can be performed using generalized Fourier Transform routines, but we have found this is already as fast as possible.
!Additionally, using only those k-points of the PW calculation prevents unecessary calculations into many k-points,
!  followed by extrapolation to the necessary ones.
SUBROUTINE real_to_bloch(real_mat,bloch_mat,l_index,k_index)
  IMPLICIT NONE

  REAL*8,DIMENSION(:,:,:),INTENT(IN)      ::  real_mat
  COMPLEX*16,DIMENSION(:,:,:),INTENT(OUT) ::  bloch_mat
  INTEGER,DIMENSION(:,:),INTENT(IN)       ::  l_index
  REAL*8,DIMENSION(:,:),INTENT(IN)        ::  k_index

  INTEGER   ::  l_half
  INTEGER   ::  nk,nl
  INTEGER   ::  ik,il,j

  COMPLEX*16  ::  arg
  !REAL*8, PARAMETER        ::  pi=3.141592653589793238462d0
  !COMPLEX*16, PARAMETER    ::  sqrt_neg_one=(0.d0,1.d0)


  !Start by making sure all dimensionalities match
  !The size of each real_mat and bloch_mat component matrix must be the same size
  !There must be an l_vector for each real matrix and a k-vector for each bloch matrix
  IF( SIZE(real_mat,1) .NE. SIZE(bloch_mat,1) .OR. SIZE(real_mat,2) .NE. SIZE(bloch_mat,2) )STOP 'real_mat and bloch_mat passed to real_to_bloch are different sizes'
  IF( SIZE(real_mat,3) .NE. SIZE(l_index,2) )STOP 'real_mat and l_index passed to real_to_bloch do not match in dimension'
  IF( SIZE(bloch_mat,3) .NE. SIZE(k_index,2) )STOP 'bloch_mat and k_index passed to real_to_bloch do not match in dimension'

  nk = SIZE(bloch_mat,3)
  nl = SIZE(real_mat,3)
  l_half = (nl+1)/2

  !The symmetry of l_vectors is once again used  arg(l_vector) = CONJG(arg(-l_vector))
  !The k-points from the PW calculation are used, otherwise extrapolation to those points would be necessary.
  bloch_mat = 0.d0

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ik,il,arg)
!$OMP DO SCHEDULE(STATIC)
  DO ik=1,nk
     !Use symmetry of l-vectors
     DO il=1,l_half-1
        arg = EXP(sqrt_neg_one*2.d0*pi*DOT_PRODUCT(k_index(:,ik),l_index(:,il)))
        bloch_mat(:,:,ik) = bloch_mat(:,:,ik) + real_mat(:,:,il)*arg
        bloch_mat(:,:,ik) = bloch_mat(:,:,ik) + real_mat(:,:,(nl+1-il))*CONJG(arg)
     ENDDO
     !Account for central unit cell (0,0,0), which is its own symmetric unit cell (don't double count)
     bloch_mat(:,:,ik) = bloch_mat(:,:,ik) + real_mat(:,:,l_half)
  ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL



END SUBROUTINE real_to_bloch



!This is a VERY similar procedure that above, simply the inverse process
!The main difference in the procedure comes from using the inversion symmetry of the Brillouin zone
!For ALL quantities it is true that A(k) = CONJG(A(-k)) for any point k in the Brillouin zone  (Time Reversal Symmetry)  SOURCE!?!?!?!?!?
!Thus the PW calculation is only done at half of the k-points (except for k=0)
!These -k components must still be included for the reverse Bloch-transform to be complete
SUBROUTINE bloch_to_real(real_mat,bloch_mat,l_index,k_index)
  IMPLICIT NONE

  REAL*8,DIMENSION(:,:,:),INTENT(OUT)      ::  real_mat
  COMPLEX*16,DIMENSION(:,:,:),INTENT(IN) ::  bloch_mat
  INTEGER,DIMENSION(:,:),INTENT(IN)       ::  l_index
  REAL*8,DIMENSION(:,:),INTENT(IN)        ::  k_index

  INTEGER   ::  l_half
  INTEGER   ::  nk,nl
  INTEGER   ::  ik,il

  COMPLEX*16  ::  arg
  REAL*8, PARAMETER        ::  pi=3.141592653589793238462d0
  COMPLEX*16, PARAMETER    ::  sqrt_neg_one=(0.d0,1.d0)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!WRITE THIS SUBROUTINE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  WRITE(6,*)'bloch_to_real subroutine needs to be written'
  real_mat = 0.d0
  STOP


END SUBROUTINE



!==============================================
!Utility Math Functions
!==============================================


!Simple double factorial function.  Using an array serves only as a minimal speed up so I will just keep the function form.
INTEGER FUNCTION dble_factor(x)

  INTEGER  :: x,n

  IF(x .LT. -1)THEN
     WRITE(6,*)'Can not campute the double factorial of a negative number, below -1'
     STOP
  ENDIF

  dble_factor = 1
  n=x

  DO
     IF( n .LE. 1 )EXIT
     dble_factor = dble_factor * n
     n = n - 2
  ENDDO


END FUNCTION dble_factor


!Binomial coefficient function
INTEGER FUNCTION binom_coeff(n,m) !'Selecting m from n'

  INTEGER  ::  n,m

  IF( m > n )THEN
     WRITE(6,*)'impossible to take this binomial coefficient, second entry can not be greater than first'
     STOP
  ENDIF

  binom_coeff = factor(n) / (factor(m) * factor(n-m))


END FUNCTION binom_coeff

!Simple factorial function.  Using an array serves only as a minimal speed up so I will just keep the function form.
INTEGER FUNCTION factor(x)

  INTEGER  :: x,n

  IF(x .LT. 0 ) STOP 'Can not campute the factorial of a negative number'

  factor = 1
  n=x

  DO 
     IF( n .LE. 1 )EXIT
     factor = factor * n       
     n = n - 1                 
  ENDDO

  
END FUNCTION factor




END MODULE projection_shared
