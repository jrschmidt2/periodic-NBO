MODULE shared
  IMPLICIT NONE

  TYPE AO_function
     REAL*8, ALLOCATABLE, DIMENSION(:) :: norm, alpha, coeff
     INTEGER                :: ngauss
     !INTEGER, DIMENSION(3)  :: lmn
     INTEGER                :: atom       
     REAL*8, DIMENSION(3)   :: pos        
     INTEGER                :: level  !Keeps track of what basis functions share the same sets of exponents, on the same atom.
     INTEGER                :: l,m     !l- and m-quantum numbers of basis function.  m actually runs from 0,2*l+1 and is onyl and index.

     INTEGER                :: ncart
     INTEGER,ALLOCATABLE    :: cart_mat(:,:)
     REAL*8,ALLOCATABLE     :: cart_coeff(:)
  
  END TYPE AO_function

  REAL*8, PARAMETER  ::  bohr = 0.529177249d0
  REAL*8, PARAMETER  ::  hartree = 27.21138386d0

  !TYPE nbo_output
  !   INTEGER, DIMENSION(:), ALLOCATABLE       ::  iatnum
  !   CHARACTER(2), DIMENSION(:), ALLOCATABLE  ::  symbols
  !   INTEGER, DIMENSION(:), ALLOCATABLE       ::  ishellmap
  !   INTEGER, DIMENSION(:), ALLOCATABLE       ::  ibasismap
  !   INTEGER, DIMENSION(:), ALLOCATABLE       ::  ilmap
  !END TYPE nbo_output


  INTERFACE temp_read
     MODULE PROCEDURE temp_read_int
     MODULE PROCEDURE temp_read_real
     MODULE PROCEDURE temp_read_comp
  END INTERFACE




  CONTAINS


!For a GTO with alphas=expon and cartesian powers=lmn, the norm is--> norm_vec
!The following subroutine comes from: Clementi and Davis, Journal of Comp. Physics; 2,223-244(1967)
SUBROUTINE norm_calc(expon, norm_vec, l) 
  IMPLICIT NONE
  
  REAL*8, DIMENSION(:), INTENT(IN)  :: expon
  REAL*8, DIMENSION(:), INTENT(OUT) :: norm_vec
  INTEGER, INTENT(IN)               :: l
  
  REAL*8, PARAMETER  ::  pi=3.14159265358979323846264338d0

  
  norm_vec = ( expon**(2*l+3) * 2.d0**(4*l+3) / pi**3 )**(0.25d0)


END SUBROUTINE norm_calc



!Based on the l- and m- indices of the basis function 'AO' the approproate cartesian component is assigned 
!The number of cartesian components as well as each's coefficient and exponents are all assigned
!There is not much fancy about this, just chucking the indices and assigning value from there
!Wherever orbitals could be grouped together I have (same ncart's, same coeff's, etc.)
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
     STOP 'can only interpret up to d-type orbitals in interp_cart'
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


!The following 3 subroutines are for reading in an array which has arbitrary formatting in a file
!They are interfaced for each particular type of variable (integer, real, complex)
!The general structure is to read all the necessary information in then assign that properly into the array.

!First is for integer type data.
SUBROUTINE temp_read_int(fn,array)
  IMPLICIT NONE

  INTEGER,INTENT(IN)  ::  fn
  INTEGER,DIMENSION(:,:),INTENT(OUT)  ::  array

  INTEGER,DIMENSION(:),ALLOCATABLE  ::  temp
  INTEGER  ::  i,j,k

  ALLOCATE(temp(SIZE(array,1)*SIZE(array,2)))

  READ(fn,*)temp

  k=0
  DO i=1,SIZE(array,1)
     DO j=1,SIZE(array,2)
        k=k+1
        array(i,j)=temp(k)
     ENDDO
  ENDDO

END SUBROUTINE temp_read_int

!Second is for real data
SUBROUTINE temp_read_real(fn,array)
  IMPLICIT NONE

  INTEGER,INTENT(IN)  ::  fn
  REAL*8,DIMENSION(:,:),INTENT(OUT)  ::  array

  REAL*8,DIMENSION(:),ALLOCATABLE  ::  temp
  INTEGER  ::  i,j,k

  ALLOCATE(temp(SIZE(array,1)*SIZE(array,2)))

  READ(fn,*)temp

  k=0
  DO i=1,SIZE(array,1)
     DO j=1,SIZE(array,2)
        k=k+1
        array(i,j)=temp(k)
     ENDDO
  ENDDO



END SUBROUTINE temp_read_real


!Finally is for complex type data.
!Note there is a third input, 'flag'.
!This is to distinguish whether the imaginary components of array will be written.
!So if the data structure 'array' is complex, but all its values are real, there will not in general be zeroes to be read in to 'temp.'
SUBROUTINE temp_read_comp(fn,array,flag)
  IMPLICIT NONE

  INTEGER,INTENT(IN)  ::  fn
  COMPLEX*16,DIMENSION(:,:),INTENT(OUT)  ::  array
  INTEGER,INTENT(IN)  ::  flag   !1-only readign in real comp;2-reading in both

  REAL*8,DIMENSION(:),ALLOCATABLE  ::  temp
  INTEGER  ::  i,j,k

  COMPLEX*16,PARAMETER  ::  imag_i = (0.d0,1.d0)

  ALLOCATE(temp(flag*SIZE(array,1)*SIZE(array,2)))

  READ(fn,*)temp
  k=0
  DO i=1,SIZE(array,1)
     DO j=1,SIZE(array,2)
        k=k+1
        array(i,j)=temp(k)
        IF( flag .EQ. 2)THEN
           k=k+1
           array(i,j)=array(i,j)+temp(k)*imag_i
        ENDIF
     ENDDO
  ENDDO

END SUBROUTINE temp_read_comp



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
  INTEGER   ::  ik,il

  COMPLEX*16  ::  arg
  REAL*8, PARAMETER        ::  pi=3.141592653589793238462d0
  COMPLEX*16, PARAMETER    ::  sqrt_neg_one=(0.d0,1.d0)


  !Start by making sure all dimensionalities match
  !The size of each real_mat and bloch_mat component matrix must be the same size
  !There must be an l_vector for each real matrix and a k-vector for each bloch matrix
  IF( SIZE(real_mat,1) .NE. SIZE(bloch_mat,1) .OR. SIZE(real_mat,2) .NE. SIZE(bloch_mat,2) )STOP 'real_mat and bloch_mat passed to real_to_bloch are different sizes'
  IF( SIZE(real_mat,3) .NE. SIZE(l_index,2) )STOP 'real_mat and l_index passed to real_to_bloch do not match in dimension'
  WRITE(6,*)SIZE(bloch_mat,3),SIZE(k_index,2)
  IF( SIZE(bloch_mat,3) .NE. SIZE(k_index,2) )STOP 'bloch_mat and k_index passed to real_to_bloch do not match in dimension'

  nk = SIZE(bloch_mat,3)
  nl = SIZE(real_mat,3)
  l_half = (nl+1)/2

  !The symmetry of l_vectors is once again used  arg(l_vector) = CONJG(arg(-l_vector))
  !The k-points from the PW calculation are used, otherwise extrapolation to those points would be necessary.
  bloch_mat = 0.d0

!!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ik,il,arg)
!!$OMP DO SCHEDULE(STATIC)
  DO ik=1,nk
     DO il=1,l_half
        arg = EXP(sqrt_neg_one*2.d0*pi*DOT_PRODUCT(k_index(:,ik),l_index(:,il)))
        bloch_mat(:,:,ik) = bloch_mat(:,:,ik) + real_mat(:,:,il)*arg
        IF( il /= l_half )THEN
           bloch_mat(:,:,ik) = bloch_mat(:,:,ik) + real_mat(:,:,(nl+1-il))*CONJG(arg)
        ENDIF
     ENDDO
  ENDDO
!!$OMP END DO NOWAIT
!!$OMP END PARALLEL



END SUBROUTINE real_to_bloch

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
  !This floor used is still far below full (1) so occupied and unoccupied will be differentiated
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
  !WRITE(6,*)'eigvalues of symm weight matrix'
  DO iband=1,nbands
     sym_weight(iband,:) = CONJG(eigvec(:,iband)) / SQRT(ABS(eigval(iband)))  !Sym_weight is now just a dummy matrix for use in the BLAS routine
     !WRITE(6,*)eigval(iband)
  ENDDO
  !WRITE(6,*)
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

  !WRITE(6,*)'orthogonalized band'
  !DO iband=1,nbands
  !   WRITE(6,'(28F6.3)')REAL(band_coeff(:,iband))
  !   !WRITE(6,*)
  !ENDDO
  !WRITE(6,*)


END SUBROUTINE do_OWSO




!This subroutine writes out the contents of a Bloch space matrix
!At each k-point, the matrices are Hermitian so only the lower triangular part is written out.
!An unformatted output file is used for memory considerations.
SUBROUTINE write_sym_matrices(fn,matrix)
  IMPLICIT NONE

  INTEGER, INTENT(IN)     ::   fn   !File number index of file that needs to be written to
  COMPLEX*16, INTENT(IN)  ::   matrix(:,:,:)  !periodic matrix to write out
  INTEGER                 ::   nk   !Total number of k vectors to write out a matrix for
  INTEGER                 ::   ndim  !Dimension of matrices
  INTEGER                 ::  i,j,ik


  IF( SIZE(matrix,1) .NE. SIZE(matrix,2) )STOP 'A non-square matrix has been passed to write_sym_matrices'
  ndim = SIZE(matrix,1)

  nk = SIZE(matrix,3)

  DO ik=1,nk
     DO i=1,ndim
       WRITE(fn)matrix(i,1:i,ik)
     ENDDO
  ENDDO

END SUBROUTINE write_sym_matrices





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





END MODULE shared
