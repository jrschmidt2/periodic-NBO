MODULE nbo_shared

  TYPE nbo_input
     INTEGER :: natom,nbasis,ng,nk,nspins
     !The mapping of the basis functions to atoms (i.e. first basis function for each atom)
     !as well as the angular momentum quantum number (l) for each basis function
     INTEGER,DIMENSION(:),POINTER :: ibasismap, ishellmap,ilmap,immap,iatnum
     REAL*8,DIMENSION(:),POINTER  :: iatval
     !The mapping of the ng matrices into spatial cells
     INTEGER,DIMENSION(:,:),POINTER :: indexg
     !The k-points used, as coefficients of reciprocal lattice vectors
     REAL*8,DIMENSION(:,:),POINTER  ::  kpt
     !Weight of the k-point in the original PW calculation
     REAL*8,DIMENSION(:),POINTER    ::  kpt_wt
     !A listing of the symbol of each atom
     CHARACTER*2,DIMENSION(:),POINTER :: symbols
     !The real-space density and overlap matrices, rho^{0g}_{mu nu}
     REAL*8,DIMENSION(:,:,:,:),POINTER :: rho0, fock0
     REAL*8,DIMENSION(:,:,:),POINTER :: s0
     !Bloch-space matrices. Density and Fock matrices can have separate spin components
     COMPLEX*16,DIMENSION(:,:,:,:),POINTER :: rhok,fockk
     COMPLEX*16,DIMENSION(:,:,:),POINTER :: sk
     REAL*8,DIMENSION(:,:,:),POINTER :: transform !transformation matrix between AO and NAO basis
  END TYPE nbo_input

  !The maximum angular momentum
  INTEGER :: lmax
  PARAMETER (lmax=3)

  TYPE nbo_output_info
     REAL*8,ALLOCATABLE          ::  occ(:)
     REAL*8,ALLOCATABLE          ::  coeff(:,:,:)
     CHARACTER(128),ALLOCATABLE  ::  label(:)
  END TYPE


  !REAL*8,ALLOCATABLE   ::  bond_coeff(:,:,:),bond_occ(:)  !Coefficients of NBO's in original basis
  !REAL*8,ALLOCATABLE   ::  lp_coeff(:,:,:),lp_occ(:)
  !REAL*8,ALLOCATABLE   ::  ryd_coeff(:,:,:),ryd_occ(:)
  !CHARACTER(128),ALLOCATABLE :: bond_lbl(:),lp_lbl(:),ryd_lbl(:)

  LOGICAL   ::  real_init
  INTEGER   ::  nkx,nky,nkz

  !Defined types that are used for visualization
  !Need to be in shared module since such variables are used in main program and visual subroutine

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

  TYPE vis_cont_struct
     INTEGER               ::  vis_start,vis_end
     INTEGER,DIMENSION(3)  ::  mesh
     INTEGER,DIMENSION(3)  ::  box_int
     REAL*8                ::  origin_fact
     LOGICAL               ::  density
  END TYPE vis_cont_struct



  !Various thresholds and options used in the NBO analysis
  REAL*8 :: core_thresh,nbo_1c_thresh,nbo_2c_thresh,polarization_thresh
  !PARAMETER (core_thresh=0.995d0,nbo_1c_thresh=0.80d0,nbo_2c_thresh=0.925d0,polarization_thresh=0.99d0)
  PARAMETER (core_thresh=0.995d0,polarization_thresh=0.99d0)

  !The threshold for whether a located hybrid is unique enough, i.e. what fraction cannot be
  !expressed in terms of previously located hybrids
  REAL*8 :: prjexp_thresh
  PARAMETER (prjexp_thresh=0.2d0)
  !Whether or not to force orthogonalization of the hybrids -- this should be OFF for
  !hypervalent ionic-like structures
  LOGICAL :: orthogonalize_hybrids
  PARAMETER (orthogonalize_hybrids=.TRUE.)




  INTERFACE remap_matrix
     MODULE PROCEDURE remap_matrix_real
     MODULE PROCEDURE remap_matrix_complex
  END INTERFACE

  CONTAINS

  !
  !Read an triangular matrix out of the file, and convert to a dense one
  !
  SUBROUTINE read_bloch_triangular(unit,nbasis,matlow)
    IMPLICIT NONE
    INTEGER ::  unit, nbasis
    COMPLEX*16 :: matlow(nbasis,nbasis)
    
    COMPLEX*16 :: temp(nbasis*(nbasis+1)/2)
    INTEGER :: i,j,k
    COMPLEX*16 ::  sqrt_neg_one

    sqrt_neg_one = (0.d0,1.d0)

    !READ(unit) temp
    k=1
    DO i=1,nbasis
       READ(unit)matlow(i,1:i)
       DO j=1,i
          !matlow(i,j)=temp(k)
          matlow(j,i)=CONJG(matlow(i,j))
          k=k+1
       ENDDO
    ENDDO
  END SUBROUTINE read_bloch_triangular


  SUBROUTINE read_real_triangular(unit,nbasis,matlow,mathigh)
    IMPLICIT NONE
    INTEGER  :: unit,nbasis
    REAL*8   :: matlow(nbasis,nbasis),mathigh(nbasis,nbasis)
     
    REAL*8 :: temp(nbasis*(nbasis+1)/2)
    INTEGER  ::  i,j,k


   !WRITE(6,*)'Test with in read_real_tri'

    READ(unit,*) (temp(i),i=1,nbasis*(nbasis+1)/2)
    k=1
    DO i=1,nbasis
       DO j=1,i
          !WRITE(6,*)k,temp(k)
          matlow(i,j)=temp(k)
          mathigh(j,i)=temp(k)
          k=k+1
       ENDDO
    ENDDO
  END SUBROUTINE read_real_triangular




  !
  !Calculates the number of electrons in the unit cell
  !
  SUBROUTINE calc_nelec(inp,is_orthog)
    USE BLAS95
    USE matutil
    IMPLICIT NONE
    TYPE(nbo_input) :: inp
    LOGICAL :: is_orthog !true if we are already in an orthogonal NAO basis (i.e. we read the NAOs from a checkpoint file)

    COMPLEX*16,DIMENSION(SIZE(inp%rhok,1),SIZE(inp%rhok,1))  ::  nelec_dummy
    INTEGER :: ik, nk, ispin
    REAL*8 :: nelec, nenergy
    
    nk=inp%nk
    nelec=0.d0
    nenergy=0.d0

    DO ispin=1,inp%nspins
       DO ik=1,nk
          CALL ZGEMM_F95(inp%rhok(:,:,ik,ispin),inp%sk(:,:,ik),nelec_dummy,'N','C',(1.d0,0.d0),(0.d0,0.d0))
          nelec=nelec+inp%kpt_wt(ik)*mattrace(nelec_dummy) 
          CALL ZGEMM_F95(inp%rhok(:,:,ik,ispin),inp%fockk(:,:,ik,ispin),nelec_dummy,'N','C',(1.d0,0.d0),(0.d0,0.d0))
          nenergy=nenergy+inp%kpt_wt(ik)*mattrace(nelec_dummy)
       ENDDO
    ENDDO
    PRINT *, 'Total number of elec. from input Bloch space matrices: ', nelec
    PRINT *, 'Total energy average from input Bloch space matrices:  ', nenergy
    PRINT *

  END SUBROUTINE calc_nelec

  !
  !Gives a remapping of the matrices based on the given mask, first
  !those which satisfy the maks, then those that do not.
  !
  SUBROUTINE partition_basis(mask,nbasis,ntrue,nfalse,ibasisremap)
    IMPLICIT NONE
    LOGICAL :: mask(:)
    INTEGER :: ntrue,nfalse,nbasis,ibasisremap(:),ibasis
    
    ntrue=0
    DO ibasis=1,nbasis
       IF (mask(ibasis)) THEN
          ntrue=ntrue+1
          ibasisremap(ibasis)=ntrue
       ENDIF
    ENDDO
    nfalse=0
    DO ibasis=1,nbasis
       IF (.NOT.mask(ibasis)) THEN
          nfalse=nfalse+1
          ibasisremap(ibasis)=ntrue+nfalse
       ENDIF
    ENDDO
  END SUBROUTINE partition_basis


  FUNCTION nmb_count(inp,iatom) RESULT (imbcount)
    IMPLICIT NONE
    TYPE(nbo_input)  ::  inp
    INTEGER :: iatom
    INTEGER :: iatnum, imbcount

    iatnum = inp%iatnum(iatom)

    IF (iatnum.GE.87) THEN
       STOP 'Principal quanutm number n > 6 to implemented'
    ELSEIF (iatnum.GE.55)THEN
       imbcount=43
    ELSEIF (iatnum.GE.37) THEN
       imbcount=27
    ELSEIF (iatnum.GE.19) THEN
       imbcount=18
    ELSEIF (iatnum.GE.11) THEN
       imbcount=9
    ELSEIF (iatnum.GE.3) THEN
       imbcount=5
    ELSE
       imbcount=1
    ENDIF

    !WRITE(6,*)imbcount
    imbcount = imbcount - (iatnum - inp%iatval(iatom))/2.d0
    !WRITE(6,*)(iatnum - inp%iatval(iatom))/2

    !WRITE(6,*)imbcount
    !WRITE(6,*)

  END FUNCTION NMB_COUNT

  !Reorder a given matrix / vector given a remapping vector; invertable via inv
  FUNCTION remap_matrix_complex(A,iremap,inv) RESULT (C)
    IMPLICIT NONE
    COMPLEX*16,DIMENSION(:,:) :: A
    INTEGER,DIMENSION(:) :: iremap
    LOGICAL :: inv
    COMPLEX*16,DIMENSION(SIZE(A,1),SIZE(A,2)) :: C
    INTEGER :: i,j,nbasis
    
    nbasis=SIZE(a,1)
    DO i=1,nbasis
       DO j=1,nbasis
          IF (inv) THEN
             C(i,j)=A(iremap(i),iremap(j))
          ELSE
             C(iremap(i),iremap(j))=A(i,j)
          ENDIF
       ENDDO
    ENDDO
  END FUNCTION remap_matrix_complex
  
  FUNCTION remap_matrix_real(A,iremap,inv) RESULT (C)
    IMPLICIT NONE
    REAL*8,DIMENSION(:,:) :: A
    INTEGER,DIMENSION(:) :: iremap
    LOGICAL :: inv
    REAL*8,DIMENSION(SIZE(A,1),SIZE(A,2)) :: C
    INTEGER :: i,j,nbasis
    
    nbasis=SIZE(a,1)
    DO i=1,nbasis
       DO j=1,nbasis
          IF (inv) THEN
             C(i,j)=A(iremap(i),iremap(j))
          ELSE
             C(iremap(i),iremap(j))=A(i,j)
          ENDIF
       ENDDO
    ENDDO
  END FUNCTION remap_matrix_real

  FUNCTION remap_vector(A,iremap,inv) RESULT (C)
    IMPLICIT NONE
    REAL*8,DIMENSION(:) :: A
    INTEGER,DIMENSION(:) :: iremap
    LOGICAL :: inv
    REAL*8,DIMENSION(SIZE(A,1)) :: C
    INTEGER :: i,nbasis
    
    nbasis=SIZE(a,1)
    DO i=1,nbasis
       IF (inv) THEN
          C(i)=A(iremap(i))
       ELSE
          C(iremap(i))=A(i)
       ENDIF
    ENDDO
  END FUNCTION remap_vector


  !If a set of kpts and g vectors was read in from the input file, the forward and reverse bloch transform can be done more efficiently doing it manually.
  !This will be possible for any system that that was projected from VASP results by BDD code. 
  !These subroutines can not be called for Gaussian results, since a set k-point grid is not a part of the results.

  !This subroutine transforms a bloch space matrix into a real space periodically resolved matrix
  SUBROUTINE bloch_to_real(inp,bloch_mat,real_mat,kpt,gvec)
    IMPLICIT NONE
    TYPE(nbo_input)  ::  inp
    COMPLEX*16,DIMENSION(:,:,:),INTENT(IN)     ::  bloch_mat
    REAL*8,DIMENSION(:,:,:),INTENT(OUT)        ::  real_mat
    REAL*8,DIMENSION(:,:),INTENT(IN)           ::  kpt
    INTEGER,DIMENSION(:,:),INTENT(IN)          ::  gvec

    COMPLEX*16      ::  arg
    INTEGER         ::  ng, nk
    INTEGER         ::  ig, ik
    REAL*8,PARAMETER           ::  pi=4.d0*ATAN(1.d0)
    COMPLEX*16,PARAMETER       ::  sqrt_neg_one=(0.d0, 1.d0)

    IF( SIZE(bloch_mat,1) /= SIZE(real_mat,1) )STOP 'Improper matching of the sizes of real and bloch space matrices'

    nk = SIZE(bloch_mat,3)
    IF( nk /= SIZE(kpt,2) )STOP 'Improper matching of bloch_mat and kpt array in bloch_to_real sub'

    ng = SIZE(real_mat,3)
    IF( ng /= SIZE(gvec,2) )STOP 'Improper matching of real_mat and gvec array in bloch_to_real sub'

    real_mat = 0.d0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ik,ig,arg)
!$OMP DO SCHEDULE(STATIC)
    DO ig=1,ng
       DO ik=1,nk
          arg = EXP(-sqrt_neg_one*2.d0*pi*DOT_PRODUCT(kpt(:,ik), gvec(:,ig)))
          real_mat(:,:,ig) = real_mat(:,:,ig) + bloch_mat(:,:,ik)*arg
          IF( inp%kpt_wt(ik) /= inp%kpt_wt(1) )THEN  !For non gamma point matrices use inversion symmetry Rho(k) = Conjg{Rho(-k)}
             real_mat(:,:,ig) = real_mat(:,:,ig) + CONJG(bloch_mat(:,:,ik)*arg)
          ENDIF
       ENDDO
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
    real_mat = real_mat * inp%kpt_wt(1)

  END SUBROUTINE bloch_to_real



  !This subroutine takes a real space periodic matrix and converts it to one in Bloch space at set k-points
  SUBROUTINE real_to_bloch(bloch_mat,real_mat,kpt,gvec)
    IMPLICIT NONE
    COMPLEX*16,DIMENSION(:,:,:),INTENT(OUT)    ::  bloch_mat
    REAL*8,DIMENSION(:,:,:),INTENT(IN)         ::  real_mat
    REAL*8,DIMENSION(:,:),INTENT(IN)           ::  kpt
    INTEGER,DIMENSION(:,:),INTENT(IN)          ::  gvec
    COMPLEX*16      ::  arg
    INTEGER         ::  ng, nk
    INTEGER         ::  ig, ik
    REAL*8,PARAMETER           ::  pi=4.d0*ATAN(1.d0)
    COMPLEX*16,PARAMETER       ::  sqrt_neg_one=(0.d0, 1.d0)
    IF( SIZE(bloch_mat,1) /= SIZE(real_mat,1) )STOP 'Improper matching of the sizes of real and bloch space matrices'
    nk = SIZE(bloch_mat,3)
    IF( nk /= SIZE(kpt,2) )STOP 'Improper matching of bloch_mat and kpt array in bloch_to_real sub'
    ng = SIZE(real_mat,3)
    IF( ng /= SIZE(gvec,2) )STOP 'Improper matching of real_mat and gvec array in bloch_to_real sub'
    bloch_mat = 0.d0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ik,ig,arg)
!$OMP DO SCHEDULE(STATIC)
    DO ik=1,nk
       DO ig=1,ng
          arg = EXP(sqrt_neg_one*2.d0*pi*DOT_PRODUCT(kpt(:,ik), gvec(:,ig)))
          bloch_mat(:,:,ik) = bloch_mat(:,:,ik) + real_mat(:,:,ig)*arg
       ENDDO
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
  END SUBROUTINE real_to_bloch



END MODULE nbo_shared
