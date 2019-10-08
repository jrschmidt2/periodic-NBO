MODULE periodic_matutil
  IMPLICIT NONE

  PRIVATE ng, indexg,  map_g

  !A module for manipulating "periodic" matrices, that those that can be represented as A_{mu nu}^{g g'}, i.e. block matrices.
  !To be periodic, we require that A_{mu nu}^{g g'} = A_{mu nu}^{0  g'-g}, that is they are translationally invarient.
  !Given this periodicity, we need only store A_{mu nu}^{0 g} for a finite number of g-vectors, assuming that the matrices fall
  !off with distance.

  !The number of blocks in our periodic matrix
  INTEGER :: ng
  !The mapping of the ng matrices into spatial cells [ i.e. (dx,dy,dz) for each of the ng blocks ]
  INTEGER,ALLOCATABLE :: indexg(:,:)
  
CONTAINS

  !Initialize the periodic matrices
  SUBROUTINE periodic_matinit(the_ng,the_indexg)
    INTEGER :: the_ng
    INTEGER,DIMENSION(3,the_ng) :: the_indexg

    ng=the_ng
    IF (ALLOCATED(indexg)) DEALLOCATE(indexg)
    ALLOCATE(indexg(3,ng))
    indexg=the_indexg
  END SUBROUTINE periodic_matinit

  FUNCTION periodic_matiden(n) RESULT(c)
    USE matutil
    IMPLICIT NONE
    REAL*8,DIMENSION(n,n,ng) :: c
    INTEGER :: n, ig

    DO ig=1,ng
       c(:,:,ig)=matiden(n)
    ENDDO
  END FUNCTION periodic_matiden


  !Multiples two periodic matrices.
  !
  FUNCTION periodic_matmul(a,b) RESULT(c)
    USE BLAS95
    IMPLICIT NONE
    REAL*8,DIMENSION(:,:,:), INTENT(in) :: a,b
    REAL*8,DIMENSION(SIZE(a,1),SIZE(a,2),SIZE(a,3)) :: c
    
    INTEGER :: ig, igp, idg
    INTEGER :: idxg, idyg, idzg, idxgp, idygp, idzgp, idx, idy, idz
    INTEGER :: n

    REAL*8 norma(SIZE(a,3)), normb(SIZE(a,3)), normab,thresh
    PARAMETER (thresh=1d-16)

    LOGICAL :: debug
    PARAMETER (debug=.FALSE.)
    INTEGER :: imul, imultot

    n=SIZE(a,1)
    !IF (ng.NE.SIZE(a,3)) STOP 'Inconsistent periodic matrix dimensions in periodic_matmul; did you call periodic_matinit?'

    !Calculate the Frobenius norm of each matrix for screening the matrix products
    DO ig=1,ng
       norma(ig)=SQRT(SUM(a(:,:,ig)**2))
       normb(ig)=SQRT(SUM(b(:,:,ig)**2))
    ENDDO
    
    !(A*B)_{mu nu}^{0 g} = \sum_g'{ A^{0 g'} * B^{g' g} }
    c=0.d0
    imul=0
    imultot=0

    !loop over all subblocks
!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(imul,imultot,a,b,c,norma,normb,ng,indexg)
!$OMP DO SCHEDULE(DYNAMIC,10)
    DO ig=1,ng
       idxg=indexg(1,ig)
       idyg=indexg(2,ig)
       idzg=indexg(3,ig)

       !loop over all possible g' vectors
       DO igp=1,ng
          idxgp=indexg(1,igp)
          idygp=indexg(2,igp)
          idzgp=indexg(3,igp)

          !get the index of the subblock corresponding to g - g'
          idx=idxg-idxgp
          idy=idyg-idygp
          idz=idzg-idzgp
          idg=find_g(idx,idy,idz)

          !Apply the Cauchy-Schwarz type inquality on the matrix norm to
          !screen for negligable contributions
          IF (idg.GT.0) THEN
             normab=norma(igp)*normb(idg)
          ENDIF
!$OMP CRITICAL
          imultot=imultot+1
!$OMP END CRITICAL
          IF (idg.GT.0.AND.normab.GT.thresh) THEN
!             c(:,:,ig)=c(:,:,ig)+MATMUL(a(:,:,igp),b(:,:,idg))
             CALL DGEMM_F95(a(:,:,igp),b(:,:,idg),c(:,:,ig),'N','N',1.d0,1.d0)
!$OMP CRITICAL
             imul=imul+1
!$OMP END CRITICAL
          ENDIF

       ENDDO
         
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    IF (debug) PRINT *, 'Total blocks multiplied: ', imul, 'of', imultot

  END FUNCTION periodic_matmul

  FUNCTION periodic_matvecmul(a,b) RESULT(c)
    USE BLAS95
    IMPLICIT NONE
    REAL*8,DIMENSION(:,:,:), INTENT(in) :: a
    REAL*8,DIMENSION(:,:), INTENT(IN)   :: b
    REAL*8,DIMENSION(SIZE(b,1),SIZE(b,2)) :: c
          
    INTEGER :: ig, igp, idg
    INTEGER :: idxg, idyg, idzg, idxgp, idygp, idzgp, idx, idy, idz
    INTEGER :: n
          
    !REAL*8 norma(SIZE(a,3)), normb(SIZE(a,3)), normab,thresh
    !PARAMETER (thresh=1d-16)

    LOGICAL :: debug
    PARAMETER (debug=.FALSE.)
    INTEGER :: imul, imultot

    !(A*B)_{mu nu}^{0 g} = \sum_g'{ A^{0 g'} * B^{g' g} }
    c=0.d0
    imul=0
    imultot=0

    !loop over all subblocks
    DO ig=1,ng
       idxg=indexg(1,ig)
       idyg=indexg(2,ig)
       idzg=indexg(3,ig)

       !WRITE(6,*)   'output of vector',ig,indexg(:,ig)

       !loop over all possible g' vectors
       DO igp=1,ng
          !WRITE(6,*)'test index     ',igp,indexg(:,igp)
          idxgp=indexg(1,igp)
          idygp=indexg(2,igp)
          idzgp=indexg(3,igp)

          !get the index of the subblock corresponding to g - g'
          idx=idxg-idxgp
          idy=idyg-idygp
          idz=idzg-idzgp
          idg=find_g(idx,idy,idz)

          imultot=imultot+1
          IF (idg.GT.0) THEN
             !WRITE(6,*)'non zero matrix',idg,indexg(:,idg)
             c(:,ig)=c(:,ig)+MATMUL(a(:,:,idg),b(:,igp))
             !CALL DGEMM_F95(a(:,:,igp),b(:,:,idg),c(:,:,ig),'N','N',1.d0,1.d0)
             imul=imul+1
          ENDIF

          !WRITE(6,*)

       ENDDO

    ENDDO

  END FUNCTION periodic_matvecmul



  !
  !Takes the transpose of a periodic matrix.
  FUNCTION periodic_transpose(a) RESULT(c)
    REAL*8,DIMENSION(:,:,:), INTENT(in) :: a
    REAL*8,DIMENSION(SIZE(a,1),SIZE(a,2),SIZE(a,3)) :: c
    INTEGER :: ig, igp

    IF (ng.NE.SIZE(a,3)) STOP 'Inconsistent periodic matrix dimensions in periodic_transpose; did you call periodic_matinit?'

    DO ig=1,ng
       igp=find_g(-indexg(1,ig),-indexg(2,ig),-indexg(3,ig))
       IF (igp.LE.0) STOP 'Malformed periodic matrix in periodic_transpose'
       c(:,:,ig)=TRANSPOSE(a(:,:,igp))
    ENDDO

  END FUNCTION periodic_transpose


!  !Takes the square root of a periodic matrix.  We accomplish this by:
!  !1) First transforming to the Block basis.
!  !2) In the Bloch basis, the matrix becomes block diagonal
!  !3) Take the square root of each (finite) block
!  !4) Transforming back from the Block basis
!  FUNCTION periodic_matsqrt(a)  RESULT(c)
!    USE matutil
!    IMPLICIT NONE
!    REAL*8,DIMENSION(:,:,:), INTENT(in) :: a
!    REAL*8 ,DIMENSION(SIZE(a,1),SIZE(a,2),SIZE(a,3)):: c
!    COMPLEX*16 , ALLOCATABLE, DIMENSION(:,:,:):: ak
!    INTEGER :: nkx, nky, nkz, nk, ik
!
!    CALL get_nk(nkx,nky,nkz)
!    nk=nkx*nky*nkz
!
!    ALLOCATE(ak(SIZE(a,1),SIZE(a,2),nk))
!
!    Ak=periodic_matbloch(A,nkx,nky,nkz)
!    
!    !Now take the square root of each of the Sk
!    !OpenMP causes problems here since the MKL LAPACK95 interface is not thread safe; I recompiled a thread safe version
!    !by adding -openmp to the makefile
!
!!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(Ak,nk)
!!$OMP DO SCHEDULE(DYNAMIC,10)
!    DO ik=1,nk
!       Ak(:,:,ik)=matsqrt(Ak(:,:,ik))
!    ENDDO
!!$OMP END DO NOWAIT
!!$OMP END PARALLEL
!
!    C=periodic_matinvbloch(Ak,nkx,nky,nkz)
!
!  END FUNCTION periodic_matsqrt
!
!
!  FUNCTION periodic_matinvsqrt(a)  RESULT(c)
!    USE matutil
!    IMPLICIT NONE
!    REAL*8,DIMENSION(:,:,:), INTENT(in) :: a
!    REAL*8 ,DIMENSION(SIZE(a,1),SIZE(a,2),SIZE(a,3)):: c
!    COMPLEX*16 , ALLOCATABLE, DIMENSION(:,:,:):: ak
!    INTEGER :: nkx, nky, nkz, nk, ik
!
!    CALL get_nk(nkx,nky,nkz)
!    nk=nkx*nky*nkz
!
!    ALLOCATE(ak(SIZE(a,1),SIZE(a,2),nk))
!
!    Ak=periodic_matbloch(A,nkx,nky,nkz)
!    
!    !Now take the square root of each of the Sk
!    !OpenMP causes problems here since the MKL LAPACK95 interface is not thread safe; I recompiled a thread safe version
!    !by adding -openmp to the makefile
!
!!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(Ak,nk)
!!$OMP DO SCHEDULE(DYNAMIC,10)
!    DO ik=1,nk
!       Ak(:,:,ik)=matinvsqrt(Ak(:,:,ik))
!    ENDDO
!!$OMP END DO NOWAIT
!!$OMP END PARALLEL
!
!    C=periodic_matinvbloch(Ak,nkx,nky,nkz)
!
!  END FUNCTION periodic_matinvsqrt

  !
  !Executes a transformation between a periodic matrix A in a localized basis and the same matrix
  !in the corresponding Bloch basis
!!$  FUNCTION periodic_matbloch(a,nkx,nky,nkz) RESULT(ak)
!!$    IMPLICIT NONE
!!$    REAL*8,DIMENSION(:,:,:), INTENT(in) :: a
!!$    COMPLEX*16 ,DIMENSION(SIZE(a,1),SIZE(a,2),nkx*nky*nkz):: ak
!!$
!!$    INTEGER :: nkx, nky, nkz, ikx, iky, ikz, ig, ik
!!$    REAL*8 :: kx, ky, kz, pi, arg
!!$    COMPLEX*16 :: sqrt_minus_one
!!$
!!$    !nk=ng
!!$    sqrt_minus_one=(0.d0, 1.d0)
!!$    pi=4.d0*ATAN(1.d0)
!!$
!!$    !Convert from A_{mu nu}^{0 g} to A_{mu _nu}^(k) via Fourier transform
!!$    Ak=0.d0
!!$    ik=0
!!$    DO ikx=1,nkx
!!$       kx=2.d0*pi*(ikx-1)/DBLE(nkx)
!!$       DO iky=1,nky
!!$          ky=2.d0*pi*(iky-1)/DBLE(nky)
!!$          DO ikz=1,nkz
!!$             kz=2.d0*pi*(ikz-1)/DBLE(nkz)
!!$             ik=ik+1
!!$
!!$             DO ig=1,ng
!!$                arg=kx*indexg(1,ig)+ky*indexg(2,ig)+kz*indexg(3,ig)
!!$                Ak(:,:,ik)=Ak(:,:,ik)+A(:,:,ig)*EXP(sqrt_minus_one*arg)
!!$             ENDDO
!!$
!!$          ENDDO
!!$       ENDDO
!!$    ENDDO
!!$
!!$  END FUNCTION periodic_matbloch

!  !aktemp and aktemp_3d really should be 'EQUIVALENCE'-ed to reduce memory
!  !usage, but this is not allowed since they are automatic arrays.  If
!  !reducing memory usage is crucial, aktemp and aktemp_3d could be passed
!  !in as dummy arguments, but pointing to the same array.  This would require
!  !a bit of work to not make use of 'SIZE' in the definition of aktemp.
!  FUNCTION periodic_matbloch(a,nkx,nky,nkz) RESULT(ak)
!    USE MKL_DFTI
!    IMPLICIT NONE
!    REAL*8,DIMENSION(:,:,:), INTENT(in) :: a
!    COMPLEX*16, DIMENSION(nkx*nky*nkz*SIZE(a,1)*SIZE(a,2)) :: aktemp
!    COMPLEX*16, DIMENSION(nkx,nky,nkz,SIZE(a,1),SIZE(a,2)) :: aktemp_3d
!    COMPLEX*16 ,DIMENSION(SIZE(a,1),SIZE(a,2),nkx*nky*nkz) :: ak
!
!    INTEGER :: nkx, nky, nkz, ig, igx, igy, igz
!    TYPE(DFTI_DESCRIPTOR), POINTER :: dfti_desc
!    INTEGER :: status, len(3)
!
!    !Convert from A_{mu nu}^{0 g} to A_{mu _nu}^(k) via Fourier transform
!
!    !Rearrange the data to make it conducive to a FFT
!    Aktemp_3d=0.d0
!    DO ig=1,ng
!       CALL map_g(ig,nkx,nky,nkz,igx,igy,igz)
!       Aktemp_3d(igx,igy,igz,:,:)=a(:,:,ig)
!    ENDDO
!    
!    Aktemp=RESHAPE(Aktemp_3d, (/SIZE(a,1)*SIZE(a,2)*nkx*nky*nkz/) )
!
!    !Do nbasis**2 3D FFTs of the data to convert the matrix to its Bloch representation
!    len=(/nkx,nky,nkz/)
!    status=DftiCreateDescriptor(dfti_desc, DFTI_DOUBLE, DFTI_COMPLEX, 3, len)
!    status=DftiSetValue(dfti_desc, DFTI_NUMBER_OF_TRANSFORMS, SIZE(a,1)*SIZE(a,2))
!    status=DftiSetValue(dfti_desc, DFTI_INPUT_DISTANCE, nkx*nky*nkz)
!    status=DftiSetValue(dfti_desc, DFTI_OUTPUT_DISTANCE, nkx*nky*nkz)
!    status=DftiCommitDescriptor(dfti_desc)
!    status=DftiComputeForward(dfti_desc, Aktemp)
!    status=DftiFreeDescriptor(dfti_desc)
!
!    !Re-rearrange the data to put it pack in a reasonable format
!    Ak=RESHAPE(Aktemp, (/SIZE(a,1), SIZE(a,2), nkx*nky*nkz/), ORDER=(/3,1,2/) )  
!
!  END FUNCTION periodic_matbloch


  !
  !Executes a transformation between a periodic matrix A in a Bloch basis and the same matrix
  !in the corresponding localized basis
!!$  FUNCTION periodic_matinvbloch(ak,nkx,nky,nkz) RESULT(a)
!!$    IMPLICIT NONE
!!$    COMPLEX*16,DIMENSION(:,:,:), INTENT(in) :: ak
!!$    REAL*8 ,DIMENSION(SIZE(ak,1),SIZE(ak,2),ng):: a
!!$
!!$    INTEGER :: nkx, nky, nkz, ikx, iky, ikz, ig, ik
!!$    REAL*8 :: kx, ky, kz, pi, arg
!!$    COMPLEX*16 :: sqrt_minus_one
!!$
!!$    sqrt_minus_one=(0.d0, 1.d0)
!!$    pi=4.d0*ATAN(1.d0)
!!$
!!$    !Convert from A_{mu nu}^{0 g} to A_{mu _nu}^(k) via Fourier transform
!!$    A=0.d0
!!$    ik=0
!!$    DO ikx=1,nkx
!!$       kx=2.d0*pi*(ikx-1)/DBLE(nkx)
!!$       DO iky=1,nky
!!$          ky=2.d0*pi*(iky-1)/DBLE(nky)
!!$          DO ikz=1,nkz
!!$             kz=2.d0*pi*(ikz-1)/DBLE(nkz)
!!$             ik=ik+1
!!$
!!$             DO ig=1,ng
!!$                arg=kx*indexg(1,ig)+ky*indexg(2,ig)+kz*indexg(3,ig)
!!$                A(:,:,ig)=A(:,:,ig)+Ak(:,:,ik)*EXP(-sqrt_minus_one*arg)
!!$             ENDDO
!!$
!!$          ENDDO
!!$       ENDDO
!!$    ENDDO
!!$    A=A/(nkx*nky*nkz)
!!$
!!$  END FUNCTION periodic_matinvbloch

!  FUNCTION periodic_matinvbloch(ak,nkx,nky,nkz) RESULT(a)
!    USE MKL_DFTI
!    IMPLICIT NONE
!    COMPLEX*16,DIMENSION(:,:,:), INTENT(in) :: ak
!    COMPLEX*16, DIMENSION(nkx,nky,nkz,SIZE(ak,1),SIZE(ak,2)) :: atemp_3d
!    COMPLEX*16, DIMENSION(nkx*nky*nkz*SIZE(ak,1)*SIZE(ak,2)) :: atemp
!    REAL*8 ,DIMENSION(SIZE(ak,1),SIZE(ak,2),nkx*nky*nkz) :: a
!
!    INTEGER :: nkx, nky, nkz, ig, igx, igy, igz
!    TYPE(DFTI_DESCRIPTOR), POINTER :: dfti_desc
!    INTEGER :: status, len(3)
!
!    !Convert from A_{mu nu}^(k) to A_{mu _nu}^{o g} via inverse Fourier transform
!    Atemp_3d=RESHAPE(ak, (/nkx,nky,nkz,SIZE(ak,1),SIZE(ak,2)/), ORDER=(/4,5,1,2,3/) )
!    Atemp=RESHAPE(Atemp_3d, (/nkx*nky*nkz*SIZE(ak,1)*SIZE(ak,2)/) )
!
!    !Do nbasis**2 3D FFTs of the data to convert the matrix to its Bloch representation
!    len=(/nkx,nky,nkz/)
!    status=DftiCreateDescriptor(dfti_desc, DFTI_DOUBLE, DFTI_COMPLEX, 3, len)
!    status=DftiSetValue(dfti_desc, DFTI_NUMBER_OF_TRANSFORMS, SIZE(ak,1)*SIZE(ak,2))
!    status=DftiSetValue(dfti_desc, DFTI_INPUT_DISTANCE, nkx*nky*nkz)
!    status=DftiSetValue(dfti_desc, DFTI_OUTPUT_DISTANCE, nkx*nky*nkz)
!    status=DftiCommitDescriptor(dfti_desc)
!    status=DftiComputeBackward(dfti_desc, Atemp)
!    status=DftiFreeDescriptor(dfti_desc)
!
!    !Re-rearrange the data to put it pack in a reasonable format
!    Atemp_3d=RESHAPE(Atemp, (/nkx,nky,nkz,SIZE(a,1),SIZE(a,2)/) )
!
!    DO ig=1,ng
!       CALL map_g(ig,nkx,nky,nkz,igx,igy,igz)
!       a(:,:,ig)=Atemp_3d(igx,igy,igz,:,:)
!    ENDDO
!
!    a=a/(nkx*nky*nkz)
!    
!  END FUNCTION periodic_matinvbloch


  !
  !Searches through the index to find the appropriate g vector
  INTEGER FUNCTION find_g(idx,idy,idz)
    IMPLICIT NONE
    INTEGER :: idx,idy,idz
    INTEGER :: j
  
    find_g = -1
    DO j=1,ng
       IF (idx.EQ.indexg(1,j).AND.idy.EQ.indexg(2,j).AND.idz.EQ.indexg(3,j)) find_g = j
    ENDDO
  END FUNCTION find_g

  !
  !Maps a g-vector into a location appropriate for Fouier transform (i.e. negative g vectors wrap around)
  SUBROUTINE map_g(ig,nkx,nky,nkz,igx,igy,igz)
    IMPLICIT NONE
    INTEGER :: ig,igx,igy,igz,nkx,nky,nkz

    igx=indexg(1,ig)
    IF (igx.lt.0) igx=nkx+igx
    igx=igx+1

    igy=indexg(2,ig)
    IF (igy.lt.0) igy=nky+igy
    igy=igy+1

    igz=indexg(3,ig)
    IF (igz.lt.0) igz=nkz+igz
    igz=igz+1  
  END SUBROUTINE map_g

  !
  !Gets the appropiate k-space size for the periodic matrix at hand
  SUBROUTINE get_nk(nkx,nky,nkz)
    IMPLICIT NONE
    INTEGER :: nkx,nky,nkz
    nkx=MAXVAL(indexg(1,:))-MINVAL(indexg(1,:))+1
    nky=MAXVAL(indexg(2,:))-MINVAL(indexg(2,:))+1
    nkz=MAXVAL(indexg(3,:))-MINVAL(indexg(3,:))+1
  END SUBROUTINE get_nk


END MODULE periodic_matutil
