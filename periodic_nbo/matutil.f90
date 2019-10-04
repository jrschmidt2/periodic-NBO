MODULE matutil
  ! Some routines require diag and ludcmp
  ! Some of these routines (but not ALL) have been convereted over to use
  ! LAPACK via the MKL LAPACK95 wrappers.  The rest should probably be
  ! converted for efficiency and consistancy.
  !
  ! Compile flags:  -I $(PATH_TO_MKL) -L $(PATH_TO_MKL) -lmkl_core -lmkl_intel_lp64 -lmkl_sequential -lmkl_blas95 -lmkl_lapack95

  INTERFACE mattrace
     MODULE PROCEDURE mattrace_real
     MODULE PROCEDURE mattrace_complex
  END INTERFACE

  INTERFACE matdiag
     MODULE PROCEDURE matdiag_real
     MODULE PROCEDURE matdiag_complex
  END INTERFACE

  INTERFACE matinv
     MODULE PROCEDURE matinv_real
     MODULE PROCEDURE matinv_complex
  END INTERFACE

  INTERFACE matsqrt
     MODULE PROCEDURE matsqrt_real
     MODULE PROCEDURE matsqrt_complex
  END INTERFACE

  INTERFACE matinvsqrt
     MODULE PROCEDURE matinvsqrt_real
     MODULE PROCEDURE matinvsqrt_complex
  END INTERFACE

  INTERFACE matunitary_trans
     MODULE PROCEDURE matunitary_trans_real
     MODULE PROCEDURE matunitary_trans_complex
     MODULE PROCEDURE matunitary_trans_complex2
  END INTERFACE

CONTAINS

  REAL*8 FUNCTION mattrace_real(dmat)
    IMPLICIT NONE
    REAL*8 , DIMENSION(:,:) :: dmat
    INTEGER :: n1, n2, i

    n1 = SIZE(dmat,1)
    n2 = SIZE(dmat,2)
    IF (n1.NE.n2) STOP 'Attempting to trace a non-square matrix'

    mattrace_real=0.d0
    DO i=1, n1
       mattrace_real=mattrace_real+dmat(i,i)
    ENDDO

    RETURN
  END FUNCTION mattrace_real

  COMPLEX*16 FUNCTION mattrace_complex(cmat)
    IMPLICIT NONE
    COMPLEX*16 , DIMENSION(:,:) :: cmat
    INTEGER :: n1, n2, i

    n1 = SIZE(cmat,1)
    n2 = SIZE(cmat,2)
    IF (n1.NE.n2) STOP 'Attempting to trace a non-square matrix'

    mattrace_complex=0.d0
    DO i=1, n1
       mattrace_complex=mattrace_complex+cmat(i,i)
    ENDDO

    RETURN
  END FUNCTION mattrace_complex

  REAL*8 FUNCTION matdet(dmat)
    IMPLICIT NONE
    REAL*8 , DIMENSION(:,:) :: dmat
    REAL*8 , DIMENSION(SIZE(dmat,1),SIZE(dmat,2)) ::  tmp
    INTEGER , DIMENSION(SIZE(dmat,1)) :: index
    INTEGER :: n1, n2, i

    n1 = SIZE(dmat,1)
    n2 = SIZE(dmat,2)
    IF (n1.NE.n2) STOP 'Attempting to take determinant of a non-square matrix'

    tmp = dmat
    CALL ludcmp(tmp,n1,n1,index,matdet)
    DO i = 1, n1
       matdet=matdet*tmp(i,i)
    ENDDO

  END FUNCTION matdet

  FUNCTION matsqrt_real(dmat)
    USE BLAS95
    IMPLICIT NONE
    REAL*8 , DIMENSION(:,:) :: dmat
    REAL*8 , DIMENSION(SIZE(dmat,1),SIZE(dmat,2)) :: matsqrt_real, vecs
    REAL*8 , DIMENSION(SIZE(dmat,1)) :: vals
    REAL*8 , DIMENSION(SIZE(dmat,1),SIZE(dmat,2)) :: dummy
    INTEGER :: n1, n2, i

    n1 = SIZE(dmat,1)
    n2 = SIZE(dmat,2)
    IF (n1.NE.n2) STOP 'Attempting to square root a non-square matrix'

    matsqrt_real = dmat
    CALL matdiag_real(matsqrt_real, vals, vecs)

    !WRITE(6,'(A,10F8.5)')' Eigenvalues in mat_sqrt',vals
    !WRITE(6,*)

    matsqrt_real=0.d0
    FORALL (i=1:n1) matsqrt_real(i,i) = SQRT(abs(vals(i)))
    !matsqrt_real = MATMUL(vecs, MATMUL(matsqrt_real, TRANSPOSE(vecs)))

    CALL DGEMM_MKL95(vecs,matsqrt_real,dummy,'N','N',1.d0,0.d0)
    CALL DGEMM_MKL95(dummy,vecs,matsqrt_real,'N','T',1.d0,0.d0)


  END FUNCTION matsqrt_real

  FUNCTION matsqrt_complex(cmat)
    USE BLAS95
    IMPLICIT NONE
    COMPLEX*16 , DIMENSION(:,:) :: cmat
    COMPLEX*16 , DIMENSION(SIZE(cmat,1),SIZE(cmat,2)) :: matsqrt_complex, vecs
    COMPLEX*16 , DIMENSION(SIZE(cmat,1),SIZE(cmat,2)) :: dummy
    REAL*8 , DIMENSION(SIZE(cmat,1)) :: vals
    INTEGER :: n1, n2, i

    n1 = SIZE(cmat,1)
    n2 = SIZE(cmat,2)
    IF (n1.NE.n2) STOP 'Attempting to square root a non-square matrix'

    matsqrt_complex = cmat
    CALL matdiag_complex(matsqrt_complex, vals, vecs)
    matsqrt_complex=0.d0
    FORALL (i=1:n1) matsqrt_complex(i,i) = SQRT(abs(vals(i)))
    !matsqrt_complex = MATMUL(vecs, MATMUL(matsqrt_complex, TRANSPOSE(CONJG(vecs))))

    CALL ZGEMM_MKL95(vecs,matsqrt_complex,dummy,'N','N',(1.d0,0.d0),(0.d0,0.d0))
    CALL ZGEMM_MKL95(dummy,vecs,matsqrt_complex,'N','C',(1.d0,0.d0),(0.d0,0.d0))

  END FUNCTION matsqrt_complex

  FUNCTION matinvsqrt_real(dmat)
    USE BLAS95
    IMPLICIT NONE
    REAL*8 , DIMENSION(:,:) :: dmat
    REAL*8 , DIMENSION(SIZE(dmat,1),SIZE(dmat,2)) :: matinvsqrt_real, vecs
    REAL*8 , DIMENSION(SIZE(dmat,1),SIZE(dmat,2)) :: dummy
    REAL*8 , DIMENSION(SIZE(dmat,1)) :: vals
    INTEGER :: n1, n2, i

    n1 = SIZE(dmat,1)
    n2 = SIZE(dmat,2)
    IF (n1.NE.n2) STOP 'Attempting to square root a non-square matrix'

    matinvsqrt_real = dmat
    CALL matdiag_real(matinvsqrt_real, vals, vecs)
    matinvsqrt_real=0.d0
    FORALL (i=1:n1) matinvsqrt_real(i,i) = 1.d0/SQRT(abs(vals(i)))
    !matinvsqrt_real = MATMUL(vecs, MATMUL(matinvsqrt_real, TRANSPOSE(vecs)))

    CALL DGEMM_MKL95(vecs,matinvsqrt_real,dummy,'N','N',1.d0,0.d0)
    CALL DGEMM_MKL95(dummy,vecs,matinvsqrt_real,'N','T',1.d0,0.d0)

  END FUNCTION matinvsqrt_real

  FUNCTION matinvsqrt_complex(cmat)
    USE BLAS95
    IMPLICIT NONE
    COMPLEX*16 , DIMENSION(:,:) :: cmat
    COMPLEX*16 , DIMENSION(SIZE(cmat,1),SIZE(cmat,2)) :: matinvsqrt_complex, vecs
    COMPLEX*16 , DIMENSION(SIZE(cmat,1),SIZE(cmat,2)) :: dummy
    REAL*8 , DIMENSION(SIZE(cmat,1)) :: vals
    INTEGER :: n1, n2, i

    n1 = SIZE(cmat,1)
    n2 = SIZE(cmat,2)
    IF (n1.NE.n2) STOP 'Attempting to square root a non-square matrix'

    matinvsqrt_complex = cmat
    CALL matdiag_complex(matinvsqrt_complex, vals, vecs)


    matinvsqrt_complex=0.d0
    FORALL (i=1:n1) matinvsqrt_complex(i,i) = 1.d0/SQRT(abs(vals(i)))

    !DO i=1,n1
    !   IF( vals(i) >= 1.d-4 )THEN
    !      matinvsqrt_complex(i,i) = 1.d0/SQRT(abs(vals(i)))
    !   ELSE
    !      matinvsqrt_complex(i,i) = 0.d0
    !   ENDIF
    !ENDDO

    !matinvsqrt_complex = MATMUL(vecs, MATMUL(matinvsqrt_complex, TRANSPOSE(CONJG(vecs))))

    CALL ZGEMM_MKL95(vecs,matinvsqrt_complex,dummy,'N','N',(1.d0,0.d0),(0.d0,0.d0))
    CALL ZGEMM_MKL95(dummy,vecs,matinvsqrt_complex,'N','C',(1.d0,0.d0),(0.d0,0.d0))

  END FUNCTION matinvsqrt_complex

  FUNCTION matexp(dmat)
    USE BLAS95
    IMPLICIT NONE
    REAL*8 , DIMENSION(:,:) :: dmat
    REAL*8 , DIMENSION(SIZE(dmat,1),SIZE(dmat,2)) :: matexp, vecs
    REAL*8 , DIMENSION(SIZE(dmat,1)) :: vals
    REAL*8 , DIMENSION(SIZE(dmat,1),SIZE(dmat,2)) :: dummy
    INTEGER :: n1, n2, i

    n1 = SIZE(dmat,1)
    n2 = SIZE(dmat,2)
    IF (n1.NE.n2) STOP 'Attempting to exponentiate a non-square matrix'

    matexp = dmat
    CALL diag(matexp, n1, n1, vals, vecs)
    matexp = 0.d0
    FORALL (i=1:n1) matexp(i,i) = EXP(vals(i))
    !matexp = MATMUL(vecs, MATMUL(matexp, TRANSPOSE(vecs)))
    CALL DGEMM_MKL95(vecs,matexp,dummy,'N','N',1.d0,0.d0)
    CALL DGEMM_MKL95(dummy,vecs,matexp,'N','T',1.d0,0.d0)

  END FUNCTION matexp

  ! returns the minor of the matrix, removing row i and column j
  FUNCTION matminor(dmat,i,j)
    IMPLICIT NONE
    REAL*8 , DIMENSION(:,:) :: dmat
    REAL*8 , DIMENSION(SIZE(dmat,1)-1,SIZE(dmat,2)-1) :: matminor
    INTEGER :: n1, n2, i, j

    n1 = SIZE(dmat,1)
    n2 = SIZE(dmat,2)

    matminor(1:i-1,1:j-1)=dmat(1:i-1,1:j-1)
    matminor(1:i-1,j:n2-1)=dmat(1:i-1,j+1:n2)
    matminor(i:n1-1,1:j-1)=dmat(i+1:n1,1:j-1)
    matminor(i:n1-1,j:n2-1)=dmat(i+1:n1,j+1:n2)

  END FUNCTION matminor

  ! solves the generalized eigenvalue equations D c=S c*lambda
  SUBROUTINE diag_general(dmat,smat,dvals,dvecs)
    USE BLAS95
    USE LAPACK95
    IMPLICIT NONE
    REAL*8 , DIMENSION(:,:) :: dmat, smat, dvecs
    REAL*8 , DIMENSION(:)   :: dvals
    REAL*8 , DIMENSION(SIZE(dmat,1),SIZE(dmat,2)) :: shalf, shalfi
    REAL*8 , DIMENSION(SIZE(dmat,1),SIZE(dmat,2)) :: dummy
    INTEGER n1,n2,nrot

    INTEGER  :: j

    n1 = SIZE(dmat,1)
    n2 = SIZE(dmat,2)

    IF (n1.NE.n2) STOP 'Attempting to diagonalize a non-square matrix'

    ! calculate the square root of S, and its inverse
    shalf = matsqrt(smat)
    shalfi = matinv(shalf)
    !call gaussj(shalfi,n1,n1,0,0,0)


    dummy = MATMUL(shalfi,MATMUL(smat,shalfi))
    !WRITE(6,*)'inverse square root test'
    !DO j=1,n1
    !   WRITE(6,'(10F10.5)')dummy(j,:)
    !ENDDO
    !WRITE(6,*)

    ! transform to the coordinates of the normal eigensystem
    !dmat = MATMUL(TRANSPOSE(shalfi),MATMUL(dmat,shalfi))
    CALL DGEMM_MKL95(shalfi,dmat,dummy,'T','N',1.d0,0.d0)
    CALL DGEMM_MKL95(dummy,shalfi,dmat,'N','N',1.d0,0.d0)
    ! get the normal eigenvalues and vectors
    call jacobi(dmat,n1,n1,dvals,dvecs,nrot)
    ! transform the eigenvectors BACK to the original coordinate system
    !dvecs = MATMUL(shalfi, dvecs)
    dummy = dvecs
    CALL DGEMM_MKL95(shalfi,dummy,dvecs,'N','N',1.d0,0.d0)
    CALL eigsrt(dvals,dvecs,n1,n1)

    !WRITE(6,*)'Eigenvalues'
    !WRITE(6,'(5F10.5)')dvals
    !WRITE(6,*)
    
  END SUBROUTINE diag_general

  FUNCTION matinv_real(dmat)
    USE LAPACK95
    IMPLICIT NONE
    REAL*8 , DIMENSION(:,:) :: dmat
    REAL*8 , DIMENSION(SIZE(dmat,1),SIZE(dmat,2)) :: matinv_real
    INTEGER, DIMENSION(SIZE(dmat,1)) :: ipiv
    INTEGER :: n1, n2, info

    n1 = SIZE(dmat,1)
    n2 = SIZE(dmat,2)
    IF (n1.NE.n2) STOP 'Attempting to invert a non-square matrix'

    matinv_real=dmat
    CALL DGETRF_MKL95(matinv_real,ipiv,info)
    CALL DGETRI_MKL95(matinv_real,ipiv,info)
  END FUNCTION matinv_real

  FUNCTION matinv_complex(cmat)
    USE LAPACK95
    IMPLICIT NONE
    COMPLEX*16 , DIMENSION(:,:) :: cmat
    COMPLEX*16 , DIMENSION(SIZE(cmat,1),SIZE(cmat,2)) :: matinv_complex
    INTEGER, DIMENSION(SIZE(cmat,1)) :: ipiv
    INTEGER :: n1, n2, info

    COMPLEX*16 , DIMENSION(SIZE(cmat,1),SIZE(cmat,2)) :: evecs
    REAL*8 , DIMENSION(SIZE(cmat,1))                  :: evals


    n1 = SIZE(cmat,1)
    n2 = SIZE(cmat,2)
    IF (n1.NE.n2) STOP 'Attempting to invert a non-square matrix'


    !evecs = cmat
    !CALL ZHEEV_MKL95(evecs,evals,'V','U',INFO)
    !DO info=1,SIZE(evals,1)
    !   WRITE(6,'(F10.5)')evals(info)
    !ENDDO
    !WRITE(6,*)



    matinv_complex=cmat
    CALL ZGETRF_MKL95(matinv_complex,ipiv,info)
    CALL ZGETRI_MKL95(matinv_complex,ipiv,info)
  END FUNCTION matinv_complex

  FUNCTION outer_product(x,y)
    IMPLICIT NONE
    REAL*8 , DIMENSION(:):: x, y
    REAL*8 :: outer_product(SIZE(x),SIZE(y))

    outer_product=SPREAD(x,2,SIZE(x))*SPREAD(y,1,SIZE(y))
  END FUNCTION outer_product

!Fortran 90 wrappers for diag routines
  SUBROUTINE matdiag_real(dmat,evals,evecs)
    USE LAPACK95
    IMPLICIT NONE
    REAL*8 , DIMENSION(:,:) :: dmat
    REAL*8 , DIMENSION(SIZE(dmat,1)) :: evals
    REAL*8 , DIMENSION(SIZE(dmat,1),SIZE(dmat,2)) :: evecs
    INTEGER :: n1, n2, info

    n1 = SIZE(dmat,1)
    n2 = SIZE(dmat,2)
    IF (n1.NE.n2) STOP 'Attempting to diagonalize a non-square matrix'

    evecs=dmat
    !CALL dsyev_mkl95(evecs,evals,'V','U',info)
    CALL diag(dmat,n1,n1,evals,evecs)

    !WRITE(6,*)'using real symmetrizing routine'

    !IF (info.NE.0) STOP 'Error in matdiag_real'

  END SUBROUTINE matdiag_real

  SUBROUTINE matdiag_complex(cmat,evals,evecs)
    USE LAPACK95
    IMPLICIT NONE
    COMPLEX*16 , DIMENSION(:,:) :: cmat
    REAL*8 , DIMENSION(SIZE(cmat,1)) :: evals
    COMPLEX*16 , DIMENSION(SIZE(cmat,1),SIZE(cmat,2)) :: evecs
    INTEGER :: info

    COMPLEX*16, DIMENSION(SIZE(cmat,1),SIZE(cmat,2))  ::  u_matrix,v_matrix,dummy
    REAL*8, DIMENSION(SIZE(cmat,1))        ::  sing_value

    INTEGER   ::  j

    IF (SIZE(cmat,1).NE.SIZE(cmat,2)) STOP 'Attempting to diagonalize a non-square matrix'

    evecs=cmat

    CALL zheev_mkl95(evecs,evals,'V','U',info)
    IF (info.NE.0) STOP 'Error in matdiag_complex'

    !dummy=cmat
    !CALL ZGESVD_MKL95(dummy,sing_value,u_matrix,v_matrix)


    !IF( evals(1) .LE. -1.d-13 )WRITE(6,*)'WHOOOOOOOOOAAAAAAAAAAAAAAA!!!   Negative eigenvalues'

    !evecs=u_matrix
    !evals=sing_value

    !DO j=1,SIZE(cmat,1)
    !   !WRITE(6,*)'vector',j
    !   IF( ABS(sing_value(SIZE(cmat,1)+1-j) - evals(j)) .GT. 1.d-9 )THEN
    !       !WRITE(6,*)'eig and single value diff',j
    !       !WRITE(6,*)sing_value(SIZE(cmat,1)+1-j),evals(j)
    !   ENDIF
    !   !WRITE(6,*)'eigenvalue     ',evals(j)
    !ENDDO
    !WRITE(6,*)

    !WRITE(6,'(A,2D13.5)')' Smallest eigennvalue found in diag for matinvsqrt',evals(1),sing_value(SIZE(cmat,1))


    !WRITE(6,*)'compare eigenvec with u_vec of matching eigen values'
    !WRITE(6,'(10D13.5)')REAL(evecs(:,1))
    !WRITE(6,*)
    !WRITE(6,'(10D13.5)')REAL(u_matrix(:,SIZE(cmat,1)+1-15))
    !WRITE(6,*)



  END SUBROUTINE matdiag_complex

  FUNCTION matiden(n)
    IMPLICIT NONE
    INTEGER :: n,i
    REAL*8,DIMENSION(n,n) :: matiden
    
    matiden=0.d0
    DO i=1,n
       matiden(i,i)=1.d0
    ENDDO
  END FUNCTION matiden

  FUNCTION matunitary_trans_real(A,U) RESULT(Ap)
    USE BLAS95
    IMPLICIT NONE
    REAL*8,DIMENSION(:,:) :: A,U
    REAL*8,DIMENSION(SIZE(U,1),SIZE(U,1)) :: Ap
    REAL*8,DIMENSION(SIZE(A,1),SIZE(U,1)) :: dummy

    !Ap=MATMUL(U,MATMUL(A,TRANSPOSE(U)))
    CALL DGEMM_MKL95(A,U,dummy,'N','T',1.d0,0.d0)
    CALL DGEMM_MKL95(U,dummy,Ap,'N','N',1.d0,0.d0)
  END FUNCTION matunitary_trans_real

  FUNCTION matunitary_trans_complex(A,U) RESULT(Ap)
    USE BLAS95
    IMPLICIT NONE
    COMPLEX*16,DIMENSION(:,:) :: A,U
    COMPLEX*16,DIMENSION(SIZE(U,1),SIZE(U,1)) :: Ap
    COMPLEX*16,DIMENSION(SIZE(A,1),SIZE(U,1)) :: dummy

    !Ap=MATMUL(U,MATMUL(A,CONJG(TRANSPOSE(U))))
    CALL ZGEMM_MKL95(A,U,dummy,'N','C',(1.d0,0.d0),(0.d0,0.d0))
    CALL ZGEMM_MKL95(U,dummy,Ap,'N','N',(1.d0,0.d0),(0.d0,0.d0))
  END FUNCTION matunitary_trans_complex

  FUNCTION matunitary_trans_complex2(A,U) RESULT(Ap)
    USE BLAS95
    IMPLICIT NONE
    COMPLEX*16,DIMENSION(:,:) :: A
    REAL*8,DIMENSION(:,:) :: U
    COMPLEX*16,DIMENSION(SIZE(U,1),SIZE(U,1)) :: Ap
    COMPLEX*16,DIMENSION(SIZE(U,1),SIZE(U,2)) :: U_comp
    COMPLEX*16,DIMENSION(SIZE(A,1),SIZE(U,1)) :: dummy

    U_comp = U
    !Ap=MATMUL(U,MATMUL(A,TRANSPOSE(U)))
    CALL ZGEMM_MKL95(A,U_comp,dummy,'N','C',(1.d0,0.d0),(0.d0,0.d0))
    CALL ZGEMM_MKL95(U_comp,dummy,Ap,'N','N',(1.d0,0.d0),(0.d0,0.d0))
  END FUNCTION matunitary_trans_complex2

  
END MODULE
