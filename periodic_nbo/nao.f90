MODULE nao
  PRIVATE
  PUBLIC :: do_nao !, do_core_projection

CONTAINS

  !
  !Actually carries out the NAO analysis and outputs the results
  !
  SUBROUTINE do_nao(inp)
    USE matutil
    USE nbo_shared
    IMPLICIT NONE
    
    TYPE(nbo_input) :: inp

    REAL*8,DIMENSION(SIZE(inp%rho0,1),SIZE(inp%rho0,2))  ::  rho_NAO

    INTEGER :: iatom,ifirst,ilast,ibasis,il,ispin
    INTEGER :: nl(0:lmax)
    CHARACTER*10 ::type
    CHARACTER*20 :: label
    REAL*8 :: occ



    rho_NAO=0.d0
    DO ispin=1,inp%nspins
       rho_NAO = rho_NAO + inp%rho0(:,:,1,ispin)
    ENDDO

    !CALL sort_nao(inp)

    PRINT *, '*****************************'
    PRINT *, '**** NATURAL POPULATIONS ****'
    PRINT *, '*****************************'
    PRINT *
    !Output the populations
    DO iatom=1,inp%natom
       ifirst=inp%ibasismap(iatom)
       ilast=inp%ibasismap(iatom+1)-1
              
       occ=mattrace(rho_NAO(ifirst:ilast,ifirst:ilast))
       WRITE(6,'(A,A3,I3,A,F11.7,A,E)')' Atom: ', inp%symbols(iatom), iatom, '  occ:', SNGL(occ), '   charge:', SNGL(inp%iatval(iatom)-occ)
    ENDDO
    PRINT *

    PRINT *, '*********************************'
    PRINT *, '**** NATURAL ATOMIC ORBITALS ****'
    PRINT *, '*********************************'
    PRINT *

    DO iatom=1,inp%natom
       ifirst=inp%ibasismap(iatom)
       ilast=inp%ibasismap(iatom+1)-1
       !Initialize the count of the shells of each angular momentum type
       DO il=0,lmax
          nl(il)=il
       ENDDO

       DO ibasis=ifirst,ilast
          !Create the labels for this shell
          IF (inp%ilmap(ibasis).EQ.0) THEN
             type='s'
             nl(0)=nl(0)+1
             WRITE(label,*) nl(0),'S' !Write to a STRING, cool
             label=ADJUSTL(label)
          ELSEIF (inp%ilmap(ibasis).EQ.1) THEN
             IF (type.EQ.'px') THEN
                type='py'
             ELSEIF (type.EQ.'py') THEN
                type='pz'
             ELSE
                type='px'
                nl(1)=nl(1)+1
                WRITE(label,*) nl(1),'P'
                label=ADJUSTL(label)
             ENDIF
          ELSEIF (inp%ilmap(ibasis).EQ.2) THEN
             IF (type.EQ.'dyz') THEN
                type='dxz'
             ELSEIF (type.EQ.'dxz') THEN
                type='dxy'
             ELSEIF (type.EQ.'dxy') THEN
                type='dx2y2'
             ELSEIF (type.EQ.'dx2y2') THEN
                type='dz2'
             ELSE
                type='dyz'
                nl(2)=nl(2)+1
                WRITE(label,*) nl(2),'D'
                label=ADJUSTL(label)
             ENDIF
          ELSEIF (inp%ilmap(ibasis).EQ.3) THEN
             IF (type.EQ.'fxyz') THEN
                type='fz(x2-y2)'
             ELSEIF (type.EQ.'fz(x2-y2)') THEN
                type='fy(3x2-y2)'
             ELSEIF (type.EQ.'fy(3x2-y2)') THEN
                type='fx(x2-3y2)'
             ELSEIF (type.EQ.'fx(x2-3y2)') THEN
                type='fxz2'
             ELSEIF (type.EQ.'fxz2') THEN
                type='fyz2'
             ELSEIF (type.EQ.'fyz2') THEN
                type='fz3'
             ELSE
                type='fxyz'
                nl(3)=nl(3)+1
                WRITE(label,*) nl(3),'F'
                label=ADJUSTL(label)
             ENDIF

          ENDIF
          !Get the occupation (just the diagonal matrix element of the density matrix,
          !which is now block diagonal by atom and angular momentum)
          occ=rho_NAO(ibasis,ibasis)

          PRINT *, ibasis, inp%symbols(iatom), ' ', type, SNGL(occ), label
       ENDDO

       PRINT *
    ENDDO
    PRINT *

  END SUBROUTINE do_nao



!!$  !
!!$  !Removes contributions from doubly occupied core orbitals from the density matrix
!!$  !
!!$  SUBROUTINE do_core_projection(inp)
!!$    USE nbo_shared
!!$    USE matutil
!!$    IMPLICIT NONE
!!$
!!$    TYPE(nbo_input) :: inp
!!$    INTEGER :: ibasis,ifirst,ilast,iatom
!!$    INTEGER :: ncore
!!$
!!$    REAL*8,DIMENSION(SIZE(inp%rho0,1),SIZE(inp%rho0,1))  ::  P,Pi
!!$
!!$
!!$    !project out any lone pairs from the density (eigenvalues ~ 2 in NAOs)
!!$    P=matiden(SIZE(P,1))
!!$    DO iatom=1,inp%natom
!!$       ifirst=inp%ibasismap(iatom)
!!$       ilast=inp%ibasismap(iatom+1)-1
!!$       ncore=0
!!$       DO ibasis=ifirst,ilast
!!$          IF (inp%rho0(ibasis,ibasis,1).GT.core_thresh) THEN
!!$             !inp%rho0(ibasis,ibasis,1)=0.d0
!!$             ncore=ncore+1
!!$             Pi=matiden(SIZE(Pi,1))
!!$             Pi(ibasis,ibasis)=0.d0
!!$             P=MATMUL(Pi,P)
!!$          ENDIF
!!$       ENDDO
!!$       IF (ncore.GT.0) PRINT *, 'Removed', ncore, 'core orbitals from atom', iatom
!!$    ENDDO
!!$    PRINT *
!!$
!!$    inp%rho0(:,:,1)=MATMUL(P,MATMUL(inp%rho0(:,:,1),P))
!!$
!!$  END SUBROUTINE do_core_projection
!!$
!!$
  !
  !Sorts the NAO by increasing occupancy within a given atom / angular momentum sub block
  !
  SUBROUTINE sort_nao(inp)
    USE sortutil
    USE matutil
    USE nbo_shared
    IMPLICIT NONE

    TYPE(nbo_input) :: inp
    REAL*8 :: occ(inp%nbasis)
    INTEGER :: iatom,ibasis,il,im,ispin
    INTEGER :: ibasisremap(inp%nbasis),isortremap(inp%nbasis)
    LOGICAL :: mask(inp%nbasis)
    INTEGER :: ig,ifirst,ilast,isize,nblock,ndum

    !Initialize the array of occupancies with the NEGATIVE of the occupancy, so we can easily
    !sort in order of decreasing occupancy
    occ = 0.d0
    DO ispin=1,inp%nspins
       DO ibasis=1,inp%nbasis
          occ(ibasis)=occ(ibasis)-inp%rho0(ibasis,ibasis,1,ispin)
       ENDDO
    ENDDO

    DO iatom=1,inp%natom
       ifirst=inp%ibasismap(iatom)
       ilast=inp%ibasismap(iatom+1)-1
       isize=ilast-ifirst+1
       DO il=0,lmax
          DO im=0,2*il
             !Group togoether all the bf of a given angular momentum on a given atom
             mask=.FALSE.
             mask(ifirst:ilast)=inp%ilmap(ifirst:ilast).EQ.il.AND.inp%immap(ifirst:ilast).EQ.im
             CALL partition_basis(mask,inp%nbasis,nblock,ndum,ibasisremap)
             occ=remap_vector(occ,ibasisremap,.FALSE.)
             DO ispin=1,inp%nspins
                DO ig=1,inp%ng
                   inp%rho0(:,:,ig,ispin)=remap_matrix(inp%rho0(:,:,ig,ispin),ibasisremap,.FALSE.)
                   !inp%transform(:,:,ig)=remap_matrix(inp%transform(:,:,ig),ibasisremap,.FALSE.)
                ENDDO
             ENDDO

             !Now sort those bf in order of decreasing occupancy
             DO ibasis=1,inp%nbasis
                isortremap(ibasis)=ibasis
             ENDDO
             CALL quick_sort(occ(1:nblock),isortremap)
             DO ispin=1,inp%nspins
                DO ig=1,inp%ng
                   inp%rho0(:,:,ig,ispin)=remap_matrix(inp%rho0(:,:,ig,ispin),isortremap,.TRUE.)
                   !inp%transform(:,:,ig)=remap_matrix(inp%transform(:,:,ig),isortremap,.TRUE.)
                ENDDO
             ENDDO

             !And re-separate the bf
             occ=remap_vector(occ,ibasisremap,.TRUE.)
             DO ispin=1,inp%nspins
                DO ig=1,inp%ng
                   inp%rho0(:,:,ig,ispin)=remap_matrix(inp%rho0(:,:,ig,ispin),ibasisremap,.TRUE.)
                   !inp%transform(:,:,ig)=remap_matrix(inp%transform(:,:,ig),ibasisremap,.TRUE.)
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO
  END SUBROUTINE sort_nao

END MODULE nao
