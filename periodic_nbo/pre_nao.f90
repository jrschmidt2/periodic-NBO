MODULE pre_nao
  USE BLAS95
  PRIVATE
  PUBLIC :: do_pre_nao

  !The thresholds to be considered a "significant" Rydberg orbital
  !Note: the original NBO paper states that w_thresh=1d-4, but the Gaussian
  !implementation uses 1d-2 (from the code!)
  REAL*8 :: w_thresh
  PARAMETER (w_thresh=1d-2)
 
  CONTAINS

  !
  !See the appendix of JCP 83 735 (1084) for details
  !
  SUBROUTINE do_pre_nao(inp)
    USE matutil
    USE periodic_matutil
    USE nbo_shared
    USE LAPACK95
    IMPLICIT NONE

    TYPE(nbo_input) :: inp

    !"W" is the weights of the pre-NAO transformation (the occupancies of the pre-NAOs)
    REAL*8, DIMENSION(SIZE(inp%rho0,1)) :: W    
    !sk,rhok are the Bloch space overlap and density matrix
    !between the pre-NAO and NAO basis
    COMPLEX*16, DIMENSION(SIZE(inp%s0,1),SIZE(inp%s0,2),inp%nk) ::  sk, rhok
    COMPLEX*16, DIMENSION(SIZE(inp%s0,1),SIZE(inp%s0,2),inp%nk) ::  ktransform
    COMPLEX*16, DIMENSION(SIZE(inp%s0,1),SIZE(inp%s0,2))        ::  energy_test,elec_test
    REAL*8, DIMENSION(SIZE(inp%s0,1),SIZE(inp%s0,2))            ::  real_energy_test
    COMPLEX*16                                                  ::  energy_sum,num_elec

    COMPLEX*16, DIMENSION(SIZE(inp%s0,1),SIZE(inp%s0,2),inp%nk) ::  svd_overlap
    COMPLEX*16, DIMENSION(SIZE(inp%s0,1),SIZE(inp%s0,2))        ::  u_matrix,v_matrix
    REAL*8, DIMENSION(SIZE(inp%s0,1))                           ::  sing_value
    COMPLEX*16, DIMENSION(SIZE(inp%s0,1),SIZE(inp%s0,2))        ::  inv_dummy
    INTEGER, DIMENSION(SIZE(inp%s0,1)) :: ipiv
    INTEGER :: info

    COMPLEX*16, DIMENSION(SIZE(inp%s0,1),SIZE(inp%s0,2))        ::  evecs
    REAL*8, DIMENSION(SIZE(inp%s0,1))                           ::  evals
    COMPLEX*16,DIMENSION(SIZE(inp%s0,1),SIZE(inp%s0,1))         ::  P,Pj,Pdummy


    INTEGER :: nk,ik,ig,j,k
    INTEGER :: ispin,nspins

    COMPLEX*16                 ::  arg
    REAL*8,PARAMETER           ::  pi=4.d0*ATAN(1.d0)
    COMPLEX*16,PARAMETER       ::  sqrt_neg_one=(0.d0,1.d0)

    REAL*8  :: t1,t2,t3,t4

    PRINT *, 'Obtaining NAO basis...'

    !Since all matrices are already in k-space, the only initialization necessary is of the transform matrix, which starts out simply as a identity matrix
    nk = inp%nk
    !WRITE(6,*)'number of k points used in NBO analysis',nk
    sk = inp%sk
    DO ik=1,nk
       ktransform(:,:,ik) = matiden(SIZE(inp%s0,1))
    ENDDO
    nspins = inp%nspins

    CALL CPU_TIME(t1)
    !Step 1:  Get the REAL density matrix (multiply both sides by overlap matrix)
    rhok = (0.d0,0.d0)
    !This step is done to all spin density matrices independently
    !Then they are combined into an averaged spin density matrix on which the NAO transformation is performed
    DO ispin=1,nspins
       CALL pre_nao_step_1(inp,inp%rhok(:,:,:,ispin),sk,nk)
       rhok = rhok + inp%rhok(:,:,:,ispin)
    ENDDO
    CALL CPU_TIME(t2)
    WRITE(6,*)'step 1 time', SNGL(t2-t1)


    !The remainder of the NAO procedure is done spin independently 
    !The spin averaged density in (rhok) is used.  Thus the same NAOs are generated for all spin types

    !Step 2:  Form the pre-NAOs, returns the pre-NAO occupancies/weights in "w"
    !This step requires the real space g=(/0,0,0/) overlap and density matrices.
    CALL real_space_matrix(inp%rho0(:,:,1,1),rhok,nk,inp)
    CALL real_space_matrix(inp%s0(:,:,1),sk,nk,inp)
    CALL pre_nao_step_2(inp,rhok,sk,ktransform,nk,w)
    CALL CPU_TIME(t1)
    WRITE(6,*)'step 2 time',SNGL(t1-t2)


    !Step 3:  Orthoganalization of the valence and Rydberg states to each other
    CALL pre_nao_step_3(inp,rhok,sk,ktransform,nk,w)
    CALL CPU_TIME(t4)
    WRITE(6,*)'step 3 time',SNGL(t4-t1)


    !Step 4a:  Do an occupancy weighted symmetric orthogonalization on the entire NMB + NRB set
    CALL pre_nao_step_4a(inp,rhok,sk,ktransform,nk,w)
    CALL CPU_TIME(t2)
    WRITE(6,*)'step 4a time',SNGL(t2-t4)


    !Step 4b:  Now go to the NAO basis, given as the eigenfunctions of each Alm block
    inp%s0 = 0.d0
    inp%s0(:,:,1) = matiden(inp%nbasis)
    !CALL real_space_matrix(inp%s0(:,:,1),sk,nk,inp)
    CALL real_space_matrix(inp%rho0(:,:,1,1),rhok,nk,inp)
    CALL pre_nao_step_2(inp,rhok,sk,ktransform,nk,w)
    CALL CPU_TIME(t1)
    WRITE(6,*)'step 4b time',SNGL(t1-t2)
    WRITE(6,*)

    WRITE(6,*)'total time for obtaining NAO basis',SNGL(t1-t3)
    WRITE(6,*)


    !Convert spin dependent density and Fock matrices from AO basis set to NAO
    !The transform matrix as it stands does this
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ik,ispin)
!$OMP DO SCHEDULE(STATIC)
    DO ik=1,nk
       DO ispin=1,nspins
          inp%rhok(:,:,ik,ispin) = matunitary_trans(inp%rhok(:,:,ik,ispin),ktransform(:,:,ik))
          inp%fockk(:,:,ik,ispin) = matunitary_trans(inp%fockk(:,:,ik,ispin),ktransform(:,:,ik))
       ENDDO
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL


    !The observables will also be tested (number of electrons, total energy) to make sure the transformation has been correct
    energy_sum = 0.d0
    num_elec = 0.d0
    DO ispin=1,nspins
       DO ik=1,nk
          CALL ZGEMM_MKL95(inp%rhok(:,:,ik,ispin),inp%fockk(:,:,ik,ispin),energy_test,'N','C',(1.d0,0.d0),(0.d0,0.d0))
          energy_sum = energy_sum + inp%kpt_wt(ik)*mattrace(energy_test)
          !CALL ZGEMM_MKL95(inp%rhok(:,:,ik,ispin),sk(:,:,ik),elec_test,'N','C',(1.d0,0.d0),(0.d0,0.d0))
          num_elec = num_elec + inp%kpt_wt(ik)*mattrace(inp%rhok(:,:,ik,ispin))  !sk should be the identity matrix now
       ENDDO
    ENDDO
    WRITE(6,*)'Energy sum in Bloch space, NAO:       ',REAL(energy_sum)
    WRITE(6,*)'Number of elec. in Bloch space, NAO: ',REAL(num_elec)


    !Now convert the transformation matrix into one that will transform the NAO basis back to the AO basis, Transpose of its current form
    DO ik=1,nk
       ktransform(:,:,ik) = TRANSPOSE(ktransform(:,:,ik))
    ENDDO

    !Finally go back to real space; note that s0 is now, by definition, the identity matrix and was already set up in step 4b
    IF( .NOT. real_init )THEN
       DO ispin=1,nspins
          CALL bloch_to_real(inp,inp%rhok(:,:,:,ispin),inp%rho0(:,:,:,ispin),inp%kpt,inp%indexg)
          CALL bloch_to_real(inp,inp%fockk(:,:,:,ispin),inp%fock0(:,:,:,ispin),inp%kpt,inp%indexg)
       ENDDO
       CALL bloch_to_real(inp,ktransform,inp%transform,inp%kpt,inp%indexg)
       !CALL bloch_to_real(inp,sk,inp%s0,inp%kpt,inp%indexg)
    ELSE
       STOP 'Need to setup taking NAO matrices back to real space for real initial calculations'
       !inp%rho0 = periodic_matinvbloch(rhok,nkx,nky,nkz)
       !!inp%s0 = periodic_matinvbloch(sk,nkx,nky,nkz)
       !inp%transform = periodic_matinvbloch(ktransform,nkx,nky,nkz)
       !inp%fock0 = periodic_matinvbloch(fockk,nkx,nky,nkz)
    ENDIF

    !Only check that the number of electrons in the real space NAO is correct
    !Energy would require real space periodic matrix multiplication, which is not worth it
    num_elec = 0.d0
    DO ispin=1,nspins
       num_elec = num_elec + mattrace(inp%rho0(:,:,1,ispin))
    ENDDO
    PRINT *, 'Number of elec. in real space, NAO:  ', REAL(num_elec)
    PRINT *  

    END SUBROUTINE do_pre_nao

    !***********************************************************
    !Subroutines for each of the sub-steps in the transformation
    !***********************************************************
    SUBROUTINE pre_nao_step_1(inp,rhok,sk,nk)
      USE nbo_shared
      IMPLICIT NONE
      COMPLEX*16, DIMENSION(:,:,:) :: rhok,sk
      TYPE(nbo_input) :: inp
      COMPLEX*16, DIMENSION(SIZE(inp%rho0,1),SIZE(inp%rho0,2))   :: dummy
      INTEGER :: ik,nk

      !Step 1:  Get the REAL density matrix
      !Transform the "bond order matrix" D to the density matrix rho  (rho = S*D*S)
!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(nk,rhok,sk)
!$OMP DO SCHEDULE(STATIC)
      DO ik=1,nk
         CALL ZGEMM_MKL95(rhok(:,:,ik),sk(:,:,ik),dummy,'N','N',(1.d0,0.d0),(0.d0,0.d0))
         CALL zGEMM_MKL95(sk(:,:,ik),dummy,rhok(:,:,ik),'N','N',(1.d0,0.d0),(0.d0,0.d0))
      ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    END SUBROUTINE pre_nao_step_1


    SUBROUTINE pre_nao_step_2(inp,rhok,sk,idenk,nk,w)
      USE periodic_matutil
      USE matutil
      USE nbo_shared
      IMPLICIT NONE

      TYPE(nbo_input) :: inp
      COMPLEX*16, DIMENSION(:,:,:) :: rhok,sk,idenk
      INTEGER :: ik,nk,im,nm
      INTEGER :: iatom,il,ishell,jshell,nshell,ifirst,jfirst
      INTEGER, DIMENSION(SIZE(inp%rho0,1)) :: ibasisremap

      !"N" is the pre-NAO transformation, "W" the associated weightes / occupancies
      COMPLEX*16, DIMENSION(SIZE(rhok,1),SIZE(rhok,2)) :: N  !N is actually filled with real valued quantities, it is made complex for use in BLAS subroutines
      REAL*8,DIMENSION(:) :: W
      !Blocks of the density,overlap,and pre-NAO transformation matrices for atom/angular momentum blocks
      REAL*8, DIMENSION(SIZE(inp%rho0,1),SIZE(inp%rho0,2)) :: rhoat,sat,Nat
      REAL*8, DIMENSION(SIZE(inp%rho0,1)) :: Wat

      COMPLEX*16, DIMENSION(SIZE(inp%rho0,1),SIZE(inp%rho0,2)) :: dummy  ! For use in matrix multiplication for idenk


      N=0.d0
      W=0.d0
      DO iatom=1,inp%natom
         DO il=0,lmax

          !Step 2a,b:  Transformation from Cartesian to spherical (omitted) and symmetry averaging
          CALL load_average(inp,iatom,il,rhoat,sat,nshell,ibasisremap)

          IF (nshell.GT.0) THEN

             !sat(nshell,:) = 0.d0
             !sat(:,nshell) = 0.d0
             !sat(nshell,nshell) = 1.d0
             !WRITE(6,*)'load averaged overlap matrix for iatom,il',iatom,il
             !DO ishell=1,nshell
             !   WRITE(6,'(10F10.5)')sat(ishell,1:nshell)
             !ENDDO

             !rhoat(nshell,:) = 0.d0
             !rhoat(:,nshell) = 0.d0
             !rhoat(nshell,nshell) = 1.d-8
             !WRITE(6,*)'load averaged density matrix for iatom,il',iatom,il
             !DO ishell=1,nshell
             !   WRITE(6,'(10F10.5)')rhoat(ishell,1:nshell)
             !ENDDO
             !WRITE(6,*)


             !Step 2c:  Calculate the pre-NAOs and their occupancies
             CALL diag_general(rhoat(1:nshell,1:nshell),sat(1:nshell,1:nshell),Wat(1:nshell),Nat(1:nshell,1:nshell))
             !CALL sort_pre_nao(Nat(1:nshell,1:nshell),Wat(1:nshell),sat(1:nshell,1:nshell))

             !WRITE(6,*)'Pre-NAO eigenvectors for iatom.il',iatom,il
             !DO ishell=1,nshell
             !   WRITE(6,*)Wat(ishell)
             !   WRITE(6,'(10F10.5)')Nat(1:nshell,ishell)
             !ENDDO
             !WRITE(6,*)

             !Now construct the pre-NAO transformation matrix "N" and the associated weights / occ. "W"
             !using the spherically averaged occupancy / transformation from each shell
             nm=2*il+1
             DO ishell=1,nshell
                DO jshell=1,nshell
                   ifirst=ibasisremap(ishell)
                   jfirst=ibasisremap(jshell)
                   DO im=1,nm
                      W(ifirst+im-1)=Wat(ishell)
                      N(ifirst+im-1,jfirst+im-1)=Nat(ishell,jshell)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
             
         ENDDO
      ENDDO

      !WRITE(6,*)'W at end of step 2'
      !OPEN(23,file='W_test')
      !WRITE(6,*)W
      !CLOSE(23)




      !Transform the density and overlap matrix to the pre-NAO basis
      !Do this in the Bloch basis to prevent numerical error from continued transofmations between Bloch and
      !localized bases
!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(nk,rhok,sk,idenk,N)
!$OMP DO SCHEDULE(STATIC)
      DO ik=1,nk
         sk(:,:,ik)=matunitary_trans(sk(:,:,ik),TRANSPOSE(N))
         rhok(:,:,ik)=matunitary_trans(rhok(:,:,ik),TRANSPOSE(N))
         dummy = idenk(:,:,ik)
         CALL ZGEMM_MKL95(N,dummy,idenk(:,:,ik),'C','N',(1.d0,0.d0),(0.d0,0.d0))
      ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

       !WRITE(6,*)'out of step 2'

    END SUBROUTINE pre_nao_step_2

    SUBROUTINE pre_nao_step_3(inp,rhok,sk,idenk,nk,w)
      USE matutil
      USE nbo_shared
      IMPLICIT NONE

      TYPE(nbo_input) :: inp
      COMPLEX*16, DIMENSION(:,:,:) :: rhok,sk,idenk
      REAL*8,DIMENSION(:) :: W
      INTEGER :: ik,nk,j
      
      INTEGER :: nvalence,nrydberg,nbasis
      INTEGER, DIMENSION(SIZE(inp%rho0,1)) :: ibasisremap
      !"Os" is the Schmidt orthogonalization matrix
      COMPLEX*16, DIMENSION(SIZE(inp%rho0,1),SIZE(inp%rho0,2)) :: Os
      !"Ow" is the population weighted symmetric orthogonalization matrix, with weights W
      COMPLEX*16, DIMENSION(SIZE(inp%rho0,1),SIZE(inp%rho0,2)) :: Ow
      COMPLEX*16, DIMENSION(SIZE(inp%rho0,1),SIZE(inp%rho0,2)) :: dummy
      LOGICAL, DIMENSION(SIZE(inp%rho0,1)) :: nmb_mask


      !WRITE(6,*)'mask for step 3'

      !Get a mask of which basis functions should be included in the NMB
      CALL create_nmb_mask(inp,w,nmb_mask)
 
      !Step 3a,b:  Schmidt orthogonalize the weakly occupied (Rydberg) orbitals to the strongly occupied
      !(natural minimal basis).
      CALL partition_basis(nmb_mask,inp%nbasis,nvalence,nrydberg,ibasisremap)


      w=remap_vector(w,ibasisremap,.FALSE.)

      !WRITE(6,*)'W vector after remapping'
      !WRITE(6,*)W

      nbasis = inp%nbasis

      !WRITE(6,*)'nvalence, nrydberg, inp%nbasis from step 3',nvalence,nrydberg,nbasis

      !WRITE(6,*)'some gram-schmidt orthogonalizations'

!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(nk,rhok,sk,idenk,nvalence,ibasisremap,nbasis,w)
!$OMP DO SCHEDULE(STATIC)

      DO ik=1,nk

         !WRITE(6,*)'kpt',ik

         !First re-arrange the matrices to separate out the valence and Rydberg states
         sk(:,:,ik)=remap_matrix(sk(:,:,ik),ibasisremap,.FALSE.)
         rhok(:,:,ik)=remap_matrix(rhok(:,:,ik),ibasisremap,.FALSE.)
         idenk(:,:,ik)=remap_matrix(idenk(:,:,ik),ibasisremap,.FALSE.)

         !Now do the actual Gram-Schmidt orthogonalization of the Rydberg to the valence states
         Os=schmidt_orthog_matrix(sk(:,:,ik),nvalence)
         sk(:,:,ik)=matunitary_trans(sk(:,:,ik),Os)
         rhok(:,:,ik)=matunitary_trans(rhok(:,:,ik),Os)
         dummy = idenk(:,:,ik)
         CALL ZGEMM_MKL95(Os,dummy,idenk(:,:,ik),'N','N',(1.d0,0.d0),(0.d0,0.d0))

      !Step 3c:  Restoration of natural character of NRB
      !(Omitted; doesn't seem important at all)

         !Now restore the expected order of the basis functions
         sk(:,:,ik)=remap_matrix(sk(:,:,ik),ibasisremap,.TRUE.)
         rhok(:,:,ik)=remap_matrix(rhok(:,:,ik),ibasisremap,.TRUE.)
         idenk(:,:,ik)=remap_matrix(idenk(:,:,ik),ibasisremap,.TRUE.)
      ENDDO

!$OMP END DO NOWAIT
!$OMP END PARALLEL

      w=remap_vector(w,ibasisremap,.TRUE.)


    END SUBROUTINE pre_nao_step_3

    SUBROUTINE pre_nao_step_4a(inp,rhok,sk,idenk,nk,w)
      USE matutil
      USE nbo_shared
      IMPLICIT NONE
      
      TYPE(nbo_input) :: inp
      COMPLEX*16, DIMENSION(:,:,:) :: rhok,sk,idenk
      REAL*8,DIMENSION(:) :: W
      INTEGER :: ik,nk,j

      !"Ow" is the population weighted symmetric orthogonalization matrix, with weights W
      COMPLEX*16, DIMENSION(SIZE(inp%rho0,1),SIZE(inp%rho0,2)) :: Ow
      !"Os" is the Schmidt orthogonalization matrix
      COMPLEX*16, DIMENSION(SIZE(inp%rho0,1),SIZE(inp%rho0,2)) :: Os

      INTEGER :: nvalence,nrydberg,nhiryd,nloryd
      INTEGER, DIMENSION(SIZE(inp%rho0,1)) :: ibasisremap, irydbasisremap
      INTEGER :: ifirst,ilast
      LOGICAL, DIMENSION(SIZE(inp%rho0,1)) :: nmb_mask

      COMPLEX*16, DIMENSION(SIZE(inp%rho0,1),SIZE(inp%rho0,2)) :: dummy
      INTEGER   :: nbasis
 
      !Get a mask of which basis functions should be included in the NMB
      CALL create_nmb_mask(inp,w,nmb_mask)
 
      !Step 4:  Do an occupancy weighted symmetric orthogonalization on the entire NMB + NRB set
      !Since we just orthogonalized the two subspaces, in principle we can do this separately for both the NMB and the NRB
      !First parition the NRB into low occupancy (<1d-4) and higher occupancy
      !CALL partition_basis(w.GT.nmb_thresh,inp%nbasis,nvalence,nrydberg,ibasisremap)
      CALL partition_basis(nmb_mask,inp%nbasis,nvalence,nrydberg,ibasisremap)
      w=remap_vector(w,ibasisremap,.FALSE.)
      CALL partition_basis(w.GT.w_thresh,inp%nbasis,nhiryd,nloryd,irydbasisremap)
      w=remap_vector(w,irydbasisremap,.FALSE.)
      w(nhiryd+1:inp%nbasis)=1.d0

      !WRITE(6,*)W
      !WRITE(6,*)'Inside of step 4'
      !WRITE(6,*)'nvalence,rydberg,nhiryd,nloryd',nvalence,nrydberg,nhiryd,nloryd

      nbasis = inp%nbasis

      !WRITE(6,*)'inside pre nao step 3'


!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(nk,rhok,sk,idenk,w,nvalence,nhiryd,nbasis,irydbasisremap,ibasisremap)
!$OMP DO SCHEDULE(STATIC)

      DO ik=1,nk

         !WRITE(6,*)'for k-point',ik

         sk(:,:,ik)=remap_matrix(sk(:,:,ik),ibasisremap,.FALSE.)
         rhok(:,:,ik)=remap_matrix(rhok(:,:,ik),ibasisremap,.FALSE.)
         idenk(:,:,ik)=remap_matrix(idenk(:,:,ik),ibasisremap,.FALSE.)
         sk(:,:,ik)=remap_matrix(sk(:,:,ik),irydbasisremap,.FALSE.)
         rhok(:,:,ik)=remap_matrix(rhok(:,:,ik),irydbasisremap,.FALSE.)
         idenk(:,:,ik)=remap_matrix(idenk(:,:,ik),irydbasisremap,.FALSE.)

         !WRITE(6,*)'first NMB orthog'

         !First deal with the NMB
         Ow=matiden(nbasis)
         Ow(1:nvalence,1:nvalence)=symmetric_orthog_matrix(sk(1:nvalence,1:nvalence,ik),w(1:nvalence))
         sk(:,:,ik)=matunitary_trans(sk(:,:,ik),Ow)
         rhok(:,:,ik)=matunitary_trans(rhok(:,:,ik),Ow)
         dummy = idenk(:,:,ik)
         CALL ZGEMM_MKL95(Ow,dummy,idenk(:,:,ik),'N','N',(1.d0,0.d0),(0.d0,0.d0))

         !WRITE(6,*)'hi occ ryd orthog'

         !Now deal with the high occupancy NRB
         Ow=matiden(nbasis)
         Ow(nvalence+1:nhiryd,nvalence+1:nhiryd)=symmetric_orthog_matrix(sk(nvalence+1:nhiryd,nvalence+1:nhiryd,ik),w(nvalence+1:nhiryd))
         sk(:,:,ik)=matunitary_trans(sk(:,:,ik),Ow)
         rhok(:,:,ik)=matunitary_trans(rhok(:,:,ik),Ow)
         dummy = idenk(:,:,ik)
         CALL ZGEMM_MKL95(Ow,dummy,idenk(:,:,ik),'N','N',(1.d0,0.d0),(0.d0,0.d0))

         !Finally with the low occupancy NRB
         !Here we first do a Gram-Schmidt orthogonalization of the low occupancy states to the remainder
         Os=schmidt_orthog_matrix(sk(:,:,ik),nhiryd)
         sk(:,:,ik)=matunitary_trans(sk(:,:,ik),Os)
         rhok(:,:,ik)=matunitary_trans(rhok(:,:,ik),Os)
         dummy = idenk(:,:,ik)
         CALL ZGEMM_MKL95(Os,dummy,idenk(:,:,ik),'N','N',(1.d0,0.d0),(0.d0,0.d0))

         !WRITE(6,*)'lo occ ryd orthog'

         !Then do a symmetric (not occupancy weighted) orthogonalization amongst the low occupancy states
         !Even though the symmetric_orthog_matrix subroutine is still used, since all weights were set to one of the low occ rydbergs, the weighting has no effect
         Ow=matiden(nbasis)
         Ow(nhiryd+1:nbasis,nhiryd+1:nbasis)=symmetric_orthog_matrix(sk(nhiryd+1:nbasis,nhiryd+1:nbasis,ik),w(nhiryd+1:nbasis))
         sk(:,:,ik)=matunitary_trans(sk(:,:,ik),Ow)
         rhok(:,:,ik)=matunitary_trans(rhok(:,:,ik),Ow)
         dummy = idenk(:,:,ik)
         CALL ZGEMM_MKL95(Ow,dummy,idenk(:,:,ik),'N','N',(1.d0,0.d0),(0.d0,0.d0))

         !Now restore the expected order of the basis functions
         sk(:,:,ik)=remap_matrix(sk(:,:,ik),irydbasisremap,.TRUE.)
         sk(:,:,ik)=remap_matrix(sk(:,:,ik),ibasisremap,.TRUE.)
         rhok(:,:,ik)=remap_matrix(rhok(:,:,ik),irydbasisremap,.TRUE.)
         rhok(:,:,ik)=remap_matrix(rhok(:,:,ik),ibasisremap,.TRUE.)
         idenk(:,:,ik)=remap_matrix(idenk(:,:,ik),irydbasisremap,.TRUE.)
         idenk(:,:,ik)=remap_matrix(idenk(:,:,ik),ibasisremap,.TRUE.)

         !WRITE(6,*)

      ENDDO

!$OMP END DO NOWAIT
!$OMP END PARALLEL

      w=remap_vector(w,irydbasisremap,.TRUE.)
      w=remap_vector(w,ibasisremap,.TRUE.)

    END SUBROUTINE pre_nao_step_4a


    !**********************************
    !Various associated helper routines
    !**********************************


    !At some points in the procedure the real space g=(/0,0,0/) matrix is needed, even though the general procedure is done in k-space (pre_nao_step2)
    !This subroutine calculates the matrix as simply the average of all the k-space matrices
    !The phase factor for all of these matrices will be 1, since DOT(g,k) is always equal to zero.
    SUBROUTINE real_space_matrix(real_matrix, k_matrix, nk, inp)
      USE nbo_shared
      IMPLICIT NONE

      TYPE(nbo_input),INTENT(IN)                ::  inp
      REAL*8, DIMENSION(:,:),INTENT(OUT)        ::  real_matrix
      COMPLEX*16, DIMENSION(:,:,:), INTENT(IN)  ::  k_matrix
      INTEGER, INTENT(IN)    :: nk

      INTEGER    ::    ik

      real_matrix = 0.d0
      DO ik=1,nk
         real_matrix = real_matrix + k_matrix(:,:,ik)
         IF( inp%kpt_wt(ik) .NE. inp%kpt_wt(1) )THEN
            real_matrix = real_matrix + CONJG(k_matrix(:,:,ik))
         ENDIF
      ENDDO
      real_matrix = real_matrix * inp%kpt_wt(1)

    END SUBROUTINE real_space_matrix



    !Returns a mask array specifying whether the corresponding basis function should be included in the NMB
    !Do this by first SORTING the occupancies, and finding a threshold corresponding to the n-th most
    !occupied pre-NAO, where n is the number of basis functions in a minimal basis.
    SUBROUTINE create_nmb_mask(inp,w,mask)
      USE sortutil
      USE nbo_shared
      IMPLICIT NONE
      TYPE(nbo_input) :: inp
      REAL*8,DIMENSION(:) :: W
      LOGICAL, DIMENSION(SIZE(W,1)) :: mask

      REAL*8,DIMENSION(SIZE(W,1)) :: Wsort
      REAL*8 :: thresh
      INTEGER :: iatom,ifirst,ilast,imbcount

      !WRITE(6,*)'inside of program creat_nmb_mask'
    
      Wsort=W
      DO iatom=1,inp%natom
         !WRITE(6,*)'for atom',iatom

         imbcount=nmb_count(inp,iatom)
         ifirst=inp%ibasismap(iatom)
         ilast=inp%ibasismap(iatom+1)-1
         !WRITE(6,*)'imbcount,ifirst,ilast',imbcount,ifirst,ilast

         CALL heap_sort(Wsort(ifirst:ilast))
         thresh=Wsort(ilast-imbcount+1)
         !WRITE(6,*)'threshhold',thresh
         mask(ifirst:ilast)=w(ifirst:ilast).GE.thresh
      ENDDO
      !WRITE(6,*)'mask',mask
      !WRITE(6,*)

    END SUBROUTINE create_nmb_mask

    !
    !Takes a set of basis functions, locates those of on a given atom (iatom) and a specified
    !angular momentum (il), and constructs a spherically averaged subblock rhoat^{Al} and sat^{Al};
    !ibasisremap stores the basis functions which contributed to each averaged matrix element,
    !i.e. the first basis function in each shell
    SUBROUTINE load_average(inp,iatom,il,rhoat,sat,nshell,ibasisremap)
      USE nbo_shared
      IMPLICIT NONE
      TYPE(nbo_input) :: inp
      !Defines the desired atom / angular momentum block of interest
      INTEGER :: iatom, il
      !Blocks of the density,overlap,and pre-NAO transformation matrices for atom/angular momentum blocks
      REAL*8, DIMENSION(:,:) :: rhoat,sat
      !The number of shells of the desired angular momentum on the atom of interest
      INTEGER :: nshell
      !The re-mapping of the basis functions from the P^{lm} to normal
      INTEGER, DIMENSION(:) :: ibasisremap

      INTEGER :: ifirst,ilast,ibasis,jfirst
      INTEGER :: icurshell,ishell,jshell
      INTEGER :: nm,im

      !keep track of the first basis function of the desired angular momentum in each shell
      ibasisremap=0
      icurshell=0
      nshell=0

      ifirst=inp%ibasismap(iatom)
      ilast=inp%ibasismap(iatom+1)-1
      DO ibasis=ifirst,ilast
         IF (inp%ilmap(ibasis).EQ.il) THEN !retain this basis function
            IF (inp%ishellmap(ibasis).NE.icurshell) THEN !we are in a new shell
               nshell=nshell+1
               ibasisremap(nshell)=ibasis
               icurshell=inp%ishellmap(ibasis)
            ENDIF
         ENDIF
      ENDDO

      !WRITE(6,*)'atom, ang mom, nshell',iatom,il,nshell

      rhoat=0.d0
      sat=0.d0
      !Step 2:  Use this mapping to construct the rho^{Al} and s^{Al} spherically averaged blocks
      nm=2*il+1
      DO ishell=1,nshell
         DO jshell=1,nshell
            ifirst=ibasisremap(ishell)
            jfirst=ibasisremap(jshell)
            DO im=1,nm
               rhoat(ishell,jshell)=rhoat(ishell,jshell)+inp%rho0(ifirst+im-1,jfirst+im-1,1,1)
               sat(ishell,jshell)=sat(ishell,jshell)+inp%s0(ifirst+im-1,jfirst+im-1,1)
            ENDDO
            !sat(ishell,jshell)=inp%s0(ifirst+nm-1,jfirst+nm-1,1)
         ENDDO
      ENDDO
      DO ishell=1,nshell
         DO jshell=1,nshell
            rhoat(ishell,jshell)=rhoat(ishell,jshell)/nm
            sat(ishell,jshell)=sat(ishell,jshell)/nm
         ENDDO
      ENDDO

    END SUBROUTINE load_average

    !
    !Returns the transformation matix which does a Gram-Schmidt
    !orthogonalization of the remaining basis vectors to the first n,
    !given an overlap matrix S
    FUNCTION schmidt_orthog_matrix(S,n) RESULT(Os)
      USE matutil
      IMPLICIT NONE
      COMPLEX*16,DIMENSION(:,:) :: S
      COMPLEX*16,DIMENSION(n,n) :: Sinv
      COMPLEX*16,DIMENSION(n,SIZE(S,2)-n) :: proj
      COMPLEX*16,DIMENSION(SIZE(S,1),SIZE(S,2)) :: Os
      INTEGER :: n,i,j,nbasis
      
      nbasis=SIZE(S,1)
     
      !This is tricky; since we are in a non-orthogonal basis set and doing only a PARTIAL gs
      !orthogonalization, we need to find the orthogonal projection of each basis function 1:n
      !onto the remaining basis function n+1:nbasis
      Sinv=matinv(S(1:n,1:n))
      !proj=MATMUL(Sinv,S(1:n,n+1:nbasis))
      CALL ZGEMM_MKL95(Sinv,S(1:n,n+1:nbasis),proj,'N','N',(1.d0,0.d0),(0.d0,0.d0))

      !Now do gs orthogonalization using those projections
      Os=0.d0
      DO i=1,nbasis
         Os(i,i)=1.d0
      ENDDO
      DO i=n+1,nbasis
         DO j=1,n
            Os(j,i)=-proj(j,i-n)
         ENDDO
      ENDDO
      Os=CONJG(TRANSPOSE(Os))
    END FUNCTION schmidt_orthog_matrix

    !
    !Return the transformation matrix for the OCCUPANCY WEIGHTED (w)
    !orthogonalization matrix
    FUNCTION symmetric_orthog_matrix(sk,w) RESULT(Ow)
      USE matutil
      IMPLICIT NONE
      COMPLEX*16 , DIMENSION(:,:):: sk
      COMPLEX*16 , DIMENSION(SIZE(sk,1),SIZE(sk,2)):: Ow
      REAL*8 , DIMENSION(:) :: w

       !Ow=W(WSW)^-1/2
       Ow=sk
       Ow=Ow*outer_product(w,w) !Ow(i,j)=Ow(i,j)*w(i)*w(j)
       Ow=matinvsqrt(Ow)
       Ow=Ow*SPREAD(w,1,SIZE(w))!Ow(i,:)=Ow(i,:)*w(i)
    END FUNCTION symmetric_orthog_matrix

    !
    !Sorts the eigenvectors of the pre-NAO transformation matrix
    !such that they maximimally preserve the order
    !of the corresponding basis function
    !
    SUBROUTINE sort_pre_nao(N,W, S)
      USE matutil
      IMPLICIT NONE
      REAL*8, DIMENSION(:,:) :: N
      REAL*8, OPTIONAL, DIMENSION(:,:) :: S
      REAL*8, DIMENSION(SIZE(N,1),SIZE(N,2)) :: overlap
      REAL*8, DIMENSION(:) :: W

      REAL*8, DIMENSION(SIZE(N,1)) :: swap
      INTEGER :: nbasis,i,j

      nbasis=SIZE(N,1)
      IF (PRESENT(S)) THEN
         !overlap=MATMUL(matsqrt(S),N)
         CALL DGEMM_MKL95(matsqrt(S),N,overlap,'N','N',1.d0,0.d0)
      ELSE
         overlap=N
      ENDIF

      DO i=1,nbasis
         j=MAXLOC(overlap(i,:)**2,1)

         swap=overlap(:,i)
         overlap(:,i)=overlap(:,j)
         overlap(:,j)=swap

         swap=N(:,i)
         N(:,i)=N(:,j)
         N(:,j)=swap

         swap(1)=W(i)
         W(i)=W(j)
         W(j)=swap(1)
      ENDDO

    END SUBROUTINE sort_pre_nao

END MODULE pre_nao
