PROGRAM projection_main
  USE rd_wavefunction
  USE projection_shared
  USE rd_basis
  USE bloch_overlap
  USE PAW
  USE blas95
  IMPLICIT NONE

  CHARACTER(128)  ::  basis_fn, VASP_fn, NBO_fn
  CHARACTER(64)   ::  mat_fn

  TYPE(AO_function), ALLOCATABLE :: AO_basis(:)

  INTEGER, ALLOCATABLE       ::  index_l(:,:)       !For each l_vector gives the number of unitcells displaces in each direction for the l vector

  COMPLEX*16,ALLOCATABLE   ::  PAW_overlap(:,:,:)
  COMPLEX*16,ALLOCATABLE   ::  bloch_band_coeff(:,:,:,:)
  !COMPLEX*16,ALLOCATABLE   ::  AO_PW_overlap(:,:,:)
  COMPLEX*16,ALLOCATABLE   ::  AO_PW_overlap(:,:)
  COMPLEX*16               ::  gdotrnu
  !COMPLEX*16,ALLOCATABLE   ::  proj_matrix(:,:,:)
  COMPLEX*16,ALLOCATABLE   ::  proj_matrix(:,:,:,:)


  COMPLEX*16,ALLOCATABLE   ::  proj_overlap(:,:,:,:)

  !COMPLEX*16,ALLOCATABLE   ::  bloch_density(:,:,:),bloch_fock(:,:,:)

  COMPLEX*16,ALLOCATABLE   ::  fock_dummy(:,:),coeff_dummy(:,:)

  COMPLEX*16,ALLOCATABLE   ::  n_elec_dummy(:,:)
  COMPLEX*16               ::  num_elec,energy_sum

  INTEGER                  ::  nu,ik,ig,iband,ispin
  !REAL*8, PARAMETER        ::  pi=3.141592653589793238462d0
  !COMPLEX*16, PARAMETER    ::  sqrt_neg_one=(0.d0,1.d0)

  REAL*8   :: ti,tf,t1,t2

  CALL CPU_TIME(ti)

  !Start by getting the names of the files to use as input (basis and wavefunction) and output (NBO)
  !If no filenames are supplied the generics will be used instead
  IF( IARGC() .LT. 1 )THEN 
    basis_fn = 'basis.inp'
  ELSE
    CALL GETARG(1, basis_fn)
  ENDIF
  IF( IARGC() .LT. 2 )THEN 
    VASP_fn = 'wavefunction.dat'
  ELSE
    CALL GETARG(2, VASP_fn)
  ENDIF
  IF( IARGC() .LT. 3 )THEN 
    NBO_fn = 'NBO.out'     
  ELSE
    CALL GETARG(3, NBO_fn) 
  ENDIF

  mat_fn = 'NBO_mat.out'
  
  WRITE(6,'(2A)')' Input file for AO basis  ',basis_fn
  WRITE(6,'(2A)')' Input file for PW calc   ',VASP_fn
  WRITE(6,'(2A)')' Output file for NBO      ',NBO_fn
  WRITE(6,'(2A)')' Matrix file for NBO      ',mat_fn
  WRITE(6,*)

  !Read in information about the underlying PW calculation
  CALL read_vasp_wavefunction(VASP_fn)

  !Read in input file for atomic orbital basis set
  CALL read_basis(basis_fn,AO_basis)
  !This subroutine has been made obsolete by a change to the data output to NBO.out
  !Basically all the unique data contained in here is now passed to the nbo program,
  !That program then generates a single file which can be read in by an auxillary fortran code to produce cube files.
  !!!Then output information on all basis functions in a readable format for a visualization program.
  !!!CALL basis_set_info_out(AO_basis)

  WRITE(6,*)'****************************************'
  WRITE(6,*)'*** Beginning Projection Calculation ***'
  WRITE(6,*)'****************************************'
  WRITE(6,*)

  !Calcualtes the overlap matrix inverse at each k-point, which is required for the projector operator.
  CALL bloch_space_overlap(AO_basis,index_l)

!  !Now ready to perform actual projection.
!  CALL CPU_TIME(t1)

  !If the PW calculation used PAW pseudopotentials, their contribution to the projection (Bloch-space overlap with AO basis) must be calculated
  !PAW_pseudo = .FALSE.
  IF( PAW_pseudo ) CALL PAW_proj_setup(AO_basis,index_l,PAW_overlap)

  !ALLOCATE(bloch_band_coeff(s_dim,nbands,nkpts,nspins),proj_matrix(s_dim,nbands,nkpts),AO_PW_overlap(s_dim,nplmax,nkpts))
  ALLOCATE(bloch_band_coeff(s_dim,nbands,nkpts,nspins),proj_matrix(s_dim,nbands,nkpts,nspins),AO_PW_overlap(nplmax,nkpts))
  AO_PW_overlap=0.d0

  !Now the k-point dependent projections will be done.
  !For each kpt the projection is done independently
  !WRITE(6,*)'Beginning of projections for each k-point'

  !Now ready to perform actual projection.
  !CALL CPU_TIME(t1)


!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ik,nu,ig,ispin,gdotrnu)
!!!!!SHARED(nkpts,s_dim,npl,nspins,AO_basis,gk,AO_PW_overlap,PAW_overlap,pw_coeff,PAW_coeff,proj_matrix,bloch_s_inv,bloch_band_coeff,PAW_pseudo)
!$OMP DO SCHEDULE(STATIC)
  DO ik=1,nkpts

     !WRITE(6,*)'kpt',ik,npl(ik)

     !!!!!!!!!!!!
     !!!!NOTE!!!!
     !!!!!!!!!!!!
     !I previously calculated the AO_PW overlap for all basis function in the top loop, 
     !this was then passed to the bottom loop, where matrix multiplication was to treat all AOs at the same time.  
     !This was faster, but had huge memore requirements.  
     !That array also had a k-point dimension, for parallelization, which only made things worse
     !Switching to one AO function at a time treated via matrix-vector multiplication is a bit slower but requires an order of magnitude lower memory requirements.
     !
     !I think a happy middle ground could be reached by treating blocks of AO's using matrix mulitplication
     !The blocks could be the number of bands, so that only the amount of memory used for the plane wave coefficients would be needed.
     !Then a few large matrix multiplications would be needed, instead of the many vector-matrix multiplications I am currently doing.
     !For another day
     !

     DO nu=1,s_dim
        !First the overlap of each nu basis function with each planewave must be calculated
        !WRITE(6,*)nu
        DO ig=1,npl(ik) 
           !WRITE(6,*)ig
           !CALL nu_g_overlap(AO_basis(nu),gk(:,ig,ik),AO_PW_overlap(nu,ig,ik))      !Overlap of AO with PW 
           !gdotrnu = EXP(sqrt_neg_one*DOT_PRODUCT(gk(:,ig,ik),AO_basis(nu)%pos)) !Phase factor accounting for AO's position off the origin
           !AO_PW_overlap(nu,ig,ik) = AO_PW_overlap(nu,ig,ik) * gdotrnu
           !!WRITE(6,*)AO_PW_overlap(nu,ig,ik)

           CALL nu_g_overlap(AO_basis(nu),gk(:,ig,ik),AO_PW_overlap(ig,ik))      !Overlap of AO with PW
           gdotrnu = EXP(sqrt_neg_one*DOT_PRODUCT(gk(:,ig,ik),AO_basis(nu)%pos)) !Phase factor accounting for AO's position off the origin
           AO_PW_overlap(ig,ik) = AO_PW_overlap(ig,ik) * gdotrnu
           !WRITE(6,*)AO_PW_overlap(ig,ik)
        ENDDO !End loop over plane wave basis

        !Then for each basis function, nu, the contribution of the planewaves is calculated for all bands at one using vector-matrix multiplication
        DO ispin=1,nspins
           !each computation is SUM(nu) S-1{mu,nu} * (SUM(g) <nu|g>*c(a,g))
           CALL ZGEMV_F95(pw_coeff(:,:,ik,ispin),AO_PW_overlap(:,ik),proj_matrix(nu,:,ik,ispin),(1.d0,0.d0),(0.d0,0.d0),'T')
        ENDDO
     ENDDO

     DO ispin=1,nspins
!        !Calculation of planewave contributions to band coeff's are done at once using matrix multiplication
!        !each computation is SUM(nu) S-1{mu,nu} * (SUM(g) <nu|g>*c(a,g))
!        CALL ZGEMM_F95(AO_PW_overlap(:,:,ik),pw_coeff(:,:,ik,ispin),proj_matrix(:,:,ik),'N','N',(1.d0,0.d0),(0.d0,0.d0))
!
!        !Since band is simply summation of PW AND PAW, PAW contirbution is simply added to PW's
!        IF( PAW_pseudo )  CALL ZGEMM_F95(PAW_overlap(:,:,ik),PAW_coeff(:,:,ik,ispin),proj_matrix(:,:,ik),'N','N',(1.d0,0.d0),(1.d0,0.d0))
!
!        !Band overlaps are finally multiplied by inverse of overlap matrix to complete projection
!        !The use of the inverse of the overlap matrix is necessary due to the non-orthogonality of the AO-basis
!        CALL ZGEMM_F95(bloch_s_inv(:,:,ik),proj_matrix(:,:,ik),bloch_band_coeff(:,:,ik,ispin),'N','N',(1.d0,0.d0),(0.d0,0.d0))

        !Since band is simply summation of PW AND PAW, PAW contirbution is simply added to PW's
        IF( PAW_pseudo )  CALL ZGEMM_F95(PAW_overlap(:,:,ik),PAW_coeff(:,:,ik,ispin),proj_matrix(:,:,ik,ispin),'N','N',(1.d0,0.d0),(1.d0,0.d0))

        !Band overlaps are finally multiplied by inverse of overlap matrix to complete projection
        !The use of the inverse of the overlap matrix is necessary due to the non-orthogonality of the AO-basis
        CALL ZGEMM_F95(bloch_s_inv(:,:,ik),proj_matrix(:,:,ik,ispin),bloch_band_coeff(:,:,ik,ispin),'N','N',(1.d0,0.d0),(0.d0,0.d0))

     ENDDO

  ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL


  DEALLOCATE(proj_matrix,AO_PW_overlap)
  DEALLOCATE(pw_coeff)
  DEALLOCATE(bloch_s_inv)
  !PAUSE


  !CALL CPU_TIME(t2)
  !WRITE(6,*)'Completed projection',SNGL(t2-t1)
  !WRITE(6,*)


  !Calculate ovelrap matrices for projected bands and use these to quantify completeness of projection.
  CALL calc_spillover(proj_overlap,bloch_band_coeff)

  !Calculate density and fock matrices in projected AO basis
  ALLOCATE(coeff_dummy(s_dim,nbands),fock_dummy(nbands,s_dim))
  ALLOCATE(bloch_density(s_dim,s_dim,nkpts,nspins),bloch_fock(s_dim,s_dim,nkpts,nspins))
  DO ispin=1,nspins
     DO ik=1,nkpts
        !Start by orthogonalizeing the bands in the AO basis, using an Occupnacy Weighted Symmetric Orthogonalization
        !Even though the PW bands are orthogonal, since the projections are not exact the projected bands lose this quality
        !This orthogonality greatly eases the calculation of both density and Fock matrices
        !Occupancy weighted symmetric orthogonaliztion is used, since it maximally maintains the character of the projected bands
        !Sacrificing the character of unoccupied bands is OKAY for us, since they do not affect the density matrix.
        CALL do_OWSO(proj_overlap(:,:,ik,ispin),weight(:,ik,ispin),bloch_band_coeff(:,:,ik,ispin))

        !DENSITY MATRIX
        ! P{mu,nu} = occ{iband} * band_coeff{mu,iband} * CONJG(band_coeff{nu,iband})
        !First scale the band coeff's by the weight of that band.
        DO iband=1,nbands
           coeff_dummy(:,iband) = bloch_band_coeff(:,iband,ik,ispin) * weight(iband,ik,ispin)
        ENDDO
        !Matrix muliplication is used to calculate the density matrix
        CALL ZGEMM_F95(coeff_dummy,bloch_band_coeff(:,:,ik,ispin),bloch_density(:,:,ik,ispin),'N','C',(1.d0,0.d0),(0.d0,0.d0))

        !FOCK MATRIX
        !Matrix multiplication is used to convert the Fock matrix into the projected basis
        ! F{mu,nu} = SUM(fock_coeff{nu,iband}*F{iband,iband}*CONJG(fock_coeff{nu,iband}))
        !First the nonorthogonality of the bloch orbtials must be addressed since the MO->AO transform is not unitary 
        !coeff's must be multiplied by Overlap matrix
        CALL ZGEMM_F95(bloch_s_mat(:,:,ik),bloch_band_coeff(:,:,ik,ispin),coeff_dummy,'N','N',(1.d0,0.d0),(0.d0,0.d0))
        !Then the trnsformation is simply a unitary transform using the new coefficients
        DO iband=1,nbands
           fock_dummy(iband,:) = eig(iband,ik,ispin)*CONJG(coeff_dummy(:,iband))
        ENDDO
        CALL ZGEMM_F95(coeff_dummy,fock_dummy,bloch_fock(:,:,ik,ispin),'N','N',(1.d0,0.d0),(0.d0,0.d0))
     ENDDO
  ENDDO


  !Now we can check the projection by comparing observables.
  !For the density matrix we will calculate the corersponding number of electrons
  !For the Fock matrix the sum of the energies of the occupied states is calculated
  !These should match the original values from the PW results
  !Q = Trace[P * ADJ(S)]
  ALLOCATE(n_elec_dummy(s_dim,s_dim))
  num_elec = 0.d0
  energy_sum = 0.d0
  DO ispin=1,nspins
     DO ik=1,nkpts
        CALL ZGEMM_F95(bloch_density(:,:,ik,ispin),bloch_s_mat(:,:,ik),n_elec_dummy,'N','N',(1.d0,0.d0),(0.d0,0.d0))
        n_elec_dummy = n_elec_dummy*kpt_wt(ik)
        DO nu=1,s_dim
           num_elec = num_elec + n_elec_dummy(nu,nu)
        ENDDO
        CALL ZGEMM_F95(bloch_density(:,:,ik,ispin),bloch_fock(:,:,ik,ispin),n_elec_dummy,'N','N',(1.d0,0.d0),(0.d0,0.d0))
        n_elec_dummy = n_elec_dummy*kpt_wt(ik)
        DO nu=1,s_dim
           energy_sum = energy_sum + n_elec_dummy(nu,nu)
        ENDDO
     ENDDO
  ENDDO

  WRITE(6,*)'Testing of projected density and Fock matrices'
  WRITE(6,*)'Bloch-space energy sum        ',energy_sum*DBLE(3-nspins)
  WRITE(6,*)'Bloch-space valence electrons ',num_elec*DBLE(3-nspins)
  WRITE(6,*)

  CALL CPU_TIME(t1)
  !WRITE(6,*)'Time for density and Fock matrices',SNGL(t1-t2)
  !WRITE(6,*)

  !Write out all information needed for NBO analysis
  !Formatting is specifc to NBO code of JRS and BDD
  CALL write_NBO_output(NBO_fn,mat_fn,AO_basis,index_l)


  CALL CPU_TIME(tf)
  WRITE(6,*)'Total processor time for projection program',SNGL(tf-ti)

CONTAINS


!This subrotuine computes the overlap matrices of the projected bands in the AO basis.
!This information is also used to calculate the spillover (quantify how much of band density is lost in projection
SUBROUTINE calc_spillover(r_mat,band_coeff)
  USE blas95
  USE rd_wavefunction
  USE projection_shared
  IMPLICIT NONE

  COMPLEX*16,ALLOCATABLE,INTENT(OUT)  ::  r_mat(:,:,:,:)  !This will contain the band overlap matrices on exit from the program.  
  COMPLEX*16,DIMENSION(:,:,:,:),INTENT(IN) :: band_coeff  !This contains the coefficient of each band in the AO-basis at each k-point

  COMPLEX*16,DIMENSION(s_dim,nbands)   ::  r_mat_dummy  !Used for BLAS subroutines
  REAL*8,DIMENSION(nspins)      ::  spillover 
  REAL*8,DIMENSION(nspins)      ::  spread
  COMPLEX*16  ::    chi_overlap
  REAL*8      ::    ind_spill
  INTEGER,DIMENSION(nspins)     ::  norm_tally

  REAL*8,DIMENSION(nspins)    ::  band_spillover
  REAL*8,DIMENSION(nspins)    ::  weight_spillover,weight_spread
  REAL*8,DIMENSION(nspins)    ::  band_weight_spillover,band_weight_spread
  INTEGER,DIMENSION(nspins)   ::  weight_count
  INTEGER,DIMENSION(nspins)   ::  tot_weight

  REAL*8,DIMENSION(nkpts,2)   ::  max_spread
  REAL*8                      ::  coeff_sum

  REAL*8,DIMENSION(n_atom,nspins)    ::  atom_spillover,atom_spread,atom_norm
  REAL*8,DIMENSION(n_atom)           ::  atom_sum
  INTEGER                     ::  iatom

  INTEGER   ::  ik,j,iband,ispin

  !All spillover analysis is writen out to this 
  OPEN(7,file='band_spillover.out')
  spillover = 0.d0
  spread = 0.d0
  weight_spread = 0.d0
  weight_spillover = 0.d0
  tot_weight = 0

  norm_tally = 0
  max_spread = 0

  atom_spillover = 0.d0
  atom_spread = 0.d0
  atom_norm = 0

  !CO_count = 0
  !CO_spillover = 0.d0
  !CO_spread = 0.d0
  !factor_test = 0.d0

  ALLOCATE(r_mat(nbands,nbands,nkpts,nspins))
  r_mat = (0.d0,0.d0)
  r_mat_dummy = (0.d0,0.d0)


  DO ik=1,nkpts

     band_spillover = 0.d0
     band_weight_spillover = 0.d0 
     band_weight_spread = 0.d0
     !weight_spillover = 0.d0
     weight_count = 0

     !WRITE(6,*)'band coefficients for kpt',ik
     !DO iband=1,nbands
     !   WRITE(6,*)'coeff for band',iband
     !   WRITE(6,*)band_coeff(:,iband,ik,1)
     !   WRITE(6,*)
     !ENDDO

     DO ispin=1,nspins

     !Overlap of projected bands is computed using matrix mulitplication
     !Each overlap is equal to SUM(mu) SUM(nu) c*S{mu,nu}c
     CALL ZGEMM_F95(bloch_s_mat(:,:,ik),band_coeff(:,:,ik,ispin),r_mat_dummy,'N','N',(1.d0,0.d0),(0.d0,0.d0))
     CALL ZGEMM_F95(band_coeff(:,:,ik,ispin),r_mat_dummy,r_mat(:,:,ik,ispin),'C','N',(1.d0,0.d0),(0.d0,0.d0))


     !Then the norm of each band at each k-point is analyzed to make sure the norm is appropriate
     !And spillover is calculated for all bands as well as occupied bands at each k-point, and overall (average over k-points)
     DO iband=1,nbands

        !IF( weight(iband,ik,ispin) > 1.d-5 )THEN
        !   coeff_sum = SUM(CONJG(band_coeff(:,iband,ik,ispin)) * band_coeff(:,iband,ik,ispin))
        !   DO iatom=1,n_atom
        !      atom_sum(iatom) = SUM(CONJG(band_coeff(atom_basis(iatom):atom_basis(iatom+1)-1,iband,ik,ispin)) &
        !                           &* band_coeff(atom_basis(iatom):atom_basis(iatom+1)-1,iband,ik,ispin))
        !      atom_sum(iatom) = atom_sum(iatom) / coeff_sum
        !   ENDDO
        !ENDIF

        chi_overlap = r_mat(iband,iband,ik,ispin)
        ind_spill = 1.d0 - REAL(chi_overlap)

        IF( REAL(chi_overlap) > 1.01d0 .AND. weight(iband,ik,ispin) > 1.d-5 )THEN
           WRITE(6,*)'The norm of band',iband,' with spin',ispin,'is larger than 1.01 for kpt',ik, SNGL(REAL(chi_overlap))
           norm_tally(ispin) = norm_tally(ispin) + 1
        ENDIF

        WRITE(7,*)'normalization for band',iband,'with spin',ispin,'kpt',ik
        WRITE(7,*)SNGL(REAL(chi_overlap))
        WRITE(7,*)

        band_spillover(ispin) = band_spillover(ispin) +  ind_spill
        spread(ispin) = spread(ispin) + ABS(ind_spill)
        IF( weight(iband,ik,ispin) > 1.d-5 )THEN
           band_weight_spillover(ispin) = band_weight_spillover(ispin) + ind_spill
           band_weight_spread(ispin) = band_weight_spread(ispin) + ABS(ind_spill)
           weight_count(ispin) = weight_count(ispin) + 1
           IF( ABS(ind_spill) .GT. ABS(max_spread(ik,1)) )THEN
              max_spread(ik,1) = ind_spill
              max_spread(ik,2) = iband
           ENDIF

           !Atomic weighted spread and spillover
           coeff_sum = SUM(CONJG(band_coeff(:,iband,ik,ispin)) * band_coeff(:,iband,ik,ispin))
           DO iatom=1,n_atom
              !Calculate atomic weighting factor
              atom_sum(iatom) = SUM(CONJG(band_coeff(atom_basis(iatom):atom_basis(iatom+1)-1,iband,ik,ispin)) &
                                   &* band_coeff(atom_basis(iatom):atom_basis(iatom+1)-1,iband,ik,ispin))
              atom_sum(iatom) = atom_sum(iatom) / coeff_sum

              !Sum up spill and spread, as well as a normalization factor
              atom_spillover(iatom,ispin) = atom_spillover(iatom,ispin) + atom_sum(iatom)*ind_spill
              atom_spread(iatom,ispin) = atom_spread(iatom,ispin) + atom_sum(iatom)*ABS(ind_spill)
              atom_norm(iatom,ispin) = atom_norm(iatom,ispin) + atom_sum(iatom)
           ENDDO

        ENDIF

     ENDDO
     spillover(ispin) = spillover(ispin) + band_spillover(ispin)
     weight_spillover(ispin) = weight_spillover(ispin) + band_weight_spillover(ispin) !/ DBLE(weight_count(ispin)))
     weight_spread(ispin) = weight_spread(ispin) + band_weight_spread(ispin) !/ DBLE(weight_count(ispin)))
     tot_weight(ispin) = tot_weight(ispin) + weight_count(ispin)

     WRITE(7,*)'For spin',ispin,'kpt',ik,' spillover is'
     WRITE(7,*)SNGL(band_spillover(ispin) / DBLE(nbands))
     WRITE(7,*)
     WRITE(7,*)'For spin',ispin,'kpt',ik,'occupied bands',weight_count(ispin),'spillover'
     WRITE(7,*)SNGL(band_weight_spillover(ispin) / DBLE(weight_count(ispin)))
     WRITE(7,*)
     WRITE(7,*)

     ENDDO


  ENDDO

  spillover = spillover / DBLE(nbands*nkpts)
  spread = spread / DBLE(nbands*nkpts)
  weight_spillover = weight_spillover / tot_weight
  weight_spread = weight_spread / tot_weight

  !!Only total spillover and number of improper norms are written to output stream
  !DO ispin=1,nspins
  !   DO j=6,7
  !      WRITE(j,*)'Total Spillover from projection for spin',ispin
  !      WRITE(j,*)SNGL(spillover(ispin))
  !      WRITE(j,*)
  !      WRITE(j,*)'Total Spread from projection for spin',ispin
  !      WRITE(j,*)SNGL(spread(ispin))
  !      WRITE(j,*)
  !      WRITE(j,*)'Total Spillover for occupied bands for spin',ispin
  !      WRITE(j,*)SNGL(weight_spillover(ispin))
  !      WRITE(j,*)
  !      WRITE(j,*)'Total spread for occupied bands for spin',ispin
  !      WRITE(j,*)SNGL(weight_spread(ispin))
  !      WRITE(j,*)
  !      WRITE(j,*)'There were',norm_tally(ispin),'improper norms out of',nbands*nkpts
  !      WRITE(j,*)
  !   ENDDO
  !ENDDO

  CLOSE(7)


  OPEN(65,file='spillover.out')

  WRITE(65,*)' ################################################### '
  WRITE(65,*)' #######   Projection Quality Measurements   ####### '
  WRITE(65,*)' ################################################### '
  WRITE(65,*)

  WRITE(65,*)' ###  Overall System  ### '
  DO ispin=1,nspins
     IF( nspins .GT. 1)WRITE(65,*)'For spin type, ',ispin
     WRITE(65,*)'Total Spillover'
     WRITE(65,*)SNGL(spillover(ispin))
     WRITE(65,*)
     WRITE(65,*)'Total Spread'
     WRITE(65,*)SNGL(spread(ispin))
     WRITE(65,*)
     WRITE(65,*)'Total Spillover for occupied bands'
     WRITE(65,*)SNGL(weight_spillover(ispin))
     WRITE(65,*)
     WRITE(65,*)'Total Spread for occupied bands'
     WRITE(65,*)SNGL(weight_spread(ispin))
     WRITE(65,*)
     WRITE(65,'(A,I4,A,I8,A)')' There were',norm_tally(ispin),' improper norms (>1.01) out of',tot_weight,' occupied bands'
     WRITE(65,*)
  ENDDO
  WRITE(65,*)
  WRITE(65,*)


  WRITE(65,*)' ###  Atomic Weighted Spillover and Spread  ### '
  DO ispin=1,nspins
     IF( nspins .GT. 1)WRITE(65,*)'For spin type,   ',ispin
     WRITE(65,*)'Atom        Spillover           Spread'
     DO iatom=1,n_atom
        atom_spillover(iatom,ispin) = atom_spillover(iatom,ispin) / atom_norm(iatom,ispin)
        atom_spread(iatom,ispin) = atom_spread(iatom,ispin) / atom_norm(iatom,ispin)
        WRITE(65,'(I5,2F20.10)')iatom,atom_spillover(iatom,ispin),atom_spread(iatom,ispin)
     ENDDO
     WRITE(65,*)
  ENDDO
  WRITE(65,*)
  WRITE(65,*)


  WRITE(65,*)' ###  Worst Bands  ### '
  WRITE(65,*)'Worst spread, averaged across all k-points'
  WRITE(65,*)SUM(max_spread(:,1))/DBLE(nkpts)
  WRITE(65,*)
  WRITE(65,*)'Worst spread at each k-point'
  DO ik=1,nkpts
     WRITE(65,'(A8,I3,A5,I4,A4,F)')'For kpt',ik,'band',INT(max_spread(ik,2)),'at',max_spread(ik,1)
  ENDDO
  WRITE(65,*)


  CLOSE(65)




END SUBROUTINE calc_spillover







END PROGRAM projection_main
