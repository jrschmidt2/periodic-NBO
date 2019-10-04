PROGRAM crystal_nbo
  USE shared
  USE bloch_overlap
  IMPLICIT NONE

  !!!!From parsed standard CRYSTAL output

  !Basis set
  TYPE(AO_function),DIMENSION(:),ALLOCATABLE  ::  AO_basis

  !System Parameters
  INTEGER  ::  nk,nbasis,natom,ibz_nk

  !Output Data for NBO.out header
  INTEGER, DIMENSION(:), ALLOCATABLE       ::  iatnum
  CHARACTER(2), DIMENSION(:), ALLOCATABLE  ::  symbols
  INTEGER, DIMENSION(:), ALLOCATABLE       ::  ishellmap
  INTEGER, DIMENSION(:), ALLOCATABLE       ::  ibasismap
  INTEGER, DIMENSION(:), ALLOCATABLE       ::  ilmap
  REAL*8, DIMENSION(:), ALLOCATABLE        ::  iatval

  !Geometric system information
  REAL*8,DIMENSION(3,3)  ::  latt_vec
  REAL*8,DIMENSION(:,:),ALLOCATABLE  ::  atom_pos


  !!!!From properties CRYAPI_OUT KRED output (KRED.dat)

  !KPoint information
  INTEGER, DIMENSION(:,:), ALLOCATABLE  ::  k_grid
  REAL*8, DIMENSION(:,:), ALLOCATABLE   ::  kpt
  REAL*8, DIMENSION(:), ALLOCATABLE     ::  sym_weight,k_weight
  INTEGER    ::  k_test, k_mesh(3)
  INTEGER, DIMENSION(:), ALLOCATABLE    ::  k_type, sym_type
  INTEGER, DIMENSION(:), ALLOCATABLE    ::  sym_count
  INTEGER, DIMENSION(:,:), ALLOCATABLE  ::  bz_k

  !Band information
  REAL*8, DIMENSION(:,:), ALLOCATABLE  ::  band_eig, sym_eig
  REAL*8, DIMENSION(:,:), ALLOCATABLE  ::  band_occ, sym_occ
  COMPLEX*16, DIMENSION(:,:,:), ALLOCATABLE  ::  band_coeff
  COMPLEX*16, DIMENSION(:,:), ALLOCATABLE  ::  band_dummy

  INTEGER,DIMENSION(:),ALLOCATABLE  ::  rearr_index

  INTEGER,DIMENSION(3)  ::  grid_test, grid_neg
  LOGICAL  :: k_fill  !All k-points in BZ are written out make sure I am not get a redundant kpoint

  INTEGER  ::  skip

  !This is variables to actually do things with our collected data

  !Overlap matrices
  COMPLEX*16, DIMENSION(:,:,:), ALLOCATABLE ::  bloch_s
  COMPLEX*16, DIMENSION(:,:), ALLOCATABLE  ::  r_mat  !overlap matrix in bands basis for testing. Used as dummy throughout

  !Density matrix
  COMPLEX*16, DIMENSION(:,:,:), ALLOCATABLE ::  bloch_rho
  !Fock matrix
  COMPLEX*16, DIMENSION(:,:,:), ALLOCATABLE ::  bloch_fock

  !Testing density matrix
  COMPLEX*16, DIMENSION(:,:), ALLOCATABLE  ::  nelec_mat
  COMPLEX*16                               ::  nelec
  REAL*8                                   ::  CRYS_nelec  !For output from CRYSTAL

  !This array holds the indices of the unit cells to be used in calculating real space overlap. 
  !The indices are counts of unit cells moved away from the central.
  INTEGER, ALLOCATABLE  ::  index_l(:,:)
  INTEGER               ::  l_half

  CHARACTER*128  ::  buffer

  !Counters
  INTEGER  ::  ibasis,igauss,icart,iatom,ik,jk,isym
  INTEGER  ::  i,j,k

  !Get the name of the file created by the python script
  CALL GETARG(1,buffer)
  OPEN(12,file=buffer)


  READ(12,*)ibz_nk
  READ(12,*)natom
  READ(12,*)nbasis

  WRITE(6,*)nbasis,'nbasis'
  WRITE(6,*)natom,'natom'
  WRITE(6,*)ibz_nk,'nk'
  WRITE(6,*)


  ALLOCATE(AO_basis(nbasis))


  DO ibasis=1,nbasis
     READ(12,*)AO_basis(ibasis)%ngauss
     ALLOCATE(AO_basis(ibasis)%norm(AO_basis(ibasis)%ngauss))
     ALLOCATE(AO_basis(ibasis)%coeff(AO_basis(ibasis)%ngauss))
     ALLOCATE(AO_basis(ibasis)%alpha(AO_basis(ibasis)%ngauss))

     READ(12,*)AO_basis(ibasis)%alpha
     READ(12,*)AO_basis(ibasis)%coeff

     READ(12,*)AO_basis(ibasis)%pos
     AO_basis(ibasis)%pos = AO_basis(ibasis)%pos / bohr
     READ(12,*)AO_basis(ibasis)%atom

     READ(12,*)AO_basis(ibasis)%l,AO_basis(ibasis)%m

     CALL norm_calc(AO_basis(ibasis)%alpha,AO_basis(ibasis)%norm,AO_basis(ibasis)%l)

     CALL interp_cart(AO_basis(ibasis))

     !WRITE(6,*)'For basis function',ibasis
     !WRITE(6,*)'Atomic center',AO_basis(ibasis)%atom
     !WRITE(6,'(A,3F12.8)')' Atom Position ',AO_basis(ibasis)%pos
     !WRITE(6,*)'number gauss',AO_basis(ibasis)%ngauss
     !WRITE(6,*)'l- and m-  ',AO_basis(ibasis)%l,AO_basis(ibasis)%m
     !WRITE(6,*)' Exponent       Coefficient    Norm'
     !DO igauss=1,AO_basis(ibasis)%ngauss
     !   WRITE(6,'(3E15.7)')AO_basis(ibasis)%alpha(igauss),AO_basis(ibasis)%coeff(igauss),AO_basis(ibasis)%norm(igauss)
     !ENDDO
     !WRITE(6,*)'Num. Cartesian',AO_basis(ibasis)%ncart
     !WRITE(6,*)'Cart. Coeff and Exponents'
     !DO icart=1,AO_basis(ibasis)%ncart
     !   WRITE(6,*)AO_basis(ibasis)%cart_coeff(icart),AO_basis(ibasis)%cart_mat(icart,:)
     !ENDDO
     !WRITE(6,*)
     !WRITE(6,*)

  ENDDO


  !ALLOCATE(bloch_s(nbasis,nbasis,nk))
  !CALL bloch_space_overlap(AO_basis,index_l,bloch_s,latt_vec)

  CALL renorm_AO(AO_basis)


  !Arrays for information about system to printed in NBO.out header
  ALLOCATE(iatnum(natom))
  ALLOCATE(symbols(natom))
  ALLOCATE(ibasismap(natom+1))
  ALLOCATE(ishellmap(nbasis))
  ALLOCATE(ilmap(nbasis))
  ALLOCATE(iatval(natom))

  READ(12,*)iatnum
  READ(12,*)symbols
  READ(12,*)ibasismap
  READ(12,*)ishellmap

  !Assume full electron calculation
  DO iatom=1,natom
     iatval(iatom) = REAL(iatnum(iatom))
  ENDDO
  !WRITE(6,*)iatval,'iatval'

  DO ibasis=1,nbasis
     ilmap(ibasis) = AO_basis(ibasis)%l
  ENDDO
  !WRITE(6,*)ilmap,'ilmap'

  !Read the geometric information of the system
  DO i=1,3
     READ(12,*)latt_vec(i,:)
  ENDDO
  latt_vec = latt_vec/bohr

  DO i=1,3
     WRITE(6,*)latt_vec(i,:)
  ENDDO
  WRITE(6,*)

  ALLOCATE(atom_pos(natom,3))
  DO iatom=1,natom
     READ(12,*)atom_pos(iatom,:)
  ENDDO
  atom_pos = atom_pos/bohr

  DO iatom=1,natom
     WRITE(6,*)atom_pos(iatom,:)
  ENDDO
  WRITE(6,*)


  !End of parsed output from standard CRYSTAL output
  CLOSE(12)


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Start analyzing output from CRYAPI_OUT call in properties
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  CALL GETARG(2,buffer)
  OPEN(23,file=buffer)

  READ(23,*)k_mesh,k_test

  IF( k_test .NE. ibz_nk )THEN
     WRITE(6,*)'The number of kpoints in the two input files does not match',ibz_nk,k_test
     STOP
  ENDIF

  !Calculate how many k-points there are in half of BZ
  !I am assuming that an odd number of k-points has been used.
  nk = (PRODUCT(k_mesh)+1)/2  
  WRITE(6,*)'nk',nk

  !IF( (PRODUCT(k_mesh)+1)/2 .NE. nk )THEN
  !   WRITE(6,*)'More than inversion symm was used for IBZ. Code can not handle this, yet.'
  !   WRITE(6,*)(PRODUCT(k_mesh)+1)/2,nk
  !   STOP
  !ENDIF


  !These are reciprocal lattice vectors.  I dont think I need them but need to still skip them
  DO i=1,3
     READ(23,*)
  ENDDO

  ALLOCATE(k_grid(ibz_nk,3))
  CALL temp_read(23,k_grid)


  !WRITE(6,*)'grid k points'
  !DO ik=1,ibz_nk
  !   WRITE(6,*)k_grid(ik,:)
  !ENDDO
  !WRITE(6,*)



  ALLOCATE(k_type(nk),sym_type(ibz_nk))
  READ(23,*)sym_type
  !I think if an odd mesh is used all k-points will be complex besides gamma
  !I should check it and stop if this is not the case.
  DO ik=2,ibz_nk
     IF( sym_type(ik).EQ.1 )STOP 'Code assumes an odd k-point grid in each direction. Change ISP.'
  ENDDO


  !Skip over symmetry opperators
  !Note these are 48 3X3 matrices written in rows of 8 => 54 rows
  DO i=1,54
     READ(23,*)
  ENDDO

  ALLOCATE(sym_weight(ibz_nk),k_weight(nk),sym_count(ibz_nk))
  READ(23,*)sym_weight

  !WRITE(6,*)'number of symm equivalent'
  k_test = 1
  sym_count(1) = 1
  IF( nk.GT.1 )THEN
     DO ik=2,ibz_nk
        sym_count(ik)=ANINT(sym_weight(ik)/sym_weight(1)/2.d0)
        k_test = k_test + sym_count(ik)
        !WRITE(6,*)ik,sym_weight(ik),sym_count(ik)
     ENDDO
  ENDIF



  IF( k_test.NE.nk )STOP 'Sym_count does not account for the correct number of symmertry redundant k-points'


  ALLOCATE(sym_eig(ibz_nk,nbasis),band_eig(nk,nbasis))
  CALL temp_read(23,sym_eig)
  sym_eig = sym_eig*hartree

  !!Now assign eigenvalues for all the symmetry equivalent bands
  !jk=0
  !DO ik=1,ibz_nk
  !   DO isym=1,sym_count(ik)
  !      jk=jk+1
  !      band_eig(jk,:)=sym_eig(ik,:)
  !   ENDDO
  !ENDDO


  ALLOCATE(sym_occ(ibz_nk,nbasis),band_occ(nk,nbasis))
  CALL temp_read(23,sym_occ)

  CRYS_nelec = 0.d0
  DO ik=1,ibz_nk
     CRYS_nelec = CRYS_nelec + SUM(sym_occ(ik,:))
     sym_occ(ik,:) = sym_occ(ik,:) / sym_weight(ik)
     !WRITE(6,*)ik
     !WRITE(6,*)sym_eig(ik,:)
     !WRITE(6,*)sym_occ(ik,:)
     !WRITE(6,*)
  ENDDO
  !WRITE(6,*)'CRYSTAL output number of electrons',CRYS_nelec










!  !Dump values for manual DOS calculation
!  OPEN(12,file='sym_eigs')
!
!  WRITE(12,'(I5,A)')ibz_nk, '  #nk'
!  WRITE(12,'(I5,A)')natom, '  #natom'
!  WRITE(12,'(I5,A)')nbasis,'  #nband'
!  WRITE(12,*)
!     
!  DO ik=1,ibz_nk
!     WRITE(12,'(D16.9,F6.1)')sym_weight(ik),ANINT(sym_weight(ik)/sym_weight(1))
!  ENDDO
!  WRITE(12,*)  
!
!  DO ik=1,ibz_nk
!     DO ibasis=1,nbasis
!        WRITE(12,*)sym_eig(ik,ibasis),2.d0*sym_occ(ik,ibasis)
!     ENDDO
!     WRITE(12,*)
!  ENDDO         
!  CLOSE(12)












  !Now transfer eigen values and occupancies to all symmetrically equivalent bands
  jk=0
  DO ik=1,ibz_nk
     DO isym=1,sym_count(ik)
        jk=jk+1
        band_eig(jk,:)=sym_eig(ik,:)
        band_occ(jk,:)=sym_occ(ik,:)
        k_type(jk) = sym_type(ik)
     ENDDO
  ENDDO


  !Now need to read in the coefficients.
  !This is the belly of the beast but I think I have set it up to be fairly easy.


  ALLOCATE(band_coeff(nk,nbasis,nbasis))
  band_coeff = 0.d0

  ALLOCATE(bz_k(nk,3))

  jk = 0
  DO j=1,PRODUCT(k_mesh)  !CRYSTAL write out coefficients for ALL kpoints in BZ, regardless of symmetry.  I'll decide which ones to take.

     READ(23,*)grid_test
     !WRITE(6,*)'grid_test',grid_test

     !Here I am checking to see if I should read in the eigenvector (dont want inverse kpoints) 
     !This code is general to any sort of symmetry of the IBZ
     k_fill = .TRUE.
     DO i=1,3
       grid_neg(i) = MODULO(-grid_test(i),k_mesh(i)) !This is the same as negative and is used since CRYSTAL does BZ from 0->1
     ENDDO
     !WRITE(6,*)'grid_neg',grid_neg

     !IF( SUM(ABS(grid_neg)) .NE. 0 )THEN  !Gamma point's negative will still be there but we want to calcualte it
     IF( jk .NE. 0 )THEN
        DO ik=1,jk  !Loop over all read in kpoints.   !!!!Loop over all previously read in BZ eigenvectors 
           DO i=1,3
              IF( grid_neg(i) .NE. bz_k(ik,i) )EXIT !Check if negative of the kpoint is in list 
              IF( i .EQ. 3 )k_fill = .FALSE.
           ENDDO
           IF( .NOT. k_fill )EXIT
        ENDDO
     ENDIF

     IF( k_fill )THEN  !Actually read in the coefficients if this is an eigenvector I care about
        !WRITE(6,*)'will read in these coeff',grid_test
        jk = jk + 1
        bz_k(jk,:)=grid_test

        !WRITE(6,*)'jk   ',jk

        IF( k_type(jk) .EQ. 1 )THEN  !Real type read in
           !WRITE(6,*)'real type read in'
           CALL temp_read(23,band_coeff(jk,:,:),1)
        ELSEIF( k_type(jk) .EQ. 0 )THEN  !Complex type read in
           !WRITE(6,*)'comp type read in'
           CALL temp_read(23,band_coeff(jk,:,:),2)
        ELSE
           WRITE(6,*)'improper k_type when reading eigenvectors',jk,k_type(jk)
           STOP
        ENDIF

        !WRITE(6,*)'completed reading in ',jk

     ELSE  !Even if I do not want this eigenvector, I still need to skip the lines in the file
        !WRITE(6,*)'not reading in these coefficients',grid_test

        skip = nbasis*nbasis
        IF( k_type(jk) .EQ. 0 )skip = skip*2
        IF( MODULO(skip,4) .NE. 0 )skip = skip + MODULO(-skip,4)
        skip = skip / 4

        !WRITE(6,*)skip

        DO i=1,skip
          READ(23,*)
        ENDDO

     ENDIF

     !WRITE(6,*)


  ENDDO

  ALLOCATE(kpt(nk,3))
  DO ik=1,nk
     DO i=1,3
        kpt(ik,i) = bz_k(ik,i) / DBLE(k_mesh(i))
        IF( kpt(ik,i) .GT. 0.5d0 )kpt(ik,i) = kpt(ik,i) - 1.d0
     ENDDO
     !WRITE(6,*)ik,bz_k(ik,:)
     !WRITE(6,*)SNGL(kpt(ik,:))
     !WRITE(6,*)
  ENDDO


  !Now that I have read in the coefficients I will need to 
  ! 1. rearrange them for proper ordering - definitely need to 
  ! 2. Maybe rescale for normalization - Just test band overlap matrix first

  !Rearrange some coefficients

  ALLOCATE(r_mat(nbasis,nbasis))

  !d_rearr = 0
  !DO i=1,5
  !   IF( MODULO(i,2) .EQ. 0 )THEN
  !      d_rearr(i) = i
  !   ELSE
  !      d_rearr(i) = MODULO(i+2,6)
  !   ENDIF
  !ENDDO
  !WRITE(6,*)d_rearr

  !rearr_index is used to rearrange the d- and f- coeffs from CRYSTAL ordering to my ordering
  !The index of the rearr_index is my ordering, 
  !the contents of the array is CRYSTAL ordering.
  ALLOCATE(rearr_index(nbasis))
  DO ibasis=1,nbasis

     IF( AO_basis(ibasis)%l .EQ. 2 )THEN

        IF( MODULO(AO_basis(ibasis)%m,2) .EQ. 0 )THEN
           rearr_index(ibasis) = ibasis
        ELSEIF( AO_basis(ibasis)%m .LT. 5 )THEN
           rearr_index(ibasis) = ibasis + 2
        ELSE
           rearr_index(ibasis) = ibasis - 4
        ENDIF

     ELSEIF( AO_basis(ibasis)%l .EQ. 3)THEN 

        IF( AO_basis(ibasis)%m .LE. 4 )THEN
           IF( MODULO(AO_basis(ibasis)%m,2) .EQ. 1 )THEN
              rearr_index(ibasis) = ibasis + 4
           ELSE
              rearr_index(ibasis) = ibasis + 2
           ENDIF
        ELSEIF( AO_basis(ibasis)%m .LE. 6 )THEN
           rearr_index(ibasis) = ibasis - 3
        ELSE
           rearr_index(ibasis) = ibasis - 6
        ENDIF

     ELSE
        rearr_index(ibasis) = ibasis
     ENDIF

  ENDDO

  !WRITE(6,*)rearr_index

  !Now, actually  rearrange the coefficients.
  DO ik=1,nk
     r_mat = band_coeff(ik,:,:) !Temporarily store coeffs
     DO ibasis=1,nbasis
        IF( rearr_index(ibasis) .NE. ibasis )THEN  !Don't need to reassign every set of coeff's
           !WRITE(6,*)ibasis,rearr_index(ibasis)
           band_coeff(ik,:,ibasis) = r_mat(:,rearr_index(ibasis))
        ENDIF
     ENDDO
     !WRITE(6,*)
  ENDDO


  !DO ik=1,nk
  !   WRITE(6,*)'for kpoint',ik
  !   DO ibasis=1,nbasis
  !      WRITE(6,*)SNGL(REAL(band_coeff(ik,ibasis,:)))
  !   ENDDO
  !   WRITE(6,*)
  !   DO ibasis=1,nbasis
  !      WRITE(6,*)SNGL(AIMAG(band_coeff(ik,ibasis,:)))
  !   ENDDO
  !   WRITE(6,*)
  !ENDDO

  IF( jk .NE. (PRODUCT(k_mesh)+1)/2 )STOP 'Incorrect number of eigenvectors read for only time reversal symm in IBZ'

  


  !Calculate bloch space overlap matrices
  ALLOCATE(bloch_s(nbasis,nbasis,nk))
  CALL bloch_space_overlap(AO_basis,index_l,bloch_s,latt_vec,kpt)


  k_weight = 2.d0 / PRODUCT(k_mesh)
  k_weight(1) = 1.d0 / PRODUCT(k_mesh)






!  !Dump values for manual DOS calculation
!  OPEN(12,file='bz_eigs')
!
!  WRITE(12,'(I5,A)')nk, '  #nk'
!  WRITE(12,'(I5,A)')natom, '  #natom'
!  WRITE(12,'(I5,A)')nbasis,'  #nband'
!  WRITE(12,*)
!
!  DO ik=1,nk
!     WRITE(12,'(D16.9,F6.1)')k_weight(ik),ANINT(k_weight(ik)/k_weight(1))
!  ENDDO
!  WRITE(12,*)
!
!  DO ik=1,nk
!     DO ibasis=1,nbasis
!        WRITE(12,*)band_eig(ik,ibasis),2.d0*band_occ(ik,ibasis)
!     ENDDO
!     WRITE(12,*)
!  ENDDO
!
!  CLOSE(12)






  !DO ik=1,nk
  !   WRITE(6,*)'overlap for kpt',ik,SNGL(kpt(ik,:))
  !   DO ibasis=1,nbasis
  !      WRITE(6,'(28F6.3)')REAL(bloch_s(ibasis,:,ik))
  !   ENDDO
  !   WRITE(6,*)
  !   DO ibasis=1,nbasis
  !      WRITE(6,'(28F6.3)')AIMAG(bloch_s(ibasis,:,ik))
  !   ENDDO
  !   WRITE(6,*)
  !ENDDO

  ALLOCATE(band_dummy(nbasis,nbasis))

  !Check orthonormality of bands
  DO ik=1,nk

     !WRITE(6,*)'starting bands',ik
     !DO ibasis=1,nbasis
     !   WRITE(6,'(28F6.3)')REAL(band_coeff(ik,ibasis,:))
     !ENDDO

     r_mat = MATMUL(CONJG(band_coeff(ik,:,:)),MATMUL(bloch_s(:,:,ik),TRANSPOSE(band_coeff(ik,:,:))))
     band_dummy = TRANSPOSE(band_coeff(ik,:,:))
     CALL do_OWSO(r_mat,band_occ(ik,:),band_dummy(:,:))

     band_coeff(ik,:,:) = TRANSPOSE(band_dummy)

     !WRITE(6,*)'orthog bands',ik
     !DO ibasis=1,nbasis
     !   WRITE(6,'(28F6.3)')REAL(band_coeff(ik,ibasis,:))
     !ENDDO
     !WRITE(6,*)

     !r_mat = MATMUL(CONJG(band_coeff(ik,:,:)),MATMUL(bloch_s(:,:,ik),TRANSPOSE(band_coeff(ik,:,:))))
     !WRITE(6,*)'band overlap for kpoint',ik
     !DO ibasis=1,nbasis
     !   WRITE(6,'(28F6.3)')REAL(r_mat(ibasis,:))
     !ENDDO
     !WRITE(6,*)
  ENDDO


  !Now lets calculate some density matrices
  ALLOCATE(bloch_rho(nbasis,nbasis,nk))
  ALLOCATE(bloch_fock(nbasis,nbasis,nk))


  ALLOCATE(nelec_mat(nbasis,nbasis))

  nelec=0.d0
  DO ik=1,nk
     !Set up an occupancy matrix
     r_mat = 0.d0
     DO ibasis=1,nbasis
        r_mat(ibasis,ibasis) = band_occ(ik,ibasis)
     ENDDO
     !Use matrix multiplication
     bloch_rho(:,:,ik) = MATMUL(TRANSPOSE(band_coeff(ik,:,:)),MATMUL(r_mat,CONJG(band_coeff(ik,:,:))))


     !WRITE(6,*)'density matrix for kpoint',ik
     !DO ibasis=1,nbasis
     !   WRITE(6,'(28F6.3)')REAL(bloch_rho(ibasis,:,ik))
     !ENDDO
     !WRITE(6,*)

     !Tests number of electrons
     nelec_mat = MATMUL(bloch_rho(:,:,ik),bloch_s(:,:,ik))
     nelec_mat = nelec_mat*k_weight(ik)
     DO ibasis=1,nbasis
        nelec = nelec + nelec_mat(ibasis,ibasis)
     ENDDO
  ENDDO
  WRITE(6,*)'total number of electrons'
  WRITE(6,*)2.d0*nelec,2.d0*CRYS_nelec


  !Now calculate the fock matrices
  !block_fock= 0.d0
  nelec=0.d0
  DO ik=1,nk
     !Set up diagonal energy matrix
     r_mat = 0.d0
     DO ibasis=1,nbasis
        r_mat(ibasis,ibasis) = band_eig(ik,ibasis)
     ENDDO

     nelec_mat = MATMUL(bloch_s(:,:,ik),TRANSPOSE(band_coeff(ik,:,:)))

     bloch_fock(:,:,ik) = MATMUL(nelec_mat,MATMUL(r_mat,TRANSPOSE(CONJG(nelec_mat))))

     !Test going back to the band basis
     nelec_mat = MATMUL(CONJG(band_coeff(ik,:,:)),MATMUL(bloch_fock(:,:,ik),TRANSPOSE(band_coeff(ik,:,:))))
     !WRITE(6,*)'reformed band basis fock kpoint',ik
     DO ibasis=1,nbasis
        !WRITE(6,'(28F6.3)')AIMAG(nelec_mat(ibasis,:))
        IF( ABS(REAL(nelec_mat(ibasis,ibasis)) - band_eig(ik,ibasis)) .GT. 1.d-6 )THEN
           WRITE(6,*)'bad reformed band energy',ik,ibasis
           WRITE(6,'(3F15.8)')REAL(nelec_mat(ibasis,ibasis)),band_eig(ik,ibasis)
        ENDIF
     ENDDO
     !WRITE(6,*)

  ENDDO

  !Now write out information to be read in for NBO analysis
  OPEN (55,file='NBO.out')
  !Begin by writing out some general info about the system
  WRITE(55,*)'comment line, output from processing CRYSTAL output'
  WRITE(55,'(I8,A7)')natom,'#natom'
  WRITE(55,'(I8,A8)')nbasis,'#nbasis'
  WRITE(55,'(I8,A8)')1,'#nspins'

  !Not all l-vectors are used, only those corresponding to neighboring unit  cells 
  WRITE(55,'(I8,A4)')125,'#ng'
  WRITE(55,'(I8,A4)')nk,'#nk'

  WRITE(55,*)
  
  !Now write out some info about basis function ordering and info about each basis function (ie. set of exponents and angular momentum)
  !The index of the first basis function on each atom.
  !There is an additional number written out for the next index AFTER the end of the basis set.
  DO j=1,natom+1
     WRITE(55,'(I5)', ADVANCE='no')ibasismap(j)
  ENDDO
  WRITE(55,'(A12)',ADVANCE='NO')'#ibasismap'
  WRITE(55,*)

  !Indices are the same for atoms with the same exponents.
  !This counts up over different atoms.
  DO j=1,nbasis
     WRITE(55,'(I5)',ADVANCE='NO')ishellmap(j)
  ENDDO
  WRITE(55,'(A11)',ADVANCE='NO')'#ishellmap'
  WRITE(55,*)

  !l-quantum number of each basis set
  DO j=1,nbasis
     WRITE(55,'(I2)',ADVANCE='NO')ilmap(j)
  ENDDO
  WRITE(55,'(A7)',ADVANCE='NO')'#ilmap'
  WRITE(55,*)

  !Then some info about the atoms present in the system (symbol and atomic number)
  !Atomic symbol of each atom.
  DO j=1,natom
     WRITE(55,'(A3)',ADVANCE='NO')symbols(j)
  ENDDO
  WRITE(55,'(A8)',ADVANCE='NO')'#symbols'
  WRITE(55,*)

  !Atomic number of each atom.
  DO j=1,natom
     WRITE(55,'(I3)',ADVANCE='NO')iatnum(j)
  ENDDO
  WRITE(55,'(A8)',ADVANCE='NO')'#iatnum'
  WRITE(55,*)

  !Number of electrons fully represented by the pseudopotential for each atom.
  DO j=1,natom
     WRITE(55,'(F6.1)',ADVANCE='NO')iatval(j)
  ENDDO
  WRITE(55,'(A8)',ADVANCE='NO')'#iatval'
  WRITE(55,*)

  WRITE(55,*)

  !Write out structural information for the system

  !First write out the lattice vectors
  DO j=1,3
     WRITE(55,'(3F20.9)')latt_vec(j,:)
  ENDDO
  WRITE(55,*)

  !Next is the atomic positions
  DO j=1,natom
     WRITE(55,'(3F20.9)')atom_pos(j,:)
  ENDDO
  WRITE(55,*)


  !Write out the l-vectors
  !The ordering of l_vectors expected by the NBO code is different than that used here, so they will be written out non-linearly
  !The first must be (/0,0,0/)
  !Then they can proceed in any order but must be printed in pairs of negatives.
  !Only those of neighboring cells will be included as those are all that will be used in the NBO analysis in searching for spanning bonds
  l_half = (SIZE(index_l,2)+1)/2
  WRITE(55,'(3I3)')index_l(:,l_half)
  !IF( .NOT. gamma_point )THEN
     DO j=1,l_half-1
        !IF( surface .AND. index_l(3,l_half-j) .NE. 0 )GOTO 55  !I am assuming the c-lattice vector is the surface norm
        !IF( MAXVAL(ABS(index_l(:,l_half-j))) < 2 )THEN
        IF( MAXVAL(ABS(index_l(:,l_half-j))) < 3 )THEN
           WRITE(55,'(3I3)')index_l(:,l_half - j)
           WRITE(55,'(3I3)')index_l(:,l_half + j)
        ENDIF
55      CONTINUE
     ENDDO
  !ENDIF

  WRITE(55,*)

  !Write out the k-vectors
  !The ordering will be that output from VASP.  Gamma point is always first.
  !The density and overlap matrices will be written out in this same order.
  !The weight of each k-point is also included. Unless hihger symmetry is included, this is 1/nkpts for the Gamma point and 2/nkpts for all other.
  DO ik=1,nk
     WRITE(55,'(4F18.10)')kpt(ik,:),k_weight(ik)
  ENDDO
  WRITE(55,*)

  WRITE(55,*)'NBO_mat.out'
  WRITE(55,*)

  !Basis set information
  WRITE(55,*)'Basis set information. Used to produce visualtion output file.'
  WRITE(55,*)
  DO ibasis=1,nbasis
     WRITE(55,'(I)')AO_basis(ibasis)%ngauss
     WRITE(55,'(10D17.9)')AO_basis(ibasis)%alpha
     WRITE(55,'(10D17.9)')AO_basis(ibasis)%coeff
     WRITE(55,'(10D17.9)')AO_basis(ibasis)%norm
     WRITE(55,'(3D16.8)')AO_basis(ibasis)%pos
     WRITE(55,'(I)')AO_basis(ibasis)%ncart
     DO j=1,AO_basis(ibasis)%ncart
        WRITE(55,'(D17.9,3I)')AO_basis(ibasis)%cart_coeff(j),AO_basis(ibasis)%cart_mat(j,:)
     ENDDO

     WRITE(55,*)
  ENDDO



  CLOSE(55)

  OPEN(66,file='NBO_mat.out',FORM='UNFORMATTED')

  !Varios Bloch-space matrices are now written out for the bands in the projected AO basis.
  !All of these are written out in a lower triangular form
  !The ordering 1)overlap, 2)density, 3)fock is expected by the NBO code
  CALL write_sym_matrices(66,bloch_s)

  CALL write_sym_matrices(66,bloch_rho) !,ispin))

  CALL write_sym_matrices(66,bloch_fock) !,ispin))


  CLOSE(66)






END PROGRAM crystal_nbo
