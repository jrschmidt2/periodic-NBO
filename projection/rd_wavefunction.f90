MODULE rd_wavefunction
  IMPLICIT NONE

  !Information about the system: number of ions, number of ion types, and type of each ion
  !Also ion positions
  INTEGER :: n_atom, ntyp
  INTEGER,ALLOCATABLE :: atom_types(:),itypes(:)
  REAL*8,ALLOCATABLE :: atoms(:,:)
  !Maximum number of kpts, bands, plane waves, PAW augmenters, and spins
  INTEGER :: nkpts, nbands, nplmax, npromax, nspins
  !The energy cutoff
  REAL*8 :: enmax
  !k-points (scaled), eigenvalues, and Fermi weights
  REAL*8,ALLOCATABLE :: kpt(:,:),eig(:,:,:),weight(:,:,:)
  REAL*8,ALLOCATABLE :: kpt_wt(:)
  !The x,y,z index of each PW component; the actual g-vector
  !is calculated below
  INTEGER,ALLOCATABLE :: igx(:,:),igy(:,:),igz(:,:)
  !The plane wave coefficients
  !COMPLEX*8,ALLOCATABLE :: cw(:,:,:), cproj(:,:,:)
  COMPLEX*16,ALLOCATABLE ::  pw_coeff(:,:,:,:),PAW_coeff(:,:,:,:)
  !REAL*8, ALLOCATABLE :: g(:,:,:) ! Array containing each g vector for each plane wave at each k-point
  REAL*8, ALLOCATABLE :: gk(:,:,:) ! Array containing each g vector for each plane wave at each k-point, has the k point added in
  INTEGER,ALLOCATABLE :: npl(:) !The number of plane waves at a give k-point
  LOGICAL             :: gamma_point  !tells whether or not the vasp output is from a gamma point only calculation

  !*************************************************************************************
  !Information about the pseudo-potential (for PAW calculations)
  TYPE pseudo
     !This includes the number of l channels, and the number used for the PAW on site terms
     INTEGER :: ldim, lmaxpaw
     !and their associated l quantum numbers
     INTEGER,ALLOCATABLE :: lps(:)
     !And also information about the radial grid, including number of points, radial values
     !and integration weights
     INTEGER :: nmax
     REAL*8,ALLOCATABLE :: r(:),si(:)
     REAL*8 :: rend
     !Finall, the (valence) all-electron and pseudo-wavefunctions, and differences
     REAL*8,ALLOCATABLE :: wae(:,:),wps(:,:),wdiff(:,:)
     !An array the will hold the value of the each wdiff channel at r=0. This value is extrapolated using the 1st and 2nd radial points
     REAL*8,ALLOCATABLE ::  center_value(:,:)
     !Number of valence electrons that are treated completely by the pseudopotential
     REAL*8  ::  val_corr, val_corr_test  !Since VASP has two variables to hold this quantity I will read both of them in and make sure they are the same
  END TYPE pseudo
  TYPE(pseudo),ALLOCATABLE :: P(:)

  LOGICAL   ::  PAW_pseudo

  REAL*8,ALLOCATABLE   ::  PAW_pos(:,:)

  CHARACTER(128)   ::  wavefunction_fn

  REAL*8    ::  a(3,3),b(3,3)   ! Direct and Reciprocal lattice vectors 
  REAL*8    ::  pwnorm   ! Normalization factor for plane waves (inverse of square root of volume)
  REAL*8    ::  unitvol  ! volume of the unit cell, to use for plane wave norm

  REAL*8, ALLOCATABLE  ::  pseudo_norm(:)    !Norm of each of band only taking into account the planewave portion
  REAL*8, ALLOCATABLE  ::  aug_norm(:)

  INTEGER  ::  kdim(3)  !Dimensionality of k-point mesh. Global so it can be used when writing the NBO.out file.

  REAL*8,PARAMETER ::  eV = 27.211652d0


  CONTAINS

    SUBROUTINE read_vasp_wavefunction(VASP_fn)

      CHARACTER(128) ::  VASP_fn

      COMPLEX*8,ALLOCATABLE :: cw(:,:,:,:), cproj(:,:,:,:)

      INTEGER ::  unity(3)=1.d0 ! Vector of ones; 'mesh' for total lattice volume
      INTEGER ::  PAW_tally, pro_per_atom
      REAL*8  ::  energy_sum, occ_sum
      INTEGER ::  iion, ityp, ikpt, iband, ipl, nplread, ipro, ispin, j, k, l 
      !INTEGER :: l,ll,m

      REAL*8  ::  kgrid(3) !For determing info on the k-point mesh

      REAL*8, PARAMETER  ::   bohr = 0.529177249d0, pi=4.d0*ATAN(1.d0), hartree=27.21138386

      OPEN(10,FILE=VASP_fn,FORM='UNFORMATTED')

      PRINT *, '*** Reading system ***'
      PRINT *
      READ(10) n_atom, ntyp
      ALLOCATE(atom_types(ntyp),atoms(3,n_atom),P(ntyp),itypes(n_atom))
      READ(10) itypes
      READ(10) atoms

      atom_types = 0
      DO j=1,ntyp
         DO k=1,n_atom
            IF( itypes(k) == j )THEN
               atom_types(j) = atom_types(j) + 1
            ENDIF
         ENDDO
         IF( atom_types(j) == 0 )THEN
            WRITE(6,*)'atom_types was not allocated propetly'
            STOP
         ENDIF
      ENDDO

      PRINT *, 'n_atom: ', n_atom
      PRINT *, 'ntyp: ', ntyp
      PRINT *, 'atom_types: ', atom_types
      PRINT *

      PRINT *, '*** Reading wavefunction *** '
      PRINT *

      READ(10) nkpts,nbands,enmax,nplmax,npromax,nspins
      PRINT *, 'nkpts:  ', nkpts
      PRINT *, 'nbands: ', nbands
      PRINT *, 'enmax:  ', enmax, 'eV'
      PRINT *, 'nplmax: ', nplmax
      PRINT *, 'npromax:', npromax
      PRINT *, 'nspins: ', nspins
      PRINT *

      gamma_point = .FALSE.
      IF( nkpts == 1 )gamma_point = .TRUE.

      ALLOCATE(kpt(3,nkpts),eig(nbands,nkpts,nspins),weight(nbands,nkpts,nspins))
      ALLOCATE(kpt_wt(nkpts))
      ALLOCATE(igx(nplmax,nkpts),igy(nplmax,nkpts),igz(nplmax,nkpts))
      ALLOCATE(cw(nplmax,nbands,nkpts,nspins),cproj(npromax,nbands,nkpts,nspins))
      ALLOCATE(pw_coeff(nplmax,nbands,nkpts,nspins),PAW_coeff(npromax,nbands,nkpts,nspins))
      ALLOCATE(npl(nkpts))
      ALLOCATE(gk(3, nplmax, nkpts))

      READ(10) a
      !VASP uses angstroms for everything, but the projection code uses bohr, so the unit cell vectors are covnerted here
      !All other position quantities are calculated from these vectors, so now all distance will be in bohr
      a = a / bohr

      PRINT *, 'Direct lattice vectors:'
      PRINT *, SNGL(a(1:3,1))
      PRINT *, SNGL(a(1:3,2))
      PRINT *, SNGL(a(1:3,3))
      PRINT *

      !Convert ion positions to absolute coordinates from fractional
      DO iion=1,n_atom
         atoms(:,iion)=atoms(1,iion)*a(:,1)+atoms(2,iion)*a(:,2)+atoms(3,iion)*a(:,3)
      ENDDO

      CALL calc_recip_lattice(a,b)
      PRINT *, 'Reciprical lattice vectors:'
      PRINT *, SNGL(b(1:3,1))
      PRINT *, SNGL(b(1:3,2))
      PRINT *, SNGL(b(1:3,3))
      PRINT *

      !Information on the k-point mesh and plane-wave basis at each k-point are then read in
      DO ikpt=1,nkpts
         READ(10) npl(ikpt), kpt(1:3,ikpt), kpt_wt(ikpt)
         READ(10) (igx(ipl,ikpt),igy(ipl,ikpt),igz(ipl,ikpt),ipl=1,npl(ikpt))
         !Now the g-vector of each plane-wave is calculated
         DO ipl=1,npl(ikpt)
            gk(:,ipl,ikpt) = igx(ipl,ikpt)*b(:,1)+igy(ipl,ikpt)*b(:,2)+igz(ipl,ikpt)*b(:,3)
            gk(:,ipl,ikpt) = gk(:,ipl,ikpt) +kpt(1,ikpt)*b(:,1)+kpt(2,ikpt)*b(:,2)+kpt(3,ikpt)*b(:,3)
         ENDDO
         !WRITE(6,*)'k-point',ikpt,SNGL(kpt(:,ikpt))
         !WRITE(6,*)'k-point weight',SNGL(kpt_wt(ikpt))
         !WRITE(6,*)
      ENDDO

      !Now analyze the k-point mesh as obtained out of VASP
      !We want to see if there are any dimensions with only a single k-point
      !This will be used for adjusting the unit cells written out to the NBO.out file
      kgrid = 1
      DO ikpt=1,nkpts
         DO j=1,3
            IF( ABS(kpt(j,ikpt)).GT.1.d-13 )THEN
               !WRITE(6,*)j,(1.d0/ABS(kpt(j,ik)))
               IF( (1.d0/ABS(kpt(j,ikpt))).GT.kgrid(j) )THEN
                  kgrid(j) = (1.d0/ABS(kpt(j,ikpt)))
               ENDIF
            ENDIF
         ENDDO
      ENDDO
      DO j=1,3
         !WRITE(6,*)kgrid(j)-NINT(kgrid(j))
         IF( ABS(kgrid(j)-NINT(kgrid(j))).GT.1.d-10 )THEN
            WRITE(6,*)'kgrid does not contain a near integer in the direction',j
            WRITE(6,*)kgrid(j)
            STOP
          ENDIF
         kdim(j) = NINT(kgrid(j))
      ENDDO

      WRITE(6,*)'Dimemsionality of k-point mesh'
      WRITE(6,'(3I8)')kdim
      WRITE(6,*)

      !Spin dependent information is now read in
      !This includes band energies and occupancies for each band at each k-point 
      !As well as coefficients for plane-waves and PAW augmenters for each band at each k-point
      DO ispin=1,nspins
         DO ikpt=1,nkpts
            READ(10) (eig(iband,ikpt,ispin), weight(iband,ikpt,ispin),iband=1,nbands)
            DO iband=1,nbands
               READ(10) (cw(ipl,iband,ikpt,ispin),ipl=1,npl(ikpt))
               READ(10) (cproj(ipro,iband,ikpt,ispin),ipro=1,npromax)
            ENDDO
         ENDDO
      ENDDO

      !The energy values in VASP are in eV, but we want them in Hartree, so they are converted here
      eig = eig / hartree

      !The number of electrons and energy expectation value for the plane-wave bands is calculated.
      !This is done to compare projected results
      energy_sum = 0.d0
      occ_sum = 0.d0
      DO ispin=1,nspins
      DO ikpt=1,nkpts
         DO iband=1,nbands
            occ_sum = occ_sum + weight(iband,ikpt,ispin)*kpt_wt(ikpt)
            energy_sum = energy_sum + eig(iband,ikpt,ispin)*weight(iband,ikpt,ispin)*kpt_wt(ikpt)
         ENDDO
      ENDDO
      ENDDO
      WRITE(6,*)'Number of electrons in VASP bands ',occ_sum * DBLE(3-nspins)
      WRITE(6,*)'Energy sum for VASP bands         ',energy_sum * DBLE(3-nspins)
      WRITE(6,*)

      !Initialize logical variable assuming this is NOT a PAW pseudopotential
      PAW_pseudo = .FALSE.

      PRINT *, '*** Reading pseudo-potential ***'
      PRINT *
      DO iion=1,ntyp

         !WRITE(6,*)'Pseudotential for ion type:',iion
         READ(10) P(iion)%ldim, P(iion)%lmaxpaw
         ALLOCATE(P(iion)%lps(P(iion)%ldim))
         READ(10) P(iion)%lps
         READ(10) P(iion)%val_corr,P(iion)%val_corr_test
         READ(10) P(iion)%nmax

         !WRITE(6,*)'ldim',P(iion)%ldim
         !WRITE(6,*)'lmaxpaw', P(iion)%lmaxpaw 
         !WRITE(6,*)'lps', P(iion)%lps
         IF( P(iion)%val_corr .NE. P(iion)%val_corr_test )THEN
            WRITE(6,*)'PVALF and PVALF_ORIG are not the same for ion type',iion
            WRITE(6,*)P(iion)%val_corr, P(iion)%val_corr_test
            STOP
         ELSE
            !WRITE(6,*)'val_corr',P(iion)%val_corr
         ENDIF
         !WRITE(6,*)'nmax', P(iion)%nmax

         ALLOCATE(P(iion)%r(P(iion)%nmax))
         ALLOCATE(P(iion)%si(P(iion)%nmax))
         ALLOCATE(P(iion)%wae(P(iion)%nmax,P(iion)%ldim))
         ALLOCATE(P(iion)%wps(P(iion)%nmax,P(iion)%ldim))
         ALLOCATE(P(iion)%wdiff(P(iion)%nmax,P(iion)%ldim))
         ALLOCATE(P(iion)%center_value(P(iion)%ldim,0:2))

         READ(10) P(iion)%r
         READ(10) P(iion)%si
         READ(10) P(iion)%wae
         READ(10) P(iion)%wps
         IF (P(iion)%nmax > 0) P(iion)%rend=P(iion)%r(P(iion)%nmax)

         PRINT *, 'Ion type: ', iion
         PRINT *, '# l channels: ', P(iion)%ldim
         PRINT *, 'lmax PAW: ', P(iion)%lmaxpaw !Evidently this is zero of non-PAW calcs, non-zero for PAW
         IF (P(iion)%lmaxpaw.GT.0) THEN
            !PRINT *, 'PAW pseudo-potential detected.'
            IF( iion .GT. 1 .AND. .NOT. PAW_pseudo )STOP ' One atom type is not PAW but another is? This aint gonna work.'
            PAW_pseudo = .TRUE.
         ELSE
            PRINT *, 'Found norm-conserving pseduo-potential.'
            IF( PAW_pseudo )STOP ' One atom type is PAW but another is not? This aint gonna work.'
            PAW_pseudo = .FALSE.
         ENDIF
         PRINT *, 'l of channel: ', P(iion)%lps
         PRINT *, '# radial pts:', P(iion)%nmax
         PRINT *

         !Even though VASP stores the AE and PS radial function separately for projection purposes we want the difference, which will be stored in wdiff
         !Also VASP stores the radial functions multiplied by a factor of r, we remove this here.
         DO j=1,P(iion)%nmax
            P(iion)%wdiff(j,:) = (P(iion)%wae(j,:) - P(iion)%wps(j,:)) / P(iion)%r(j)
            P(iion)%wae(j,:) = P(iion)%wae(j,:) / P(iion)%r(j)
            P(iion)%wps(j,:) = P(iion)%wps(j,:) / P(iion)%r(j)
         ENDDO

         !The PAW radial grids must be converted from angstrom to bohr
         P(iion)%r = P(iion)%r / bohr
         P(iion)%si = P(iion)%si / bohr

         !For off-site overlaps the integral of the PAW along the radial grid must be computed
         P(iion)%center_value = 0.d0
         DO k=1,P(iion)%ldim
            DO l=0,2
               DO j=1,P(iion)%nmax
                  P(iion)%center_value(k,l) = P(iion)%center_value(k,l) + P(iion)%wdiff(j,k) * P(iion)%r(j)**(2+l) * P(iion)%si(j)
               ENDDO
               !WRITE(6,*)'Integrated PAW value for channel',k,'of order',l,P(iion)%center_value(k,l)
            ENDDO
         ENDDO
         !WRITE(6,*)

         IF( PAW_pseudo )THEN
            !Adjust the maximum radius of the core region into bohr. This parameter is only used when doing a 3D integration of the core region (density test)
            P(iion)%rend = P(iion)%rend / bohr
            !PRINT *, 'Max radial dist. (bohr): ', SNGL(P(iion)%rend)
            !PRINT *,SNGL(P(iion)%r(P(iion)%nmax))
            !PRINT *
            IF( P(iion)%r(P(iion)%nmax) /= P(iion)%rend )THEN
               WRITE(6,*)'there are problems with bohr angstrom comparison for r-end'
               STOP
            ENDIF
         ENDIF
         !WRITE(6,*)

      ENDDO

      !Now an array containing the coordinate of the atom which contains a given projector function will be filled
      !The array is only filled if the pseudopotential is a PAW
      !All atoms are looped over, how many projectors are on an atom is calcualted and that many spots in the array are filled with the location of the atom
      !PAW_tally is used to placehold in the array for already filled positions.
      IF( PAW_pseudo )THEN
         ALLOCATE(PAW_pos(npromax,3))
         PAW_tally = 0
         DO k=1,n_atom 
            ityp = itypes(k)
            pro_per_atom = SUM(2*P(ityp)%lps+1)
            !WRITE(6,*)'PAW augmenters on atom',k,pro_per_atom
            DO j=PAW_tally+1,PAW_tally+pro_per_atom
               PAW_pos(j,:) = atoms(:,k)
               !WRITE(6,*)'PAW_pos',j
               !WRITE(6,*)PAW_pos(j,:)
            ENDDO
            PAW_tally = PAW_tally + pro_per_atom
         ENDDO

         WRITE(6,*)'Total number of PAW-augmenter functions, per unit cell',PAW_tally
         IF( PAW_tally /= npromax )STOP 'the PAW_tally parameter is not correct coming out of filling PAW_pos array'
         WRITE(6,*)
      ENDIF

      CLOSE(10)

      !WRITE(6,*)'Band normalization test'
      !DO ikpt=1,nkpts
      !   WRITE(6,*)'For kpt',ikpt
      !   DO iband=1,nbands
      !      WRITE(6,*)'Norm for band',iband,SUM(CONJG(cw(:,iband,ikpt,1))*cw(:,iband,ikpt,1))
      !   ENDDO
      !   WRITE(6,*)
      !ENDDO


      !Include planewave norm -> 1/sqrt(volume of unit cell)
      CALL grid_volume(unitvol, a, unity)
      cw = cw *  1.d0/SQRT(unitvol)

      !Convert the coefficient arrays to ones with doulbe precision for use in BLAS subroutines
      PAW_coeff = cproj * bohr**(1.5) !Also need to convert this to inverse bohr from inverse angstroms
      pw_coeff  = cw



ENDSUBROUTINE read_vasp_wavefunction

  !Does a simple linear interpolation to calculat the radial wavefunction at
  !a given radial distance.  Note that VASP actually tabulates the radial
  !wfn times r, so it is necessary to divide by r to get the desired value.
SUBROUTINE get_radial_wavefunction(w,rvals,r,val)
  IMPLICIT NONE
  REAL*8 :: w(:),rvals(:),r,val
  INTEGER :: n,nmax

  nmax=SIZE(w,1)
  val=0.d0
  DO n=1,nmax-1
     IF (rvals(n) <= r .AND. rvals(n+1) > r) THEN
        val=w(n)+(w(n+1)-w(n))/(rvals(n+1)-rvals(n))*(r-rvals(n))
     ENDIF
  ENDDO
  !val=val/r

END SUBROUTINE get_radial_wavefunction


  !********************************************************************
  !********************************************************************
  ! VASP ROUTINES (somewhat modified
  !********************************************************************
  !********************************************************************

  !************************* SETYLM ************************************
  !
  ! calculate spherical harmonics for a set of grid points up to
  ! LMAX.  Note that these are missing a factor of r**(-l) from the "true"
  ! cartesian speherical harmonics!
  ! written by Georg Kresse and updated by JRS
  !*********************************************************************

SUBROUTINE SETYLM(LYDIM,INDMAX,YLM,X,Y,Z)
  IMPLICIT NONE
  INTEGER LYDIM           ! maximum L
  INTEGER INDMAX          ! number of points (X,Y,Z)
  REAL*8 YLM(:,:)        ! spherical harmonics
  REAL*8 X(:),Y(:),Z(:)  ! x,y and z coordinates

  ! local variables
  REAL*8 FAK
  INTEGER IND

  REAL*8, PARAMETER  ::   bohr = 0.529177249d0, pi=4.d0*ATAN(1.d0)

  !-----------------------------------------------------------------------
  ! runtime check of workspace
  !-----------------------------------------------------------------------
  IF ( UBOUND(YLM,2) < (LYDIM+1)**2) THEN
     WRITE(0,*)'internal ERROR: SETYLM, insufficient L workspace'
     STOP
  ENDIF

  IF ( UBOUND(YLM,1) < INDMAX) THEN
     WRITE(0,*)'internal ERROR: SETYLM, insufficient INDMAX workspace'
     STOP
  ENDIF

  FAK=1/(2.d0 * SQRT(PI))
  !-----------------------------------------------------------------------
  ! here is the code for L=0, hard coded
  !-----------------------------------------------------------------------
  IF (LYDIM <0) GOTO 100
  !DIR$ IVDEP
  !OCL NOVREC
  DO IND=1,INDMAX
     YLM(IND,1)=FAK
  ENDDO
  !-----------------------------------------------------------------------
  ! here is the code for L=1, once again hard coded
  !-----------------------------------------------------------------------
  IF (LYDIM <1) GOTO 100
  !DIR$ IVDEP
  !OCL NOVREC
  DO IND=1,INDMAX
     YLM(IND,2)  = (FAK*SQRT(3.d0))*Y(IND)
     YLM(IND,3)  = (FAK*SQRT(3.d0))*Z(IND)
     YLM(IND,4)  = (FAK*SQRT(3.d0))*X(IND)
  ENDDO
  !-----------------------------------------------------------------------
  ! code for L=2,
  !-----------------------------------------------------------------------
  IF (LYDIM <2) GOTO 100
  !DIR$ IVDEP
  !OCL NOVREC
  DO IND=1,INDMAX
     YLM(IND,5)= (FAK*SQRT(15.d0))  *X(IND)*Y(IND)
     YLM(IND,6)= (FAK*SQRT(15.d0))  *Y(IND)*Z(IND)
     YLM(IND,7)= (FAK*SQRT(5.d0)/2.d0)*(3*Z(IND)*Z(IND)-1)
     YLM(IND,8)= (FAK*SQRT(15.d0))  *X(IND)*Z(IND)
     YLM(IND,9)= (FAK*SQRT(15.d0)/2.d0)*(X(IND)*X(IND)-Y(IND)*Y(IND))
  ENDDO
  !-----------------------------------------------------------------------
  ! initialize all componentes L>2 to zero (Actually unimplemented)
  !-----------------------------------------------------------------------
  IF (LYDIM <3) GOTO 100
  STOP 'L > 2 unimplemented in SETYLM'

100 CONTINUE

END SUBROUTINE SETYLM







SUBROUTINE calc_recip_lattice(a,b)
  IMPLICIT NONE
  REAL*8,INTENT(IN) :: a(3,3)
  REAL*8,INTENT(OUT) :: b(3,3)

  REAL*8 :: norm

  REAL*8,PARAMETER :: pi = 4.d0*datan(1.d0)

  norm = 1.d0/DOT_PRODUCT(a(:,1),cross(a(:,2),a(:,3)))
  norm = norm * 2.d0*pi

  b(:,1) = cross(a(:,2),a(:,3))*norm
  b(:,2) = cross(a(:,3),a(:,1))*norm
  b(:,3) = cross(a(:,1),a(:,2))*norm

END SUBROUTINE calc_recip_lattice

FUNCTION cross(a, b) 
  REAL*8,INTENT (in) :: a(3), b(3)
  REAL*8 :: cross(3)

  cross(1) = a(2) * b(3) - a(3) * b(2)
  cross(2) = a(3) * b(1) - a(1) * b(3)
  cross(3) = a(1) * b(2) - a(2) * b(1)
END FUNCTION cross


SUBROUTINE grid_volume(gridvol, b, mesh)
  IMPLICIT NONE
  REAL*8, INTENT(IN)    ::   b(3,3)          ! Reciprocal lattice vectors
  INTEGER, INTENT(IN)   ::   mesh(3)         ! Number of grid points in each direction, inverse gives lenght of grid vector
  REAL*8, INTENT(OUT)   ::   gridvol         ! Volume of box to use for integration
  REAL*8                ::   grid_vec(3,3)    ! Vectors of gridbox in direction of recip lattice vectors
  INTEGER               ::   j
 
  DO j=1,3
  grid_vec(:,j)=(1.d0/DBLE(mesh(j)))*b(:,j)
  ENDDO
  
  gridvol = DOT_PRODUCT(grid_vec(:,1), cross(grid_vec(:,2),grid_vec(:,3)))


END SUBROUTINE grid_volume







ENDMODULE
