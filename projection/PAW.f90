MODULE PAW
  USE projection_shared
  IMPLICIT NONE
  PRIVATE
  PUBLIC PAW_proj_setup

  REAL*8,PARAMETER   ::  PAW_screen=100.d0
  !COMPLEX*16, PARAMETER    ::  sqrt_neg_one=(0.d0,1.d0)
  REAL*8, PARAMETER  ::   bohr = 0.529177249d0 !, pi=4.d0*ATAN(1.d0)
  REAL*8,PARAMETER   ::   twosqrtpi = 2.d0 * SQRT(pi)

  TYPE PAW_type
    INTEGER              ::  ncart
    REAL*8,ALLOCATABLE   ::  cart_coeff(:)
    INTEGER,ALLOCATABLE  ::  deriv_ind(:,:)
    LOGICAL              ::  second_order
  END TYPE  



CONTAINS

!This subroutine obtains the real-space overlap of the AO basis functions and PAW-augmenters over a sufficient number of unit cells
!Bloch-space overlap is needed for the projection, which can be obtained by linear transform of periodic real space overlaps
SUBROUTINE PAW_proj_setup(AO_basis,AO_l,PAW_overlap)
  USE rd_wavefunction
  IMPLICIT NONE

  TYPE(AO_function), DIMENSION(:)    ::  AO_basis
  INTEGER,DIMENSION(:,:),INTENT(IN)  ::  AO_l
  COMPLEX*16,ALLOCATABLE,INTENT(OUT)  :: PAW_overlap(:,:,:)

  INTEGER,ALLOCATABLE   :: PAW_l(:,:)  !Set of unit cells to use for calculating periodic overlap.  Subset of those used for all AO overlap 
  REAL*8,ALLOCATABLE   ::  real_PAW_overlap(:,:,:)  !Holds the overlap functions for each unit cell.  Will later be used to Bloch-space overlaps
  INTEGER    ::  nl, l_half  !Index of central unit cell in 'PAW_l'. The unit cells are again mirrored across this point within 'PAW_l', just as in 'index_l'.

  LOGICAL    ::  on_site, spher_harm_match  !For determining if PAW and AO are on the same atomic center and if so, if they have the same spherical harmonic
  REAL*8     ::  gauss  !In all overlap calculations, this variable holds the AO contribution
  !REAL*8     ::  ang_corr(4) !When on-site overlaps are calculated, the GTOs lack of separate spherical norm is accounted for.
  INTEGER    ::  l,ll,m                !Terms for keeping track of spherical harmonic of PAW functions
  INTEGER    ::  nu,ipro,iatom,itype   !Counters used in loops over AO and PAW augmenter functions
  REAL*8     ::  bohr_kpt(3)

  !COMPLEX*16  :: arg

  INTEGER    ::  j,igauss,ik,il

  REAL*8   ::  t1,t2,t3

  TYPE(PAW_type),DIMENSION(5)  ::  d_PAW


  !Determine appropriate unit cells to use in real-space PAW-AO overlap calculations, and store indices in PAW_l.
  CALL setup_PAW_l(AO_l,PAW_l)
  nl = SIZE(PAW_l,2)
  l_half = (nl+1) / 2

  ALLOCATE(real_PAW_overlap(s_dim,npromax,SIZE(PAW_l,2)), PAW_overlap(s_dim,npromax,nkpts))
  real_PAW_overlap = 0.d0

  CALL CPU_TIME(t1)

  !Set up PAW_types for d-type PAW. To be used for calculating their overlap with 2nd order Taylor series components.
  CALL setup_d_PAW(d_PAW)

  !WRITE(6,*)'Calculating real space PAW overlaps'

  !Here the overlap are finally obtained
  !Looping over the projector functions is not so straightforward (ipro=1,npro). So I use nested loops over atoms, augmentor channels, l-, m-angular quantum numbers
  !For each ipro and nu pair, first overlap within the central unit cell is calculated.
  !This can either corresponds to a pair on the same atomic center or on different
  !If the AO basis function and PAW augmenter are centered on the same atom, the overlap is performed using a numerical integral on the gird defining the PAW
  !If the pair are centered on different atoms, a Taylor expansion of the AO function around the PAW center is used.
  !For overlap between unit cell all pairs are necessarily off-site and the same approximation will be used.
  DO nu=1,s_dim
     !WRITE(6,*)'nu',nu
     ipro = 1

     DO iatom=1,n_atom
        itype = itypes(iatom)

        on_site = .FALSE.
        IF( AO_basis(nu)%atom == iatom )on_site=.TRUE. !Test to see if nu and ipro are centered on the same atom

        DO l=1,P(itype)%ldim

           ll=P(itype)%lps(l)

           DO m=1,2*ll+1
              !WRITE(6,*)'ipro',ipro
              spher_harm_match = .FALSE.

              !If the PAW and augmenter are centered on the same atom, the spherical-harmonics are checked
              !If they are the same (l- and m- quantum numbers match) the overlap is non-zero.  
              IF( on_site )THEN
                 !The l-quantum number is initally screened, then a function is called to test the m-quantum number
                 IF( AO_basis(nu)%l == ll )THEN
                    !WRITE(6,*)'found on-site and same l, testing m-quantum'
                    CALL spher_harm_test(spher_harm_match,ll,m,AO_basis(nu)%m)
                 ENDIF
              ENDIF

              IF( spher_harm_match )THEN

                 !Onsite overlap will be treated completely. 
                 !The radial component will be numerically integrated along the existing grid of the PAW
                 DO j=1,P(itype)%nmax
                    !Gaussian component is calculated at each radial point.  Note, there is no cartesian component here since it distributes out of sum over gaussians
                    gauss = 0.d0
                    DO igauss=1,AO_basis(nu)%num_gauss
                       gauss = gauss + AO_basis(nu)%coeff(igauss)*AO_basis(nu)%norm(igauss)*EXP(-AO_basis(nu)%alpha(igauss)*P(itype)%r(j)**2)
                    ENDDO

                    !The value of the GTO is multiplied by the value of the augmenter, a weight factor, and radial factor.
                    !For the radial factor the exponents come from:
                    !  2 - Spherical integration Jacobian
                    !  ll - Cartesian aspect of GTO. Specific components (i.e. x,y,z) are separated into the spherical portion as x/r or xy/r^2, leaving r^ll.
                    real_PAW_overlap(nu,ipro,l_half) = real_PAW_overlap(nu,ipro,l_half) + gauss*P(itype)%wdiff(j,l)*P(itype)%r(j)**(2+ll)*P(itype)%si(j)
                 ENDDO

                 !Correction for fact that GTO's norm is not separated into radial and angular components 
                 !Correction equals spherical harmonic norm ^ (-1/2)
                 real_PAW_overlap(nu,ipro,l_half) = real_PAW_overlap(nu,ipro,l_half) * twosqrtpi / SQRT(DBLE(dble_factor(2*ll+1)))

              ENDIF

              !Now we will calculate the off site overlaps (if AO and PAW do not share the same atomic center).
              !This will need to cover all unit cells, to properly capture the real space periodicity for the Bloch-transform
              !Off-site overlap is non-trivial since the PAW is represented by a radial grid multiplied by a spherical harmonic.
              !Thus, to deal with it we use a Taylor exapnsion of the AO-function around the center of the PAW-function
              !This makes the overlap integarl again separable in to a raidal and angular component
              !The angular component is again trivial as a Taylor series is akin to an expansion in spherical harmonics (0th order = s; 1st=p; 2nd=d; etc.)
              !The radial component is then simply a constant term from the AO-fucntion (Taylor exp. coeff.) & the integral of the PAW-radial component (constant)

              !We assume that off site overlaps of d-type PAWs are zero as their oscillatory nature would be orhtogonal to a smooth AO-function.
              !IF( ll .GE. 0 )GOTO 110

              !For s-type PAW, the AO contribution is simply the value of the AO fucntion at the center of the PAW function
              IF( ll == 0 )THEN

                 !The symmetry of the l-indices is not used since distance are nto related in a straightforward fashion
                 !However, it avoids including l_half, which is not always computed
                 DO il=1,l_half-1
                    !WRITE(6,*)'il',il
                    CALL stype_PAW_offsite(AO_basis(nu),ipro,PAW_l(:,il),P(itype)%center_value(l,:),gauss)
                    real_PAW_overlap(nu,ipro,il) = gauss  !P(itype)%center_value(l) * gauss
                    CALL stype_PAW_offsite(AO_basis(nu),ipro,PAW_l(:,nl+1-il),P(itype)%center_value(l,:),gauss)
                    real_PAW_overlap(nu,ipro,nl+1-il) =  gauss !P(itype)%center_value(l) * gauss
                 ENDDO

                 !If the PAW and AO are centered in the same unit cell, but NOT on the same atom, off-site overlap will need to be computed
                 IF( .NOT. on_site )THEN
                    CALL stype_PAW_offsite(AO_basis(nu),ipro,(/0,0,0/),P(itype)%center_value(l,:),gauss)
                    real_PAW_overlap(nu,ipro,l_half) = gauss !P(itype)%center_value(l) * gauss
                 ENDIF

              !For p-type PAW, the AO contribution is the gradient of the AO function in the PAW's m-index direction, calculated at the PAW center
              ELSEIF( ll == 1 )THEN

                 DO il=1,l_half-1
                    !WRITE(6,*)'il',il,PAW_l(:,il)
                    CALL ptype_PAW_offsite(AO_basis(nu),ipro,MODULO(m,3)+1,PAW_l(:,il),P(itype)%center_value(l,:),gauss)
                    real_PAW_overlap(nu,ipro,il) = gauss  !P(itype)%center_value(l) * gauss
                    CALL ptype_PAW_offsite(AO_basis(nu),ipro,MODULO(m,3)+1,PAW_l(:,nl+1-il),P(itype)%center_value(l,:),gauss)
                    real_PAW_overlap(nu,ipro,nl+1-il) = gauss ! P(itype)%center_value(l) * gauss
                    !WRITE(6,*)
                 ENDDO

                 !Similar to s-type, PAW and AO on different centers within the same unit cell are treated.
                 IF( .NOT. on_site )THEN
                    CALL ptype_PAW_offsite(AO_basis(nu),ipro,MODULO(m,3)+1,(/0,0,0/),P(itype)%center_value(l,:),gauss)
                    real_PAW_overlap(nu,ipro,l_half) = gauss !P(itype)%center_value(l) * gauss
                 ENDIF

              ELSEIF( ll == 2 )THEN
              
                 DO il=1,l_half-1
                    CALL dtype_PAW_offsite(AO_basis(nu),ipro,PAW_l(:,il),d_PAW(m),P(itype)%center_value(l,:),gauss)
                    real_PAW_overlap(nu,ipro,il) = gauss
                    CALL dtype_PAW_offsite(AO_basis(nu),ipro,PAW_l(:,nl+1-il),d_PAW(m),P(itype)%center_value(l,:),gauss) 
                    real_PAW_overlap(nu,ipro,nl+1-il) = gauss
                 ENDDO

                 IF( .NOT. on_site )THEN
                    CALL dtype_PAW_offsite(AO_basis(nu),ipro,(/0,0,0/),d_PAW(m),P(itype)%center_value(l,:),gauss)
                    real_PAW_overlap(nu,ipro,l_half) = gauss
                 ENDIF

              ENDIF

110           CONTINUE
              !WRITE(6,*)
              ipro=ipro+1
           ENDDO
        ENDDO
       
     ENDDO

     !WRITE(6,*)
  ENDDO

  !Perform a Bloch transform to get the necessary reciprocal space overlaps
  PAW_overlap = 0.d0
  CALL real_to_bloch(real_PAW_overlap,PAW_overlap,PAW_l,-kpt)
  !The phase factor e^ikr also needs to be incorporated into the overlap.
  !Rigorourly this term is not constant in space, but since the PAW-augmenter functions are localized, we use a constant term at the center of the PAW.
  DO ik=1,nkpts
     bohr_kpt = kpt(1,ik)*b(:,1) + kpt(2,ik)*b(:,2) + kpt(3,ik)*b(:,3)
     !WRITE(6,*)'bohr_kpt',bohr_kpt
     DO ipro=1,npromax
        PAW_overlap(:,ipro,ik) = PAW_overlap(:,ipro,ik) * EXP(sqrt_neg_one*DOT_PRODUCT(PAW_pos(ipro,:),bohr_kpt))
     ENDDO
  ENDDO

  CALL CPU_TIME(t2)
  !WRITE(6,*)'Finished PAW overlap setup',SNGL(t2-t1)
  !WRITE(6,*)

END SUBROUTINE PAW_proj_setup

  !This subroutine is used to calculate the overlap of a given basis function with an s-type PAW that is offsite
  !This requires the value of the basis function at the center of the PAW
  !This value must then be multiplied appropriately in the main program
SUBROUTINE stype_PAW_offsite(AO,ipro,PAW_l,radial_int,offsite)
  USE rd_wavefunction
  IMPLICIT NONE

  TYPE(AO_function),INTENT(IN)             ::  AO

  INTEGER, INTENT(IN)                      ::  ipro  !Index of the given basis function and augmentor
  INTEGER, DIMENSION(3), INTENT(IN)        ::  PAW_l   !Indices of the unit cell which the basis function will move around into
  REAL*8, DIMENSION(3), INTENT(IN)         ::  radial_int  !Radial component of overlap integral, is simply the integral of the PAW radial function
  REAL*8, INTENT(OUT)                      ::  offsite   !output value

  REAL*8, DIMENSION(3)      ::    lvec    !How much the position of atom 2 must be shifted to place it into the appropriate unit cell
  REAL*8, DIMENSION(3)      ::    rdiff
  REAL*8                    ::    dist, screen, cartesian
  REAL*8                    ::    kroenecker, arg1, arg2
  INTEGER, DIMENSION(3)     ::    kroen_lmn
  INTEGER                   ::    j, k, icart

  REAL*8     ::  zeroth, second


  offsite = 0.d0

  lvec = 0.d0
  DO j=1,3
     lvec = lvec + a(:,j)*PAW_l(j)
  ENDDO

  rdiff = PAW_pos(ipro,:) - (AO%pos + lvec)
  dist = DOT_PRODUCT(rdiff,rdiff)
  !WRITE(6,*)'rdiff',rdiff


  !First the overlap of the zeroth order Taylor expansion with the s-type PAW
  zeroth = 0.d0

  cartesian = 0.d0
  DO icart=1,AO%ncart
     cartesian = cartesian + AO%cart_coeff(icart)*PRODUCT(rdiff**AO%cart_mat(icart,:))
  ENDDO
  !DO k=AO%num_gauss,1,-1  !The exponents are run through in ascending order, once one is found that is too high, the loop is exited since all others will be too high
  DO k=1,AO%num_gauss
     screen = dist * AO%alpha(k)
     !IF( screen .GT. PAW_screen )EXIT 
     zeroth = zeroth + AO%coeff(k)*AO%norm(k)*EXP(-screen)
  ENDDO
  zeroth = zeroth * cartesian * radial_int(1) * twosqrtpi 

  !Then the overlap of the second odrer Taylor expansion components with the s-type PAW
  second = 0.d0

  !The second derivative in each direction must be calculated. Mixed partials do not need to be accounted for
  kroenecker = 0.d0
  arg1 = 0.d0
  arg2 = 0.d0
  DO j=1,3
     DO icart=1,AO%ncart
        IF( AO%cart_mat(icart,j) .GE. 2 )THEN
           kroen_lmn = AO%cart_mat(icart,:)
           kroen_lmn(j) = kroen_lmn(j) - 2
           kroenecker = kroenecker + AO%cart_coeff(icart) * AO%cart_mat(icart,j) * (AO%cart_mat(icart,j)-1) * PRODUCT(rdiff**kroen_lmn)
        ENDIF
        arg1 = arg1 + 4.d0 * AO%cart_coeff(icart) * rdiff(j)**2 * PRODUCT(rdiff**AO%cart_mat(icart,:))
        arg2 = arg2 - AO%cart_coeff(icart) * (2.d0 + 4.d0*AO%cart_mat(icart,j)) * PRODUCT(rdiff**AO%cart_mat(icart,:))
     ENDDO
  ENDDO

  !IF( kroenecker .NE. 0.d0 )WRITE(6,*)'nonzero kroenecker',kroenecker,AO%l,AO%m

  !WRITE(6,*)'AO_l and_m',AO%l,AO%m
  !WRITE(6,*)'kroenecker',kroenecker
  !WRITE(6,*)'arg1,arg2 ',arg1,arg2
  !WRITE(6,*)

  !DO k=AO%num_gauss,1,-1
  DO k=1,AO%num_gauss
     screen = dist * AO%alpha(k)
     !IF( screen .GT. PAW_screen)EXIT
     second = second + (kroenecker + arg1*AO%alpha(k)**2 + arg2*AO%alpha(k)) * AO%coeff(k)*AO%norm(k)*EXP(-screen)
  ENDDO
  second = second * radial_int(3) * (twosqrtpi / 6.d0)  !The 1/6 factor is: 1/3 from the spher harmonic integration and 1/2 from taylor expansion

  !Then the overlap of both components is combined to give the overall overlap of the AO with the s-type PAW
  offsite = zeroth + second

  !WRITE(6,*)'s-type offsite',offsite

ENDSUBROUTINE stype_PAW_offsite


  !This subroutine is used in calculating the overlap of given basis function with a p-type augmentor function that is offsite
  !This requires the gradient of the basis function at the center of the PAW, in the direction corresponding to the PAW's spherical harmonic (px->x)
  !offsite will hold this value on exit which is then multiplied by other factors in the main program
SUBROUTINE ptype_PAW_offsite(AO,ipro,PAW_m,PAW_l,radial_int,offsite)
  USE rd_wavefunction
  IMPLICIT NONE

  TYPE(AO_function),INTENT(IN)             ::  AO 
  INTEGER, INTENT(IN)                      ::  ipro, PAW_m  !Indices of the basis funct, augmentor, and direction of gradient required
  INTEGER, DIMENSION(3), INTENT(IN)        ::  PAW_l        !Indices of the unit cell which the basis function will move around into
  REAL*8, DIMENSION(3), INTENT(IN)         ::  radial_int   !Radial component of overlap integral.  Simple integral of radial component of PAW.
  REAL*8, INTENT(OUT)                      ::  offsite      !output  

  REAL*8, DIMENSION(3)      ::    lvec    !How much the position of atom 2 must be shifted to place it into the appropriate unit cell
  REAL*8, DIMENSION(3)      ::    rdiff 
  REAL*8                    ::    dist, screen
  REAL*8                    ::    kroenecker, arg   !Used in actual gradient calculation, determined by each basis set augmentor pair
  INTEGER, DIMENSION(3)     ::    kroen_lmn
  INTEGER                   ::    k,j,icart

  lvec = 0.d0
  DO j=1,3
     lvec = lvec + a(:,j)*PAW_l(j)
  ENDDO

  rdiff = PAW_pos(ipro,:) - (AO%pos + lvec)
  dist = DOT_PRODUCT(rdiff,rdiff)
  !WRITE(6,*)'dist',dist

  offsite = 0.d0

  kroenecker = 0.d0
  arg = 0.d0
  DO icart = 1,AO%ncart
     IF( AO%cart_mat(icart,PAW_m) .GE. 1 )THEN !Test if the cartesian part of the overlap will be non-zero.
         kroen_lmn = AO%cart_mat(icart,:)
         kroen_lmn(PAW_m) = kroen_lmn(PAW_m) - 1
         kroenecker = kroenecker + AO%cart_coeff(icart) * AO%cart_mat(icart,PAW_m) * PRODUCT(rdiff**kroen_lmn)
     ENDIF
     arg = arg - 2.d0 * AO%cart_coeff(icart) * rdiff(PAW_m) * PRODUCT(rdiff**AO%cart_mat(icart,:))
  ENDDO

  !WRITE(6,*)'PAW_m     ',PAW_m
  !WRITE(6,*)'AO_l and_m',AO%l,AO%m
  !WRITE(6,*)'kroenecker',kroenecker
  !WRITE(6,*)'arg       ',arg


  !Now that the orbital has been identified and kroenecker and arg defined, the gradient of nu at PAW-ipro's center can be computed
  !DO k=AO%num_gauss,1,-1
  DO k=1,AO%num_gauss
     screen = dist * AO%alpha(k)
     !IF( screen .GT. PAW_screen )EXIT
     offsite = offsite + (kroenecker + AO%alpha(k)*arg) * AO%coeff(k)*AO%norm(k)*EXP(-screen)
  ENDDO

  offsite = offsite * radial_int(2) * (twosqrtpi/SQRT(3.d0))

  !WRITE(6,*)'radial_int 2',radial_int(2)
  !WRITE(6,*)'p-type offsite',offsite


END SUBROUTINE ptype_PAW_offsite



SUBROUTINE dtype_PAW_offsite(AO,ipro,PAW_l,d_PAW,radial_int,offsite)
  USE rd_wavefunction
  IMPLICIT NONE

  TYPE(AO_function),INTENT(IN)             ::  AO
  INTEGER, INTENT(IN)                      ::  ipro        !Indices of the basis funct, augmentor, and direction of gradient required
  INTEGER, DIMENSION(3), INTENT(IN)        ::  PAW_l       !Indices of the unit cell which the basis function will move around into
  TYPE(PAW_type),INTENT(IN)                ::  d_PAW       !Information on the particular PAW for which overlap is being calculated
  REAL*8, DIMENSION(3), INTENT(IN)         ::  radial_int  !Radial component of overlap integral.  Simple integral of radial component of PAW.
  REAL*8, INTENT(OUT)                      ::  offsite     !output  

  REAL*8, DIMENSION(3)      ::    lvec    !How much the position of atom 2 must be shifted to place it into the appropriate unit cell
  REAL*8, DIMENSION(3)      ::    rdiff
  REAL*8                    ::    dist, screen, gaussian
  REAL*8                    ::    kroenecker, arg1, arg2   !Used in actual gradient calculation, determined by each basis set augmentor pair
  INTEGER, DIMENSION(2,3)   ::    kroen_lmn
  INTEGER, DIMENSION(3)     ::    arg_lmn
  INTEGER                   ::    k,j,icart,jcart

  lvec = 0.d0
  DO j=1,3
     lvec = lvec + a(:,j)*PAW_l(j)
  ENDDO

  rdiff = PAW_pos(ipro,:) - (AO%pos + lvec)
  dist = DOT_PRODUCT(rdiff,rdiff)
  !WRITE(6,*)'rdiff',rdiff

  offsite = 0.d0

  DO icart=1,d_PAW%ncart
     kroenecker = 0.d0
     arg1 = 0.d0
     arg2 = 0.d0

     DO jcart=1,AO%ncart

        kroen_lmn(1,:) = AO%cart_mat(jcart,:)
        IF( kroen_lmn(1,d_PAW%deriv_ind(icart,1)) .GE. 1 )THEN
           kroen_lmn(1,d_PAW%deriv_ind(icart,1)) = kroen_lmn(1,d_PAW%deriv_ind(icart,1)) - 1
           kroen_lmn(2,:) = kroen_lmn(1,:)
           IF( kroen_lmn(2,d_PAW%deriv_ind(icart,2)) .GE. 1 )THEN
              kroen_lmn(2,d_PAW%deriv_ind(icart,2)) = kroen_lmn(2,d_PAW%deriv_ind(icart,2)) - 1
              kroenecker = kroenecker + AO%cart_coeff(jcart)*AO%cart_mat(jcart,d_PAW%deriv_ind(icart,1)) &
                               &* kroen_lmn(1,d_PAW%deriv_ind(icart,2)) * PRODUCT(rdiff**kroen_lmn(2,:))
           ENDIF
        ENDIF

        arg1 = arg1 + (AO%cart_coeff(jcart) * 4.d0 * rdiff(d_PAW%deriv_ind(icart,1)) * rdiff(d_PAW%deriv_ind(icart,2)) * PRODUCT(rdiff**AO%cart_mat(jcart,:)))

        IF( d_PAW%second_order )THEN
           arg2 = arg2 - (AO%cart_coeff(jcart) * 2.d0 * PRODUCT(rdiff**AO%cart_mat(jcart,:)))
        ENDIF
        DO j=1,2
           IF( AO%cart_mat(jcart,d_PAW%deriv_ind(icart,j)) .GE. 1 )THEN
              k = MODULO(j,2)+1
              arg_lmn = AO%cart_mat(jcart,:)
              arg_lmn(d_PAW%deriv_ind(icart,j)) = arg_lmn(d_PAW%deriv_ind(icart,j)) - 1
              arg_lmn(d_PAW%deriv_ind(icart,k)) = arg_lmn(d_PAW%deriv_ind(icart,k)) + 1
 
              arg2 = arg2 - (AO%cart_coeff(jcart) * 2.d0 * AO%cart_mat(jcart,d_PAW%deriv_ind(icart,j)) * PRODUCT(rdiff**arg_lmn))
           ENDIF
        ENDDO

     ENDDO

     !WRITE(6,*)'AO_l and_m',AO%l,AO%m
     !WRITE(6,*)'kroenecker',kroenecker
     !WRITE(6,*)'arg1,arg2 ',arg1,arg2
     !WRITE(6,*)

     gaussian = 0.d0
     !DO k=AO%num_gauss,1,-1
     DO k=1,AO%num_gauss
        screen = dist * AO%alpha(k)
        !IF( screen .GT. PAW_screen)EXIT
        gaussian = gaussian + (kroenecker + arg1*AO%alpha(k)**2 + arg2*AO%alpha(k)) * AO%coeff(k)*AO%norm(k)*EXP(-screen)
     ENDDO

     offsite = offsite + d_PAW%cart_coeff(icart)*gaussian

  ENDDO

  offsite = offsite * radial_int(3) * (twosqrtpi / SQRT(15.d0))



END SUBROUTINE dtype_PAW_offsite




!Check if the AO-basis function described by 'test_lmn' has the same spherical harmonic 
! as the PAW augmenter fucntion described by ll and m angular momenta quantum numbers
!There is not much fancy about this subroutine, especially for d-orbitals and beyond, where it is simply brute checking
SUBROUTINE spher_harm_test(spher_harm_match,ll,PAW_m,AO_m)
  IMPLICIT NONE

  LOGICAL,INTENT(INOUT)            ::  spher_harm_match
  INTEGER,INTENT(IN)               ::  ll,PAW_m,AO_m

  !INTEGER  :: m_test

  !The basic structure of this subroutine is a series of IF statements
  !On entry to the program, 'spher_harm_match'=FALSE, and will only change it the program reaches a TRUE entry
  !If no TRUEs are found after exiting the IF-tree, teh program will jump to the end and 'spher_harm_match' remains FALSE
  !Each TRUE corresponds to matching spherical harmonics and leads to the program jumping to changing 'spher_harm_match' as TRUE

  !The preliminary test is to make sure the l-quantum number (s-,p-,etc.) matches. If not, the m-quantum number does not matter
  IF( ll .EQ. 0 )THEN  !For s-type, if l- matches, m- necessarily does as well
     GOTO 20
  ELSEIF( ll .EQ. 1 )THEN  !For p-type, m-must be checked.  VASP ordering is y,z,x; which is just one place shifted from x,y,z used for AO p-types.
     IF( MODULO(PAW_m,3)+1 == AO_m )THEN
        GOTO 20
     ENDIF
  ELSEIF( ll .EQ. 2 )THEN  !For d-type, VASP ordering is more different form mine, so I just manually check every indices.

     SELECT CASE (PAW_m)
        CASE (1)
           IF( AO_m == 3 )GOTO 20       !xy PAW/orbital
        CASE (2)
           IF( AO_m == 1 )GOTO 20       !yz PAW/orbital
        CASE (3)
           IF( AO_m == 5 )GOTO 20       !z^2 PAW/orbital
        CASE (4)
           IF( AO_m == 2 )GOTO 20       !xz PAW/orbital
        CASE (5)
           IF( AO_m == 4 )GOTO 20       !x^2-y^2 PAW/orbital
        CASE DEFAULT
           STOP 'Improper PAW_m value passed in d-testing for spher_harm_check'
     END SELECT

  ELSEIF( ll .EQ. 3 )THEN  !f-type ordering also done using brute checking

     SELECT CASE (PAW_m)
        CASE (1)
           IF( AO_m == 3 )GOTO 20       !y(3x^2-y) PAW/orbital
        CASE (2)
           IF( AO_m == 1 )GOTO 20       !xyz PAW/orbital
        CASE (3)
           IF( AO_m == 6 )GOTO 20       !yz^2 PAW/orbital
        CASE (4)
           IF( AO_m == 7 )GOTO 20       !z^3 PAW/orbital
        CASE (5)
           IF( AO_m == 5 )GOTO 20       !xz^2 PAW/orbital
        CASE (6)
           IF( AO_m == 2 )GOTO 20       !(x^2-y^2)z PAW/orbital
        CASE (7)
           IF( AO_m == 4 )GOTO 20       !x(x^2-3y^2) PAW/orbital
        CASE DEFAULT
           STOP 'Improper PAW_m value passed in f-testing for spher_harm_check'
     END SELECT

  ELSE
     STOP 'Can not check spherical harmonics beyond f-type'
  ENDIF

  GOTO 40

20 CONTINUE

  !WRITE(6,*)'glory DAY!!! matching spherical harmonics'
  spher_harm_match = .TRUE.

40 CONTINUE


END SUBROUTINE




SUBROUTINE  setup_d_PAW(d_PAW)
  IMPLICIT NONE

  TYPE(PAW_type),DIMENSION(5),INTENT(INOUT)   ::  d_PAW

  INTEGER  ::  j,icart

  DO j=1,5

     IF( j .EQ. 3 .OR. j .EQ. 5 )GOTO 25

     d_PAW(j)%ncart = 1
     d_PAW(j)%second_order = .FALSE.
     CALL array_allocation(d_PAW(j))
     d_PAW(j)%cart_coeff = 1.d0

     IF( j .LT. 3 )THEN
        d_PAW(j)%deriv_ind(1,1) = j
        d_PAW(j)%deriv_ind(1,2) = j+1
     ELSE
        d_PAW(j)%deriv_ind(1,1) = 1
        d_PAW(j)%deriv_ind(1,2) = 3
     ENDIF
     GOTO 35

25   CONTINUE

     d_PAW(j)%second_order = .TRUE.

     IF( j .EQ. 5 )THEN
        d_PAW(j)%ncart = 2
        CALL array_allocation(d_PAW(j))
        d_PAW(j)%cart_coeff(1) = 0.5d0
        d_PAW(j)%cart_coeff(2) = -0.5d0
     ELSE
        d_PAW(j)%ncart = 3
        CALL array_allocation(d_PAW(j))
        d_PAW(j)%cart_coeff = -0.5d0 / SQRT(3.d0)
        d_PAW(j)%cart_coeff(3) = 1.d0 / SQRT(3.d0)
     ENDIF

     DO icart=1,d_PAW(j)%ncart
        d_PAW(j)%deriv_ind(icart,:) = icart
     ENDDO

35   CONTINUE

     !WRITE(6,*)'For PAW              ',j
     !WRITE(6,*)'number of cartesian  ',d_PAW(j)%ncart      
     !DO icart=1,d_PAW(j)%ncart
     !   WRITE(6,*)'Cartesian component  ',icart
     !   WRITE(6,*)d_PAW(j)%cart_coeff(icart),d_PAW(j)%deriv_ind(icart,:)
     !ENDDO
     !WRITE(6,*)

  ENDDO


  CONTAINS

    !Simple subroutine that allocates the cartesian arrays for each individual PAW.
    !This is utilized to clean up setup_d_PAW
    SUBROUTINE array_allocation(ind_PAW)
      IMPLICIT NONE

      TYPE(PAW_type)    ::  ind_PAW

      ALLOCATE(ind_PAW%cart_coeff(ind_PAW%ncart))
      ALLOCATE(ind_PAW%deriv_ind(ind_PAW%ncart,2))

    END SUBROUTINE

END SUBROUTINE



!The overlap of the AO and augmenter functions is needed in Bloch-space for the projection, therefore the periodic real space overlaps are needed.
!The set of unit cells used for AO overlap is overkill as the augmenters are non-zero outside the core region.
!This subroutine generates a new set of l-vector indices to use for calculating PAW overlaps.
SUBROUTINE  setup_PAW_l(AO_l,PAW_l)
  IMPLICIT NONE

  INTEGER,DIMENSION(:,:),INTENT(IN)  ::  AO_l       !Set of l_vectors used in AO periodic overlap calculations.  
  INTEGER,ALLOCATABLE,INTENT(OUT)    ::  PAW_l(:,:) !Set of trimmed lvectors to be used for PAW-AO overlap calculations, trimmed version of AO_l

  LOGICAL,DIMENSION(SIZE(AO_l,2))    ::  PAW_keep   !Used in screening process

  INTEGER  ::  shell_max,move_screen  !Cutoff criteria for inclusion of unit cells into PAW overlap calculations
  INTEGER  ::  nl,shell_tot           !Total number of unit cells in initial AO_l and number that will be used in PAW overlap calculations
  INTEGER  ::  il,jl

  !The real space distance cutoff will be exactly half that used in the AO overlap (1 AO function is taken into account instead of 2)
  !'shell_max' is the number of unit cells necessary to span this cutoff, in a straight line
  shell_max = CEILING(MAXVAL(AO_l) / 2.d0)
  !The cube of unit cells represented by shell_max, will be trimmed into a spherical shape, based on how many "moves" the unit cell is away from the central position
  move_screen = MAX(shell_max+2, 3) 


  nl=SIZE(AO_l,2)
  PAW_keep = .FALSE.
  shell_tot = 0

  !Screening is done similary to that in lvec_init (bloch_overlap MODULE), in terms of unit cell translations away from the central cell.
  !The additional condition is necessary here because of the +1 used to define 'move_screen'
  !This extends the possible moves beyond shell_max if all moves were made in the same direction off the central unit cell
  !Not taking this into account would introduce 6 unnecessary unit cell pair overlap calculations (beyond shell_max overlap is zero)
  DO il=1,nl
    !IF( SUM(ABS(AO_l(:,il))) <= move_screen .AND. MAXVAL(ABS(AO_l(:,il))) /= move_screen )THEN
        shell_tot = shell_tot + 1
        PAW_keep(il) = .TRUE.
    !ENDIF
  ENDDO

  ALLOCATE(PAW_l(3,shell_tot))
  jl = 0
  DO il=1,shell_tot
     DO       
       jl = jl + 1
       IF( PAW_keep(jl) )THEN
          PAW_l(:,il) = AO_l(:,jl)
          EXIT
       ENDIF
     ENDDO
  ENDDO

  !WRITE(6,*)shell_tot,'unit cells out of',nl,'will be used in PAW offsite'
  !WRITE(6,*)



END SUBROUTINE setup_PAW_l






END MODULE PAW
