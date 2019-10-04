MODULE bloch_overlap
  USE shared
  IMPLICIT NONE
  PRIVATE
  PUBLIC bloch_space_overlap, renorm_AO

CONTAINS



!This subroutine calculates the inverse overlap matrices at each k-point used in the planewave calculation
!The first step is to calculate the real space overlap over a sufficient number of unit cells 
!Linear combinations of these real space overlaps (Bloch transforms) then yield the k-space overlap
!LAPACK subroutines are then used to calculate the inverse overlap matrices, which are required for the projection

!The process for calculating the Bloch-space matrices comes from and uses the notation of:
! Kudin and Scuseria, Chem. Phys. Lett.; 289, 611-616 (1998)

SUBROUTINE bloch_space_overlap(AO_basis,index_l,bloch_s,latt_vec,kpt)
  !USE mkl95_lapack
  !USE mkl95_blas
  !USE shared
  IMPLICIT NONE

  TYPE(AO_function), DIMENSION(:), INTENT(IN) :: AO_basis

  INTEGER  ::  s_dim

  REAL*8, DIMENSION(:,:,:), ALLOCATABLE ::  per_s_mat

  COMPLEX*16, DIMENSION(:,:,:), INTENT(OUT)  ::  bloch_s

  REAL*8, DIMENSION(:,:), INTENT(IN)   ::  kpt

  !This array holds the indices of the unit cells to be used in calculating real space overlap. 
  !The indices are counts of unit cells moved away from the central.
  INTEGER, ALLOCATABLE, INTENT(OUT)  ::  index_l(:,:)  
  INTEGER               ::  num_l,l_half  !Total number of unit cells, and the index of the central unit cell
  REAL*8, DIMENSION(3)  ::  lvec   !Actual distance vector from origin of central unit cell to origin of unit cell of interest

  REAL*8, DIMENSION(3,3), INTENT(IN)  :: latt_vec

  INTEGER   ::  il,ik,i,j

  REAL*8    ::  t1,t2

  s_dim = SIZE(AO_basis)


  WRITE(6,*)'*** Bloch Overlap Calculation ***'

  !The first step is to determine all real space unit cells that must be included in calculating overlaps
  !The indices of these unit cells (relative to a central 0,0,0) are then stored in index_l
  CALL lvec_init(index_l,num_l,AO_basis,latt_vec)

  !In general the overlap matrix between each unit cell pair must be kept, and then used in the Bloch-transform later.
  !However for a gamma-point calculation this is not true as the Bloch-transform is simply the summation of all real space matrices
  !Thus a smaller (REAL*8 type) matrix will be used for that particular case.
  !This attempts to save memory which is more likely an issue for cases where gamma point calculations can be used (large unit cell --> lots of atoms)

  ALLOCATE(per_s_mat(s_dim,s_dim,num_l))

  !Certain symmetries among l-vectors will be utilized below, that require the need for the l_half variable
  !Additionally based on how index_l has been constructed, l_half corresponds to the central unit cell (0,0,0)
  l_half = (num_l + 1)/2

  !The real space overlaps will now be calculated for all relevant unit cells (l-vectors)
  WRITE(6,*)'Calcualting Real Space Overlap Matrices'
  CALL CPU_TIME(t1)


  !Only half of the matrices acutally need to be computed since symmetry can be used
  !S{mu,nu},{0,l} = TRANSPOSE( S{mu,nu},{0,-l} )
  !Based on how index_l was constructed, -l is the mirror of l across l_half
  !Thus, the negative of the first l-vector is the last l-vecotr and so one.
  per_s_mat = 0.d0
  DO il=l_half,1,-1
     !lvec contains the vector from the origin of unit cell (0,0,0) to the unit cell indexed by il 
     lvec = 0.d0
     DO j=1,3
        lvec = lvec + latt_vec(j,:)*index_l(j,il)
     ENDDO
     CALL real_overlap(lvec,AO_basis,per_s_mat(:,:,il))

     per_s_mat(:,:,(num_l + 1 - il)) = TRANSPOSE(per_s_mat(:,:,il))

     !IF( SUM(ABS(index_l(:,il))) .LE. 1 )THEN
     !   WRITE(6,*)'real space overlap matrix',il,index_l(:,il)
     !   DO j=1,s_dim
     !      WRITE(6,'(14F11.7)')per_s_mat(j,1:14,il)
     !      !WRITE(6,*)per_s_mat(j,j,il)
     !   ENDDO
     !   WRITE(6,*)
     !ENDIF

  ENDDO


  CALL CPU_TIME(t2)

  WRITE(6,*)'Finished real space overlap computations',SNGL(t2-t1)
  WRITE(6,*)


  !Now the real-space overlap matrices will be processed to obtain the overlap matrix inverses at each k-point from the PW calculation.
  !The inversion is more straighforward to do in Bloch-space (independent at eack k-point), so the Bloch-transform is performed first

  !Convert from S_{mu nu}^{0 l} to S_{mu nu}^(k) via Fourier transform
  bloch_s = 0.d0
  CALL real_to_bloch(per_s_mat,bloch_s,index_l,TRANSPOSE(kpt))
  CALL CPU_TIME(t1)

  WRITE(6,*)'Finished Bloch transform of overlap matrice',SNGL(t1-t2)
  WRITE(6,*)



END SUBROUTINE bloch_space_overlap


!This subroutine reinforces that every AO functions is normalized.  
!In assumbling the basis set each individual Gaussian is normalized, but there is no guarantee the linear combo into an AO is also normalized
!Since CRYSTAL assumes perfect normalization, I need to scale coefficients to account for this.

!This is just a repurposed version of the real_overlap subroutine (SEE BELOW)
!It only calculates norms of all the basis functions using the overlap helper subroutines
!The formulas for overlap all come from: Clementi and Davis, Journal of Comp. Physics; 2,223-244(1967)
SUBROUTINE renorm_AO(AO_basis)
  IMPLICIT NONE

  TYPE(AO_function),DIMENSION(:)    ::  AO_basis 

  INTEGER   ::  s_dim

  INTEGER   ::  nu,mu
  INTEGER   ::  nugauss,mugauss

  REAL*8  ::  old_norm  !, DIMENSION(:), ALLOCATABLE    ::  old_norm

  REAL*8,DIMENSION(3)  ::  PA,PB  !Vectors from each nuclei to the center of the resulting gaussian
  REAL*8               ::  a1,a2  !Shorthand for the exponents of the basis functions of interest
  REAL*8, DIMENSION(3) ::  sigma  !Vector containg the 'sigma function' described in the Clementi and Davis paper
  REAL*8               ::  integral,cartesian  !Placeholders for 1) S-type overlap 2)factor from sigma functions that account for higher order harmonics

  INTEGER   :: i,j,k

  REAL*8, PARAMETER        ::  pi=3.141592653589793238462d0


  s_dim = SIZE(AO_basis,1)

  !All basis function pairs are looped over, and the resulting overlap added into the array s_mat
  DO nu=1,s_dim

        old_norm = 0.d0

        !Overlap is defined for each Gaussian function, but each basis function is a summation of Gaussians
        !Thus, for all basis function pairs the overlap is really a sum over all possible pairs of Gaussians
        DO nugauss=1,AO_basis(nu)%ngauss
           DO mugauss=1,AO_basis(nu)%ngauss

              !The full AO_basis notation is cumbersome to use and exponents are repeatedly used.  a1 and a2 are used for clarity.
              a1 = AO_basis(nu)%alpha(nugauss)
              a2 = AO_basis(nu)%alpha(mugauss)

              !This expression is the integral of the gaussian that results from the intersection of the particular Gaussian components of nu & mu
              !This is the complete overlap expression for two S-type gaussians.  But higher order harmonics complicate the calculation
              integral = AO_basis(nu)%norm(nugauss) * AO_basis(nu)%norm(mugauss) * (pi/(a1+a2))**(1.5) !* EXP( -screen ) !dist is always zero

              !These vectors are the distance from each respective atomic center to the location of the Gaussian that results from the intersection 
              ! of the two particular Gaussians described by a1 and a2
              PA = 0.d0  !(- a2 / (a1 + a2)) * dist_vec
              PB = 0.d0  !( a1/ (a1 + a2)) * dist_vec

              !The overlap formulas used all depend on the use of Cartesian Gaussian type orbitals
              !The true spherical harmonic Gaussians are then formed as linear combinations of these.
              !Thus the overlap for a particular Spherical gaussian is a linear combination of overlaps for cartesian GTO's
              !The linear combinations used come from: Schlegel and Frisch, Intl. Journal of Quant. Chem.; 54, 83-87, (1995)
              !The 'SIGMA' function is calculated for all Cartesian-GTO's that make up each Spherical-GTO contained in nu & mu
              CALL cartesian_component(AO_basis(nu),AO_basis(nu),PA,PB,a1,a2,cartesian)

              old_norm = old_norm + AO_basis(nu)%coeff(nugauss)*AO_basis(nu)%coeff(mugauss)*cartesian*integral

           ENDDO
        ENDDO

     !WRITE(6,*)nu,old_norm
     AO_basis(nu)%coeff = AO_basis(nu)%coeff / SQRT(old_norm)

  ENDDO
  !WRITE(6,*)

  !DO nu=1,s_dim
  !   AO_basis(nu)%coeff = AO_basis(nu)%coeff / SQRT(old_norm(nu))
  !ENDDO


END SUBROUTINE renorm_AO




!This subroutine calculates the overlap between GTO basis sets in two unit cells, separated by 'lvec'
!Upon exit the variable 's_mat' contains the overlap matrix.
!It should be noted that in the case of lvec=0 this is just a traditional overlap routine

!The formulas for overlap all come from: Clementi and Davis, Journal of Comp. Physics; 2,223-244(1967)
SUBROUTINE real_overlap(lvec,AO_basis,s_mat)
  IMPLICIT NONE

  REAL*8,DIMENSION(3),INTENT(IN)    ::  lvec
  TYPE(AO_function),DIMENSION(:)    ::  AO_basis 
  REAL*8,DIMENSION(:,:),INTENT(OUT) ::  s_mat

  INTEGER   ::  s_dim

  INTEGER   ::  nu,mu,igauss,jgauss
  INTEGER   ::  nugauss,mugauss

  REAL*8,DIMENSION(3)  ::  dist_vec  !Vector between the respective centers of the basis functions of interest
  REAL*8               ::  dist      !Distance between the centers containing the basis functions of interest
  REAL*8               ::  screen
  REAL*8,PARAMETER     ::  overlap_screen = 45.d0
  REAL*8,DIMENSION(3)  ::  PA,PB  !Vectors from each nuclei to the center of the resulting gaussian
  REAL*8               ::  a1,a2  !Shorthand for the exponents of the basis functions of interest
  REAL*8, DIMENSION(3) ::  sigma  !Vector containg the 'sigma function' described in the Clementi and Davis paper
  REAL*8               ::  integral,cartesian  !Placeholders for 1) S-type overlap 2)factor from sigma functions that account for higher order harmonics

  INTEGER   :: i,j,k

  REAL*8, PARAMETER        ::  pi=3.141592653589793238462d0

  s_dim = SIZE(AO_basis,1)

  !All basis function pairs are looped over, and the resulting overlap added into the array s_mat
  !The transpose symmetry of molecular overlap calculations can in general not be used here, since nu and mu are not necessarily in the same unit cell
  !SO s_mat(mu,nu) IS NOT THE SAME AS s_mat(nu,mu)
  DO nu=1,s_dim
     DO mu=1,s_dim

        dist_vec = AO_basis(nu)%pos - (AO_basis(mu)%pos + lvec)
        dist = DOT_PRODUCT(dist_vec,dist_vec)

        !Overlap is defined for each Gaussian function, but each basis function is a summation of Gaussians
        !Thus, for all basis function pairs the overlap is really a sum over all possible pairs of Gaussians
        DO nugauss=1,AO_basis(nu)%ngauss
           DO mugauss=1,AO_basis(mu)%ngauss

              !The full AO_basis notation is cumbersome to use and exponents are repeatedly used.  a1 and a2 are used for aestheic reasons.
              a1 = AO_basis(nu)%alpha(nugauss)
              a2 = AO_basis(mu)%alpha(mugauss)

              screen = dist*(a1*a2)/(a1+a2)
              IF( screen > overlap_screen )CYCLE 

              !This expression is the integral of the gaussian that results from the intersection of the particular Gaussian components of nu & mu
              !This is the complete overlap expression for two S-type gaussians.  But higher order harmonics complicate the calculation
              integral = AO_basis(nu)%norm(nugauss) * AO_basis(mu)%norm(mugauss) * (pi/(a1+a2))**(1.5) * EXP( -screen )

              !These vectors are the distance from each respective atomic center to the location of the Gaussian that results from the intersection 
              ! of the two particular Gaussians described by a1 and a2
              PA = (- a2 / (a1 + a2)) * dist_vec
              PB = ( a1/ (a1 + a2)) * dist_vec

              !The overlap formulas used all depend on the use of Cartesian Gaussian type orbitals
              !The true spherical harmonic Gaussians are then formed as linear combinations of these.
              !Thus the overlap for a particular Spherical gaussian is a linear combination of overlaps for cartesian GTO's
              !The linear combinations used come from: Schlegel and Frisch, Intl. Journal of Quant. Chem.; 54, 83-87, (1995)
              !The 'SIGMA' function is calculated for all Cartesian-GTO's that make up each Spherical-GTO contained in nu & mu
              CALL cartesian_component(AO_basis(nu),AO_basis(mu),PA,PB,a1,a2,cartesian)

              s_mat(nu,mu) = s_mat(nu,mu) + AO_basis(nu)%coeff(nugauss)*AO_basis(mu)%coeff(mugauss)*cartesian*integral

           ENDDO
        ENDDO

        !WRITE(6,*)'real space overlap for basis functions',nu,mu
        !WRITE(6,'(D15.6)')s_mat(nu,mu)
        !WRITE(6,*)
     ENDDO
  ENDDO


END SUBROUTINE real_overlap




!The subroutine computes the SIGMA function of Clementi and Davis, Journal of Comp. Physics; 2,223-244(1967)
!This takes into account the effects of the higher order harmonics of the GTO in the overlap integral
!Since the Spherical-GTO's are linear combinations of the Cartesian-GTO's, the SIGMA function is computed for all component Cartesian-GTOs
SUBROUTINE cartesian_component(AO_nu,AO_mu,PA,PB,a1,a2,cartesian)
  IMPLICIT NONE

  TYPE(AO_function),INTENT(IN)       ::   AO_nu,AO_mu
  !INTEGER,DIMENSION(3,3),INTENT(IN)  ::   nu_lmn_mat,mu_lmn_mat
  !REAL*8, DIMENSION(3),INTENT(IN)    ::   nu_cart_coeff,mu_cart_coeff 
  REAL*8,DIMENSION(3),INTENT(IN)     ::   PA,PB
  REAL*8,INTENT(IN)                  ::   a1,a2
  REAL*8,INTENT(OUT)                 ::   cartesian

  REAL*8,DIMENSION(3)    ::   sigma_vec
  INTEGER              ::  roof
  REAL*8               ::  f_sigma     !Holds the result of the f_function calculation for each i in each sigma summatio

  INTEGER  ::  i,k
  INTEGER  :: nucart, mucart

  !The component Cartesian-GTOs are looped over for each Spherical-GTO nu & mu
  !For up to f-level, there is a maximum of three Cartesian-GTO components to each Spherical_GTO (more are needed for higher levels)
  !However, not all Spherical-GTOs have this many components so unnecessary calculations are skipped
  cartesian = 0
  DO nucart=1,AO_nu%ncart
     DO mucart=1,AO_mu%ncart
  
        sigma_vec = 0
        DO k=1,3
           roof = (AO_nu%cart_mat(nucart,k) + AO_mu%cart_mat(mucart,k)) / 2
           DO i=0,roof
              CALL f_funct(2*i, AO_nu%cart_mat(nucart,k), AO_mu%cart_mat(mucart,k), PA(k), PB(k), f_sigma)
             sigma_vec(k) = sigma_vec(k) + f_sigma * (dble_factor(2*i-1) / (2.d0*(a1+ a2))**i)
           ENDDO
        ENDDO

        cartesian = cartesian + PRODUCT(sigma_vec)*AO_nu%cart_coeff(nucart)*AO_mu%cart_coeff(mucart)
  
     ENDDO 
  ENDDO
  !WRITE(6,*)


END SUBROUTINE



!The subroutine computes the  f_function of Clementi and Davis, Journal of Comp. Physics; 2,223-244(1967)
SUBROUTINE f_funct(j, l, m, a, b, f_output)
  IMPLICIT NONE

  INTEGER,INTENT(IN)           ::  j    !Index of the f_function
  INTEGER,INTENT(IN)           ::  l,m
  REAL*8,INTENT(IN)            ::  a,b
  REAL*8,INTENT(OUT)           ::  f_output  !output from function

  INTEGER                      ::  fmax, fmin   !determines the upper and lower limit for the summation
  REAL*8                       ::  coeff1, coeff2   !resulting binomial 
  INTEGER                      ::  i

  f_output = 0.d0

  fmax = MIN(j,l)
  fmin = MAX(0,j-m)

  DO i=fmin,fmax
     f_output = f_output + binom_coeff(l,i) * binom_coeff(m,j-i) * a**(l-i) * b**(m+i-j)
  ENDDO


ENDSUBROUTINE f_funct




!This subroutine sets up an array 'index_l' that contains the indices of all real space unit cells that must be included in calculating real space overlaps
SUBROUTINE lvec_init(index_l,num_l,AO_basis,latt_vec)
  IMPLICIT NONE

  INTEGER, ALLOCATABLE, INTENT(OUT)            ::  index_l(:,:)
  TYPE(AO_function), DIMENSION(:), INTENT(IN)  ::  AO_basis
  INTEGER, INTENT(OUT)                         ::  num_l
  REAL*8, DIMENSION(3,3), INTENT(IN)           ::  latt_vec

  INTEGER,DIMENSION(3)       ::  n_lvec
  INTEGER,ALLOCATABLE        ::  l_test(:,:)       !All l-vectors possible will first be calculated, and filled here
  INTEGER                    ::  move_screen,screen_tally !Then they will each be tested for how many unitcells shifted away from 0. How many l-vectors to screen
  LOGICAL,ALLOCATABLE        ::  screen(:)         !Whether or not a particular l-vector i

  !LOGICAL                    ::  surface
  INTEGER                    ::  surface_norm

  INTEGER :: il,i,j,k


  !Test if the system is a surface, based on the k-point mesh
  !CALL surface_test(surface,surface_norm)

  !How far it takes for the most diffuse basis function to die off to an acceptable level and how many unit cells this corresponds to in each direction
  CALL cutoff(n_lvec,AO_basis,latt_vec)
  !n_lvec(3) = 5
  num_l = PRODUCT(n_lvec)
  WRITE(6,*)'Initial number of real space unit cells ',num_l

  ALLOCATE(l_test(3,num_l))
  ALLOCATE(screen(num_l))

  !Based on the total number of unitcells to move in each direction the indices of each l-vector are computed
  !The l-vectors start at the most negative and end with the most positive
  !Thus the l-vector of (0,0,0) is il=(num_l+1)/2
  !How I determine the number of vectors in each direction in lvec_init guarantees an odd number of vectors in total so all arithmetic will gauarantee integers
  il = 0
  DO i=1,n_lvec(1)
     DO j=1,n_lvec(2)
        DO k=1,n_lvec(3)
           il = il + 1
           l_test(1,il) = i - (n_lvec(1)+1)/2
           l_test(2,il) = j - (n_lvec(2)+1)/2
           l_test(3,il) = k - (n_lvec(3)+1)/2
        ENDDO
     ENDDO
  ENDDO

  !The above procedure generates a cube of possible unit cells to use in the overlap calculation, centered at the central unit cell.
  !The unit cells at the corners of this cube are actually much further away from the central unit cell, than those in the center of a face.
  !The overlap with these corner unit cells is certainly minimal and we will therefore screen them out.
  !The resutling set of unit cells will more closely resemble a sphere.

  !To actually screen unit cells, we will limit how many orthogonal shifts the particular unit cell can be from the central position. SUM(ABS(l_test(:,il)))
  !The variable 'move-screen' is set at the maximum allowable position.
  !This is set as one more than the distance from the center to the face.
  !However, this screen must always be at least 3 to include all 27 nearest neighbor unit cells
  move_screen = MAX(3, (MAXVAL(n_lvec)+3)/2)
  !move_screen = 0.d0
  !DO k=1,3
  !   move_screen = move_screen + CEILING(((n_lvec(k)-1)/2)/SQRT(3.d0))
  !ENDDO

  !WRITE(6,*)'move_screen',move_screen
  !WRITE(6,*)

  !Then all possible unit cells (l_test) are screened to exclude removed 'corner-type' unit cells
  !Additionally if the system is a surface, we exclude any unit cells out of the plane of the surface direction.
  screen_tally=0
  DO il=1,num_l
     !IF( SUM(ABS(l_test(:,il))) > move_screen )THEN
     !   screen_tally = screen_tally + 1
     !   screen(il) = .TRUE.
     !ELSEIF( surface .AND. ABS(l_test(surface_norm,il)) .EQ. (n_lvec(surface_norm)-1)/2 )THEN
     !   screen_tally = screen_tally + 1
     !   screen(il) = .TRUE.
     !ELSE
     screen(il) = .FALSE.
     !ENDIF
  ENDDO

  !After screening is complete, the total number of unit cells is adjusted ('num_l') and 'index_l' is filled
  num_l = num_l - screen_tally
  WRITE(6,*)'Screened number of real space unit cells',num_l
  ALLOCATE(index_l(3,num_l))

  screen_tally = 0
  DO il=1,num_l
     DO
       screen_tally = screen_tally + 1
       IF( .NOT. screen(screen_tally) )EXIT
     ENDDO
     index_l(:,il) = l_test(:,screen_tally)
     !WRITE(6,*)index_l(:,il)
  ENDDO

  WRITE(6,*)


END SUBROUTINE lvec_init




!This subroutine calculates the number of real space unit cells that must be included in each direction to ensure a certain overlap minimum
!This is done by first finding the GTO, with the smallest exponent.
!The necessary distance for this functiuon to die off to a desired level is then converted into unit cell amounts.
SUBROUTINE cutoff(n_lvec,AO_basis,latt_vec)
  IMPLICIT NONE

  INTEGER,DIMENSION(3),INTENT(OUT)           ::  n_lvec  !Contains the total number of unit cells that must be included in all directions on exit.
  TYPE(AO_function),DIMENSION(:),INTENT(IN)  ::  AO_basis
  REAL*8, DIMENSION(3,3),INTENT(IN)          ::  latt_vec

  INTEGER    ::  s_dim

  INTEGER,DIMENSION(2)  :: nu_min     !index of basis function with minimum exponent, second index is that of the minimum exponent
  REAL*8                :: alpha_min  !placeholder used for finding minimum exponent
  REAL*8  ::  l_cutoff           !Actual distance 
  INTEGER ::  nu,igauss,k        !Counters
  REAL*8  ::  thold              !Threshold for value of basis function to determine "far enough" away 1E-10 is used here
  REAL*8  ::  test_value         !Used in testing l_cutoff, compared to thold

  s_dim = SIZE(AO_basis,1)

  alpha_min = AO_basis(1)%alpha(1)     !This is just seeding for testing alpha minima.

  !First the smallest exponent of any single gaussian is found among the whole basis set.
  !Since the exponential portion dominates in the long range, this is the sole factor in determining the most diffuse function
  DO nu=1,s_dim
     DO igauss=1,AO_basis(nu)%ngauss
        IF( AO_basis(nu)%alpha(igauss) .LE. alpha_min )THEN
            alpha_min = AO_basis(nu)%alpha(igauss)
            nu_min = (/nu,igauss/)
            !WRITE(6,*)'new min alpha',alpha_min,nu,igauss
        ENDIF
     ENDDO     
  ENDDO

  !In the case of SP formatting, both have the same exponent values, but the s-type will be found above since it is listed first.
  !However, the p-type are more diffuse and should actually be tested.
  !Here we test for SP-formatting and adjust nu_min and sum_lmn in accordance.
  !IF( AO_basis(nu_min(1)+1)%l == AO_basis(nu_min(1))%l + 1)THEN
  !    IF( AO_basis(nu_min(1))%alpha(nu_min(2)) == AO_basis(nu_min(1)+1)%alpha(nu_min(2)))THEN
  !        nu_min(1) = nu_min(1) + 1
  !    ENDIF
  !ENDIF

  WRITE(6,*)'found lowest exponent'
  WRITE(6,*)'nu_min',nu_min
  WRITE(6,*)'exponent is ',AO_basis(nu_min(1))%alpha(nu_min(2))
  WRITE(6,*)

  !Now l_cutoff is calculated.
  !An open loop is used, counting up.
  !Once a distance is found where the value of the GTO is small enough, the loop is exited.
  thold = 1.d-12
  l_cutoff = 0.d0
  DO 
     l_cutoff = l_cutoff + 1.d0
     test_value = l_cutoff**AO_basis(nu_min(1))%l * AO_basis(nu_min(1))%coeff(nu_min(2)) * AO_basis(nu_min(1))%norm(nu_min(2)) &
                   & * EXP(-AO_basis(nu_min(1))%alpha(nu_min(2))*l_cutoff**2)
     IF( ABS(test_value) < thold )THEN
        EXIT
     ENDIF
  ENDDO 

  WRITE(6,*)'Real space distance cutoff, bohr',l_cutoff

  !Using l_cutoff the number of unit cells necessary in each direction is calculated in n_lvec(:)
  !No matter what at least 3 unit cells in each direction are used which corresponds to a cube of only nearest neighbors.
  !This shell is always necessary to capture edge effects
  DO k=1,3
     n_lvec(k) = CEILING(2.d0*l_cutoff / SQRT(SUM(latt_vec(k,:)*latt_vec(k,:)))) + 1
     n_lvec(k) = 2*n_lvec(k) - 1  !Include unit cells in the negative direction
  ENDDO

  WRITE(6,'(A,3I)')' Numer of unit cells needed in each direction',n_lvec

  WRITE(6,*)



END SUBROUTINE cutoff


END MODULE bloch_overlap
