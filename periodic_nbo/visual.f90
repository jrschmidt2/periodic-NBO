MODULE visual

CONTAINS

SUBROUTINE NBO_visualization(AO_basis,ao_coeff,atom_pos,indexg,a,atomic_number,vis_control)
  USE nbo_shared  !Where the defined types are stored
  IMPLICIT NONE

  TYPE(AO_function),INTENT(IN)  ::  AO_basis(:)  !Basis functions
  TYPE(vis_cont_struct),INTENT(IN)  :: vis_control  !Info from config file. Type limits the size of the this subroutine call. 
  REAL*8,INTENT(IN)   ::  ao_coeff(:,:,:,:)   !Coefficients for each NBO in AO basis
  REAL*8,INTENT(IN)   ::  atom_pos(:,:)  !Positions of atoms in central unit cell
  INTEGER,INTENT(IN)  ::  indexg(:,:)  !Indices of unit cells to use
  REAL*8,INTENT(IN)   ::  a(3,3)  !Lattice vectors. First index is cartesian, second is lattice vector. 
  INTEGER,INTENT(IN)  ::  atomic_number(:)  !Atomic number. Makes the cube file more complete

  !System parameters
  INTEGER                ::  nbasis,nnbo,ng,nspin,natom

  !Used for constructing chunk of atoms from bulk to use
  REAL*8, DIMENSION(3)   ::  image_pos
  REAL*8, DIMENSION(3,2) ::  image_screen
  INTEGER                ::  image_count, image_tot
  LOGICAL, ALLOCATABLE   ::  image_disp(:)


  !For standardized file names
  !I make all possible file names and create them only when the file is called
  CHARACTER(20),ALLOCATABLE ::  file_names(:,:)
  INTEGER                 ::  file_num


  !Parameters for gridding
  INTEGER, DIMENSION(3)   ::  mesh  !mesh size to use
  REAL*8, DIMENSION(3,3)  ::  box   !box size to grid
  REAL*8, DIMENSION(3)    ::  r,origin!     !radius used in griding and origin for the cube file
  REAL*8                  ::  gridvol


  !For calculating AO values
  REAL*8                  ::  gaussian, gauss_coeff
  REAL*8                  ::  rsqr
  REAL*8                  ::  cartesian
  REAL*8,DIMENSION(3)     ::  r_pos
  REAL*8                  ::  screen


  !Values to plot and control for calculations
  REAL*8                  ::  density
  REAL*8                  ::  wave_funct
  INTEGER                 ::  dens_exp


  !Counters
  INTEGER                :: ig,inbo,nu,ispin
  INTEGER                 :: ix,iy,iz, i,j,k

  !Start by determining some system parameters from array dimensions
  natom = SIZE(atom_pos,2)
  nnbo = SIZE(ao_coeff,2)
  nbasis = SIZE(ao_coeff,1)
  IF( nbasis.NE.nnbo )STOP 'There are a different number of AOs and NBOs in visualization subroutine. Something is wrong'
  ng=SIZE(ao_coeff,3)
  nspin=SIZE(ao_coeff,4)

  !WRITE(6,*)'Range of NBOs to visualize'
  !WRITE(6,*)vis_control%vis_start,vis_control%vis_end
  !WRITE(6,*)

  !The density is potentaily desired outside of the central unit cell, for instance a bond that spans unit cells
  !Thus, we will numerically integrate a grid over the volume 'box' which will contain the central unit cell and surrounding space.
  !The value of box is the absolute value of the lengths of the box
  !'origin' then account for the shift of the grid origin off center

  !WRITE(6,*)'latt_vec'
  !WRITE(6,'(3F11.5)')a


  DO i=1,3
     box(:,i) = a(:,i)*DBLE(vis_control%box_int(i))
  ENDDO
  origin = vis_control%origin_fact*(a(:,1)+a(:,2)+a(:,3))
  !The mesh parameter determines how many grid points are used in each direction for the numerical integration
  !For non-cubic boxes, different components of this vector should be assigned so that resoultion is similar in all directions.
  mesh = vis_control%mesh

  WRITE(6,*)'Box dimensions in bohr'
  WRITE(6,'(3F11.5)')box
  WRITE(6,*)
  WRITE(6,*)'Origin in bohr'
  WRITE(6,'(3F11.5)')origin
  WRITE(6,*)
  WRITE(6,*)'Grid resoultion for each vector'
  WRITE(6,'(3I11)')mesh
  WRITE(6,*)


  !Final control parameter is whether density or wave functions will be plotted
  !The exponent will be applied to the wavefunction when calculating the density. 
  IF( vis_control%density )THEN
     dens_exp = 2
  ELSE
     dens_exp = 1
  ENDIF


  !Prepare the file names of all the .cube files for each NBO of each type of spin
  !If there is only one spin type, it is unecessary to give different names based on spoin type
  ALLOCATE(file_names(nnbo,nspin))
  IF( nspin .EQ. 1 )THEN
     DO inbo=1,nnbo
        IF( inbo < 10 )THEN
           WRITE(file_names(inbo,1),'(A5,I1,A5)')'nbo_',inbo,'.cube'
        ELSEIF( inbo < 100 )THEN
           WRITE(file_names(inbo,1),'(A5,I2,A5)')'nbo_',inbo,'.cube'
        ELSEIF( inbo < 1000 )THEN
           WRITE(file_names(inbo,1),'(A5,I3,A5)')'nbo_',inbo,'.cube'
        ELSE
           WRITE(6,*)'the code is not set up to write out the denisty for more than 999 NBOs'
           STOP
        ENDIF
     ENDDO
  ELSE  !For spin polarized calculations call the first spin type alpha and the second type beta NBO's
     DO inbo=1,nnbo
        IF( inbo < 10 )THEN
           WRITE(file_names(inbo,1),'(A11,I1,A5)')'alpha_nbo_',inbo,'.cube'
           WRITE(file_names(inbo,2),'(A10,I1,A5)')'beta_nbo_',inbo,'.cube'
        ELSEIF( inbo < 100 )THEN
           WRITE(file_names(inbo,1),'(A11,I2,A5)')'alpha_nbo_',inbo,'.cube'
           WRITE(file_names(inbo,2),'(A10,I2,A5)')'beta_nbo_',inbo,'.cube'
        ELSEIF( inbo < 1000 )THEN
           WRITE(file_names(inbo,1),'(A11,I3,A5)')'alpha_nbo_',inbo,'.cube'
           WRITE(file_names(inbo,2),'(A10,I3,A5)')'beta_nbo_',inbo,'.cube'
        ELSE
           WRITE(6,*)'the code is not set up to write out the denisty for more than 999 NBOs'
           STOP
        ENDIF
     ENDDO
  ENDIF

  !Calculates the grid volume for the specified mesh dimensions
  !This is used in numerically integrating the density
  CALL grid_volume(gridvol, box, mesh)
  !WRITE(6,'(A,F10.5)')' integration grid volume',gridvol
  !WRITE(6,*)

  !To aid in visualization we will want more than just the atoms in the unit cell.
  !To do this we will scan over all atom positions in all unit cells, and store those 'close' to the central unit cell
  !Start by establishing cutoffs for determining 'closeness'
  image_screen = 0.d0
  DO j=1,3
    image_screen(j,1) = MIN( origin(j), MINVAL(0.5d0*box(j,:)))
    image_screen(j,2) = MAXVAL(box(j,:)) + origin(j)
  ENDDO

  !Then every position is scanned over (each atom in all unit cells characterized by an l-vector)
  !Those within the cutoffs determined above are stored
  ALLOCATE(image_disp(natom*ng))
  image_disp = .FALSE.
  image_count = 0
  image_tot = 0
  DO j=1,natom
     !WRITE(6,*)'atom',j,atom_pos(:,j)
     DO ig=1,ng
        image_tot = image_tot + 1
        image_pos = atom_pos(:,j)
        DO k=1,3
           image_pos = image_pos + indexg(k,ig)*a(:,k)
        ENDDO
        IF( image_pos(1) > image_screen(1,1) .AND. image_pos(1) < image_screen(1,2))THEN
          IF( image_pos(2) > image_screen(2,1) .AND. image_pos(2) < image_screen(2,2))THEN
            IF( image_pos(3) > image_screen(3,1) .AND. image_pos(3) < image_screen(3,2))THEN
              !WRITE(6,*)index_l(il,:)
              !WRITE(6,*)image_pos
              image_disp(image_tot) = .TRUE.
              image_count = image_count + 1
            ENDIF
          ENDIF
        ENDIF
     ENDDO
     !WRITE(6,*)
  ENDDO
  !WRITE(6,*)'image count',image_count
  !WRITE(6,*)

  !Then a file containing the dimensions of the central unit cell as well as the positions of all atoms in the central and surrounding unit cells
  OPEN(10, file='lattice_vec.cube')
  WRITE(10, *) 'Cube file generated by write_cube subroutine'
  WRITE(10, *) 'Density'
  WRITE(10, '(I5,3F12.6)' ) image_count, 0.d0, 0.d0, 0.d0
  DO i=1,3
     WRITE(10, '(I5,3F12.6)' ) (mesh(i)/vis_control%box_int(i)), a(1,i)/(mesh(i)/vis_control%box_int(i)), a(2,i)/(mesh(i)/vis_control%box_int(i)), a(3,i)/(mesh(i)/vis_control%box_int(i))
  ENDDO
  !WRITE(10, '(I5,3F12.6)' ) mesh(2), a(1,2)/mesh(2), a(2,2)/mesh(2), a(3,2)/mesh(2)
  !WRITE(10, '(I5,3F12.6)' ) mesh(3), a(1,3)/mesh(3), a(2,3)/mesh(3), a(3,3)/mesh(3)
  image_tot = 0
  DO j=1,natom
     DO ig=1,ng
        image_tot = image_tot + 1
        IF( image_disp(image_tot) )THEN
           image_pos = atom_pos(:,j)
           DO k=1,3
              image_pos = image_pos + indexg(k,ig)*a(:,k)
           ENDDO
           WRITE(10, '(I5,4F12.6)' ) atomic_number(j), 1.0, image_pos(1), image_pos(2), image_pos(3)
        ENDIF
     ENDDO 
  ENDDO

  !The remainder of the lattice_vec.cube file is filled with appropriate dimension and formatted zero-valued density
  !This prevents VMD from having trouble opening the file and displaying the appropriate unit cell vectors
  DO j=1,PRODUCT(DBLE(mesh/vis_control%box_int))
     WRITE(10,"(E13.5)",ADVANCE='NO')0.d0
     IF( MODULO(MODULO(j,mesh(3)/vis_control%box_int(3)),6) == 0 )WRITE(10,*)
  ENDDO

  CLOSE(10)

  WRITE(6,*)'Now creating .cube file for desired orbitals'

  !Finally a .cube file for each nbo, for each spin is filled with appropriate density on the 'mesh' grid
  !Note: this grid can extend beyond the unit cell based on the box variable
  DO ispin=1,nspin
     IF( nspin.GT.1 )WRITE(6,*)'For spin type',ispin

     DO inbo=vis_control%vis_start,vis_control%vis_end
        WRITE(6,*)'Orbital Number: ',inbo

        !Start by opening the appropriate file (names created above) and writing the cube file header 
        file_num = inbo+6
        OPEN(file_num, file=file_names(inbo,ispin))
        write(file_num, *) 'Cube file generated by write_cube subroutine'
        write(file_num, *) 'Density'
        write(file_num, '(I5,3F12.6)' ) natom, origin(1), origin(2), origin(3)
        write(file_num, '(I5,3F12.6)' ) mesh(1), box(1,1)/mesh(1), box(2,1)/mesh(1), box(3,1)/mesh(1)
        write(file_num, '(I5,3F12.6)' ) mesh(2), box(1,2)/mesh(2), box(2,2)/mesh(2), box(3,2)/mesh(2)
        write(file_num, '(I5,3F12.6)' ) mesh(3), box(1,3)/mesh(3), box(2,3)/mesh(3), box(3,3)/mesh(3)
        DO j=1,natom
           WRITE(file_num, '(I5,4F12.6)' ) atomic_number(j), 1.0, atom_pos(1,j), atom_pos(2,j), atom_pos(3,j)
        ENDDO

        !Then the grid is looped over for the total box and density tabulated
        DO ix=1,mesh(1)
           DO iy=1,mesh(2)
              DO iz=1,mesh(3)
                 density = 0.d0

                 r = ((ix-1)/DBLE(mesh(1)))*box(:,1) + ((iy-1)/DBLE(mesh(2)))*box(:,2) + ((iz-1)/DBLE(mesh(3)))*box(:,3)
                 r = r + origin

                 wave_funct = 0.d0
                 !The effect of each basis function in each unit cell must be calculated for each grid point
                 DO nu=1,nbasis
                    DO ig=1,ng

                       !In many unit cells, the coefficient is zero, so there is no need to calculate a function's value
                       IF( ao_coeff(nu,inbo,ig,ispin) /= 0.d0 )THEN
                          r_pos = r - (AO_basis(nu)%pos + indexg(1,ig)*a(:,1) + indexg(2,ig)*a(:,2) + indexg(3,ig)*a(:,3))
                          rsqr = DOT_PRODUCT(r_pos,r_pos)

                          !The guassian component of nu at r is now calculated. Screening is used
                          gaussian = 0.d0
                          DO j=AO_basis(nu)%num_gauss,1,-1  !In going backwards I assume alpha's are stored in descending order
                             screen = AO_basis(nu)%alpha(j) * rsqr
                             IF( screen .GT. 45.d0 )EXIT    !In the loop each gaussian is less diffuse than the last, so all following gaussians will also be screened
                             gaussian = gaussian + AO_basis(nu)%coeff(j)*AO_basis(nu)%norm(j)*EXP(-screen)
                          ENDDO

                          !The cartesian component and multiplication with the coefficient of the AO in the NBO is only done, it the gaussian component is non-zero
                          !Since many gaussians are screened out because of length, this screens alot of unnecessary calculations
                          IF( gaussian .NE. 0.d0 )THEN
                              cartesian = 0.d0
                              DO j=1,AO_basis(nu)%ncart
                                 cartesian = cartesian + AO_basis(nu)%cart_coeff(j)*PRODUCT(r_pos**AO_basis(nu)%cart_mat(j,:))
                              ENDDO
                              wave_funct = wave_funct + gaussian*cartesian*ao_coeff(nu,inbo,ig,ispin)
                          ENDIF
                       ENDIF

                    ENDDO
                 ENDDO

                 !The wave function value is then squared to obtain a density
                 IF( wave_funct .NE. 0.d0 )THEN
                    density = wave_funct**dens_exp
                    IF( ABS(density) .LT. 1.d-30 )density = 0.d0  !Too low of values in the cube file mess up VMD.  Since they are small anyways I round down to zero
                 ENDIF

                 !The density is then written in the appropriate format of a .cube file
                 WRITE(file_num, "(E13.5)", advance="no") density
                 IF ( MOD(iz,6) == 0 )THEN
                    WRITE(file_num,*)
                 ENDIF

              ENDDO  !End of loop in mesh through fastest direction (z)

              WRITE(file_num,*)

           ENDDO
        ENDDO  !End of loop over mesh

        CLOSE(file_num)


     ENDDO  !End of loop over band

  ENDDO


END SUBROUTINE NBO_visualization

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

FUNCTION cross(a, b)
  REAL*8,INTENT (in) :: a(3), b(3)
  REAL*8 :: cross(3)

  cross(1) = a(2) * b(3) - a(3) * b(2)
  cross(2) = a(3) * b(1) - a(1) * b(3)
  cross(3) = a(1) * b(2) - a(2) * b(1)
END FUNCTION cross



END MODULE visual
