!This is an NBO analysis code appropriate for periodic systems

PROGRAM nbo_main
  USE periodic_matutil
  USE nbo_shared
  USE matutil
  USE pre_nao
  USE nao
  USE nbo
  USE visual

  IMPLICIT NONE

  !Command line arguments
  CHARACTER*128 :: buffer

  !Define total number of atoms, basis functions, and (relevant) real-space matrices (i.e. overlap and density)
  !See PRB, 61, 16440, "Linear-scaling DFT with Gaussian orbitals...", Kudin and Scuseria for details of the
  !real-space implementation of periodic SCF and definitions of matrices S^{og}_{mu nu} and P^{og}_{mu nu}
  INTEGER :: natom,nbasis,ng,nk,nspins
  INTEGER,ALLOCATABLE ::  nnbo(:)
  !The mapping of the basis functions to atoms (i.e. first basis function for each atom)
  !as well as the angular momentum quantum number (l) for each basis function
  INTEGER,DIMENSION(:),ALLOCATABLE,TARGET :: ibasismap,ishellmap,ilmap,immap,iatnum
  REAL*8,DIMENSION(:),ALLOCATABLE,TARGET  :: iatval
  !The mapping of the ng matrices into spatial cells
  INTEGER,DIMENSION(:,:),ALLOCATABLE,TARGET ::  indexg
  !The k-points to used, as coefficients of reciprocal lattice vectors
  REAL*8,DIMENSION(:,:),ALLOCATABLE,TARGET  ::  kpt
  !The weight of each k-point as used in the plane wave calculation.  From space group of the Brillouin zone
  REAL*8,DIMENSION(:),ALLOCATABLE,TARGET    ::  kpt_wt
  !A listing of the symbol of each atom
  CHARACTER*2,DIMENSION(:),ALLOCATABLE,TARGET :: symbols
  !Will hold the calculated real-space density and overlap matrices, rho^{0g}_{mu nu}
  REAL*8,DIMENSION(:,:,:,:),ALLOCATABLE,TARGET :: rho0, fock0
  REAL*8,DIMENSION(:,:,:),ALLOCATABLE,TARGET :: s0
  !The initial k-space density and overlap matrices
  COMPLEX*16,DIMENSION(:,:,:,:),ALLOCATABLE,TARGET :: rhok, fockk
  COMPLEX*16,DIMENSION(:,:,:),ALLOCATABLE,TARGET :: sk
  REAL*8,DIMENSION(:,:,:),ALLOCATABLE,TARGET :: transform  !Will hold the real space transformation matrix to go from NHO's back to AO basis

  REAL*8,ALLOCATABLE      ::   output_coeff(:,:,:)  !Hold the coefficients of each orbital in the AO basis, can then be used for visualization
  COMPLEX*16,ALLOCATABLE  ::   bloch_coeff(:,:,:,:)
  COMPLEX*16,ALLOCATABLE  ::   fock_nbo(:,:,:,:), rho_nbo(:,:,:,:)
  REAL*8,ALLOCATABLE      ::   rho_dummy(:,:,:,:), real_fock_nbo(:,:,:,:)

  !A structure which encapsulates all above information
  TYPE(nbo_input) :: inp


  REAL*8,ALLOCATABLE    ::  output_occ(:)
  REAL*8                ::  energy_diff, occ_diff,perturb,occ_factor

  CHARACTER(32),ALLOCATABLE  ::   nbo_lbl(:)

  COMPLEX*16,ALLOCATABLE    ::   energy_test(:,:)
  COMPLEX*16                ::   energy_sum, occ_sum  !For testing bloch-space matrices in NBO basis

  REAL*8     :: ti,tf


  !The following variables are only passed into the NBO code for the purpose of being output for visualization.
  !Therefore I will not make them targets.
  !Any modification to utilize the atom positions should make them into targets and include a pointer in the 'inp' type contained in nbo_shared.f90
!  TYPE AO_function
!     REAL*8, ALLOCATABLE, DIMENSION(:) :: norm, alpha, coeff
!     INTEGER                :: num_gauss
!     !INTEGER, DIMENSION(3)  :: lmn
!     INTEGER                :: atom
!     REAL*8, DIMENSION(3)   :: pos
!     INTEGER                :: level  !Keeps track of what basis functions share the same sets of exponents, on the same atom.
!     INTEGER                :: l,m     !l- and m-quantum numbers of basis function.  m actually runs from 0,2*l+1 and is onyl and index.
!
!     INTEGER                :: ncart
!     INTEGER,ALLOCATABLE    :: cart_mat(:,:)
!     REAL*8,ALLOCATABLE     :: cart_coeff(:)
!
!  END TYPE AO_function

  TYPE(AO_function),ALLOCATABLE  ::  AO_basis(:)  !Info on all basis functions

  REAL*8, ALLOCATABLE    ::  atom_pos(:,:)   !Atomic positions within the central unit cell (in bohr?)
  REAL*8,DIMENSION(3,3)  ::  latt_vec        !Real space unit cell vectors (in bohr?)


!  TYPE vis_cont_struct
!     INTEGER               ::  vis_start,vis_end
!     INTEGER,DIMENSION(3)  ::  mesh
!     INTEGER,DIMENSION(3)  ::  box_int
!     REAL*8                ::  origin_fact
!     LOGICAL               ::  density
!  END TYPE vis_cont_struct

  TYPE(vis_cont_struct)  :: vis_control

  REAL*8,ALLOCATABLE  ::  ao_coeff(:,:,:,:)

  !End of visualization variables


  !Checkpoint file information
  LOGICAL :: checkpoint_exists, write_checkpoint
  !Configuration file control
  LOGICAL :: write_config

  !Variables used for visualization output
  LOGICAL  ::  visualize  !Control output
!  INTEGER  ::  vis_start,vis_end
!  INTEGER  ::  mesh(3),box_int(3)
!  REAL*8   ::  origin_fact
!  LOGICAL  ::  write_density !Control whether density or wavefunctions will be written



  !Temporary variables
  INTEGER :: ig,inbo,ik,i,j,k,ispin,nu
  COMPLEX*16  ::  arg

  INTEGER,ALLOCATABLE  :: nbond(:),nlp(:),nryd(:)

  REAL*8,PARAMETER           ::  pi=4.d0*ATAN(1.d0)
  COMPLEX*16,PARAMETER       ::  sqrt_minus_one=(0.d0, 1.d0)
  REAL*8,PARAMETER           ::  aukcal=627.509469d0


  CALL CPU_TIME(ti)

  checkpoint_exists=.FALSE.
  write_checkpoint=.FALSE.

  !Do IO; first read the input file
  CALL GETARG(1,buffer)
  OPEN(10, FILE=buffer, STATUS='old')
  CALL read_input_file
  CLOSE(10)

  !See if a checkpoint file was specified
  IF (IARGC().GT.1) THEN
     CALL GETARG(2,buffer)
     OPEN(10, FILE=buffer, STATUS='old', FORM='unformatted', ERR=20)

     PRINT *, 'Reading density matrix in NAO basis from checkpoint file: ', buffer
     READ(10) inp%rho0
     READ(10) inp%transform
     READ(10) inp%fock0
     inp%s0=periodic_matiden(inp%nbasis)
     checkpoint_exists=.TRUE.
     GOTO 30

20   write_checkpoint=.TRUE.

30 ENDIF

  !Check number of electrons and do transform to an orthogonal basis
  CALL calc_nelec(inp,checkpoint_exists)

  IF (.NOT.checkpoint_exists) THEN
     !Convert to a symmetrically orthoganlized basis, weighted by the occupations of the pre-NAOs.
     !This is the NAO basis.
     CALL do_pre_nao(inp)
  ENDIF
  IF (write_checkpoint) THEN
     PRINT *, 'Writing density matrix in NAO basis to checkpoint file: ', buffer
     OPEN(10, FILE=buffer, STATUS='new', FORM='unformatted')
     WRITE(10) inp%rho0
     WRITE(10) inp%transform
     WRITE(10) inp%fock0
  ENDIF
  CLOSE(10)

  !Now get down to buisiness; first output the NAOs (this is easy, since we are already in the
  !NAO basis due to the prior transformations
  CALL do_nao(inp)
  
  !***This has been removed so that all one center orbitals are treated as lone pairs.
  !***The lone pairs were not actually projected, but acutally depleted. Less rigorous
  !***It was actually not even depletion (which would be equivalent in this case) just zeroing diagonals
  !Project off the core electrons (with occupancy ~2) from the density matrix
  !CALL do_core_projection(inp)
  

  !Now we will set up some of the parameters for the NBO search
  !I am going to hard code in some values, but they can also be determined by the user via the nbo.config file
  !The code will look for this file and read it if it exists
  !If the file does not exist upon run, defaults will be used and an nbo.config file will be generated.
  write_config = .FALSE.
  OPEN(40, FILE='nbo.config', STATUS='old', ERR=50)

  !PRINT *, 'Reading in parameters from existing nbo.config file. '
  READ(40,*)buffer
  READ(40,*)buffer
  READ(40,*)buffer
  READ(40,*)buffer
  READ(40,*)nbo_1c_thresh
  READ(40,*)nbo_2c_thresh
  READ(40,*)buffer
  READ(40,*)visualize
  READ(40,*)vis_control%density
  READ(40,*)vis_control%vis_start,vis_control%vis_end
  READ(40,*)vis_control%mesh
  READ(40,*)vis_control%box_int
  READ(40,*)vis_control%origin_fact
  !Here you could read in more parameters if we want to add some functionality
  CLOSE(40)

!  WRITE(6,*)buffer
!  WRITE(6,*)nbo_1c_thresh
!  WRITE(6,*)nbo_2c_thresh
!  WRITE(6,*)visualize
!  WRITE(6,*)vis_control%density
!  WRITE(6,*)vis_control%vis_start,vis_control%vis_end
!  WRITE(6,*)vis_control%mesh
!  WRITE(6,*)vis_control%box_int
!  WRITE(6,*)vis_control%origin_fact

  !STOP

  !Now should make sure those values all make sense

  !Make sure cutoffs aren't negative
  IF( nbo_1c_thresh.LT.0.d0 .OR. nbo_2c_thresh.LT.0.d0 )THEN
     WRITE(6,'(A40)')'Found a threshhold cutoff below 0 e.'
     WRITE(6,'(A13,F6.2)')'One center:  ',nbo_1c_thresh
     WRITE(6,'(A13,F6.2)')'Two center:  ',nbo_2c_thresh
     WRITE(6,'(A65)')'This would be impossible. Please check nbo.config'
     STOP
  ENDIF


  !Make sure the cutoffs are not too large (defintion of 'too' is spin dependent)
  IF( nspins.EQ.2 )THEN

     IF( nbo_1c_thresh.GT.1.d0 .OR. nbo_2c_thresh.GT.1.d0 )THEN
        WRITE(6,'(A70)')'Found a threshhold cutoff above 1 e for a spin polarized calculation'
        WRITE(6,'(A13,F6.2)')'One center:  ',nbo_1c_thresh
        WRITE(6,'(A13,F6.2)')'Two center:  ',nbo_2c_thresh
        WRITE(6,'(A65)')'This would be impossible. Please check nbo.config'
        STOP
     ENDIF

  ELSE

     IF( nbo_1c_thresh.GT.2.d0 .OR. nbo_2c_thresh.GT.2.d0 )THEN
        WRITE(6,'(A75)')'Found a threshhold cutoff above 2 e for a non-spin polarized calculation'
        WRITE(6,'(A13,F6.2)')'One center:  ',nbo_1c_thresh
        WRITE(6,'(A13,F6.2)')'Two center:  ',nbo_2c_thresh
        WRITE(6,'(A65)')'This would be impossible. Please check nbo.config'
        STOP
     ENDIF

     !Also make sure they are not likely cutoffs from a spin calculation
     IF( nbo_1c_thresh.LT.1.d0 .OR. nbo_2c_thresh.LT.1.d0 )THEN
        WRITE(6,'(A75)')'Found a threshhold cutoff below 1 e for a non-spin polarized calculation'
        WRITE(6,'(A13,F6.2)')'One center:  ',nbo_1c_thresh
        WRITE(6,'(A13,F6.2)')'Two center:  ',nbo_2c_thresh
        WRITE(6,'(A65)')'While possible, that is very low.  Please check nbo.config'
        STOP
     ENDIF

     !Adjust the stored value of the threshhold 
     !I use one threshhold evaluation for all spins and multiply the stored value by (3-nspins)
     !So I need the read in values to be divided by two, for when its multiplied by two later...
     nbo_1c_thresh = nbo_1c_thresh / 2.d0
     nbo_2c_thresh = nbo_2c_thresh / 2.d0

  ENDIF


  !Make sure visualization parameters make sense
  IF( visualize )THEN

     IF( vis_control%vis_start.LT.1 )THEN
        WRITE(6,*)'Lower bound for NBOs to visualize must be positive.'
        WRITE(6,*)'Vis_start: ',vis_control%vis_start
        WRITE(6,*)'Please check nbo.config'
        STOP
     ENDIF

     IF( vis_control%vis_end.GT.nbasis )THEN
        WRITE(6,*)'Upper bound for NBOs to visualize is larger than total number of AOs.'
        WRITE(6,*)'Vis_end:   ',vis_control%vis_end
        WRITE(6,*)'# of NBOs: ',nbasis
        WRITE(6,*)'Please check nbo.config'
        STOP
     ENDIF

     IF( vis_control%vis_end.LT.vis_control%vis_start )THEN
        WRITE(6,*)'Start of visualization range is larger than end point'
        WRITE(6,*)'Vis_start: ',vis_control%vis_start
        WRITE(6,*)'Vis_end:   ',vis_control%vis_end
        WRITE(6,*)'Please check nbo.config'
     ENDIF

     IF( MINVAL(vis_control%mesh).LT.1 )THEN
        WRITE(6,*)'All values of mesh must be positive'
        WRITE(6,*)'Mesh: ',vis_control%mesh
        WRITE(6,*)'Please check nbo.config'
        STOP
     ENDIF


  ENDIF


  GOTO 60

50  write_config = .TRUE.

  !Default values if nothing has been read in 

  !Occupancy cutoffs
  nbo_1c_thresh=0.80d0
  nbo_2c_thresh=0.925d0
  !Visualization controls
  visualize=.FALSE.
  vis_control%density=.FALSE.
  vis_control%vis_start=0
  vis_control%vis_end=-1
  vis_control%mesh=(/0,0,0/)
  vis_control%box_int=(/1,1,1/)
  vis_control%origin_fact=0.d0

  OPEN(60, file='nbo.config', STATUS='new')
  WRITE(60,'(A69)')'#####################################################################'
  WRITE(60,'(A69)')'#####Configuration file for the periodic NBO code of JRS and BDD#####'
  WRITE(60,'(A69)')'#####################################################################'
  WRITE(60,*)
  WRITE(60,'(A30)')'#####NBO search parameters#####'
  WRITE(60,'(F6.2,A40)')nbo_1c_thresh*DBLE(3-nspins),'  #Occupancy cutoff for one-center NBOs'
  WRITE(60,'(F6.2,A40)')nbo_2c_thresh*DBLE(3-nspins),'  #Occupancy cutoff for two-center NBOs'
  WRITE(60,'(A48)')'#####Visualization output control parameters#####'
  WRITE(60,'(L6,A59)')visualize,'  #Control over printing of .cube files for visualization.'
  WRITE(60,'(L6,A70)')vis_control%density,'  #density - Whether density (T) or wavefunctions (F) are visualized.'
  WRITE(60,'(2I3,A70)')vis_control%vis_start,vis_control%vis_end,'  #vis_start vis_end - Start and end of NBOs to print .cube files for'
  WRITE(60,'(3I2,A76)')vis_control%mesh,'  #mesh - Number of points along each lattice vectors to use in .cube files'
  WRITE(60,'(3I2,A85)')vis_control%box_int,'  #box_int - Number of unit cell to use for .cube file. See READ_ME.txt for guidance'
  WRITE(60,'(F6.2,A82)')vis_control%origin_fact,'  #origin_fact - Shift of the origin for .cube file. See READ_ME.txt for guidance'
  WRITE(60,*)
  CLOSE(60)


60  CONTINUE


  WRITE(6,'(A)')' ******************************* '
  WRITE(6,'(A)')' ****** SEARCH PARAMETERS ****** '
  WRITE(6,'(A)')' ******************************* '
  IF( .NOT.write_config )THEN
     WRITE(6,*)'Parameters read in from existing nbo.config file. '
  ELSE
     WRITE(6,*)'Default parameters. To change, modify the nbo.config file. '
  ENDIF
  WRITE(6,'(F5.2,A40)')nbo_1c_thresh*DBLE(3-nspins),'  #Occupancy cutoff for one-center NBOs'
  WRITE(6,'(F5.2,A40)')nbo_2c_thresh*DBLE(3-nspins),'  #Occupancy cutoff for two-center NBOs'
  WRITE(6,*)
  WRITE(6,*)


  ALLOCATE(rho_dummy(nbasis,nbasis,ng,nspins))
  rho_dummy = inp%rho0

  !Now do the final NBO analysis
  !ALLOCATE(bond_out(inp%nspins),lp_out(inp%nspins),ryd_out(inp%nspins))
  ALLOCATE(output(nspins),nnbo(nspins))
  DO ispin=1,nspins
     CALL do_nbo(inp,ispin,nnbo(ispin))
  ENDDO
  CALL CPU_TIME(tf)
  WRITE(6,*)'total time for NBO analysis',SNGL(tf-ti)


  !Check to make the sure a full number of NBO's was found.
  !Basically did a bond get rejected after orthogonalization.
  !A square matrix is necessary for the unitary transforms.
  DO ispin=1,nspins
     IF( nnbo(ispin) /= nbasis )THEN
        WRITE(6,*)'Incorrect number of NBOs for spin',ispin
     ENDIF
  ENDDO




  !!!!!!!!!!!!!!!!!
  !This is STOP #1!
  !!!!!!!!!!!!!!!!!
  !STOP


  !This section of the code is for visualization purposes.
  IF( .NOT.visualize )GOTO 90



  WRITE(6,*)
  WRITE(6,*)
  WRITE(6,'(A)')' ******************************* '
  WRITE(6,'(A)')' ******** VISUALIZATION ******** '
  WRITE(6,'(A)')' ******************************* '
  WRITE(6,*)
  IF( vis_control%density )THEN
     WRITE(6,*)'Densities will be written out for NBOs in the following range'
  ELSE
     WRITE(6,*)'Wavefunctions will be written out for NBOs in the following range'
  ENDIF
  WRITE(6,'(A10,I5,A8,I5)')'Start:  ',vis_control%vis_start,'End:  ',vis_control%vis_end
  WRITE(6,*)


  ALLOCATE(ao_coeff(nbasis,nnbo(1),ng,nspins))

  !Convert the nbo orbitals with coeff in the NAO basis into coeff in the AO basis for visualization.
  DO ispin=1,nspins
     DO inbo=1,nnbo(ispin)
        !output(ispin)%coeff(:,inbo,:) = periodic_matvecmul(transform,output(ispin)%coeff(:,inbo,:))
        ao_coeff(:,inbo,:,ispin) = periodic_matvecmul(transform,output(ispin)%coeff(:,inbo,:))
     ENDDO
  ENDDO

  CALL NBO_visualization(AO_basis,ao_coeff,atom_pos,indexg,latt_vec,iatnum,vis_control)

 
  !CALL NBO_visualization
 
  STOP

  !Now actually write out a file to be read in for use in visualization
  OPEN(95,file='nbo_vis.out')
  WRITE(95,*)"Output of lattice vector and ao coeffs from periodic NBO code"
  WRITE(95,*)

  !System information
  WRITE(95,*)natom, '! number of atoms in central unit cell'
  WRITE(95,*)nbasis, '! number of basis functions per unti cell'
  WRITE(95,*)nnbo(1), '! number of popssible lonepairs and NBOs, set of coefficients'
  WRITE(95,*)ng, '! number of l_vectors, unit cell pairs'
  WRITE(95,*)nspins,'! number of unique spins for which NBOs have been obtained'
  WRITE(95,*)

  !Lattice cell vectors
  DO j=1,3
     WRITE(95,*)latt_vec(:,j)
  ENDDO
  WRITE(95,*)

  !Atomic information
  DO j=1,natom
     WRITE(95,*)iatnum(j)
     WRITE(95,*)atom_pos(:,j)
  ENDDO
  WRITE(95,*)

  !Real space unit cells
  DO ig=1,ng
     WRITE(95,*)indexg(:,ig)
  ENDDO
  WRITE(95,*)

  !Information on the basis set
  DO nu=1,nbasis
     WRITE(95,'(I)')AO_basis(nu)%num_gauss
     WRITE(95,'(10D17.9)')AO_basis(nu)%alpha
     WRITE(95,'(10D17.9)')AO_basis(nu)%coeff
     WRITE(95,'(10D17.9)')AO_basis(nu)%norm
     WRITE(95,'(3D16.8)')AO_basis(nu)%pos
     WRITE(95,'(I)')AO_basis(nu)%ncart
     DO j=1,AO_basis(nu)%ncart
        WRITE(95,'(D17.9,3I)')AO_basis(nu)%cart_coeff(j),AO_basis(nu)%cart_mat(j,:)
     ENDDO

     WRITE(95,*)
  ENDDO


  !Coefficients of each NBO in the AO basis
  DO ispin=1,nspins
     DO ig=1,ng
        DO inbo=1,nnbo(1)
           WRITE(95,*)output(ispin)%coeff(:,inbo,ig)
        ENDDO
        WRITE(95,*)
     ENDDO
  ENDDO

  CLOSE(95)

  WRITE(6,*)'done with writing vis output'



90  CONTINUE

  !!!!!!!!!!!!!!!!!
  !This is STOP #2!
  !!!!!!!!!!!!!!!!!
  STOP


  !This section of the code is for energy analysis in the NBO basis, which is currently not numerically meaningful.
  !The only useful thing is that the sum of the trace of [rho*fock] summed over k-points should still match the original VASP value.

  !Check to make the sure a full number of NBO's was found. 
  !Basically did a bond get rejected after orthogonalization.
  !A square matrix is necessary for the unitary transforms.
  DO ispin=1,nspins
     IF( nnbo(ispin) /= nbasis )THEN
        WRITE(6,*)'Incorrect number of NBOs for spin',ispin
     ENDIF
  ENDDO

  ALLOCATE(rho_nbo(nbasis,nbasis,nk,nspins),fock_nbo(nbasis,nbasis,nk,nspins),real_fock_nbo(nbasis,nbasis,ng,nspins))
  ALLOCATE(bloch_coeff(nbasis,nbasis,nk,nspins))


  inp%fockk = 0.d0
  inp%rhok = 0.d0
  bloch_coeff = 0.d0

  IF( real_init )STOP 'Real space intialized calculations are not setup for past bond search'

  !Go to Bloch space (unitary transform will then be straight forward)
  DO ispin=1,nspins
     CALL real_to_bloch(inp%fockk(:,:,:,ispin),inp%fock0(:,:,:,ispin),inp%kpt,inp%indexg)
     CALL real_to_bloch(inp%rhok(:,:,:,ispin),rho_dummy(:,:,:,ispin),inp%kpt,inp%indexg)
     CALL real_to_bloch(bloch_coeff(:,:,:,ispin),output(ispin)%coeff,inp%kpt,inp%indexg)
  ENDDO

  ALLOCATE(energy_test(nbasis,nbasis))
  energy_sum = 0.d0
  occ_sum = 0.d0

  !Perform unitary transform and multiplications to test sum of orbital energies
  DO ik=1,nk
     DO ispin=1,nspins
        fock_nbo(:,:,ik,ispin)=matunitary_trans(inp%fockk(:,:,ik,ispin),TRANSPOSE(bloch_coeff(:,:,ik,ispin)))
        rho_nbo(:,:,ik,ispin)=matunitary_trans(inp%rhok(:,:,ik,ispin),TRANSPOSE(bloch_coeff(:,:,ik,ispin)))

        energy_test = MATMUL(fock_nbo(:,:,ik,ispin),rho_nbo(:,:,ik,ispin))
        DO j=1,nnbo(ispin)
           energy_sum = energy_sum + energy_test(j,j)*inp%kpt_wt(ik)
           occ_sum = occ_Sum + rho_nbo(j,j,ik,ispin)*inp%kpt_wt(ik)
        ENDDO
     ENDDO
  ENDDO
  WRITE(6,*)'energy test in nbo basis set    ',energy_sum  
  WRITE(6,*)'occupancy test in nbo basis set ',occ_sum
  WRITE(6,*)

  !Then take matrixces in NBO basis back to real space
  DO ispin=1,nspins
     CALL bloch_to_real(inp,fock_nbo(:,:,:,ispin),real_fock_nbo(:,:,:,ispin),inp%kpt,inp%indexg)
     CALL bloch_to_real(inp,rho_nbo(:,:,:,ispin),rho_dummy(:,:,:,ispin),inp%kpt,inp%indexg)
  ENDDO

  !Dump out some information on the NBO's (diagonal amtrix elements)
  DO ispin=1,nspins
  IF( nspins .GT. 1 )WRITE(6,*)"NBOs for spin",ispin
  WRITE(6,'(A,4I5)')'fock matrix of nbos in central unit cell'
  DO inbo=1,nnbo(ispin)
     WRITE(6,*)inbo,real_fock_nbo(inbo,inbo,1,ispin)
  ENDDO
  WRITE(6,*)
  
  WRITE(6,*)'density matrix of nbos in central unit cell'
  DO inbo=1,nnbo(ispin)
     WRITE(6,*)inbo,rho_dummy(inbo,inbo,1,ispin)
  ENDDO
  WRITE(6,*)
  ENDDO


  !!!!!!!!!!!!!!!!!
  !This is STOP #3!
  !!!!!!!!!!!!!!!!!
  STOP


  !This is second order perturbation analysis.
  !This suffers from very high numerical inaccuracies based on the projection process.
  !From testing, the results seem to be what you would expect
  !i.e. back bonding from metal surfaces into pi* orbitals being the largest factor
  !But I don't trust any of it.  
  !This relies on the real space Fock matrices obtained in the above step

  DO ispin=1,nspins

  IF( nspins .GT. 1 )WRITE(6,*)'perturbation analysis for spin',ispin

  !Looking at perturbations driven donations from j -> i nbo 
  DO j=1,nnbo(ispin)

     !Make sure the orbital is capable of donation, based on a high occupancy
     IF( output(ispin)%occ(j) < 0.7d0*DBLE(3-nspins) )GOTO 9 

     DO i=1,nnbo(ispin)

        !Check that the acceptor orbital is unique and has a low enough occupancy to accept
        IF( i == j )GOTO 40
        IF( output(ispin)%occ(i) > 0.7d0*DBLE(3-nspins) )GOTO 40

        energy_diff = real_fock_nbo(i,i,1,ispin) - real_fock_nbo(j,j,1,ispin)
        occ_diff = output(ispin)%occ(j) - output(ispin)%occ(i)

        !WRITE(6,'(A,I4,A15,I4,A15)')'looking at nbos',j,output(ispin)%label(j),i,output(ispin)%label(i)
        !WRITE(6,*)'donation is possible'
        !WRITE(6,*)
        !Since donation is possible, loop over all possible unit cells to see if a big enough delocalization exists
        DO ig=1,ng

           !The occupancy used is NOT simply the occupancy of the donor orbital
           !The perturbative energy lowering comes from a mixing of the occupancy of the two orbitals, thus their occupancies should be added
           !However if it will result in an occupancy over 2, some electrons will go to the higher energy split and cancel out the lowering effect
           !Thus occ_factor represents how many electrons will be serving to lower the energy in the resulting orbital mixing
           occ_factor = output(ispin)%occ(j) + output(ispin)%occ(i)
           IF( occ_factor > DBLE(3-nspins) )occ_factor = 2.d0*DBLE(3-nspins) - occ_factor
           perturb = occ_factor * real_fock_nbo(j,i,ig,ispin)**2 / energy_diff
           perturb = perturb * aukcal

           IF( perturb > 0.25d0*DBLE(3-nspins) )THEN
               WRITE(6,'(A,I5,A15,A,I5,A15,I4)')'Relevant perturbation found ',j,output(ispin)%label(j),'to ',i,output(ispin)%label(i),ig
               WRITE(6,*)'Energy diff and coupling',SNGL(energy_diff),SNGL(real_fock_nbo(j,i,ig,ispin))
               WRITE(6,*)'pertubation energy',SNGL(perturb),'kcal/mol'
               WRITE(6,*)
           ENDIF

        ENDDO

70      CONTINUE

        !WRITE(6,*)

40   ENDDO
9 ENDDO


  ENDDO


  STOP


CONTAINS

  !
  !Reads the input file and allocates necessary memory
  !All matrices read in are in reciprocal space, NOT real space
  !
  SUBROUTINE read_input_file
    USE matutil
    IMPLICIT NONE
    INTEGER :: ibasis,im,ik,ispin
    !LOGICAL :: real_init
    !INTEGER :: nkx,nky,nkz
    REAL*8  :: num_elec

    CHARACTER(64)  :: NBO_mat_fn
    CHARACTER(128) :: comment

      real_init = .TRUE.

      READ(10,*) !read comment line
      READ(10,*) natom
      READ(10,*) nbasis
      READ(10,*) nspins
      READ(10,*) ng

      ALLOCATE(ibasismap(natom+1))
      READ(10,*,ERR=20)ibasismap
      GOTO 30

!20    WRITE(6,*)'The input file has an entry for nk, so the matrices must be in bloch space'
20    CONTINUE
      real_init = .FALSE.
      nk = ibasismap(1)


30    CONTINUE
      !WRITE(6,*)'real start boolean   ',real_init


      IF( real_init )STOP 'NBO code is no longer compatible with Gaussian output'

      ALLOCATE(ishellmap(nbasis))
      ALLOCATE(ilmap(nbasis))
      ALLOCATE(immap(nbasis))
      ALLOCATE(symbols(natom))
      ALLOCATE(iatnum(natom))
      ALLOCATE(iatval(natom))
      ALLOCATE(indexg(3,ng))

      IF( .NOT. real_init )THEN
         ALLOCATE(kpt(3,nk))
         ALLOCATE(kpt_wt(nk))
         READ(10,*) ibasismap
      ENDIF
      READ(10,*) ishellmap
      READ(10,*) ilmap
      READ(10,*) symbols
      READ(10,*) iatnum
      IF( real_init )THEN
         iatval = iatnum
      ELSE
         READ(10,*) iatval
      ENDIF

      !WRITE(6,*)'start of new read'

      !Read in structural information
      DO im=1,3 !im is just a dummy variable in this context
         READ(10,*)latt_vec(:,im)
      ENDDO

      !WRITE(6,*)'now reading positions'

      ALLOCATE(atom_pos(3,natom))
      DO im=1,natom !im is just a dummy variable in this context
         READ(10,*)atom_pos(:,im)
      ENDDO

      !WRITE(6,*)'made it to reading g-vectors'

      !Read in information about real space g-vectors.
      !These will only be for nearest neighbor cells where bonds we will be searched for
      DO ig =1,ng
         READ(10,*) indexg(1,ig), indexg(2,ig), indexg(3,ig)
      ENDDO
      CALL periodic_matinit(ng,indexg)

      IF( real_init )THEN
         CALL get_nk(nkx,nky,nkz)
         nk=nkx*nky*nkz
      ENDIF


      ALLOCATE(rho0(nbasis,nbasis,ng,nspins))
      ALLOCATE(s0(nbasis,nbasis,ng))
      ALLOCATE(fock0(nbasis,nbasis,ng,nspins))
      ALLOCATE(transform(nbasis,nbasis,ng))
      
      ALLOCATE(rhok(nbasis,nbasis,nk,nspins))
      ALLOCATE(sk(nbasis,nbasis,nk))
      ALLOCATE(fockk(nbasis,nbasis,nk,nspins))

      !WRITE(6,*)'allocated everything'

      IF( real_init )THEN

          STOP 'Need to resetup matrix input for real space initialization systems'


!         CALL read_real_triangular(10,nbasis,s0(:,:,1),s0(:,:,1))
!         DO ig=2,ng,2
!            CALL read_real_triangular(10,nbasis,s0(:,:,ig),s0(:,:,ig+1))
!            CALL read_real_triangular(10,nbasis,s0(:,:,ig+1),s0(:,:,ig))
!         ENDDO
!         CALL read_real_triangular(10,nbasis,rho0(:,:,1),rho0(:,:,1))
!         DO ig=2,ng,2
!            CALL read_real_triangular(10,nbasis,rho0(:,:,ig),rho0(:,:,ig+1))
!            CALL read_real_triangular(10,nbasis,rho0(:,:,ig+1),rho0(:,:,ig))
!         ENDDO
!
!         OPEN(12,file='NBODATA.54')
!         DO j=1,3
!            READ(12,*)
!         ENDDO
!         fock0=0.d0
!         CALL read_real_triangular(12,nbasis,fock0(:,:,1),fock0(:,:,1))
!         !WRITE(6,*)'test of read in fock matrix'
!         !DO ibasis=1,nbasis
!         !   WRITE(6,'(27F10.5)')fock0(ibasis,:,1)
!         !   WRITE(6,*)
!         !ENDDO
!         
!         CLOSE(12)
!
!         !STOP 'Fock matrix input has not been implemented for real space start'
!
!         DO ig=1,ng
!            num_elec=num_elec + mattrace(MATMUL(rho0(:,:,ig),TRANSPOSE(s0(:,:,ig))))
!         ENDDO
!         WRITE(6,*)'number of initial electrons in real space',2.d0*num_elec
!         num_elec = 0.d0
!         DO ig=1,ng
!            num_elec=num_elec + mattrace(MATMUL(rho0(:,:,ig),TRANSPOSE(fock0(:,:,ig))))
!         ENDDO
!         WRITE(6,*)'initial real space energy sum',2.d0*num_elec
!
!         sk=0.d0
!         rhok=0.d0
!         sk=periodic_matbloch(s0,nkx,nky,nkz)
!         rhok=periodic_matbloch(rho0,nkx,nky,nkz)
!         fockk=periodic_matbloch(fock0,nkx,nky,nkz)

      ELSE

         !WRITE(6,*)'reading in k-pts'

         !Read in k-vectors indexing the overlap and density matrices.
         !These are in the same order as the matrices 
         DO ik=1,nk
            READ(10,*)kpt(:,ik),kpt_wt(ik)
            !WRITE(6,*)ik,kpt_wt(ik)
         ENDDO

         READ(10,*)NBO_mat_fn

         !WRITE(6,*)'Input matrices will be read in from',NBO_mat_fn

         OPEN(66,file=NBO_mat_fn,FORM='UNFORMATTED')

         !read in the overlap and density matrices s^{k}_{mu nu}
         !At each k-point each of these matrices are hermitian
         !Therefore only half of each matrix is written out in the file, and symmetry is used to file in the remainder
         DO ik=1,nk
            CALL read_bloch_triangular(66,nbasis,sk(:,:,ik))
         ENDDO

         DO ispin=1,nspins
            DO ik=1,nk
               CALL read_bloch_triangular(66,nbasis,rhok(:,:,ik,ispin))
            ENDDO
         ENDDO

         DO ispin=1,nspins
            DO ik=1,nk
               CALL read_bloch_triangular(66,nbasis,fockk(:,:,ik,ispin))
            ENDDO
         ENDDO

         CLOSE(66)

      ENDIF

      !Now going to read in information on the basis set.
      !This is only passed in here to then output into the visualization

      ALLOCATE(AO_basis(nbasis))
      !First read in the comment about this datas purpose
      READ(10,*)comment

      !WRITE(6,*)'reading basis set info'

      DO ibasis=1,nbasis
         READ(10,*)AO_basis(ibasis)%num_gauss
         ALLOCATE(AO_basis(ibasis)%alpha(AO_basis(ibasis)%num_gauss),AO_basis(ibasis)%coeff(AO_basis(ibasis)%num_gauss),AO_basis(ibasis)%norm(AO_basis(ibasis)%num_gauss))
         READ(10,*)AO_basis(ibasis)%alpha
         READ(10,*)AO_basis(ibasis)%coeff
         READ(10,*)AO_basis(ibasis)%norm
         READ(10,*)AO_basis(ibasis)%pos
         READ(10,*)AO_basis(ibasis)%ncart
         ALLOCATE(AO_basis(ibasis)%cart_coeff(AO_basis(ibasis)%ncart),AO_basis(ibasis)%cart_mat(AO_basis(ibasis)%ncart,3))
         DO im=1,AO_basis(ibasis)%ncart
            READ(10,*)AO_basis(ibasis)%cart_coeff(im),AO_basis(ibasis)%cart_mat(im,:)
         ENDDO
         READ(10,*)
      ENDDO



      !correct for the fact we (maybe) just read the alpha density
      rhok=rhok * DBLE(3-nspins)
   
      !Finally, initalize immap to contain the m quantum number for each basis function
      !Assumes first basis function is s-type
      immap=0
      DO ibasis=2,nbasis
         IF (ilmap(ibasis).NE.ilmap(ibasis-1).OR.ishellmap(ibasis).NE.ishellmap(ibasis-1)) THEN
            DO im=0,2*ilmap(ibasis)
               immap(ibasis+im)=im
            ENDDO
         ENDIF
      ENDDO

      !Now store these values in our "system" structure to keep them all together
      inp%natom=natom
      inp%nbasis=nbasis
      inp%ng=ng
      inp%nk=nk
      inp%nspins=nspins
      inp%ibasismap=>ibasismap
      inp%ishellmap=>ishellmap
      inp%ilmap=>ilmap
      inp%immap=>immap
      inp%indexg=>indexg
      IF( .NOT. real_init )THEN
          inp%kpt=>kpt
          inp%kpt_wt=>kpt_wt
      ENDIF
      inp%symbols=>symbols
      inp%iatnum=>iatnum
      inp%iatval=>iatval
      inp%rho0=>rho0
      inp%s0=>s0
      inp%fock0=>fock0
      inp%transform=>transform
      inp%sk=>sk
      inp%rhok=>rhok
      inp%fockk=>fockk

    END SUBROUTINE read_input_file


END PROGRAM nbo_main


