MODULE nbo
  USE nbo_shared
  USE BLAS95
  PRIVATE
  PUBLIC :: do_nbo, output

  TYPE nbo_bond
     INTEGER :: iatom,jatom
     INTEGER :: ig,jg
     INTEGER :: inho,jnho
     REAL*8 :: occ
  END TYPE nbo_bond

  INTEGER :: nnbo,maxnbo,maxao,maxbondorlpperatom
  !PARAMETER (maxbondorlpperatom=24)
  PARAMETER (maxbondorlpperatom=100)

  !INTEGER  ::  num_lp

  TYPE(nbo_output_info), ALLOCATABLE   ::  output(:)


  REAL*8  ::  lp_cutoff, bond_cutoff

  CONTAINS

    !
    !Actually do the NBO search
    !
    SUBROUTINE do_nbo(inp,ispin,nnbo_out)
      USE matutil
      USE periodic_matutil
      USE nbo_shared
      IMPLICIT NONE

      TYPE(nbo_input) :: inp

      INTEGER :: ispin   !NBO'a are found independently for each spin type
      INTEGER :: nnbo_out

      INTEGER :: ig,jg,kg,iatom,jatom,ibasis,jbasis
      INTEGER :: inbo,il
      INTEGER :: ifirst,ilast,isize,jfirst,jlast,jsize
      INTEGER :: iifirst,iilast,jjfirst,jjlast

      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: rhoat,nbo_coeff,Tnho
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: nho_coeff
      REAL*8, ALLOCATABLE, DIMENSION(:) :: nbo_occ
      INTEGER :: nnho(inp%natom),nnho_init(inp%natom)

      REAL*8 :: occ,polarization,fracl_hybrid1(0:lmax),fracl_hybrid2(0:lmax)
      CHARACTER*20 :: fmt_statement

      INTEGER :: ibondorlp,nbondorlp,nbond,nlpair,nryd
      INTEGER :: inho,jnho
      TYPE(nbo_bond),ALLOCATABLE,DIMENSION(:) :: bonds
      TYPE(nbo_bond),ALLOCATABLE,DIMENSION(:) :: lpairs
      TYPE(nbo_bond),ALLOCATABLE,DIMENSION(:) :: rydberg

      REAL*8  ::  val_occ, nl_occ, ryd_occ
      REAL*8,ALLOCATABLE ::  Pi(:,:),proj_vec(:)



      INTEGER :: j

      !WRITE(6,*)'inside of the do_nbo subroutine'

      !Find the LARGEST number of basis functions per atom; twice that is the size of the rhoat
      !matrices we need for our 2 atom NBO search
      maxao=0
      DO iatom=1,inp%natom
         ifirst=inp%ibasismap(iatom)
         ilast=inp%ibasismap(iatom+1)-1
         isize=ilast-ifirst+1
         maxao=MAX(maxao,isize)
      ENDDO
      maxnbo=maxao*2

      ALLOCATE(rhoat(maxnbo,maxnbo))
      ALLOCATE(nbo_coeff(maxnbo,maxnbo))
      ALLOCATE(Tnho(maxnbo,2))
      ALLOCATE(nbo_occ(maxnbo))
      ALLOCATE(bonds(inp%natom*maxbondorlpperatom))
      ALLOCATE(lpairs(inp%natom*maxbondorlpperatom))
      ALLOCATE(rydberg(inp%natom*maxbondorlpperatom))
      ALLOCATE(nho_coeff(maxao,maxbondorlpperatom,inp%natom))

      ALLOCATE(Pi(maxnbo,maxnbo),proj_vec(maxnbo))

      !Initialize the counts of potential bonds and per-atom-hybrids to zero
      nbondorlp=0
      nnho=0
      nnho_init=0

      nbond = 0
      nlpair = 0
      nryd = 0

      nnbo_out = 0

      lp_cutoff = nbo_1c_thresh * DBLE(3-inp%nspins)
      bond_cutoff = nbo_2c_thresh * DBLE(3-inp%nspins)
      !WRITE(6,*)'lone pair and bond cutoffs',lp_cutoff,bond_cutoff
      !WRITE(6,*)


      !WRITE(6,*)'got to calling the 1 center search'
      CALL do_1c_search(inp,lpairs,nnho,nlpair,rhoat,nbo_coeff,nbo_occ,nho_coeff,ispin)
      WRITE(6,*)"Number of non-bonding NBO's",nlpair
      CALL do_2c_search(inp,bonds,nnho,nbond,rhoat,nbo_coeff,nbo_occ,nho_coeff,ispin)
      WRITE(6,*)'Number of potential bonds',nbond
      WRITE(6,*)

      !Optionally orthogonalize the hybrids -- this is done in traditional NBO analysis, but can break the symmetry
      !of the resonance structures in hypervalent / ionic structures.  See PCCP 10, 5207 (2008) for details.
      IF (orthogonalize_hybrids) CALL do_orthogonalize_hybrids(inp,nho_coeff,nnho,.TRUE.,nnho_init)

      !WRITE(6,*)'finished orthogonalizing hybrids'

      !The factor of 2 in bond_coeff is to account for anti bonds
      !ALLOCATE(bond_out(ispin)%coeff(inp%nbasis,2*bond

      ALLOCATE(output(ispin)%occ(inp%nbasis),output(ispin)%coeff(inp%nbasis,inp%nbasis,inp%ng),output(ispin)%label(inp%nbasis))
      output(ispin)%coeff = 0.d0

      !WRITE(6,*)'allocated output arrays'

      !ALLOCATE(bond_coeff(inp%nbasis,2*nbond,inp%ng),lp_coeff(inp%nbasis,nlpair,inp%ng))
      !ALLOCATE(bond_occ(2*nbond),lp_occ(nlpair))
      !ALLOCATE(bond_lbl(2*nbond),lp_lbl(nlpair))
      !bond_coeff = 0.d0
      !lp_coeff = 0.d0

      val_occ = 0.d0
      nl_occ = 0.d0

      !Now the lonepairs will be looked at again (now of orthogonalized hybrids)
      !Information will be output to the screen and coefficients stored in lp_coeff  

      IF( inp%nspins .GT. 1 )WRITE(6,'(A,I2,A)')' *** NBO Results for spin',ispin,' ***'

      IF( nlpair > 0 )THEN

      WRITE(6,'(A)')' **************************** '
      WRITE(6,'(A)')' **** NATURAL LONE PAIRS **** '
      WRITE(6,'(A)')' **************************** '
      WRITE(6,*)
      DO ibondorlp=1,nlpair
         iatom=lpairs(ibondorlp)%iatom
         ig=lpairs(ibondorlp)%ig
         inho=lpairs(ibondorlp)%inho
         occ=lpairs(ibondorlp)%occ
         ifirst=inp%ibasismap(iatom)
         ilast=inp%ibasismap(iatom+1)-1
         isize=ilast-ifirst+1
         iifirst=1
         iilast=iifirst+isize-1
    
         !Obtain the coefficients of the lone pair in the NAO basis
         output(ispin)%coeff(ifirst:ilast,ibondorlp,1)=nho_coeff(1:isize,inho,iatom)

         !Determine the hybridization by identifying, the type of each basis function and its contribution to the hybrid
         fracl_hybrid1 = 0.d0 
         DO ibasis=ifirst,ilast
            il=inp%ilmap(ibasis)
            fracl_hybrid1(il)=fracl_hybrid1(il)+output(ispin)%coeff(ibasis,ibondorlp,1)**2
         ENDDO
         fracl_hybrid1=fracl_hybrid1/SUM(fracl_hybrid1)*100

         val_occ = val_occ + occ

         nnbo_out = nnbo_out + 1

         !Output some information about the lonepair to the screen. NOTE: coefficients are already normalized since the coefficient of the NHO is 1
         IF( occ > (core_thresh+(2-inp%nspins)) )THEN
            PRINT *,'CR ',inp%symbols(iatom),iatom
            WRITE(output(ispin)%label(ibondorlp),'(A4,A3,I2)')'CR',inp%symbols(iatom),iatom
         ELSE
            PRINT *,'LP ',inp%symbols(iatom),iatom
            WRITE(output(ispin)%label(ibondorlp),'(A4,A3,I2)')'LP',inp%symbols(iatom),iatom
         ENDIF
         WRITE(6,*)'Occ:',SNGL(occ)
         WRITE(6,'(A16,A4,F15.9,A4,F15.9,A4,F15.9,A4,F15.9)') ' Hybridization ', 's%',SNGL(fracl_hybrid1(0)),'p%',SNGL(fracl_hybrid1(1)),'d%',SNGL(fracl_hybrid1(2)),'f%',SNGL(fracl_hybrid1(3))
         WRITE(fmt_statement,'("(",i4,"f10.5)")')isize
         WRITE(*,fmt_statement) SNGL(output(ispin)%coeff(ifirst:ilast,ibondorlp,1))

         output(ispin)%occ(ibondorlp) = occ

         PRINT *

      ENDDO

      WRITE(6,'(A)')'--------------------------------------------------------'
      WRITE(6,*)

      ENDIF

      !Prepare a matrix to project out the bonding orbital density, so that Rydberg states can then be found.
      !P=0.d0
      !P=matiden(SIZE(P,1))
       

      IF( nbond > 0 )THEN

      !Loop over all potential bonds located previously
      WRITE(6,'(A)')' ******************************* '
      WRITE(6,'(A)')' **** NATURAL BOND ORBITALS **** '
      WRITE(6,'(A)')' ******************************* '
      WRITE(6,*)
      DO ibondorlp=1,nbond

         iatom=bonds(ibondorlp)%iatom
         jatom=bonds(ibondorlp)%jatom
         ig=bonds(ibondorlp)%ig
         jg=bonds(ibondorlp)%jg
         inho=bonds(ibondorlp)%inho
         jnho=bonds(ibondorlp)%jnho

         IF (jatom.EQ.-1) GOTO 50

         ifirst=inp%ibasismap(iatom)
         ilast=inp%ibasismap(iatom+1)-1
         isize=ilast-ifirst+1
         iifirst=1
         iilast=iifirst+isize-1

         jfirst=inp%ibasismap(jatom)
         jlast=inp%ibasismap(jatom+1)-1
         jsize=jlast-jfirst+1
         jjfirst=iilast+1
         jjlast=jjfirst+jsize-1

         nbo_coeff=0.d0
         nbo_occ=0.d0

         rhoat=0.d0
         rhoat(iifirst:iilast,iifirst:iilast)=inp%rho0(ifirst:ilast,ifirst:ilast,ig,ispin)
         rhoat(jjfirst:jjlast,jjfirst:jjlast)=inp%rho0(jfirst:jlast,jfirst:jlast,ig,ispin)
         rhoat(iifirst:iilast,jjfirst:jjlast)=inp%rho0(ifirst:ilast,jfirst:jlast,jg,ispin)
         rhoat(jjfirst:jjlast,iifirst:iilast)=TRANSPOSE(inp%rho0(ifirst:ilast,jfirst:jlast,jg,ispin))

         !WRITE(6,'(A,4I4)')'Diatomic block of density matrix',ibondorlp,iatom,jatom,jg
         !DO j=iifirst,jjlast
         !   WRITE(6,'(18F9.4)')rhoat(j,iifirst:jjlast)
         !ENDDO
         !WRITE(6,*)


         !Re-express this density matrix in the basis of the hybrid orbitals on the atoms
         !Tnho is a jjlast x nnbo matrix whose columns contain the NHOs for iatom and jatom in the NAO basis
         nnbo=2
         Tnho=0.d0
         Tnho(iifirst:iilast,1)=nho_coeff(1:isize,inho,iatom)
         Tnho(jjfirst:jjlast,2)=nho_coeff(1:jsize,jnho,jatom)
         !Transform rhoat into the NHO basis
         rhoat(1:nnbo,1:nnbo)=matunitary_trans(rhoat(iifirst:jjlast,iifirst:jjlast),TRANSPOSE(Tnho(iifirst:jjlast,1:nnbo)))

         !PRINT *, 'off diag', rhoat(1,2)**2

         !nnbo=isize+jsize
         CALL matdiag(rhoat(1:nnbo,1:nnbo),nbo_occ(1:nnbo),nbo_coeff(1:nnbo,1:nnbo))

         !Now transform the NBO coefficients BACK from the NHO basis to the NAO basis
         nbo_coeff(iifirst:jjlast,1:nnbo)=MATMUL(Tnho(iifirst:jjlast,1:nnbo),nbo_coeff(1:nnbo,1:nnbo))

         !WRITE(6,*)'nbo_occ to test',nbo_occ

         IF (MAXVAL(nbo_occ).GT.bond_cutoff) THEN

            DO inbo=1,nnbo
               occ=nbo_occ(inbo)
               polarization=SUM(nbo_coeff(iifirst:iilast,inbo)**2)

               fracl_hybrid1=0.d0
               fracl_hybrid2=0.d0
               DO ibasis=iifirst,iilast
                  il=inp%ilmap(ifirst+ibasis-iifirst)
                  fracl_hybrid1(il)=fracl_hybrid1(il)+nbo_coeff(ibasis,inbo)**2
               ENDDO
               fracl_hybrid1=fracl_hybrid1/SUM(fracl_hybrid1)*100
               DO jbasis=jjfirst,jjlast
                  il=inp%ilmap(jfirst+jbasis-jjfirst)
                  fracl_hybrid2(il)=fracl_hybrid2(il)+nbo_coeff(jbasis,inbo)**2
               ENDDO
               fracl_hybrid2=fracl_hybrid2/SUM(fracl_hybrid2)*100

               IF (MAX(polarization,1.d0-polarization).LE.polarization_thresh) THEN
                  IF (occ.GT.bond_cutoff)THEN
                     IF (jg.NE.1) THEN
                        PRINT *, 'BD ', inp%symbols(iatom), iatom, inp%symbols(jatom), jatom, "' (", &
                             & inp%indexg(1,jg), inp%indexg(2,jg), inp%indexg(3,jg), ")"
                     ELSE
                        PRINT *, 'BD ', inp%symbols(iatom), iatom, inp%symbols(jatom), jatom
                     ENDIF
                     val_occ = val_occ + occ
                     WRITE(output(ispin)%label(nlpair+2*ibondorlp+1-inbo),'(A4,A3,I2,A3,I2)')'BD',inp%symbols(iatom), iatom, inp%symbols(jatom), jatom
                  ELSE
                     IF (jg.NE.1) THEN
                        PRINT *, 'BD* ', inp%symbols(iatom), iatom, inp%symbols(jatom), jatom, "' (", &
                             & inp%indexg(1,jg), inp%indexg(2,jg), inp%indexg(3,jg), ")"
                     ELSE
                        PRINT *, 'BD* ', inp%symbols(iatom), iatom, inp%symbols(jatom), jatom
                     ENDIF
                     nl_occ = nl_occ + occ
                     WRITE(output(ispin)%label(nlpair+2*ibondorlp+1-inbo),'(A4,A3,I2,A3,I2)')'BD*',inp%symbols(iatom), iatom, inp%symbols(jatom), jatom
                  ENDIF
                  
                  PRINT *, 'Occ: ', SNGL(occ), 'Polarization:', SNGL(polarization*100),'%',SNGL((1.d0-polarization)*100),'%'
                  WRITE(6,'(A13,A4,F15.9,A4,F15.9,A4,F15.9,A4,F15.9)')' Hybrid 1: ', 's%',SNGL(fracl_hybrid1(0)),'p%',SNGL(fracl_hybrid1(1)),'d%',SNGL(fracl_hybrid1(2)),'f%',SNGL(fracl_hybrid1(3))
                  WRITE(fmt_statement,'("(",i4,"f10.5)")')isize
                  WRITE(*,fmt_statement) SNGL(nbo_coeff(iifirst:iilast,inbo)/SQRT(SUM(nbo_coeff(iifirst:iilast,inbo)**2)))
                  output(ispin)%coeff(ifirst:ilast,nlpair+2*ibondorlp+1-inbo,ig) = nbo_coeff(iifirst:iilast,inbo)
 
                  WRITE(6,'(A13,A4,F15.9,A4,F15.9,A4,F15.9,A4,F15.9)') ' Hybrid 2: ', 's%',SNGL(fracl_hybrid2(0)),'p%',SNGL(fracl_hybrid2(1)),'d%',SNGL(fracl_hybrid2(2)),'f%',SNGL(fracl_hybrid2(3))
                  WRITE(fmt_statement,'("(",i4,"f10.5)")')jsize
                  WRITE(*,fmt_statement) SNGL(nbo_coeff(jjfirst:jjlast,inbo)/SQRT(SUM(nbo_coeff(jjfirst:jjlast,inbo)**2)))
                  output(ispin)%coeff(jfirst:jlast,nlpair+2*ibondorlp+1-inbo,jg) = nbo_coeff(jjfirst:jjlast,inbo)

                  output(ispin)%occ(nlpair+2*ibondorlp+1-inbo) = occ

                  nnbo_out = nnbo_out + 1

                  PRINT *

                  !A projection operator can be constructed as the product of projections of each orbital
                  !Pi corresponds to that for the single orbital, and P is for the collection of orbitals
                  !Even though the bonds can span unit cells, the projection assumes all hybirds are  within the central unit cell
                  !This is the only matrix that will be diagonalized for Rydberg states, and is the only one that must be depleted
                  Pi=0.d0
                  Pi(iifirst:jjlast,iifirst:jjlast)=matiden(jjlast)
                  !proj_vec(iifirst:iilast) = nbo_coeff(iifirst:iilast,inbo)
                  !Prettyse sure I am just using depletion now
                  !For orthogonal hybrids, its the same as projection
                  Pi(iifirst:jjlast,iifirst:jjlast)=Pi(iifirst:jjlast,iifirst:jjlast) - outer_product(nbo_coeff(iifirst:jjlast,inbo),nbo_coeff(iifirst:jjlast,inbo))
                  !IF( jg /= ig )Pi = Pi - outer_product(bond_coeff(:,2*ibondorlp+1-inbo,jg),bond_coeff(:,2*ibondorlp+1-inbo,jg))
                  !proj_dummy = P
                  !P = MATMUL(Pi,P)
                  !CALL DGEMM_F95(Pi,proj_dummy,P,'N','N',1.d0,0.d0)

                  rhoat=0.d0
                  rhoat(iifirst:iilast,iifirst:iilast)=inp%rho0(ifirst:ilast,ifirst:ilast,ig,ispin)
                  rhoat(jjfirst:jjlast,jjfirst:jjlast)=inp%rho0(jfirst:jlast,jfirst:jlast,ig,ispin)
                  rhoat(iifirst:iilast,jjfirst:jjlast)=inp%rho0(ifirst:ilast,jfirst:jlast,jg,ispin)
                  rhoat(jjfirst:jjlast,iifirst:iilast)=TRANSPOSE(inp%rho0(ifirst:ilast,jfirst:jlast,jg,ispin))

                  rhoat(iifirst:jjlast,iifirst:jjlast)=matunitary_trans(rhoat(iifirst:jjlast,iifirst:jjlast),Pi(iifirst:jjlast,iifirst:jjlast))
 
                  !WRITE(6,'(A,4I5)')'Diatomic block of density matrix after projection',ibondorlp,iatom,jatom,jg
                  !DO j=iifirst,jjlast
                  !   WRITE(6,'(18F9.4)')rhoat(j,iifirst:jjlast)
                  !ENDDO
                  !WRITE(6,*)
                  !WRITE(6,*)

                  inp%rho0(ifirst:ilast,ifirst:ilast,ig,ispin)=rhoat(iifirst:iilast,iifirst:iilast)
                  inp%rho0(jfirst:jlast,jfirst:jlast,ig,ispin)=rhoat(jjfirst:jjlast,jjfirst:jjlast)
                  inp%rho0(ifirst:ilast,jfirst:jlast,jg,ispin)=rhoat(iifirst:iilast,jjfirst:jjlast)

                  kg = find_g(-inp%indexg(1,jg),-inp%indexg(2,jg),-inp%indexg(3,jg))
                  IF( kg == -1 )THEN
                      WRITE(6,*)'Error in looking for kg index in removing nbo density from rho0',ibondorlp,jg
                      STOP
                  ENDIF
                  !WRITE(6,*)'kg index for density removal',kg,jg
                  !WRITE(6,*)

                  inp%rho0(jfirst:jlast,ifirst:ilast,kg,ispin)=rhoat(jjfirst:jjlast,iifirst:iilast)

               ELSE
                  WRITE(6,*)'nbo excluded beacuse polarization is below threshold'
               ENDIF
            ENDDO
         ENDIF

50    ENDDO

      ENDIF

      !Now project out all of the bonding and antibonding density, to leave solely the dentisy associated with Rydberg states
      !inp%rho0(:,:,1) = MATMUL(P,(MATMUL(inp%rho0(:,:,1),P)))
      !CALL DGEMM_F95(inp%rho0(:,:,1),P,proj_dummy,'N','N',1.d0,0.d0)
      !CALL DGEMM_F95(P,proj_dummy,inp%rho0(:,:,1),'N','N',1.d0,0.d0)


      !PRINT *, 'made it to doing rydberg search'
      CALL do_ryd_search(inp,rydberg,nnho,nryd,rhoat,nbo_coeff,nbo_occ,nho_coeff,nnho_init,ispin)
      !PRINT *,'made it past rydberg seach to orthogonalization test'
      IF (orthogonalize_hybrids) CALL do_orthogonalize_hybrids(inp,nho_coeff,nnho,.FALSE.,nnho_init)

      !ALLOCATE(ryd_coeff(inp%nbasis,nryd,inp%ng),ryd_occ(nryd))
      !ALLOCATE(ryd_lbl(nryd))
      !ryd_coeff = 0.d0
      ryd_occ = 0.d0
      !ryd_lbl = 'RYD'

      WRITE(6,'(A)')' ******************************** '
      WRITE(6,'(A)')' **** NATURAL RYDBERG STATES **** '
      WRITE(6,'(A)')' ******************************** '
      WRITE(6,*)

      DO ibondorlp=1,nryd
         iatom=rydberg(ibondorlp)%iatom
         ig=rydberg(ibondorlp)%ig
         inho=rydberg(ibondorlp)%inho
         occ=rydberg(ibondorlp)%occ
         ifirst=inp%ibasismap(iatom)
         ilast=inp%ibasismap(iatom+1)-1
         isize=ilast-ifirst+1
         iifirst=1
         iilast=iifirst+isize-1

         !Obtain the coefficients of the rydberg state in the NAO basis
         output(ispin)%coeff(ifirst:ilast,nlpair+2*nbond+ibondorlp,1)=nho_coeff(1:isize,inho,iatom)

         !Determine the hybridization by identifying, the type of each basis function and its contribution to the hybrid
         fracl_hybrid1 = 0.d0
         DO ibasis=ifirst,ilast
            il=inp%ilmap(ibasis)
            fracl_hybrid1(il)=fracl_hybrid1(il)+output(ispin)%coeff(ibasis,nlpair+2*nbond+ibondorlp,1)**2
         ENDDO
         fracl_hybrid1=fracl_hybrid1/SUM(fracl_hybrid1)*100

         output(ispin)%occ(nlpair+2*nbond+ibondorlp) = occ
         WRITE(output(ispin)%label(nlpair+2*nbond+ibondorlp),'(A4,A3,I2)')'RYD',inp%symbols(iatom),iatom

         ryd_occ = ryd_occ + occ

         nnbo_out = nnbo_out + 1

         !Output some information about the rydberg orbital to the screen. NOTE: coefficients are already normalized since the coefficient of the NBO is 1
         !This is only output for rydbergs with some occupancy
         IF( occ .GT. 1.d-3 )THEN
            PRINT *,'RYD ',inp%symbols(iatom),iatom
            WRITE(6,*)'Occ:',SNGL(occ)
            WRITE(6,'(A16,A4,F15.9,A4,F15.9,A4,F15.9,A4,F15.9)') ' Hybridization ', 's%',SNGL(fracl_hybrid1(0)),'p%',SNGL(fracl_hybrid1(1)),'d%',SNGL(fracl_hybrid1(2)),'f%',SNGL(fracl_hybrid1(3))
            WRITE(fmt_statement,'("(",i4,"f10.5)")')isize
            WRITE(*,fmt_statement) SNGL(output(ispin)%coeff(ifirst:ilast,nlpair+2*nbond+ibondorlp,1))
            PRINT *
         ENDIF

      ENDDO


      WRITE(6,*)'Total occupancy of Lewis type valence NBOs     ',SNGL(val_occ)
      WRITE(6,*)'Total occupancy of Non-Lewis type valence NBOs ',SNGL(nl_occ)
      WRITE(6,*)'Total occupancy of Rydberg NBOs                ',SNGL(ryd_occ)
      WRITE(6,*)'Total valence occupancy                        ',SNGL(val_occ+nl_occ+ryd_occ)
      WRITE(6,*)

      !STOP

    END SUBROUTINE do_nbo


    !This subroutine does a 1c search for lone pairs, and then depletes the density matrix of any that are
    !located.
    SUBROUTINE do_1c_search(inp,lpairs,nnho,nbondorlp,rhoat,nbo_coeff,nbo_occ,nho_coeff,ispin)
      USE nbo_shared
      USE matutil
      IMPLICIT NONE

      TYPE(nbo_input) :: inp
      TYPE(nbo_bond),DIMENSION(:) :: lpairs

      REAL*8, DIMENSION(:,:) :: rhoat,nbo_coeff
      REAL*8, DIMENSION(:,:,:) :: nho_coeff
      REAL*8, DIMENSION(:) :: nbo_occ
      INTEGER, DIMENSION(:) :: nnho
      INTEGER :: nbondorlp

      INTEGER :: ig,iatom,ifirst,ilast,isize
      INTEGER :: iifirst,iilast
      INTEGER :: inbo,nnbo,inho
      REAL*8 :: lpcontrib(inp%nbasis), occ
      INTEGER :: ibondorlp

      INTEGER ::   nnho_init

      INTEGER ::  ispin  !Which density matrix is currently being examined for NBO's

      REAL*8,DIMENSION(SIZE(inp%rho0,1),SIZE(inp%rho0,1))  :: P,Pi,Pdummy

      LOGICAL  ::  rydberg
      LOGICAL,ALLOCATABLE  ::  non_proj(:)

      !First do a search for "one-center" bonds (lone pairs)
      ig=1
      DO iatom=1,inp%natom
         !WRITE(6,*)'loop over atomic center',iatom
         ifirst=inp%ibasismap(iatom)
         ilast=inp%ibasismap(iatom+1)-1
         isize=ilast-ifirst+1
         iifirst=1
         iilast=iifirst+isize-1

         nnbo=isize
         nbo_coeff=0.d0
         nbo_occ=0.d0

         rhoat(iifirst:iilast,iifirst:iilast)=inp%rho0(ifirst:ilast,ifirst:ilast,ig,ispin)
         CALL matdiag(rhoat(iifirst:iilast,iifirst:iilast),nbo_occ(iifirst:iilast),nbo_coeff(iifirst:iilast,iifirst:iilast))

         !WRITE(6,*)'ready to store information on a hybrid'
         !Store all lone-pair hybrid orbitals
         DO inbo=1,nnbo
            !WRITE(6,*)'testing nbo',inbo,'of atom',iatom
            occ=nbo_occ(inbo)
            !WRITE(6,*)'Occ to test',occ
            IF (occ.GT.lp_cutoff) THEN
               !WRITE(6,*)'NBO found'
               !store this bond for future reference
               nnho(iatom)=nnho(iatom)+1
               nho_coeff(1:isize,nnho(iatom),iatom)=nbo_coeff(iifirst:iilast,inbo)

               nbondorlp=nbondorlp+1
               lpairs(nbondorlp)%iatom=iatom
               lpairs(nbondorlp)%jatom=-1
               lpairs(nbondorlp)%ig=ig
               lpairs(nbondorlp)%jg=-1
               lpairs(nbondorlp)%inho=nnho(iatom)
               lpairs(nbondorlp)%jnho=-1
               lpairs(nbondorlp)%occ=occ
            ENDIF
         ENDDO


      ENDDO


      !WRITE(6,*)'made it to deleting the lone pair density'
      !Next project the "one-center bond" contributions from the density matrix
      P=matiden(SIZE(inp%rho0,1))
      DO ibondorlp=1,nbondorlp
         iatom=lpairs(ibondorlp)%iatom
         ig=lpairs(ibondorlp)%ig
         inho=lpairs(ibondorlp)%inho
         occ=lpairs(ibondorlp)%occ
         ifirst=inp%ibasismap(iatom)
         ilast=inp%ibasismap(iatom+1)-1
         isize=ilast-ifirst+1
         iifirst=1
         iilast=iifirst+isize-1

         !Construct the contribution of the lone-pair denisty in lpcontrib
         lpcontrib=0.d0
         lpcontrib(ifirst:ilast)=nho_coeff(1:isize,inho,iatom)
         
         !WRITE(6,*)lpcontrib
         !Multiply the projection operator for this obrital to the previously constructed product of projection operators
         Pi=matiden(SIZE(inp%rho0,1))
         Pi=Pi-outer_product(lpcontrib,lpcontrib)
         Pdummy=P
         !P=MATMUL(Pi,P)
         CALL DGEMM_F95(Pi,Pdummy,P,'N','N',1.d0,0.d0)

         !IF( occ > core_thresh*DBLE(3-inp%nspins) )THEN
         !   PRINT *, 'Removing CR from atom',iatom
         !ELSE
         !   PRINT *, 'Removing LP from atom',iatom
         !ENDIF
      ENDDO
      !PRINT *

      !inp%rho0(:,:,1) = MATMUL(P,(MATMUL(inp%rho0(:,:,1),P)))
      CALL DGEMM_F95(inp%rho0(:,:,1,ispin),P,Pdummy,'N','N',1.d0,0.d0)
      CALL DGEMM_F95(P,Pdummy,inp%rho0(:,:,1,ispin),'N','N',1.d0,0.d0)


    END SUBROUTINE do_1c_search

    !This subroutine determines the Rydbegr states on each atomic center.  
    !These functions serve to finish out the basis set of NBO's so that is spans the same space as the original AO set.
    !They are all expected to have a low occupancy but can also serve an important role in delocaliztion effects of the Lewis-type NBO's
    !This subroutine is very similar to the lone pair search, except:
    ! -A different diagonalization routine is used that works to minimize the mixing of degenerate eigenvectors(hybrids). 
    !    This maintains a pure s,p, or d characater of the Rydbergs hybrids
    !    This will not create a set of perfectly orthogonal hybrids so some orthogonalization will need to follow this subroutine
    ! -Since the diagonalization will return the same number of hybrids as NAO's, some of these will overlap with the previously found bonding hybrids.
    !    Therefore the unique Rydberg hybrids must be seperated out.
    ! -There is no projection.  These are the last NBO's possible and represent the completion of the NBO basis.
    !    To remove them would leave a density matrix composed solely of numerical noise.
    SUBROUTINE do_ryd_search(inp,lpairs,nnho,nbondorlp,rhoat,nbo_coeff,nbo_occ,nho_coeff,nnho_init,ispin)
      USE nbo_shared
      USE matutil
      IMPLICIT NONE

      TYPE(nbo_input) :: inp
      TYPE(nbo_bond),DIMENSION(:) :: lpairs

      REAL*8, DIMENSION(:,:) :: rhoat,nbo_coeff
      REAL*8, DIMENSION(:,:,:) :: nho_coeff
      REAL*8, DIMENSION(:) :: nbo_occ
      INTEGER, DIMENSION(:) :: nnho
      INTEGER :: nbondorlp

      INTEGER :: ig,iatom,ifirst,ilast,isize
      INTEGER :: iifirst,iilast
      INTEGER :: inbo,nnbo,inho
      REAL*8 :: lpcontrib(inp%nbasis), occ
      INTEGER :: ibondorlp

      INTEGER,DIMENSION(:) ::   nnho_init

      INTEGER  ::  ispin   !Which density matrix is currently being analyzed for NBO's

      LOGICAL,DIMENSION(SIZE(nbo_occ,1))  ::  non_ryd
      INTEGER,DIMENSION(1) ::  non_loc

      !INTEGER,DIMENSION(SIZE(inp%rho0,1))  :: il_map

      !INTEGER  :: il, il_count, tot_count
      !INTEGER,DIMENSION(0:lmax) :: lfirst,llast,llfirst,lllast
      INTEGER  :: i,j,k !,lj,lk

      INTEGER  ::  nrot

      !First do a search for "one-center" bonds (lone pairs)
      ig=1
      DO iatom=1,inp%natom
         !WRITE(6,*)'loop over atomic center',iatom
         ifirst=inp%ibasismap(iatom)
         ilast=inp%ibasismap(iatom+1)-1
         isize=ilast-ifirst+1
         iifirst=1
         iilast=iifirst+isize-1

         nnbo=isize
         nbo_coeff=0.d0
         nbo_occ=0.d0

         nnho_init(iatom) = nnho(iatom)
         !WRITE(6,*)'nnho_init',nnho_init(iatom)

         !WRITE(6,*)'atomic block density matrix'
         !DO j=ifirst,ilast
         !   WRITE(6,'(18D9.2)')inp%rho0(j,ifirst:ilast,1)
         !ENDDO

         !Perform diagonalization using Jacobian rotations
         rhoat(iifirst:iilast,iifirst:iilast)=inp%rho0(ifirst:ilast,ifirst:ilast,ig,ispin)
         CALL matdiag(rhoat(iifirst:iilast,iifirst:iilast),nbo_occ(iifirst:iilast),nbo_coeff(iifirst:iilast,iifirst:iilast))

         !Now that the hybrids have been found, we need to screen out those that correspond to previously found (and projected out) hybrids
         CALL find_ryd(nnho_init(iatom),nbo_occ,non_ryd,isize)

         !non_ryd = .TRUE.
         !Store all rydberg hybrid orbitals
         DO inbo=1,nnbo
            !WRITE(6,*)'testing nbo',inbo,'of atom',iatom
            occ=nbo_occ(inbo)
            !WRITE(6,*)'Occ to test',occ
            IF (non_ryd(inbo)) THEN
               !WRITE(6,*)'NBO found'
               !store this bond for future reference
               nnho(iatom)=nnho(iatom)+1
               nho_coeff(1:isize,nnho(iatom),iatom)=nbo_coeff(iifirst:iilast,inbo)

               nbondorlp=nbondorlp+1
               lpairs(nbondorlp)%iatom=iatom
               lpairs(nbondorlp)%jatom=-1
               lpairs(nbondorlp)%ig=ig
               lpairs(nbondorlp)%jg=-1
               lpairs(nbondorlp)%inho=nnho(iatom)
               lpairs(nbondorlp)%jnho=-1
               lpairs(nbondorlp)%occ=occ
            ENDIF
         ENDDO

         !WRITE(6,*)

      ENDDO


    END SUBROUTINE do_ryd_search



    !This subroutine identifies which of the found Rydberg hybrids is actually redundant with the previously projected hybrids.
    !Since that density was projected out their eigenvalues should have a very small magnitude.
    !Thus, we assume the eigenvectors with the smallest absolute eigenvalues are redundant.  This is not a fool proof plan, but has worked so far.
    !NOTE: Even though the density matrix is still 'hermitian', 
    !      the values are so low that numerical noise starts to dominate and small negative eigenvalues occur.
    !A more rigorous test would be to calculate the overlap of all found hybrids with the previous hybrids and remove those with highest overlap.
    SUBROUTINE find_ryd(nnho_proj,occ,ryd_bool,nbasis)
      IMPLICIT NONE

      INTEGER,INTENT(IN)               ::  nnho_proj
      REAL*8,DIMENSION(:),INTENT(IN)   ::  occ
      LOGICAL,DIMENSION(:),INTENT(OUT) ::  ryd_bool
      INTEGER,INTENT(IN)               ::  nbasis

      INTEGER   ::  norb
      INTEGER,DIMENSION(1)  ::  proj_loc
      INTEGER   ::  i,j,k

      ryd_bool = .TRUE.

      !First just make sure the back end of the ryd_bool vector contains falses. 
      !The vector is bigger than the number of hybrids on some atoms, so just be sure and make sure these blank entries are not selected.
      !Otherwise the occupancy is exactly zero and these vectors will be selected.
      DO j=1,SIZE(ryd_bool,1)
         IF( j .GT. nbasis )ryd_bool(j)=.FALSE.
      ENDDO

      DO j=1,nnho_proj
         proj_loc=MINLOC(ABS(occ),ryd_bool)
         !WRITE(6,*)'location of lowest occupancy eigenvalue',proj_loc
         ryd_bool(proj_loc(1)) = .FALSE.
      ENDDO


    END SUBROUTINE find_ryd


    !This subroutine does a 2c search for BONDS.
    SUBROUTINE do_2c_search(inp,bonds,nnho,nbondorlp,rhoat,nbo_coeff,nbo_occ,nho_coeff,ispin)
      USE nbo_shared
      USE matutil
      IMPLICIT NONE

      TYPE(nbo_input) :: inp
      TYPE(nbo_bond),DIMENSION(:) :: bonds

      REAL*8, DIMENSION(:,:) :: rhoat,nbo_coeff
      REAL*8, DIMENSION(:,:,:) :: nho_coeff
      REAL*8, DIMENSION(:) :: nbo_occ
      INTEGER, DIMENSION(:) :: nnho
      INTEGER :: nbondorlp

      INTEGER :: ig,jg,iatom,jatom, j
      INTEGER :: ifirst,ilast,isize,iifirst,iilast
      INTEGER :: jfirst,jlast,jsize,jjfirst,jjlast
      INTEGER :: inbo,nnbo,inho,jnho
      REAL*8 :: occ

      INTEGER  :: ispin   ! Which density matrix is being analyzed for NBOs

      REAL*8 :: polarization,projection
      LOGICAL :: new_bond_found

      LOGICAL :: debug_search
      PARAMETER (debug_search=.FALSE.)

      !Next do a search for "two-center" bonds (normal covalent bonds)
      !First do an NBO search within the central unit cell, and between the central cell and the adjacent ones
      !in order to locate potential hybrid orbitals
      ig=1
      DO iatom=1,inp%natom
         ifirst=inp%ibasismap(iatom)
         ilast=inp%ibasismap(iatom+1)-1
         isize=ilast-ifirst+1
         iifirst=1
         iilast=iifirst+isize-1

         DO jg=1,inp%ng
            !Avoid looking past cells which do not neighbor the central cell, since the FFT introduces numerical issues in distant cells
            !due to truncation, etc.
            IF (MAXVAL(ABS(inp%indexg(:,jg))).GT.1) GOTO 30

            DO jatom=1,inp%natom
               jfirst=inp%ibasismap(jatom)
               jlast=inp%ibasismap(jatom+1)-1
               jsize=jlast-jfirst+1
               jjfirst=iilast+1
               jjlast=jjfirst+jsize-1

               !Avoid double counting bonds in the unit cell
               IF (jg.EQ.1.AND.iatom.GE.jatom) GOTO 40 !why does FORTRAN not have a real continue statement!

               nnbo=isize+jsize
               nbo_coeff=0.d0
               nbo_occ=0.d0

               rhoat=0.d0
               rhoat(iifirst:iilast,iifirst:iilast)=inp%rho0(ifirst:ilast,ifirst:ilast,ig,ispin)
               rhoat(jjfirst:jjlast,jjfirst:jjlast)=inp%rho0(jfirst:jlast,jfirst:jlast,ig,ispin)
               rhoat(iifirst:iilast,jjfirst:jjlast)=inp%rho0(ifirst:ilast,jfirst:jlast,jg,ispin)
               rhoat(jjfirst:jjlast,iifirst:iilast)=TRANSPOSE(inp%rho0(ifirst:ilast,jfirst:jlast,jg,ispin))

               !WRITE(6,*)'Diatomic block of density matrix',iatom,jatom,jg
               !DO j=iifirst,jjlast
               !   WRITE(6,'(18F9.4)')rhoat(j,iifirst:jjlast)
               !ENDDO

               CALL matdiag(rhoat(iifirst:jjlast,iifirst:jjlast),nbo_occ(iifirst:jjlast),nbo_coeff(iifirst:jjlast,iifirst:jjlast))
               !Remove trivial degeneracies due to essentially zero overlap
               !CALL remove_degenerate_pairs(nbo_occ(iifirst:jjlast),nbo_coeff(iifirst:jjlast,iifirst:jjlast),iifirst,iilast,jjfirst,jjlast)

               IF (debug_search) THEN
                  PRINT *, 'Searching: ', iatom, jatom, jg, MAXVAL(nbo_occ(iifirst:jjlast))
               ENDIF

               !Store all potential hybrid orbitals
               DO inbo=1, nnbo
                  occ=nbo_occ(inbo)
                  new_bond_found=.TRUE.
                  polarization=SUM(nbo_coeff(iifirst:iilast,inbo)**2)

                  IF (debug_search) THEN
                     PRINT *, inbo, occ, polarization
                  ENDIF

                  IF (occ.GT.bond_cutoff.AND.MAX(polarization,1.d0-polarization).LE.polarization_thresh) THEN
                     !First normalize the hybrids
                     nbo_coeff(iifirst:iilast,inbo)=nbo_coeff(iifirst:iilast,inbo)/SQRT(polarization)
                     nbo_coeff(jjfirst:jjlast,inbo)=nbo_coeff(jjfirst:jjlast,inbo)/SQRT(1.d0-polarization)

                     !Check whether the hybrid is unique; in case where we are not orthogonalize the hybrids, keep them ALL
                     !see PCCP 10, 5207 (2008) for details and some discussion of this point.
                     projection=prjexp(nho_coeff(1:isize,1:nnho(iatom),iatom),nnho(iatom),nbo_coeff(iifirst:iilast,inbo))
                     IF (projection.GT.prjexp_thresh.OR..NOT.orthogonalize_hybrids) THEN
                        !Store the hybrid
                        nnho(iatom)=nnho(iatom)+1
                        IF (debug_search) PRINT *, 'Adding hybrid', iatom, jatom, jg, occ, polarization
                        IF ((nnho(iatom).GT.isize.AND.orthogonalize_hybrids).OR.nnho(iatom).GT.maxbondorlpperatom) STOP 'Too many hybrids on atom'
                        nho_coeff(1:isize,nnho(iatom),iatom)=nbo_coeff(iifirst:iilast,inbo)
                        inho=nnho(iatom)
                     ELSE
                        new_bond_found=.FALSE.
                        IF (debug_search) PRINT *, 'Omitting hybrid', iatom, jatom, jg, projection
                     ENDIF
                     
                     !Check whether the hybrid is unique
                     projection=prjexp(nho_coeff(1:jsize,1:nnho(jatom),jatom),nnho(jatom),nbo_coeff(jjfirst:jjlast,inbo))
                     IF (projection.GT.prjexp_thresh.OR..NOT.orthogonalize_hybrids) THEN
                        !Store the hybrid
                        nnho(jatom)=nnho(jatom)+1
                        IF (debug_search) PRINT *, 'Adding hybrid', jatom, iatom, jg, occ, polarization
                        IF ((nnho(jatom).GT.jsize.AND.orthogonalize_hybrids).OR.nnho(jatom).GT.maxbondorlpperatom) STOP 'Too many hybrids on atom'
                        nho_coeff(1:jsize,nnho(jatom),jatom)=nbo_coeff(jjfirst:jjlast,inbo)
                        jnho=nnho(jatom)
                     ELSE
                        new_bond_found=.FALSE.
                        IF (debug_search) PRINT *, 'Omitting hybrid', jatom, iatom, jg, projection
                     ENDIF
                     
                     !store this bond for future reference
                     IF (new_bond_found) THEN
                        nbondorlp=nbondorlp+1
                        IF (nbondorlp.GT.SIZE(bonds)) STOP 'Exceeded max # of bonds'

                        bonds(nbondorlp)%iatom=iatom
                        bonds(nbondorlp)%jatom=jatom
                        bonds(nbondorlp)%ig=ig
                        bonds(nbondorlp)%jg=jg
                        bonds(nbondorlp)%inho=inho
                        bonds(nbondorlp)%jnho=jnho
                        bonds(nbondorlp)%occ=occ
                     ENDIF
                     
                  ENDIF
               ENDDO
40          ENDDO
30       ENDDO
      ENDDO

    END SUBROUTINE do_2c_search




    !This subroutine removes TRUE degenerate pairs of eigenvalues / eigenvectors from the NBOs.
    !These originate because two distant symmetrically equivalent atoms without any coupling
    !may still mix (due the the eigensolver), thus giving rise to a "bond".  The signature of this
    !will be a degenerate pair of eigenvalues, with mirrored polarization.
    SUBROUTINE remove_degenerate_pairs(nbo_occ,nbo_coeff,iifirst,iilast,jjfirst,jjlast)
      USE nbo_shared
      IMPLICIT NONE
      REAL*8,  DIMENSION(:,:) :: nbo_coeff
      REAL*8,  DIMENSION(:) :: nbo_occ
      INTEGER :: iifirst,iilast,jjfirst,jjlast

      INTEGER :: inbo,jnbo,nnbo
      REAL*8 :: poli,polj,delta,overlap1,overlap2

      REAL*8 :: degen_thresh, overlap_thresh
      PARAMETER (degen_thresh=0.4d0,overlap_thresh=0.2d0)

      nnbo=SIZE(nbo_coeff,1)
      DO inbo=1,nnbo
         DO jnbo=inbo+1,nnbo
            poli=SUM(nbo_coeff(iifirst:iilast,inbo)**2)
            polj=SUM(nbo_coeff(iifirst:iilast,jnbo)**2)
            overlap1=ABS(DOT_PRODUCT(nbo_coeff(iifirst:iilast,inbo),nbo_coeff(iifirst:iilast,jnbo)))/SQRT(poli*polj)
            overlap2=ABS(DOT_PRODUCT(nbo_coeff(jjfirst:jjlast,inbo),nbo_coeff(jjfirst:jjlast,jnbo)))/SQRT((1.d0-poli)*(1.d0-polj))
            delta=ABS(nbo_occ(inbo)-nbo_occ(jnbo))

            !If we find two degenerate NBOs, then a true trivial degeneracy should be due to corresponding
            !degenerate values in the density matrix
            IF (delta.LT.degen_thresh.AND.overlap1.GT.overlap_thresh.AND.overlap2.GT.overlap_thresh) THEN
               nbo_occ(inbo)=0.d0
               nbo_occ(jnbo)=0.d0
            ENDIF

         ENDDO

      ENDDO

    END SUBROUTINE remove_degenerate_pairs


    !The hybrids on a given center are typically not exactly orthoganal, which means that we are double-using
    !some of the NAOs in constructing these hybrids.  Therefore do a symmetric orthogonalization on each set of
    !hybrids
    SUBROUTINE do_orthogonalize_hybrids(inp,nho_coeff,nnho,sym,nnho_init)
      USE nbo_shared
      USE matutil
      IMPLICIT NONE
      
      TYPE(nbo_input) :: inp
      REAL*8,DIMENSION(:,:,:) :: nho_coeff
      !"s" is the overlap matrix of the nho's
      REAL*8,DIMENSION(SIZE(nho_coeff,1),SIZE(nho_coeff,2)) :: s, sisqrt, Os
      INTEGER :: nnho(:),nnho_init(:)

      LOGICAL :: sym

      INTEGER :: iatom
      INTEGER :: ifirst,ilast,isize

      INTEGER :: j

      DO iatom = 1, inp%natom
         ifirst=inp%ibasismap(iatom)
         ilast=inp%ibasismap(iatom+1)-1
         isize=ilast-ifirst+1

 
         s(1:nnho(iatom),1:nnho(iatom))=MATMUL(TRANSPOSE(nho_coeff(1:isize,1:nnho(iatom),iatom)),nho_coeff(1:isize,1:nnho(iatom),iatom))
         !WRITE(6,*)'s of nhos of atom',iatom
         !DO j=1,nnho(iatom)
         !   WRITE(6,'(30F5.2)')s(j,1:nnho(iatom))
         !ENDDO
         !WRITE(6,*)

         IF( sym )THEN

            sisqrt(1:nnho(iatom),1:nnho(iatom))=matinvsqrt(s(1:nnho(iatom),1:nnho(iatom)))
            nho_coeff(1:isize,1:nnho(iatom),iatom)=MATMUL(nho_coeff(1:isize,1:nnho(iatom),iatom),sisqrt(1:nnho(iatom),1:nnho(iatom)))

         ELSE

            IF( nnho_init(iatom) /= 0 )THEN

            Os(1:nnho(iatom),1:nnho(iatom))=schmidt_orthog_matrix(s(1:nnho(iatom),1:nnho(iatom)),nnho_init(iatom))
            nho_coeff(1:isize,1:nnho(iatom),iatom)=MATMUL(nho_coeff(1:isize,1:nnho(iatom),iatom),Os(1:nnho(iatom),1:nnho(iatom)))
            s(nnho_init(iatom)+1:nnho(iatom),nnho_init(iatom)+1:nnho(iatom))=MATMUL(TRANSPOSE(nho_coeff(1:isize,nnho_init(iatom)+1:nnho(iatom),iatom)),nho_coeff(1:isize,nnho_init(iatom)+1:nnho(iatom),iatom))
            sisqrt(nnho_init(iatom)+1:nnho(iatom),nnho_init(iatom)+1:nnho(iatom))=matinvsqrt(s(nnho_init(iatom)+1:nnho(iatom),nnho_init(iatom)+1:nnho(iatom)))
            nho_coeff(1:isize,nnho_init(iatom)+1:nnho(iatom),iatom)=MATMUL(nho_coeff(1:isize,nnho_init(iatom)+1:nnho(iatom),iatom),sisqrt(nnho_init(iatom)+1:nnho(iatom),nnho_init(iatom)+1:nnho(iatom)))

            ENDIF

         ENDIF

         s(1:nnho(iatom),1:nnho(iatom))=MATMUL(TRANSPOSE(nho_coeff(1:isize,1:nnho(iatom),iatom)),nho_coeff(1:isize,1:nnho(iatom),iatom))
         !WRITE(6,*)'s of orthogonalized nhos of atom',iatom
         !DO j=1,nnho(iatom)
         !   WRITE(6,'(30F5.2)')s(j,1:nnho(iatom))
         !ENDDO
         !WRITE(6,*)

      ENDDO

    END SUBROUTINE do_orthogonalize_hybrids

    FUNCTION schmidt_orthog_matrix(S,n) RESULT(Os)
      USE matutil
      USE BLAS95
      IMPLICIT NONE
      REAL*8,DIMENSION(:,:) :: S
      REAL*8,DIMENSION(n,n) :: Sinv
      REAL*8,DIMENSION(n,SIZE(S,2)-n) :: proj
      REAL*8,DIMENSION(SIZE(S,1),SIZE(S,2)) :: Os
      INTEGER :: n,i,j,nbasis
     
      nbasis=SIZE(S,1)
      
      !This is tricky; since we are in a non-orthogonal basis set and doing only a PARTIAL gs
      !orthogonalization, we need to find the orthogonal projection of each basis function 1:n
      !onto the remaining basis function n+1:nbasis
      Sinv=matinv(S(1:n,1:n))
      !proj=MATMUL(Sinv,S(1:n,n+1:nbasis))
      CALL DGEMM_F95(Sinv,S(1:n,n+1:nbasis),proj,'N','N',1.d0,0.d0)

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
      !Os=TRANSPOSE(Os)
    END FUNCTION schmidt_orthog_matrix




    !This subroutine calculates the expectation value of 1-P with a proposed hybrid, where P is the projector into the span on the
    !previously located hybrids.  In othe words, this tells us how much of the proposed hybrid is unique.
    FUNCTION prjexp(nho_coeff,nnho,proposed_nho)
      USE matutil
      IMPLICIT NONE
      REAL*8 :: prjexp

      REAL*8,DIMENSION(:,:) :: nho_coeff
      INTEGER :: nnho
      REAL*8 :: proposed_nho(:)

      INTEGER :: inho
      REAL*8,DIMENSION(SIZE(nho_coeff,1),SIZE(nho_coeff,1)) :: P,Pi

      P=matiden(SIZE(P,1))
      DO inho=1,nnho
         Pi=matiden(SIZE(Pi,1))
         Pi=Pi-outer_product(nho_coeff(:,inho),nho_coeff(:,inho))
         P=MATMUL(Pi,P)
      ENDDO
      
      prjexp=DOT_PRODUCT(proposed_nho,MATMUL(P,proposed_nho))
    END FUNCTION prjexp


END MODULE nbo
