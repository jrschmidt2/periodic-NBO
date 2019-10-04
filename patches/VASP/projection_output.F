#include "symbol.inc"


      MODULE projection_output

      CONTAINS



      !This subroutine was written by JRS and BDD to print information from VASP to a standardized format
      !This is intended to be used with the projection algortihm described in: JCTC 2012, 8, 1902-1911.
      !The information is output to an unformatted file called wavefunction.dat
      !As written, this is not compatible with VASP runs that do not use the default value for NPAR

      SUBROUTINE output_NBO(WDES,W,LATT_INI,T_INFO,P,CQIJ)

      USE prec
      USE wave
      USE lattice
      USE poscar
      USE pseudo
      IMPLICIT REAL(q) (A-H,O-Z)

      TYPE (latt)     LATT_INI
      TYPE (wavedes)  WDES
      TYPE (wavedes1) WDES1
      TYPE (wavespin) W
      TYPE (type_info) T_INFO
      TYPE (potcar)   P(WDES%NTYP)
      OVERLAP :: CQIJ(:,:,:,:)
! local work arrays
      COMPLEX(q),ALLOCATABLE :: CW(:),EIG(:)
#ifdef WAVECAR_double
      COMPLEX(q), ALLOCATABLE :: CRD(:)
#else
      COMPLEX(qs),  ALLOCATABLE :: CRD(:)
#endif
      GDEF,ALLOCATABLE :: CPROJ(:),CPROJ_LOC(:)

      NPL_TOT = MAXVAL(WDES%NPLWKP_TOT)
      NPRO_TOT=WDES%NPRO_TOT
      ALLOCATE(CW(NPL_TOT),CRD(NPL_TOT),EIG(WDES%NB_TOT),CPROJ(WDES%NPRO_TOT),CPROJ_LOC(WDES%NPRO_TOT))

      NODE_ME=0
      IONODE=0
#ifdef MPI
      NODE_ME=WDES%COMM%NODE_ME
      IONODE =WDES%COMM%IONODE
#endif

      io_begin
      OPEN(99,FILE='wavefunction.dat',FORM='UNFORMATTED',STATUS='UNKNOWN')

      WRITE(6,*)'Writing out the wavefunction.dat file'

      !First write out some information about the system in general (number of ion, number of ion types, and type of each ion)
      !Also the POSTION of each ion
      !WRITE(99) WDES%NIONS, WDES%NTYP
      WRITE(99) T_INFO%NIONS, T_INFO%NTYP
      WRITE(99) T_INFO%ITYP
      WRITE(99) T_INFO%POSION

      !Write out the number of k points, number of bands, energy cutoff of basis set, number of planewaves and projectors, number of spins
      WRITE(99) WDES%NKPTS, WDES%NB_TOT, WDES%ENMAX, NPL_TOT, NPRO_TOT, WDES%ISPIN
      WRITE(99) ((LATT_INI%A(I,J),I=1,3),J=1,3)
      io_end

      !Write out information about each k-point that is spin independent
      !For each k point, write out the number of plane wave coefficients, the k-vector, and the symmetry weight
      !Also write out the index of each PW component for each k point (i.e. corresponding g vector written as a
      !linear combination of the reciprical lattice vectors, e.g. g=(igx*kx+igy*ky+igz*kz)+(vkptx*kx+vkpty*ky+vkptz*kz))
      DO K=1,WDES%NKPTS
         CALL SETWDES(WDES,WDES1,K)
         io_begin
         NPL=WDES%NPLWKP_TOT(K)
         WRITE(99) NPL, WDES%VKPT(1:3,K),WDES%WTKPT(K)
         WRITE(99) (WDES%IGX(J,K),WDES%IGY(J,K),WDES%IGZ(J,K),J=1,NPL)
         io_end
      ENDDO


      !For each k point, write out the eigenvalues & occupations (for each band)
      DO ISP=1,WDES%ISPIN
      DO K=1,WDES%NKPTS
         CALL SETWDES(WDES,WDES1,K)
         io_begin
         NPL=WDES%NPLWKP_TOT(K)
         WRITE(99) (DBLE(REAL(W%CELTOT(J,K,ISP),KIND=q)),W%FERTOT(J,K,ISP),J=1,WDES%NB_TOT)
         io_end
         !Then write out the plane wave coefficients for the band
         DO J=1,WDES%NB_TOT
            CALL MRG_PW_BAND(WDES1, J, CW, W%CW(1,1,K,ISP))
#ifdef MPI
            IND=0
            DO NC=1,WDES1%NCOL
               DO I=1,WDES1%PL_COL(NC)
                  CRD(IND+I)=CW(WDES1%PL_INDEX(NC)+I)
               ENDDO
               IND=IND+WDES1%PL_COL(NC)
            ENDDO
#else
            CRD(1:NPL)=W%CW(1:NPL,J,K,ISP)
#endif

!     Write out the pseudo-wavefunction (which is the same as the AE wavefunciton for a norm
!     conserving pseudo-potential)
            do_io WRITE(99) (CRD(I),I=1,NPL)

!     Write out the PAW projected part of the wavefunction
#ifdef MPI
            !Get the local storage index for this band, if bands are distributed amongst processors
            NB_LOC=NB_LOCAL(J,WDES1)
            CPROJ=0.d0
            IF (NB_LOC.GT.0) THEN
               !If this one band is split amongst several procs, merge the data
               CALL MRG_PROJ(WDES1, CPROJ(1), W%CPROJ(1,NB_LOC,K,ISP))
            ENDIF
            !Now get the desired projector coefficients for the desired by by summation (all zero but one)
            CALLMPI( M_sum_g(WDES1%COMM_INTER,CPROJ, WDES1%NPRO_TOT))
#else
            CPROJ(1:NPRO_TOT)=W%CPROJ(1:NPRO_TOT,NB_LOC,K,ISP)
#endif

            CRD(1:NPRO_TOT)=CPROJ(1:NPRO_TOT)
            do_io WRITE(99) (CRD(I),I=1,NPRO_TOT)

         ENDDO
      ENDDO
      ENDDO



      ! Now output the pseudo-potential information for each ion type (this is necessary for PAW to reconstruct that AE wavefunction)
      io_begin
      DO I=1,T_INFO%NTYP
         !Write out the number of l channels, and the max l used for the PAW calculation
         !WRITE(99) P(I)%LDIM, P(I)%LMAX_CALC
         WRITE(99) P(I)%LMAX, P(I)%LMAX_CALC
         !Write out their associated l quantum numbers
         WRITE(99) P(I)%LPS
         !Write out information about the valence of the original atoma s well as the pseudopotential
         WRITE(99) P(I)%ZVALF, P(I)%ZVALF_ORIG
         !Write out the radial grid information, number of points and r values, and integration weights
         WRITE(99) P(I)%R%NMAX
         WRITE(99) P(I)%R%R
         WRITE(99) P(I)%R%SI
         !Write out the (valence) all-electron and pseudo wavefunctions on the r grid
         WRITE(99) P(I)%WAE
         WRITE(99) P(I)%WPS
      ENDDO

      CLOSE(99)
      io_end

      DEALLOCATE(CW,CRD,EIG)

      RETURN
      END SUBROUTINE


      END MODULE
