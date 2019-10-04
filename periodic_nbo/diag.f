C     Computes all eigenvalues and eigenvectors of a real symmetric matrix a, which is of size n
C     by n, stored in a physical np by np array. On output, elements of a above the diagonal are 
C     destroyed. d returns the eigenvalues  of a in its first n elements. v is a matrix with the 
C     same logical and physical dimensions as a, whose columns contain, on output, the normalized
C     eigenvectors of a.
      SUBROUTINE diag(a,n,np,d,v)
      IMPLICIT NONE
      INTEGER n,np,nrot
      REAL*8 a(n,n),d(n),v(n,n)
      
      call jacobi(a,n,np,d,v,nrot)
      call eigsrt(d,v,n,np)
      return
      end    

      SUBROUTINE jacobi(a,n,np,d,v,nrot)
      IMPLICIT NONE
      INTEGER n,np,nrot,NMAX
      REAL*8 a(np,np),d(np),v(np,np)
      PARAMETER (NMAX=2000)
C     Computes all eigenvalues and eigenvectors of a real symmetric matrix a, which is of size n
C     by n, stored in a physical np by np array. On output, elements of a above the diagonal are
C     destroyed. d returns the eigenvalues of a in its first n elements. v is a matrix with the same
C     logical and physical dimensions as a, whose columns contain, on output, the normalized
C     eigenvectors of a. nrot returns the number of Jacobi rotations that were required.
      INTEGER i,ip,iq,j
      REAL*8 c,g,h,s,sm,t,tau,theta,tresh,b(NMAX),z(NMAX)
      if (n.gt.NMAX) STOP 'n > nmax in jacobi'
      do ip=1,n                 !Initialize to the identity matrix.
         do iq=1,n
            v(ip,iq)=0.
         enddo
         v(ip,ip)=1.
      enddo
      do ip=1,n
         b(ip)=a(ip,ip)         !Initialize b and d to the diagonal of a.
         d(ip)=b(ip)
         z(ip)=0.               !This vector will accumulate terms of the form tapq
                                !as in equation (11.1.14). 
      enddo
      nrot=0
      do i=1,50
         sm=0.
         do ip=1,n-1         !Sum o.-diagonal elements.
            do iq=ip+1,n
               sm=sm+dabs(a(ip,iq))
            enddo
         enddo
         if(sm.eq.0.) return    !The normal return, which relies on quadratic convergence
                                !to machine underflow
         if(i.lt.4)then
            tresh=0.2*sm/n**2   !on the first three sweeps.
         else
            tresh=0.            !thereafter.
         endif
         do ip=1,n-1
            do iq=ip+1,n
               g=100.*dabs(a(ip,iq))
                                !After four sweeps, skip the rotation if the off-diagonal element is small.
               if((i.gt.4).and.(dabs(d(ip))+g.eq.dabs(d(ip)))
     *              .and.(dabs(d(iq))+g.eq.dabs(d(iq))))then
                  a(ip,iq)=0.
               else if(dabs(a(ip,iq)).gt.tresh)then
                  h=d(iq)-d(ip)
                  if(dabs(h)+g.eq.dabs(h))then
                     t=a(ip,iq)/h
                  else
                     theta=0.5*h/a(ip,iq) !Equation (11.1.10).
                     t=1./(dabs(theta)+dsqrt(1.+theta**2))
                     if(theta.lt.0.)t=-t
                  endif
                  c=1./dsqrt(1+t**2)
                  s=t*c
                  tau=s/(1.+c)
                  h=t*a(ip,iq)
                  z(ip)=z(ip)-h
                  z(iq)=z(iq)+h
                  d(ip)=d(ip)-h
                  d(iq)=d(iq)+h
                  a(ip,iq)=0.
                  do j=1,ip-1   !Case of rotations
                     g=a(j,ip)
                     h=a(j,iq)
                     a(j,ip)=g-s*(h+g*tau)
                     a(j,iq)=h+s*(g-h*tau)
                  enddo
                  do j=ip+1,iq-1 !Case of rotations p < j < q.
                     g=a(ip,j)
                     h=a(j,iq)
                     a(ip,j)=g-s*(h+g*tau)
                     a(j,iq)=h+s*(g-h*tau)
                  enddo
                  do j=iq+1,n   !Case of rotations q < j
                     g=a(ip,j)
                     h=a(iq,j)
                     a(ip,j)=g-s*(h+g*tau)
                     a(iq,j)=h+s*(g-h*tau)
                  enddo
                  do j=1,n
                     g=v(j,ip)
                     h=v(j,iq)
                     v(j,ip)=g-s*(h+g*tau)
                     v(j,iq)=h+s*(g-h*tau)
                  enddo
                  nrot=nrot+1
               endif
            enddo
         enddo
         do ip=1,n
            b(ip)=b(ip)+z(ip)
            d(ip)=b(ip)         !Update d with the sum of tapq
            z(ip)=0.            !and reinitialize z.
         enddo
      enddo
      pause 'too many iterations in jacobi'
      return
      END

      SUBROUTINE eigsrt(d,v,n,np)
      IMPLICIT NONE
      INTEGER n,np
      REAL*8 d(np),v(np,np)
C     Given the eigenvalues d and eigenvectors v as output from jacobi (§11.1) or tqli (§11.3),
C     this routine sorts the eigenvalues into accending order, and rearranges the columns of v
C     correspondingly. The method is straight insertion.
      INTEGER i,j,k
      REAL p*8
      do i=1,n-1
         k=i
         p=d(i)
         do j=i+1,n
            if(d(j).le.p)then ! changed from ge to le for accending order
               k=j
               p=d(j)
            endif
         enddo
         if(k.ne.i)then
            d(k)=d(i)
            d(i)=p
            do j=1,n
               p=v(j,i)
               v(j,i)=v(j,k)
               v(j,k)=p
            enddo
         endif
      enddo
      return
      END

      SUBROUTINE ludcmp(a,n,np,indx,d)
      INTEGER n,np,indx(n),NMAX
      REAL*8 d,a(np,np),TINY
      PARAMETER (NMAX=500,TINY=1.0e-20) 
C     Largest expected n, and a small number.
C     Given a matrix a(1:n,1:n), with physical dimension np by np, this routine replaces it by
C     the LU decomposition of a rowwise permutation of itself. a and n are input. a is output,
C     arranged as in equation (2.3.14) above; indx(1:n) is an output vector that records the
C     row permutation e.ected by the partial pivoting; d is output as ±1 depending on whether
C     the number of row interchanges was even or odd, respectively. This routine is used in
C     combination with lubksb to solve linear equations or invert a matrix.
      INTEGER i,imax,j,k
      REAL*8 aamax,dum,sum,vv(NMAX) !vv stores the implicit scaling of each row.
      d=1.                      !No row interchanges yet.
      do i=1,n               !Loop over rows to get the implicit scaling information.
         aamax=0.
         do j=1,n
            if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
         enddo
         if (aamax.eq.0.) pause 'singular matrix in ludcmp' !No nonzero largest element.
         vv(i)=1./aamax         !Save the scaling.
      enddo
      do j=1,n               !This is the loop over columns of Crout’s method.
         do i=1,j-1          !This is equation (2.3.12) except for i = j.
            sum=a(i,j)
            do k=1,i-1
               sum=sum-a(i,k)*a(k,j)
            enddo
            a(i,j)=sum
         enddo
         aamax=0.               !Initialize for the search for largest pivot element.
         do i=j,n               !This is i = j of equation (2.3.12) and i = j+1. . .N
                                !of equation (2.3.13). 
            sum=a(i,j)
            do k=1,j-1
               sum=sum-a(i,k)*a(k,j)
            enddo
            a(i,j)=sum
            dum=vv(i)*abs(sum)  !Figure of merit for the pivot.
            if (dum.ge.aamax) then !Is it better than the best so far?
               imax=i
               aamax=dum
            endif
         enddo
         if (j.ne.imax)then     !Do we need to interchange rows?
            do k=1,n         !Yes, do so...
               dum=a(imax,k)
               a(imax,k)=a(j,k)
               a(j,k)=dum
            enddo
            d=-d                !...and change the parity of d.
            vv(imax)=vv(j)      !Also interchange the scale factor.
         endif
         indx(j)=imax
         if(a(j,j).eq.0.)a(j,j)=TINY
C     If the pivot element is zero the matrix is singular (at least to the precision of the algorithm).
C     For some applications on singular matrices, it is desirable to substitute TINY
C     for zero.
         if(j.ne.n)then         !Now, .nally, divide by the pivot element.
            dum=1./a(j,j)
            do i=j+1,n
               a(i,j)=a(i,j)*dum
            enddo
         endif
      enddo                     !Go back for the next column in the reduction.
      return
      END

c
c Subroutine for Gauss-Jordan elimination--Numerical Recipes p.28-29.
c
      SUBROUTINE gaussj(a,n,np,b,m,mp)
      IMPLICIT REAL*8 (a-h,o-z)
      PARAMETER (nmax=100)
      DIMENSION a(np,np),b(nmax,mp),ipiv(nmax),indxr(nmax),
     *     indxc(nmax)
c
      DO 11 j=1,n
         ipiv(j)=0
 11   CONTINUE
c
      DO 22 i=1,n   !This is the main loop over the columns to be reduced.
         big=0d0
         DO 13 j=1,n   !This is the outer loop of the search for a
c                       pivot element.
            IF (ipiv(j).NE.1) THEN
               DO 12 k=1,n
                  IF (ipiv(k).EQ.0) THEN
                     IF (DABS(a(j,k)).GE.big) THEN
                        big=DABS(a(j,k))
                        irow=j
                        icol=k
                     END IF
                  ELSE IF (ipiv(k).GT.1) THEN
                     PAUSE 'Singular matrix'
                  END IF
 12            CONTINUE
            END IF
 13      CONTINUE
         ipiv(icol)=ipiv(icol)+1
c
c We now have the pivot element, so we interchange rows, if needed, to
c put the pivot element on the diagonal.  The columns are not physically
c interchanged, only relabeled: indx(i), the column of the ith pivot
c element, is the ith column that is reduced, while indxr(i) is the row
c in which that pivot element was originally located.  If indxr(i).NE.
c indxc(i) there is an implied column interchange.  With this form of
c bookkeeping, the solution B's will end up in the correct order, and
c the inverse matrix will be scrambled by columns.
c
         IF (irow.NE.icol) THEN
            DO 14 l=1,n
               dum=a(irow,l)
               a(irow,l)=a(icol,l)
               a(icol,l)=dum
 14         CONTINUE
            DO 15 l=1,m
               dum=b(irow,l)
               b(irow,l)=b(icol,l)
               b(icol,l)=dum
 15         CONTINUE
         END IF
c
c We are now ready to divide the pivot row by the pivot element, located
c at irow and icol.
c
         indxr(i)=irow
         indxc(i)=icol
         IF (a(icol,icol).EQ.0) PAUSE 'Singular matrix'
         pivinv=1/a(icol,icol)
         a(icol,icol)=1d0
         DO 16 l=1,n
            a(icol,l)=a(icol,l)*pivinv
 16      CONTINUE
         DO 17 l=1,m
            b(icol,l)=b(icol,l)*pivinv
 17      CONTINUE
c
c Next we reduce the rows--except for the pivot one of course.
c
         DO 21 ll=1,n
            IF (ll.NE.icol) THEN
               dum=a(ll,icol)
               a(ll,icol)=0d0
               DO 18 l=1,n
                  a(ll,l)=a(ll,l)-a(icol,l)*dum
 18            CONTINUE
               DO 19 l=1,m
                  b(ll,l)=b(ll,l)-b(icol,l)*dum
 19            CONTINUE
            END IF
 21      CONTINUE
 22   CONTINUE   !This ends the main loop over columns of the reduction.
c
c It only remains to unscramble the solution in view of the column
c interchanges.  We do this by interchanging pairs of columns in the
c reverse order that the permutation was built up.
c
      DO 24 l=n,1,-1
         IF (indxr(l).NE.indxc(l)) THEN
            DO 23 k=1,n
               dum=a(k,indxr(l))
               a(k,indxr(l))=a(k,indxc(l))
               a(k,indxc(l))=dum
 23         CONTINUE
         END IF
 24   CONTINUE
c
      RETURN
      END
