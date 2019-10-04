MODULE sortutil

  INTERFACE quick_sort
     MODULE PROCEDURE quick_sort1
     MODULE PROCEDURE quick_sort2
  END INTERFACE

  CONTAINS

    !In place heap sort, courtesty of NR (F77 version
    !ported to F90)
    SUBROUTINE heap_sort(ra)
      REAL*8 :: ra(:)
      INTEGER :: i,ir,j,l,n
      REAL*8 :: rra

      n=SIZE(ra,1)

      IF (n.LT.2) RETURN
      l=n/2+1
      ir=n
10    CONTINUE
        IF(l.GT.1)THEN
          l=l-1
          rra=ra(l)
        ELSE
          rra=ra(ir)
          ra(ir)=ra(1)
          ir=ir-1
          IF(ir.EQ.1)THEN
            ra(1)=rra
            RETURN
          ENDIF
        ENDIF
        i=l
        j=l+l
20      IF(j.LE.ir)THEN
          IF(j.LT.ir)THEN
            IF(ra(j).LT.ra(j+1))j=j+1
          ENDIF
          IF(rra.LT.ra(j))THEN
            ra(i)=ra(j)
            i=j
            j=j+j
          ELSE
            j=ir+1
          ENDIF
        GOTO 20
        ENDIF
        ra(i)=rra
      GOTO 10
    END SUBROUTINE heap_sort

    SUBROUTINE quick_sort1(arr)
      INTEGER :: N, M,NSTACK
      REAL*8 :: arr(:)
      PARAMETER (M=7,NSTACK=50)
      INTEGER :: i,ir,j,jstack,k,l,istack(NSTACK)
      REAL*8 a,temp

      N=SIZE(arr,1)
      jstack=0
      l=1
      ir=n
1     IF(ir-l.LT.M)THEN
        DO 12 j=l+1,ir
          a=arr(j)
          DO 11 i=j-1,l,-1
            IF(arr(i).LE.a)GOTO 2
            arr(i+1)=arr(i)
11        CONTINUE
          i=l-1
2         arr(i+1)=a
12      CONTINUE
        IF(jstack.EQ.0)RETURN
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      ELSE
        k=(l+ir)/2
        temp=arr(k)
        arr(k)=arr(l+1)
        arr(l+1)=temp
        IF(arr(l).GT.arr(ir))THEN
          temp=arr(l)
          arr(l)=arr(ir)
          arr(ir)=temp
        ENDIF
        IF(arr(l+1).GT.arr(ir))THEN
          temp=arr(l+1)
          arr(l+1)=arr(ir)
          arr(ir)=temp
        ENDIF
        IF(arr(l).GT.arr(l+1))THEN
          temp=arr(l)
          arr(l)=arr(l+1)
          arr(l+1)=temp
        ENDIF
        i=l+1
        j=ir
        a=arr(l+1)
3       CONTINUE
          i=i+1
        IF(arr(i).LT.a)GOTO 3
4       CONTINUE
          j=j-1
        IF(arr(j).GT.a)GOTO 4
        IF(j.LT.i)GOTO 5
        temp=arr(i)
        arr(i)=arr(j)
        arr(j)=temp
        GOTO 3
5       arr(l+1)=arr(j)
        arr(j)=a
        jstack=jstack+2
        IF(jstack.GT.NSTACK)PAUSE 'NSTACK too small in sort'
        IF(ir-i+1.GE.j-l)THEN
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        ELSE
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        ENDIF
      ENDIF
      GOTO 1
    END SUBROUTINE quick_sort1

    !
    !Quick sorts arr with corresponding arrangements in brr;
    !thus brr can be used to idenitfy the permutations
    SUBROUTINE quick_sort2(arr,brr)
      INTEGER :: n,M,NSTACK
      REAL*8 :: arr(:)
      INTEGER :: brr(:)
      PARAMETER (M=7,NSTACK=50)
      INTEGER :: i,ir,j,jstack,k,l,istack(NSTACK)
      REAL*8 :: a,b,temp

      N=SIZE(arr,1)
      jstack=0
      l=1
      ir=n
1     IF(ir-l.LT.M)THEN
        DO 12 j=l+1,ir
          a=arr(j)
          b=brr(j)
          DO 11 i=j-1,l,-1
            IF(arr(i).LE.a)GOTO 2
            arr(i+1)=arr(i)
            brr(i+1)=brr(i)
11        CONTINUE
          i=l-1
2         arr(i+1)=a
          brr(i+1)=b
12      CONTINUE
        IF(jstack.EQ.0)RETURN
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      ELSE
        k=(l+ir)/2
        temp=arr(k)
        arr(k)=arr(l+1)
        arr(l+1)=temp
        temp=brr(k)
        brr(k)=brr(l+1)
        brr(l+1)=temp
        IF(arr(l).GT.arr(ir))THEN
          temp=arr(l)
          arr(l)=arr(ir)
          arr(ir)=temp
          temp=brr(l)
          brr(l)=brr(ir)
          brr(ir)=temp
        ENDIF
        IF(arr(l+1).GT.arr(ir))THEN
          temp=arr(l+1)
          arr(l+1)=arr(ir)
          arr(ir)=temp
          temp=brr(l+1)
          brr(l+1)=brr(ir)
          brr(ir)=temp
        ENDIF
        IF(arr(l).GT.arr(l+1))THEN
          temp=arr(l)
          arr(l)=arr(l+1)
          arr(l+1)=temp
          temp=brr(l)
          brr(l)=brr(l+1)
          brr(l+1)=temp
        ENDIF
        i=l+1
        j=ir
        a=arr(l+1)
        b=brr(l+1)
3       CONTINUE
          i=i+1
        IF(arr(i).LT.a)GOTO 3
4       CONTINUE
          j=j-1
        IF(arr(j).GT.a)GOTO 4
        IF(j.LT.i)GOTO 5
        temp=arr(i)
        arr(i)=arr(j)
        arr(j)=temp
        temp=brr(i)
        brr(i)=brr(j)
        brr(j)=temp
        GOTO 3
5       arr(l+1)=arr(j)
        arr(j)=a
        brr(l+1)=brr(j)
        brr(j)=b
        jstack=jstack+2
        IF(jstack.GT.NSTACK)PAUSE 'NSTACK too small in sort2'
        IF(ir-i+1.GE.j-l)THEN
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        ELSE
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        ENDIF
      ENDIF
      GOTO 1
    END SUBROUTINE quick_sort2

END MODULE sortutil
