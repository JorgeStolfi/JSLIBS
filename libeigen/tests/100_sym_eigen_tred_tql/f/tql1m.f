      subroutine tql1m(n,d,e,ierr)
c
      integer n,ierr
      double precision d(n),e(n)
c
c     this subroutine is a translation of the algol procedure tql1,
c     num. math. 11, 293-306(1968) by bowdler, martin, reinsch, and
c     wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 227-240(1971).
c
c     this subroutine finds the eigenvalues of a symmetric
c     tridiagonal matrix by the ql method.
c
c     on input
c
c        n is the order of the matrix.
c
c        d contains the diagonal elements of the input matrix.
c
c        e contains the subdiagonal elements of the input matrix
c          in its last n-1 positions.  e(1) is arbitrary.
c
c      on output
c
c        d contains the eigenvalues in ascending order.  if an
c          error exit is made, the eigenvalues are correct and
c          ordered for indices 1,2,...ierr-1, but may not be
c          the smallest eigenvalues.
c
c        e has been destroyed.
c
c        ierr is set to
c          zero       for normal return,
c          j          if the j-th eigenvalue has not been
c                     determined after 30 iterations.
c
c     calls pythag for  dsqrt(a*a + b*b) .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      integer i,j,r,m,ii,r1,r2,mmr
      double precision c,c2,c3,dr1,er1,f,g,h,p,rr,s,s2,tst1,tst2,pythag

      ierr = 0
      if (n .eq. 1) go to 1001
c
      do 100 i = 2, n
  100 e(i-1) = e(i)
      e(n) = 0.0d0
      call tq1s290(n,d,e,ierr)
 1001 return
      end
c
c
c
      subroutine tq1s290(n,d,e,ierr)
      integer n,ierr
      double precision d(n),e(n)
      
      integer i,j,r,m,ii,r1,r2,mmr
      double precision c,c2,c3,dr1,er1,f,g,h,p,rr,s,s2,tst1,tst2,pythag

      f = 0.0d0
      tst1 = 0.0d0
      
      do 290 r = 1, n
         j = 0
         h = dabs(d(r)) + dabs(e(r))
         if (tst1 .lt. h) tst1 = h
c     .......... look for small sub-diagonal element ..........
         do 110 m = r, n
            tst2 = tst1 + dabs(e(m))
            if (tst2 .eq. tst1) go to 120
c     .......... e(n) is always zero, so there is no exit
c                through the bottom of the loop ..........
  110    continue
c
  120    if (m .eq. r) go to 210
  130    if (j .eq. 30) go to 1000
         j = j + 1
c     .......... form shift ..........
         r1 = r + 1
         r2 = r1 + 1
         g = d(r)
         p = (d(r1) - g) / (2.0d0 * e(r))
         rr = pythag(p,1.0d0)
         d(r) = e(r) / (p + dsign(rr,p))
         d(r1) = e(r) * (p + dsign(rr,p))
         dr1 = d(r1)
         h = g - d(r)
         if (r2 .gt. n) go to 145
c
         do 140 i = r2, n
  140    d(i) = d(i) - h
c
  145    f = f + h
         call tq1s200(n,d,e,r,r1,m,dr1)
         tst2 = tst1 + dabs(e(r))
         if (tst2 .gt. tst1) go to 130
  210    p = d(r) + f
         call tq1s230(n,d,p,r)
  290 continue
      return
c     .......... set error -- no convergence to an
c                eigenvalue after 30 iterations ..........
 1000 ierr = r
      return
      end
c
c
c
      subroutine tq1s200(n,d,e,r,r1,m,dr1)
      integer n,r,r1,m
      double precision d(n),e(n),dr1
c     .......... ql transformation ..........

      integer i,j,ii,mmr
      double precision c,c2,c3,er1,f,g,h,p,rr,s,s2,tst1,tst2,pythag

         p = d(m)
         c = 1.0d0
         c2 = c
         er1 = e(r1)
         s = 0.0d0
         mmr = m - r
c     .......... for i=m-1 step -1 until r do -- ..........
         do 200 ii = 1, mmr
            c3 = c2
            c2 = c
            s2 = s
            i = m - ii
            g = c * e(i)
            h = c * p
            rr = pythag(p,e(i))
            e(i+1) = s * rr
            s = e(i) / rr
            c = p / rr
            p = c * d(i) - s * g
            d(i+1) = h + s * (c * g + s * d(i))
  200    continue
c
         p = -s * s2 * c3 * er1 * e(r) / dr1
         e(r) = s * p
         d(r) = c * p
      return
      end
c
c
c
      subroutine tq1s230(n,d,p,r)
      integer n,r
      double precision d(n),p
c     .......... order eigenvalues ..........
      integer i,ii
      double precision c,c2,c3,dr1,er1,f,g,h,rr,s,s2,tst1,tst2,pythag

      if (r .eq. 1) go to 250
c     .......... for i=r step -1 until 2 do -- ..........
      do 230 ii = 2, r
         i = r + 2 - ii
         if (p .ge. d(i-1)) go to 270
         d(i) = d(i-1)
  230 continue
c
  250    i = 1
  270    d(i) = p
      return
      end
