      subroutine tql2m(nm,n,d,e,z,ierr)
c
      integer n,nm,ierr
      double precision d(n),e(n),z(nm,n)
c
c     this subroutine is a translation of the algol procedure tql2,
c     num. math. 11, 293-306(1968) by bowdler, martin, reinsch, and
c     wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 227-240(1971).
c
c     this subroutine finds the eigenvalues and eigenvectors
c     of a symmetric tridiagonal matrix by the ql method.
c     the eigenvectors of a full symmetric matrix can also
c     be found if  tred2  has been used to reduce this
c     full matrix to tridiagonal form.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        d contains the diagonal elements of the input matrix.
c
c        e contains the subdiagonal elements of the input matrix
c          in its last n-1 positions.  e(1) is arbitrary.
c
c        z contains the transformation matrix produced in the
c          reduction by  tred2, if performed.  if the eigenvectors
c          of the tridiagonal matrix are desired, z must contain
c          the identity matrix.
c
c      on output
c
c        d contains the eigenvalues in ascending order.  if an
c          error exit is made, the eigenvalues are correct but
c          unordered for indices 1,2,...,ierr-1.
c
c        e has been destroyed.
c
c        z contains orthonormal eigenvectors of the symmetric
c          tridiagonal (or full) matrix.  if an error exit is made,
c          z contains the eigenvectors associated with the stored
c          eigenvalues.
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

      integer i,j,k,r,m,ii,r1,r2,mmr
      double precision c,c2,c3,dr1,er1,f,g,h,p,rr,s,s2,tst1,tst2,hypot

      ierr = 0
      if (n .eq. 1) go to 1001
c
      do 100 i = 2, n
  100 e(i-1) = e(i)
      e(n) = 0.0d0
c
      call tq2s240(nm,n,d,e,z,ierr)
      if (ierr .ne. 0) go to 1001
      call tq2s300(nm,n,d,z)
      go to 1001
 1001 return
      end
c
c
c
      subroutine tq2s240(nm,n,d,e,z,ierr) 
      integer n,nm,ierr
      double precision d(n),e(n),z(nm,n)
c
      double precision dr1,f,g,h,p,rr,tst1,tst2,hypot
      integer i,j,k,r,m,r1,r2
c
      f = 0.0d0
      tst1 = 0.0d0

      do 240 r = 1, n
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
  120    continue
         if (m .eq. r) go to 220
  130    continue
         if (j .eq. 30) go to 1000
         j = j + 1
c     .......... form shift ..........
         r1 = r + 1
         r2 = r1 + 1
         g = d(r)
         p = (d(r1) - g) / (2.0d0 * e(r))
         rr = hypot(p,1.0d0)
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
         call tq2s200(nm,n,d,e,z,m,r,r1,dr1)
         tst2 = tst1 + dabs(e(r))
         if (tst2 .gt. tst1) go to 130
  220    d(r) = d(r) + f
  240 continue
      return
c     .......... set error -- no convergence to an
c                eigenvalue after 30 iterations ..........
 1000 ierr = r
      return
      end
c
c
c         
      subroutine tq2s200(nm,n,d,e,z,m,r,r1,dr1)
      integer n,nm,m,r,r1
      double precision d(n),e(n),z(nm,n),dr1
c     .......... ql transformation ..........
      
      double precision c,c2,c3,er1,g,h,p,rr,s,s2,hypot
      integer i,j,k,ii,mmr
      
         p = d(m)
         c = 1.0d0
         c2 = c
         er1 = e(r1)
         s = 0.0d0
         mmr = m - r
         print 7000, "!! m=", m, " r=", r, " r1=", r1, " dr1=", dr1
7000     format(a,i1,a,i1,a,i1,a,f20.14)
c     .......... for i=m-1 step -1 until r do -- ..........
         do 200 ii = 1, mmr
            i = m - ii
4000        format(a)
            print 7100, "!! i=", i, " p=", p, " c=", c, " s=", s
7100        format(a,i1,a,f20.14,a,f20.14,a,f20.14)
            c3 = c2
            c2 = c
            s2 = s
            g = c * e(i)
            h = c * p
            rr = hypot(p,e(i))
            print 7200, "!! i=", i, " s=", s, " e(i)=", e(i), " rr=", rr
7200        format(a,i1,a,f20.14,a,f20.14,a,f24.18)
            e(i+1) = s * rr
            s = e(i) / rr
            c = p / rr
            p = c * d(i) - s * g
            d(i+1) = h + s * (c * g + s * d(i))
            print 6000, '!! i = ', i, '  d(i+1) =  ', d(i+1)
 6000       format(a,i1,a,f20.14)            
c     .......... form vector ..........
            do 180 k = 1, n
               h = z(k,i+1)
               z(k,i+1) = s * z(k,i) + c * h
               z(k,i) = c * z(k,i) - s * h
  180       continue
c
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
      subroutine tq2s300(nm,n,d,z)
      integer n,nm
      double precision d(n),z(nm,n)
c     .......... order eigenvalues and eigenvectors ..........
      double precision p
      
      print 4100, '!! nev = ', n
 4100 format(a,i1)     
      do 666 i = 1, n
        print 4200, '!!   ', d(i)
  666 continue
 4200 format(a,f14.10)     
      do 300 ii = 2, n
         i = ii - 1
         k = i
         p = d(i)
c
         do 260 j = ii, n
            if (d(j) .ge. p) go to 260
            k = j
            p = d(j)
  260    continue
c
         if (k .eq. i) go to 300
         d(k) = d(i)
         d(i) = p
c
         do 280 j = 1, n
            p = z(j,i)
            z(j,i) = z(j,k)
            z(j,k) = p
  280    continue
c
  300 continue
      return
      end
      