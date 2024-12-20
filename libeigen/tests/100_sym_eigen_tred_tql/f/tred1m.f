      subroutine tred1m(nm,n,a,d,e,e2)
c
      integer nm,n
      double precision a(nm,n),d(n),e(n),e2(n)
c
c     this subroutine is a translation of the algol procedure tred1,
c     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
c
c     this subroutine reduces a real symmetric matrix
c     to a symmetric tridiagonal matrix using
c     orthogonal similarity transformations.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        a contains the real symmetric input matrix.  only the
c          lower triangle of the matrix need be supplied.
c
c     on output
c
c        a contains information about the orthogonal trans-
c          formations used in the reduction in its strict lower
c          triangle.  the full upper triangle of a is unaltered.
c
c        d contains the diagonal elements of the tridiagonal matrix.
c
c        e contains the subdiagonal elements of the tridiagonal
c          matrix in its last n-1 positions.  e(1) is set to zero.
c
c        e2 contains the squares of the corresponding elements of e.
c          e2 may coincide with e if the squares are not needed.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      integer i,ii

      do 100 i = 1, n
         d(i) = a(n,i)
         a(n,i) = a(i,i)
  100 continue
c     .......... for i=n step -1 until 1 do -- ..........
      do 300 ii = 1, n
        i = n + 1 - ii
        call t1s300(nm,n,a,d,e,e2,i)
  300 continue
      return
      end
c
c
c
      subroutine t1s300(nm,n,a,d,e,e2,i)
      integer nm,n,i
      double precision a(nm,n),d(n),e(n),e2(n)

      integer j,k,r
      double precision f,g,h,scale

      r = i - 1
      h = 0.0d0
      scale = 0.0d0
      if (r .lt. 1) go to 130
c        -- if (i ge 2): --
c        .......... scale row (algol tol then not needed) ..........
         do 120 k = 1, r
  120    scale = scale + dabs(d(k))
c
         if (scale .ne. 0.0d0) go to 140
c
            do 125 j = 1, r
               d(j) = a(r,j)
               a(r,j) = a(i,j)
               a(i,j) = 0.0d0
  125       continue
c
  130 e(i) = 0.0d0
      e2(i) = 0.0d0
      return
c
  140 do 150 k = 1, r
         d(k) = d(k) / scale
         h = h + d(k) * d(k)
  150 continue
c
      e2(i) = scale * scale * h
      f = d(r)
      g = -dsign(dsqrt(h),f)
      e(i) = scale * g
      h = h - f * g
      d(r) = f - g
      if (r .eq. 1) go to 285
      call t1s280(nm,n,a,d,e,e2,r,f,g,h)
c
  285 do 290 j = 1, r
         f = d(j)
         d(j) = a(r,j)
         a(r,j) = a(i,j)
         a(i,j) = f * scale
  290 continue
      return
      end
c
c
c
      subroutine t1s280(nm,n,a,d,e,e2,r,f,g,h)
      integer nm,n,r
      double precision a(nm,n),d(n),e(n),e2(n)
      double precision f,g,h

      integer j,k,jp1

c     .......... form a*u ..........
      do 170 j = 1, r
  170 e(j) = 0.0d0
c
      do 240 j = 1, r
         f = d(j)
         g = e(j) + a(j,j) * f
         jp1 = j + 1
         if (r .lt. jp1) go to 220
c
         do 200 k = jp1, r
            g = g + a(k,j) * d(k)
            e(k) = e(k) + a(k,j) * f
  200    continue
c
  220    e(j) = g
  240 continue
c     .......... form p ..........
      f = 0.0d0
c
      do 245 j = 1, r
         e(j) = e(j) / h
         f = f + e(j) * d(j)
  245 continue
c
      h = f / (h + h)
c     .......... form q ..........
      do 250 j = 1, r
  250 e(j) = e(j) - h * d(j)
c     .......... form reduced a ..........
      do 280 j = 1, r
         f = d(j)
         g = e(j)
c
         do 260 k = j, r
  260    a(k,j) = a(k,j) - f * e(k) - g * d(k)
c
  280 continue
      return
      end
      