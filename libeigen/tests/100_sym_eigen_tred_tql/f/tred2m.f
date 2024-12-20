      subroutine tred2m(nm,n,a,d,e,z)
c
      integer n,nm
      double precision a(nm,n),d(n),e(n),z(nm,n)
c
c     this subroutine is a translation of the algol procedure tred2,
c     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
c
c     this subroutine reduces a real symmetric matrix to a
c     symmetric tridiagonal matrix using and accumulating
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
c        d contains the diagonal elements of the tridiagonal matrix.
c
c        e contains the subdiagonal elements of the tridiagonal
c          matrix in its last n-1 positions.  e(1) is set to zero.
c
c        z contains the orthogonal transformation matrix
c          produced in the reduction.
c
c        a and z may coincide.  if distinct, a is unaltered.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      integer i,j,ii

      do 100 i = 1, n
         do 80 j = i, n
   80      z(j,i) = a(j,i)
         d(i) = a(n,i)
  100 continue
c
      if (n .eq. 1) go to 510
c     .......... for i=n step -1 until 2 do -- ..........
      do 300 ii = 2, n
         i = n + 2 - ii
         call t2s300(nm,n,a,d,e,z,i)
  300 continue
c     .......... accumulation of transformation matrices ..........
      do 500 i = 2, n
        call t2s500(nm,n,a,d,e,z,i)
  500 continue
c
  510 do 520 i = 1, n
         d(i) = z(n,i)
         z(n,i) = 0.0d0
  520 continue
c
      z(n,n) = 1.0d0
      e(1) = 0.0d0
      return
      end
c
c
c
      subroutine t2s300(nm,n,a,d,e,z,i)         
      integer n,nm,i
      double precision a(nm,n),d(n),e(n),z(nm,n)

      integer j,k,r,jp1
      double precision f,g,h,hh,scale

      r = i - 1
      h = 0.0d0
      scale = 0.0d0
      if (r .lt. 2) go to 130
c     .......... scale row (algol tol then not needed) ..........
      do 120 k = 1, r
  120 scale = scale + dabs(d(k))
c
      if (scale .ne. 0.0d0) go to 140
  130 e(i) = d(r)
c
      do 135 j = 1, r
         d(j) = z(r,j)
         z(i,j) = 0.0d0
         z(j,i) = 0.0d0
  135 continue
c
      go to 290
c
  140 do 150 k = 1, r
         d(k) = d(k) / scale
         h = h + d(k) * d(k)
  150 continue
c
      f = d(r)
      g = -dsign(dsqrt(h),f)
      e(i) = scale * g
      h = h - f * g
      d(r) = f - g
c     .......... form a*u ..........
      do 170 j = 1, r
  170 e(j) = 0.0d0
c
      do 240 j = 1, r
         f = d(j)
         z(j,i) = f
         g = e(j) + z(j,j) * f
         jp1 = j + 1
         if (r .lt. jp1) go to 220
c
         do 200 k = jp1, r
            g = g + z(k,j) * d(k)
            e(k) = e(k) + z(k,j) * f
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
      hh = f / (h + h)
c     .......... form q ..........
      do 250 j = 1, r
  250 e(j) = e(j) - hh * d(j)
c     .......... form reduced a ..........
      do 280 j = 1, r
         f = d(j)
         g = e(j)
c
         do 260 k = j, r
  260    z(k,j) = z(k,j) - f * e(k) - g * d(k)
c
         d(j) = z(r,j)
         z(i,j) = 0.0d0
  280 continue
c
  290 d(i) = h
      return 
      end
c
c
c
      subroutine t2s500(nm,n,a,d,e,z,i)        
      integer n,nm,i
      double precision a(nm,n),d(n),e(n),z(nm,n)

      integer j,k,r
      double precision g,h
      
      r = i - 1
      z(n,r) = z(r,r)
      z(r,r) = 1.0d0
      h = d(i)
      if (h .eq. 0.0d0) go to 380
c
      do 330 k = 1, r
  330 d(k) = z(k,i) / h
c
      do 360 j = 1, r
         g = 0.0d0
c
         do 340 k = 1, r
  340    g = g + z(k,i) * z(k,j)
c
         do 360 k = 1, r
            z(k,j) = z(k,j) - g * d(k)
  360 continue
c
  380 do 400 k = 1, r
  400 z(k,i) = 0.0d0
c
      return
      end