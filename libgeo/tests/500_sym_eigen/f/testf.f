*     Last edited on 2002-12-31 17:07:08 by stolfi */
c     test of tred1.f, tred2.f
c
      double precision a(10,7), z(10,7), d(7), e(7), e2(7), v(7)
      integer i, j, nm, n, ierr
      logical test1, test2
      n = 7
      nm = 10
      test1 = .true.
      test2 = .true.
      
      if (.not. test1) go to 100
        print 4000, ''
        print 4000, 'F - WITHOUT Z'
        print 4000, ''
        call filla(nm,n,a)
        print 4000, 'input matrix'
        call prta(nm,n,a)
        call tred1(nm,n,a,d,e,e2)
        print 4000, 'tridiagonal form'
        call prtri(n,d,e)
        call tql1(n,d,e,ierr)
        print 4000, 'eigenvalues'
        call preval(n,d,ierr)
        print 4000, ''
 100    continue
      
      if (.not. test2) go to 200
        print 4000, ''
        print 4000, 'F - WITH Z'
        print 4000, ''
        call filla(nm,n,a)
        print 4000, 'input matrix'
        call prta(nm,n,a)
        call tred2(nm,n,a,d,e,z)
        print 4000, 'tridiagonal form'
        call prtri(n,d,e)
        print 4000, 'orthogonal transf'
        call prevec(nm,n,z,0)
        call filla(nm,n,a)
        call appsim(nm,n,z,a,v)
        print 4000, 'transformed matrix'
        call prta(nm,n,a)
        call tql2(nm,n,d,e,z,ierr)
        print 4000, 'eigenvalues'
        call preval(n,d,ierr)
        print 4000, 'eigenvectors'
        call prevec(nm,n,z,ierr)
        call filla(nm,n,a)
        call appsim(nm,n,z,a,v)
        print 4000, 'transformed matrix'
        call prta(nm,n,a)
        print 4000, ''
 200    continue
 
 4000 format(a)

      stop
      end
      
c     generate a symmetric matrix a(n,n)
      subroutine filla(nm,n,a)
      integer nm, n
      double precision a(nm,n)
c
      integer i,j
      double precision s, t
      do 100 i = 1,n
        do 100 j = 1,i
          s = i - (n+1-j)
          t = i*j
          a(i,j) = 1.0d0/((dabs(s)+1.0)*dsqrt(t))
          a(j,i) = a(i,j)
 100  continue
      return
      end

c     print matrix a(n,n)
      subroutine prta(nm,n,a)
      integer nm, n
      double precision a(nm,n)
c
      integer i,j
      print 2003
      do 200 i = 1,n
        print 2001,  (a(i,j),j=1,n)
 2001   format(10f10.6)
 200  continue
      print 2003
 2003 format()
      return
      end

c     prints symmetric tridiagonal matrix
c     d(1..n) is diagonal, e(2..n) is sub-diagonal
      subroutine prtri(n,d,e)
      integer n
      double precision d(n), e(n)
c
      integer i,j
      print 3003 
      do 300 i = 1,n
        do 310 j = 1,i
 330      if (j .le. i-2) print 3001, ' ' 
 3001     format(a8,$)
 310    continue
        if (i .gt. 1) print 3002, e(i)
        print 3002, d(i)
        if (i .lt. n) print 3002, e(i+1)
 3002   format(f10.6,$)
        print 3003
 300  continue
      print 3003
 3003 format()
      return
      end

c     prints eigenvalue list d(1..n)
      subroutine preval(n,d,ierr)
      integer n
      double precision d(n)
c
      integer i, nev
      print 3003 
      nev = n
      if (ierr .ne. 0) nev = ierr-1
      do 300 i = 1,nev
        print 3002, d(i)
 3002   format(f10.6)
 300  continue
      print 3003
 3003 format()
      return
      end
      
c     prints eigenvectors z(i,*)
      subroutine prevec(nm,n,z,ierr)
      integer nm,n
      double precision z(nm,n)
      integer i, j, nev
      print 3003
      nev = n
      if (ierr .ne. 0) nev = ierr-1
      do 300 j = 1,nev
        print 3002, (z(i,j), i=1,n)
 3002   format(10f10.6)
 300  continue
      print 3003
 3003 format()
      return
      end
      
c     applies similarity transformation z to a
c     i.e. computes (z^t)*a*z, using v as temp storage
      subroutine appsim(nm,n,z,a,v)
      integer nm,n
      double precision v(n), z(nm,n), a(nm,n)
c
      double precision s
      integer i, j, k

c     set "a" to "z" transposed times "a":
      do 100 j = 1,n
c       set "v" to "z" transposed times the column "j" of "a":
        do 200 i = 1,n
          s = 0.0d0
          do 300 k = 1,n
            s = s + z(k,i)*a(k,j)
 300        continue
          v(i) = s
 200      continue
c       now store "v" into colum "j" of "a":
        do 250 i = 1,n
          a(i,j) = v(i)
 250      continue
 100    continue

c     set "a" to "a" times "z":
      do 500 i = 1,n
c       set "v" to row "i" of "a" times "z":
        do 600 j = 1,n
          s = 0.0d0
          do 700 k = 1,n
            s = s + a(i,k)*z(k,j)
 700        continue
          v(j) = s
 600      continue
c       now store "v" into row "i" of "a":
        do 650 j = 1,n
          a(i,j) = v(j)
 650      continue
 500    continue

      return
      end
    
   
