*     Last edited on 2002-12-31 17:07:08 by stolfi */
c     test of tred1.f, tred2.f
c
      double precision a(10,7), z(10,7), d(7), e(7), e2(7), v(7)
      integer i, j, nm, n, nt, it, ierr
      logical test1, test2
      n = 7
      nm = 10
      nt = 3
      test1 = .true.
      test2 = .true.
      
      do 500 it = 0,nt-1
        print 5000, it
        if (.not. test1) go to 100
          print 4000, ''
          print 4000, 'F - WITHOUT Z'
          print 4000, ''
          call filla(nm,n,a,it)
          print 4000, 'F - input matrix'
          call prta(nm,n,a)
          call tred1(nm,n,a,d,e,e2)
          print 4000, 'F - tridiagonal form'
          call prtri(n,d,e)
          call tql1(n,d,e,ierr)
          print 4000, 'F - eigenvalues'
          call preval(n,d,ierr)
          print 4000, ''
 100    continue
      
        if (.not. test2) go to 200
          print 4000, ''
          print 4000, 'F - WITH Z'
          print 4000, ''
          call filla(nm,n,a,it)
          print 4000, 'F - input matrix'
          call prta(nm,n,a)
          call tred2(nm,n,a,d,e,z)
          print 4000, 'F - tridiagonal form'
          call prtri(n,d,e)
          print 4000, 'F - orthogonal transf (tridiagonal)'
          call prevec(nm,n,z,0)
          call filla(nm,n,a,it)
          call appsim(nm,n,z,a,v)
          print 4000, 'F - transformed matrix (tridiagonal)'
          call prta(nm,n,a)
          call tql2(nm,n,d,e,z,ierr)
          print 4000, 'F - eigenvalues'
          call preval(n,d,ierr)
          print 4000, 'F - eigenvectors'
          call prevec(nm,n,z,ierr)
          call filla(nm,n,a,it)
          call appsim(nm,n,z,a,v)
          print 4000, 'F - transformed matrix (diagonal)'
          call prta(nm,n,a)
          print 4000, ''
 200    continue
 500  continue
 
 4000 format(a)
 5000 format('=== F test ',i1,' ===')

      stop
      end
      
c     generate a symmetric matrix a(n,n)
      subroutine filla(nm,n,a,it)
      integer nm, n
      double precision a(nm,n)
c
      integer i,j,im,jm
      double precision s, t, ui0, ui1, uj0, uj1
      do 100 i = 1,n
        do 110 j = 1,i
          if (it .ne. 0) go to 310
c           * Randomish: *         
            s = i + j - n - 1
            t = i*j
            a(i,j) = 0.0001d0/((dabs(s)+1.0)*dsqrt(t))
            go to 400
 310      if (it .ne. 1) go to 320
c           * Tridiagonal, with 2x2 blocks: *
            if (i .eq. j) a(i,j) = 2.0d0
            if (abs(i-j) .eq. 1) a(i,j) = 1.0d0
            if (abs(i-j) .gt. 1) a(i,j) = 0.0d0
            go to 400
 320      if (it .ne. 2) go to 330      
c           * Rank 2: *
            im = i-1
            jm = j-1
            ui0 = 1.0d0/(1.0d0+im)
            ui1 = 1.0d0/(1.0d0+im*im)
            uj0 = 1.0d0/(1.0d0+jm)
            uj1 = 1.0d0/(1.0d0+jm*jm)
            a(i,j) = ui0*uj0 + ui1*uj1
            go to 400
 330      print 4000, '** invalid {it}'            
          stop
 400      a(j,i) = a(i,j)
 110    continue
 100  continue
      return
 4000 format(a)     
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
 2001   format(10f14.10)
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
 3001     format(a14,$)
 310    continue
        if (i .gt. 1) print 3002, e(i)
        print 3002, d(i)
        if (i .lt. n) print 3002, e(i+1)
 3002   format(f14.10,$)
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
 3002   format(f14.10)
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
 3002   format(10f14.10)
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
    
   
