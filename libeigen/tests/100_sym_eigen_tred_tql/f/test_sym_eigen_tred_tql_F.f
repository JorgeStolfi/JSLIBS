c     # Last edited on 2024-12-05 17:27:48 by stolfi
c     test of tred1.f, tred2.f tql1.f, tql2.f
c     in their "m" (slightly more structured) versions.
c
      integer nt,tnum

      print 4000, 'starting ...'
      
c     number of test matrix types:
      nt = 5
            
      print 4000, 'opening files ...'
      open(unit=3, status="new", file="out/version1.txt")
      open(unit=4, status="new", file="out/version2.txt")

      do 10 tnum = 0, nt-1            
         print 4100, '  tests with tnum = ', tnum, ' ...'
         call dotests( 3, 1, tnum)
         call dotests( 5, 2, tnum)
         call dotests(10, 7, tnum)
   10 continue
      close(3)
      close(4)
      
 4000 format(a)
 4100 format(a,i1,a)
      stop
      end

      subroutine dotests(nm,n,tnum)
      integer nm, n, tnum
c     tests the two versions of the tridiagonalization and QL procedures
c     on an {n x n} array {a} of type {tnum}. the array will be allocated
c     with dimension {nm,n} where {nm} must be {n} or more. writes the
c     results to output units 3 and 4, respectively.

      double precision a(nm,n), z(nm,n), d(n), e(n), e2(n)

      print 4300, '  tests with nm = ', nm, ' n = ', n, ' tnum = ', tnum, ' ...'

      call ttr1(3, nm,n,a,d,e,e2,tnum)
      call ttr2(4, nm,n,a,d,e,z,tnum)

      return
 4300 format(a,i1,a,i1,a,i1,a)      
      end
      
c
c
c
      subroutine ttr1(unit, nm,n,a,d,e,e2,tnum)
      integer unit,nm,n,tnum
      double precision a(nm,n),d(n),e(n),e2(n)
c
c     tests tred1, tql1      
c      
      integer ierr
         
      write(unit,4000) ''
      write(unit,5000) n,tnum 
 5000 format('=== F OUT - TRED1/TQL1 (NO Z) N=',i1,' TYPE=',i1,' ===')
      write(unit,4000) ''

c     ........ convert matrix a to tridiagonal form ........      
      call filla(nm,n,a,tnum)
      write(unit,4000) 'F OUT - input matrix'
      call prta(unit,nm,n,a)
      call tred1m(nm,n,a,d,e,e2)
      write(unit,4000) 'F OUT - tridiagonal form'
      call prtri(unit,n,d,e)

c     ........ compute eigenvalues ........      
      call tql1m(n,d,e,ierr)
      write(unit,4000) 'F OUT - eigenvalues'
      call preval(unit,n,d,ierr)
      return
 4000 format(a)
      end
c
c
c
      subroutine ttr2(unit,nm,n,a,d,e,z,tnum)
      integer unit,n,nm,tnum
      double precision a(nm,n),d(n),e(n),z(nm,n)
c
c     tests tred2, tql2      
c      
      double precision v(n)
      integer ierr
          
      write(unit,4000) ''
      write(unit,5000) n,tnum 
 5000 format('=== F OUT - TRED2/TQL2 (Z.NE.A) N=',i1,' TYPE=',i1,' ===')
      write(unit,4000) ''
      
c     ........ convert matrix a to tridiagonal form ........      
      call filla(nm,n,a,tnum)
      write(unit,4000) 'F OUT - input matrix'
      call prta(unit,nm,n,a)
      call tred2m(nm,n,a,d,e,z)
      write(unit,4000) 'F OUT - tridiagonal form'
      call prtri(unit,n,d,e)
      write(unit,4000) 'F OUT - orthogonal transf (tridiagonal)'
      call prevec(unit,nm,n,z,0)

c     ........ check matrix z is tridiagonal ........  
      call filla(nm,n,a,tnum)
      call appsim(nm,n,z,a,v)
      write(unit,4000) 'F OUT - transformed matrix (tridiagonal)'
      call prta(unit,nm,n,a)
      
c     ........ compute eigenvalues and eigenvectors ........      
      call tql2m(nm,n,d,e,z,ierr)
      print 4000, 'returned from tql2m'
      print 4001, 'ierr = ', ierr
      write(unit,4000) 'F OUT - eigenvalues'
      call preval(unit,n,d,ierr)
      write(unit,4000) 'F OUT - eigenvectors'
      call prevec(unit,nm,n,z,ierr)
      
c     ........ check matrix z is diagonal ........  
      call filla(nm,n,a,tnum)
      call appsim(nm,n,z,a,v)
      write(unit,4000) 'F OUT - transformed matrix (diagonal)'
      call prta(unit,nm,n,a)
      return
 4000 format(a)
 4001 format(a,i1)
      end
 
      subroutine filla(nm,n,a,tnum)
      integer nm,n,tnum
      double precision a(nm,n)
c
c     fills {a(1..n,1..n)} with a symmetric matrix of type {tnum}
c
      integer i,j,k,r,im,jm,rm
      double precision s,t, ui0,ui1,uj0,uj1, aij,rir,rjr,dr

      print 4200, '    filling matrix n = ', n, ' tnum = ', tnum, ' ...'
 4200 format(a,i1,a,i1,a)

      do 100 i = 1,n
        do 110 j = 1,i
          if (tnum .ne. 0) go to 310
c           * Diagonal: */
             a(i,j) = 0.0d0
             if (i .eq. j) a(i,j) = i
             go to 400
 310      if (tnum .ne. 1) go to 320
c            * Tridiagonal, with 2x2 blocks: *
             k = (i+1)/2
             a(i,j) = 0.0d0
             if (abs(i-j) .gt. 1) go to 400
             if (i .eq. j) a(i,j) = 1.5d0 + 2.0d0*k
             if (modulo(i+j,4) .ne. 3) go to 400
             a(i,j) = 0.5d0
             go to 400
 320      if (tnum .ne. 2) go to 330      
c            * Randomish: *         
             s = i + j - n - 1
             t = i*j
             a(i,j) = 0.0001d0/((dabs(s)+1.0)*dsqrt(t))
             go to 400
 330      if (tnum .ne. 3) go to 340      
c            * Rank 2: *
             im = i-1
             jm = j-1
             ui0 = 1.0d0/(1.0d0+im)
             ui1 = 1.0d0/(1.0d0+im*im)
             uj0 = 1.0d0/(1.0d0+jm)
             uj1 = 1.0d0/(1.0d0+jm*jm)
             a(i,j) = ui0*uj0 + ui1*uj1
             go to 400
 340      if (tnum .ne. 4) go to 350 
             aij = 0d0
             im = i-1
             jm = j-1
             do 342 r = 1, n
                rm = r-1
                rir = cos((im*r*2.0d0/n + 0.25d0)*pi)
                rjr = cos((jm*r*2.0d0/n + 0.25d0)*pi)
                dr = r
                aij = aij + rir*dr*rjr/(2*n)
 342         continue
             a(i,j) = aij
             go to 400
 350      print 4000, '** invalid {tnum}'            
          stop
 400      a(j,i) = a(i,j)
 110    continue
 100  continue
      return
 4000 format(a)     
      end

      subroutine prta(u,nm,n,a)
      integer u,nm,n
      double precision a(nm,n)
c
c     prints matrix a(n,n) to unit u
c
      integer i,j
      write(u,2003)
      do 200 i = 1,n
        write(u,2001) (a(i,j),j=1,n)
 2001   format(10f14.10)
 200  continue
      write(u,2003)
 2003 format()
      return
      end

      subroutine prtri(u,n,d,e)
      integer u,n
      double precision d(n), e(n)
c
c     prints symmetric tridiagonal matrix to unit u
c     d(1..n) is diagonal, e(2..n) is sub-diagonal
c
      integer i,j
      write(u,3003) 
      do 300 i = 1,n
        do 310 j = 1,i
 330      if (j .le. i-2) write(u,3001) ' ' 
 3001     format(a14,$)
 310    continue
        if (i .gt. 1) write(u,3002) e(i)
        write(u,3002) d(i)
        if (i .lt. n) write(u,3002) e(i+1)
 3002   format(f14.10,$)
        write(u,3003) 
 300  continue
      write(u,3003) 
 3003 format()
      return
      end

      subroutine preval(u,n,d,ierr)
      integer u,n,ierr
      double precision d(n)
c
c     prints eigenvalue list d(1..n) to unit u
c
      integer i, nev
      write(u,3003) 
      nev = n
      if (ierr .ne. 0) nev = ierr-1
      do 300 i = 1,nev
        write(u,3002) d(i)
 3002   format(f14.10)
 300  continue
      write(u,3003) 
 3003 format()
      return
      end
      
      subroutine prevec(u,nm,n,z,ierr)
      integer u,nm,n,ierr
      double precision z(nm,n)
c
c     prints eigenvectors z(i,*) to unit u
c
      integer i, j, nev
      write(u,3003) 
      nev = n
      if (ierr .ne. 0) nev = ierr-1
      do 300 j = 1,nev
        write(u,3002) (z(i,j), i=1,n)
 3002   format(10f14.10)
 300  continue
      write(u,3003) 
 3003 format()
      return
      end
      
      subroutine appsim(nm,n,z,a,v)
      integer nm,n
      double precision v(n), z(nm,n), a(nm,n)
c
c     applies similarity transformation z to a
c     i.e. computes (z^t)*a*z, using v as temp storage
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
    
   
