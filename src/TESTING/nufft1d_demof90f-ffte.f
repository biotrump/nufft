cc Copyright (C) 2004-2009: Leslie Greengard and June-Yub Lee
cc Contact: greengard@cims.nyu.edu
cc
cc This software is being released under a FreeBSD license
cc (see license.txt in this directory).
c
      program testfft
      implicit none
c
c --- local variables
c
      integer i,ier,iflag,j,k1,mx,ms,nj
      parameter (mx=10 000)
      real*4 xj(mx), sk(mx)
      real*4 eps,pi
      real*8 err
      parameter (pi=3.141592653589793238462643383279502884197e0)
      complex*8 cj(mx),cj0(mx),cj1(mx)
      complex*8 fk0(mx),fk1(mx)
c
c     --------------------------------------------------
c     create some test data
c     --------------------------------------------------
      ms = 90
      nj = 128
      do k1 = -nj/2, (nj-1)/2
         j = k1+nj/2+1
         xj(j) = pi * cos(-pi*j/nj)
         cj(j) = cmplx( sin(pi*j/nj), cos(pi*j/nj))
c         print*,xj(j),cj(j)
      enddo
c
c     --------------------------------------------------
c     start tests
c     --------------------------------------------------
c
      iflag = 1
      print*,' Start 1D testing: ', ' nj =',nj, ' ms =',ms
      do i = 1,5
         if (i.eq.1) eps=1e-3
         if (i.eq.2) eps=1e-4
         if (i.eq.3) eps=1e-5
         if (i.eq.4) eps=1e-6
         if (i.eq.5) eps=1e-7
c extended/quad precision tests
c         if (i.eq.5) eps=1d-20
c         if (i.eq.6) eps=1d-24
c         if (i.eq.7) eps=1d-28
c         if (i.eq.8) eps=1d-32
	 print*,' '
  	 print*,' Requested precision eps =',eps
	 print*,' '
c
c     -----------------------
c     call 1D Type1 method
c     -----------------------
c
		iflag = -1
         call dirft1d1f(nj,xj,cj,iflag, ms,fk0)
         call nufft1d1ff90_ffte(nj,xj,cj,iflag,eps, ms,fk1,ier)
         call errcompf(fk0,fk1,ms,err)
         print *,' ier = ',ier
         print *,' type 1 error = ',err
c
c     -----------------------
c     call 1D Type2 method
c     -----------------------
c
         call dirft1d2f(nj,xj,cj0,iflag, ms,fk0,ier)
         call nufft1d2ff90_ffte(nj,xj,cj1,iflag, eps, ms,fk0,ier)
         call errcompf(cj0,cj1,nj,err)
         print *,' ier = ',ier
         print *,' type 2 error = ',err
c
c     -----------------------
c     call 1D Type3 method
c     -----------------------
         do k1 = 1, ms
            sk(k1) = 48*cos(k1*pi/ms)
         enddo
         call dirft1d3f(nj,xj,cj,iflag, ms,sk,fk0)
         call nufft1d3ff90_ffte(nj,xj,cj,iflag,eps, ms,sk,fk1,ier)
         call errcompf(cj0,cj1,nj,err)
         print *,' ier = ',ier
         print *,' type 3 error = ',err
      enddo
      stop
      end
c
c
c
c
c
      subroutine errcompf(fk0,fk1,n,err)
      implicit none
      integer k,n
      complex*8 fk0(n), fk1(n)
      real *8 salg,ealg,err
c
      ealg = 0e0
      salg = 0e0

      do k = 1, n
         ealg = ealg + REAL(cabs(fk1(k)-fk0(k)))**2
         salg = salg + REAL(cabs(fk0(k)))**2
C         print *,fk1(k),fk0(k)
      enddo
C      print *,'ealg,salg:',ealg,salg
      err =sqrt(ealg/salg)
      return
      end
