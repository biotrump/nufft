cc Copyright (C) 2004-2009: Leslie Greengard and June-Yub Lee
cc Contact: greengard@cims.nyu.edu
cc
cc This software is being released under a FreeBSD license
cc (see license.txt in this directory).
cc
c
c  NUFFT 1.3 release notes:
c
c  These codes are asymptotically fast (O(N log N)), but not optimized.
c
c  1) We initialize the FFT on every call.
c
c  2) We do not precompute the exponentials involved in "fast Gaussian
c  gridding".
c
c  3) We do not block structure the code so that irregularly placed points
c  are interpolated (gridded) in a cache-aware fashion.
c
c  4) We use the Netlib FFT library (www.netlib.org)
c     rather than the state of the art FFTW package (www.fftw.org).
c
c  Different applications have different needs, and we have chosen
c  to provide the simplest code as a reasonable efficient template.
C
C	!!!!!!!!!!!
C	F90 calls C :
c	http://docs.oracle.com/cd/E19059-01/stud.8/817-5066/11_cfort.html
c**********************************************************************
      subroutine nufft1d1ff90_ffte(nj,xj,cj,iflag,eps,ms,fk,ier)
C		use moudle to load global vars
      USE NUFFTModule
C	TYPE(C_PTR), i.e., "void *", variable needs "ISO_C_BINDING" to define the TYPE(C_PTR)
      USE , intrinsic ::ISO_C_BINDING
      implicit none

		interface
C		fortran always uses call-by-address, so it passes pointers, not value to c functions.
C		ffts_plan_t *ffts_init_1d(size_t n, int sign);
C 		size_t for 64 bits platform is 8 bytes, fortran uses "call-by-addr",
C		so wrong size passing corrupts the memory
		function ffts_init_1df ( n, sign ) result(ptr) bind ( c )
			use, intrinsic :: iso_c_binding
			implicit none
			TYPE(C_PTR) 	:: ptr
			integer ( c_long ) :: n
			integer ( c_int ) :: sign
		end function ffts_init_1df

C		void ffts_execute (ffts_plan_t *plan , const void *input, void *output )
		subroutine ffts_executef (plan , input, output ) bind ( c )
			use, intrinsic :: iso_c_binding
			implicit none
			TYPE(C_PTR) :: plan
			REAL(C_FLOAT) :: input(*)
			REAL(C_FLOAT) :: output(*)
		end subroutine ffts_executef

C		void ffts_free(ffts_plan_t *plan);
		subroutine ffts_freef ( plan ) bind ( c )
			use, intrinsic :: iso_c_binding
			implicit none
			TYPE(C_PTR) :: plan
		end subroutine ffts_freef

C		unsigned find_best_pow2f(size_t *n);the nearest power of 2 next to n
		function find_best_pow2f ( n ) bind ( c )
			use, intrinsic :: iso_c_binding
			implicit none
			integer ( c_long ) :: n
			integer ( c_int ) :: find_best_pow2f
		end function find_best_pow2f
		end interface
C	special vars passed to c functions
C	https://people.sc.fsu.edu/~jburkardt/f_src/f90_calls_c/f90_calls_c.html
		integer (c_int) :: iflag
		integer (c_long) :: nf1

      integer ier,istart,iw1,iwtot,iwsav
      integer j,jb1,jb1u,jb1d,k1,ms,next235,nj,nspread
      real*4 cross,cross1,diff1,eps,hx,pi,rat,r2lamb,t1,tau
      real*4 xc(-147:147),xj(nj)
      parameter (pi=3.141592653589793238462e0)
      complex*8 cj(nj),fk(-ms/2:(ms-1)/2),zz,ccj
C 		void *
      TYPE(C_PTR) :: fftp

c ----------------------------------------------------------------------
C variable size 1-D real*4 arrary: fw[]
      real*4, allocatable, target :: fw(:)
!dir$ attributes align:64
c ----------------------------------------------------------------------
c     if (iflag .ge. 0) then
c
c               1  nj
c     fk(k1) = -- SUM cj(j) exp(+i k1 xj(j))  for -ms/2 <= k1 <= (ms-1)/2
c              nj j=1
c
c     else
c
c               1  nj
c     fk(k1) = -- SUM cj(j) exp(-i k1 xj(j))  for -ms/2 <= k1 <= (ms-1)/2
c              nj j=1
c
c     References:
c
c     [DR] Fast Fourier transforms for nonequispaced data,
c          A. Dutt and V. Rokhlin, SIAM J. Sci. Comput. 14,
c          1368-1383, 1993.
c
c     [GL] Accelerating the Nonuniform Fast Fourier Transform,
c          L. Greengard and J.-Y. Lee, SIAM Review 46, 443-454 (2004).
c
c ----------------------------------------------------------------------
c     INPUT:
c
c     nj     number of sources   (integer)
c     xj     location of sources (real *4)
c
c            on interval [-pi,pi].
c
c     cj     strengths of sources (complex *8)
c     iflag  determines sign of FFT (see above)
c     eps    precision request  (between 1.0d-33 and 1.0d-1)
c               recomended value is 1d-15 for double precision calculations
c     ms     number of Fourier modes computed (-ms/2 to (ms-1)/2 )
c
c     OUTPUT:
c
c     fk     Fourier transform values (complex *16)
c     ier    error return code
c
c            ier = 0  => normal execution.
c            ier = 1  => precision eps requested is out of range.
c
c     The type 1 NUFFT proceeds in three steps (see [GL]).
c
c     1) spread data to oversampled regular mesh using convolution with
c        a Gaussian
c     2) compute FFT on uniform mesh
c     3) deconvolve each Fourier mode independently
c          (mutiplying by Fourier transform of Gaussian).
c
c ----------------------------------------------------------------------
c
c     The oversampled regular mesh is defined by
c
c     nf1 = rat*ms  points, where rat is the oversampling ratio.
c
c     For simplicity, we set
c
c         rat = 2 for eps > 1.0d-11
c         rat = 3 for eps <= 1.0d-11.
c
c     The Gaussian used for convolution is:
c
c        g(x) = exp(-x^2 / 4tau)
c
c     It can be shown [DR] that the precision eps is achieved when
c
c     nspread = int(-log(eps)/(pi*(rat-1d0)/(rat-.5d0)) + .5d0)
c     and tau is chosen as
c
c     tau = pi*lambda/(ms**2)
c     lambda = nspread/(rat(rat-0.5)).
c
c     Note that the Fourier transform of g(x) is
c
c     G(s) = exp(-s^2 tau) = exp(-pi*lambda s^2/ms^2)
c
c
c ----------------------------------------------------------------------
c     Fast Gaussian gridding is based on the following observation.
c
c     Let hx = 2*pi/nf1. In gridding data onto a regular mesh with
c     spacing nf1, we shift the source point xj by pi so
c     that it lies in [0,2*pi] to simplify the calculations.
c     Since we are viewing the function
c     as periodic, this has no effect on the result.
c
c     For source (xj+pi), let kb*hx denote the closest grid point and
c     let  kx*hx be a regular grid point within the spreading
c     distance. We can write
c
c     (xj+pi) - kx*hx = kb*hx + diff*hx - kx*hx = diff*hx - (kx-kb)*hx
c
c     where diff = (xj+pi)/hx - kb.
c
c     Let t1 = hx*hx/(4 tau) = pi/(nf1*nf1)/lambda*ms*ms
c                            = pi/lambda/(rat*rat)
c
c     exp(-( (xj+pi) -kx*hx)**2 / 4 tau)
c         = exp(-pi/lamb/rat^2 *(diff - (kx-kb))**2)
c         = exp(-t1 *(diff - (kx-kb))**2)
c         = exp(-t1*diff**2) * exp(2*t1*diff)**k * exp(-t1*k**2)
c           where k = kx-kb.
c
c************************************************************************
c
c     Precision dependent parameters
c
c     rat is oversampling parameter
c     nspread is number of neighbors to which Gaussian gridding is
c     carried out.
c -------------------------------
      ier = 0
      if ((eps.lt.1d-33).or.(eps.gt.1d-1)) then
         ier = 1
         return
      endif
      if (eps.le.1e-11) then
         rat = 3.0e0
      else
         rat = 2.0e0
      endif
      nspread = int(-log(eps)/(pi*(rat-1e0)/(rat-.5e0)) + .5e0)
      nf1 = rat*ms
      if (2*nspread.gt.nf1) then
         nf1 = next235(DBLE(2e0*nspread))
      endif

C	ffts supports only radix2
C	WRITE(6,*) ">> nf1=",nf1
	if(ffte .eq. E_FFTS) then
		nf1=find_best_pow2f(nf1)
	endif
C	WRITE(6,*) "<< nf1=",nf1

c
c     lambda (described above) = nspread/(rat*(rat-0.5d0))
c     It is more convenient to define r2lamb = rat*rat*lambda
c
      r2lamb = rat*rat * nspread / (rat*(rat-.5e0))
      hx = 2*pi/nf1
c
c     -----------------------------------
c     Compute workspace size and allocate
c     -----------------------------------
      iw1 = 2*nf1
      iwsav = iw1+nspread+1
C	iwsav should be 16 bytes aligned
      if(mod(iwsav,16).ne.0) then
		iwsav = iwsav + 16 - mod(iwsav,16)
      endif
C	if(ffte .eq. E_FFTS) then
CC	ffts has its own workspace, but it needs input/output buffer
C		iwtot = iwsav+2*nf1+15
C	else
CC	ffte needs this extra space
      iwtot = iwsav+4*nf1+15
C	endif

C allocate size "0 to iwtot", 1-D real*4 arrary: fw[0:iwtot]
      allocate ( fw(0:iwtot) )
c
c     ---------------------------------------------------------------
c     Precompute spreading constants and initialize fw
c     to hold one term needed for fast Gaussian gridding
c     ---------------------------------------------------------------
      t1 = pi/r2lamb
      do k1 = 1, nspread
         fw(iw1+k1) = exp(-t1*k1**2)
      enddo

C Thomas init : CALL ZFFT1F(fw(0), nf1, 0, fw(iwsav))
C Thomas dcfftb(nf1,fw(0),fw(iwsav)) ==> CALL ZFFT1F(fw(0), nf1, 1, fw(iwsav))
C Thomas dcfftf(nf1,fw(0),fw(iwsav)) ==> CALL ZFFT1F(fw(0), nf1, -1, 0 8zvgfkk3[r(iwsav))
C
	if(ffte .eq. E_FFTS) then
C	FFTS
C	iflag, < 0, forward, >=0 back
		fftp = ffts_init_1df(nf1, iflag)
	else if(ffte .eq. E_FFTE) then
C	FFTE
		CALL		ZFFT1F(fw(0), nf1, 0, fw(iwsav))
	else
C	netlib fftpack
		call dcffti(nf1,fw(iwsav))
	endif
	WRITE(6,*) "NUFFT type1: init fftp=",fftp, nf1,iwsav,iflag
C	CALL		DUMPF(fw(iwsav), nf1)

c
c     ---------------------------------------------------------------
c     Initialize fine grid data to zero.
c     ---------------------------------------------------------------
      do k1 = 0, 2*nf1-1
         fw(k1) = cmplx(0e0,0e0)
      enddo
c
c     ---------------------------------------------------------------
c     Loop over sources (1,...,nj)
c
c     1. find closest mesh point (with periodic wrapping if necessary)
c     2. spread source data onto nearest nspread grid points
c        using fast Gaussian gridding.
c
c     The following is a little hard to read because it takes
c     advantage of fast gridding and optimized to minimize the
c     the number of multiplies in the inner loops.
c
c    ---------------------------------------------------------------
c
      do j = 1, nj
         ccj = cj(j)/real(nj)

         jb1 = int((xj(j)+pi)/hx)
         diff1 = (xj(j)+pi)/hx - jb1
         jb1 = mod(jb1, nf1)
         if (jb1.lt.0) jb1=jb1+nf1
c
         xc(0) = exp(-t1*diff1**2)
         cross = xc(0)
         cross1 = exp(2e0*t1 * diff1)
         do k1 = 1, nspread
            cross = cross * cross1
            xc(k1) = fw(iw1+k1)*cross
         enddo
         cross = xc(0)
         cross1 = 1e0/cross1
         do k1 = 1, nspread-1
            cross = cross * cross1
            xc(-k1) = fw(iw1+k1)*cross
         enddo
c
         jb1d = min(nspread-1, jb1)
         jb1u = min(nspread, nf1-jb1-1)
         do k1 = -nspread+1, -jb1d-1
	    istart = 2*(jb1+k1+nf1)
            zz=xc(k1)*ccj
            fw(istart)=fw(istart)+real(zz)
            fw(istart+1)=fw(istart+1)+imag(zz)
         enddo
         do k1 = -jb1d, jb1u
	    istart = 2*(jb1+k1)
            zz=xc(k1)*ccj
            fw(istart)=fw(istart)+real(zz)
            fw(istart+1)=fw(istart+1)+imag(zz)
         enddo
         do k1 = jb1u+1, nspread
	    istart = 2*(jb1+k1-nf1)
            zz=xc(k1)*ccj
            fw(istart)=fw(istart)+real(zz)
            fw(istart+1)=fw(istart+1)+imag(zz)
         enddo
      enddo
c
c     ---------------------------------------------------------------
c     Compute 1D FFT and carry out deconvolution.
c
c     There is a factor of (-1)**k1 needed to account for the
c     FFT phase shift.
c     ---------------------------------------------------------------
c

C Thomas init : CALL ZFFT1F(fw(0), nf1, 0, fw(iwsav))
C Thomas backward/inverse FFT : dcfftb(nf1,fw(0),fw(iwsav)) ==> CALL ZFFT1F(fw(0), nf1, 1, fw(iwsav))
C Thomas forward FFT : dcfftf(nf1,fw(0),fw(iwsav)) ==> CALL ZFFT1F(fw(0), nf1, -1, fw(iwsav))

	if(ffte .eq. E_FFTS) then
C	output=[real, imaginary, real, imaginary ... interleaving format]
		call ffts_executef(fftp, fw, fw(iwsav) )
C	copy the output buffer back to the input buffer
C	output : complex, complex.... => real,image,real,image....
C	memcpy(fw(0), fw(iwsav), 2*nf1)
		fw(0:2*nf1-1) = fw(iwsav:iwsav+2*nf1-1)
		call ffts_freef(fftp)
	else
		if (iflag .ge. 0) then
			if(ffte .eq. E_FFTE) then
				CALL		ZFFT1F(fw(0), nf1, 1, fw(iwsav))
			else
				call dcfftb(nf1,fw(0),fw(iwsav))
			endif
		else
			if(ffte .eq. E_FFTE) then
				CALL		ZFFT1F(fw(0), nf1, -1, fw(iwsav))
			else
				call dcfftf(nf1,fw(0),fw(iwsav))
			endif
		endif
	endif
	WRITE(6,*) "NUFFT type1:",nf1,iwsav
C		CALL		DUMPF(fw(0), nf1)

c
      tau = pi * r2lamb / real(nf1)**2
      cross1 = 1e0/sqrt(r2lamb)
      zz = cmplx(fw(0),fw(1))
      fk(0) = cross1*zz
      do k1 = 1, (ms-1)/2
         cross1 = -cross1
         cross = cross1*exp(tau*real(k1)**2)
	 zz = cmplx(fw(2*k1),fw(2*k1+1))
         fk(k1) = cross*zz
	 zz = cmplx(fw(2*(nf1-k1)),fw(2*(nf1-k1)+1))
         fk(-k1) = cross*zz
      enddo
      if (ms/2*2.eq.ms) then
         cross = -cross1*exp(tau*real(ms/2)**2)
         zz = cmplx(fw(2*nf1-ms),fw(2*nf1-ms+1))
         fk(-ms/2) = cross*zz
      endif
      deallocate(fw)
      return
      end
c
c
c
c
c
************************************************************************
      subroutine nufft1d2ff90_ffte(nj,xj,cj, iflag,eps, ms,fk,ier)
C		use moudle to load global vars
      USE NUFFTModule
C	TYPE(C_PTR), i.e., "void *", variable needs "ISO_C_BINDING" to define the TYPE(C_PTR)
      USE , intrinsic ::ISO_C_BINDING
      implicit none

		interface
C		fortran always uses call-by-address, so it passes pointers, not value to c functions.
C		ffts_plan_t *ffts_init_1d(size_t n, int sign);
C 		size_t for 64 bits platform is 8 bytes, fortran uses "call-by-addr",
C		so wrong size passing corrupts the memory
		function ffts_init_1df ( n, sign ) result(ptr) bind ( c )
			use, intrinsic :: iso_c_binding
			implicit none
			TYPE(C_PTR) 	:: ptr
			integer ( c_long ) :: n
			integer ( c_int ) :: sign
		end function ffts_init_1df

C		void ffts_execute (ffts_plan_t *plan , const void *input, void *output )
		subroutine ffts_executef (plan , input, output ) bind ( c )
			use, intrinsic :: iso_c_binding
			implicit none
			TYPE(C_PTR) :: plan
			REAL(C_FLOAT) :: input(*)
			REAL(C_FLOAT) :: output(*)
		end subroutine ffts_executef

C		void ffts_free(ffts_plan_t *plan);
		subroutine ffts_freef ( plan ) bind ( c )
			use, intrinsic :: iso_c_binding
			implicit none
			TYPE(C_PTR) :: plan
		end subroutine ffts_freef

C		unsigned find_best_pow2f(size_t *n);the nearest power of 2 next to n
		function find_best_pow2f ( n ) bind ( c )
			use, intrinsic :: iso_c_binding
			implicit none
			integer ( c_long ) :: n
			integer ( c_int ) :: find_best_pow2f
		end function find_best_pow2f
		end interface
C	special vars passed to c functions
C	https://people.sc.fsu.edu/~jburkardt/f_src/f90_calls_c/f90_calls_c.html
		integer (c_int) :: iflag
		integer (c_long) :: nf1

      integer ier,iw1,iwsav,iwtot,j,jb1,jb1u,jb1d,k1
      integer ms,next235,nj,nspread,nw
      real*4 cross,cross1,diff1,eps,hx,pi,rat,r2lamb,t1
      real*4 xj(nj),xc(-147:147)
      parameter (pi=3.141592653589793238462643383279502884197e0)
      complex*8 cj(nj), fk(-ms/2:(ms-1)/2)
      complex*8 zz
C 		void *
      TYPE(C_PTR) :: fftp

c ----------------------------------------------------------------------
      real*4, allocatable :: fw(:)
!dir$ attributes align:64
c ----------------------------------------------------------------------
c     if (iflag .ge. 0) then
c
c              (ms-1)/2
c     cj(j) =    SUM      fk(k1) exp(+i k1 xj(j))  for j = 1,...,nj
c              k1= -ms/2
c
c     else
c
c              (ms-1)/2
c     cj(j) =    SUM      fk(k1) exp(-i k1 xj(j))  for j = 1,...,nj
c              k1= -ms/2
c
c ----------------------------------------------------------------------
c     INPUT:
c
c     nj     number of output values   (integer)
c     xj     location of output values (real *8 array)
c     iflag  determines sign of FFT (see above)
c     eps    precision request  (between 1.0d-33 and 1.0d-1)
c               recomended value is 1d-15 for double precision calculations
c     ms     number of Fourier modes given  [ -ms/2: (ms-1)/2 ]
c     fk     Fourier coefficient values (complex *16 array)
c
c     OUTPUT:
c
c     cj     output values (complex *16 array)
c     ier    error return code
c
c            ier = 0  => normal execution.
c            ier = 1  => precision eps requested is out of range.
c
c
c     The type 2 algorithm proceeds in three steps (see [GL]).
c
c     1) deconvolve (amplify) each Fourier mode first
c     2) compute inverse FFT on uniform fine grid
c     3) spread data to regular mesh using Gaussian
c
c
c     See subroutine nufft1d1ff90_ffte(nj,xj,cj,iflag,eps,ms,fk,ier)
c     for more comments on fast gridding and parameter selection.
c
************************************************************************
c
c     Precision dependent parameters
c
c     rat is oversampling parameter
c     nspread is number of neighbors to which Gaussian gridding is
c     carried out.
c     -------------------------------
c
c     lambda (described above) = nspread/(rat*(rat-0.5d0))
c     It is more convenient to define r2lamb = rat*rat*lambda
c
c     -------------------------------
      ier = 0
      if ((eps.lt.1d-33).or.(eps.gt.1d-1)) then
         ier = 1
         return
      endif
      if (eps.le.1d-11) then
         rat = 3.0e0
      else
         rat = 2.0e0
      endif
c
      nspread = int(-log(eps)/(pi*(rat-1e0)/(rat-.5e0)) + .5d0)
      nf1 = rat*ms
      if (2*nspread.gt.nf1) then
         nf1 = next235(DBLE(2e0*nspread))
      endif

C	ffts supports only radix2
	if(ffte .eq. E_FFTS) then
		nf1=find_best_pow2f(nf1)
	endif

c
      r2lamb = rat*rat * nspread / (rat*(rat-.5e0))
      hx = 2*pi/nf1
c
c     -----------------------------------
c     Compute workspace size and allocate
c     -----------------------------------
      iw1 = 2*nf1
      iwsav = iw1 + nspread+1
C	iwsav should be 16 bytes aligned
      if(mod(iwsav,16).ne.0) then
		iwsav = iwsav + 16 - mod(iwsav,16)
      endif
	if(ffte .eq. E_FFTS) then
C	ffts has its own workspace, but it needs input/output buffer
		iwtot = iwsav+2*nf1+15
	else
C	ffte needs this extra space
      iwtot = iwsav+4*nf1+15
	endif

      allocate ( fw(0:iwtot))
c
c     ---------------------------------------------------------------
c     Precompute spreading constants and initialize fw
c     to hold one term needed for fast Gaussian gridding
c     ---------------------------------------------------------------
      t1 = pi/r2lamb
      do k1 = 1, nspread
         fw(iw1+k1) = exp(-t1*k1**2)
      enddo

C Thomas init : CALL ZFFT1F(fw(0), nf1, 0, fw(iwsav))
C Thomas backward/inverse FFT : dcfftb(nf1,fw(0),fw(iwsav)) ==> CALL ZFFT1F(fw(0), nf1, 1, fw(iwsav))
C Thomas forward FFT : dcfftf(nf1,fw(0),fw(iwsav)) ==> CALL ZFFT1F(fw(0), nf1, -1, fw(iwsav))

	if(ffte .eq. E_FFTS) then
C	FFTS
C	iflag, < 0, forward, >=0 back
		fftp = ffts_init_1df(nf1, iflag)
	else if(ffte .eq. E_FFTE) then
C	FFTE
		CALL		ZFFT1F(fw(0), nf1, 0, fw(iwsav))
	else
C	netlib fftpack
		call dcffti(nf1,fw(iwsav))
	endif
	WRITE(6,*) "NUFFT type2: init fftp=",fftp, nf1,iwsav,iflag

c
c     ---------------------------------------------------------------
c     Deconvolve and compute inverse 1D FFT
c     (A factor of (-1)**k is needed to shift phase.)
c     ---------------------------------------------------------------
c
      t1 = pi * r2lamb / real(nf1)**2
      cross1 = 1e0/sqrt(r2lamb)
      zz = cross1*fk(0)
      fw(0) = real(zz)
      fw(1) = imag(zz)
      do k1 = 1, (ms-1)/2
         cross1 = -cross1
         cross = cross1*exp(t1*real(k1)**2)
         zz = cross*fk(k1)
         fw(2*k1) = real(zz)
         fw(2*k1+1) = imag(zz)
         zz = cross*fk(-k1)
         fw(2*(nf1-k1)) = real(zz)
         fw(2*(nf1-k1)+1) = imag(zz)
      enddo
      cross = -cross1*exp(t1*real(ms/2)**2)
      if (ms/2*2.eq.ms) then
	 zz = cross*fk(-ms/2)
         fw(2*nf1-ms) = real(zz)
         fw(2*nf1-ms+1) = imag(zz)
      endif
      do k1 = (ms+1)/2, nf1-ms/2-1
         fw(2*k1) = cmplx(0e0, 0e0)
         fw(2*k1+1) = cmplx(0e0, 0e0)
      enddo
c

C Thomas init : CALL ZFFT1F(fw(0), nf1, 0, fw(iwsav))
C Thomas backward/inverse FFT : dcfftb(nf1,fw(0),fw(iwsav)) ==> CALL ZFFT1F(fw(0), nf1, 1, fw(iwsav))
C Thomas forward FFT : dcfftf(nf1,fw(0),fw(iwsav)) ==> CALL ZFFT1F(fw(0), nf1, -1, fw(iwsav))

	if(ffte .eq. E_FFTS) then
C	output=[real, imaginary, real, imaginary ... interleaving format]
		call ffts_executef(fftp, fw(0), fw(iwsav))
C	copy the output buffer back to the input buffer
C	output : complex, complex.... => real,image,real,image....
C	memcpy(fw(0), fw(iwsav), 2*nf1)
		fw(0:2*nf1-1) = fw(iwsav:iwsav+2*nf1-1)
		call ffts_freef(fftp)
	else
		if (iflag .ge. 0) then
			if(ffte .eq. E_FFTE) then
				CALL		ZFFT1F(fw(0), nf1, 1, fw(iwsav))
			else
				call dcfftb(nf1,fw(0),fw(iwsav))
			endif
		else
			if(ffte .eq. E_FFTE) then
				CALL		ZFFT1F(fw(0), nf1, -1, fw(iwsav))
			else
				call 	dcfftf(nf1,fw(0),fw(iwsav))
			endif
		endif
	endif
	WRITE(6,*) "NUFFT type 2:",nf1,iwsav
C		CALL		DUMPF(fw(0), nf1)

c
c     ---------------------------------------------------------------
c     Loop over target points (1,...,nj)
c
c       1. find closest mesh point (with periodic wrapping if needed)
c       2. get contributions from regular fine grid to target
c          locations using Gaussian convolution.
c     ---------------------------------------------------------------
      t1 = pi/r2lamb
      do j = 1, nj
         cj(j) = cmplx(0e0,0e0)
         jb1 = int((xj(j)+pi)/hx)
         diff1 = (xj(j)+pi)/hx - jb1
         jb1 = mod(jb1, nf1)
         if (jb1.lt.0) jb1=jb1+nf1
         xc(0) = exp(-t1*diff1**2)
         cross = xc(0)
         cross1 = exp(2e0*t1 * diff1)
         do k1 = 1, nspread
            cross = cross * cross1
            xc(k1) = fw(iw1+k1)*cross
         enddo
         cross = xc(0)
         cross1 = 1e0/cross1
         do k1 = 1, nspread-1
            cross = cross * cross1
            xc(-k1) = fw(iw1+k1)*cross
         enddo
c
         jb1d = min(nspread-1, jb1)
         jb1u = min(nspread, nf1-jb1-1)
         do k1 = -nspread+1, -jb1d-1
	    zz = cmplx(fw(2*(jb1+k1+nf1)),fw(2*(jb1+k1+nf1)+1))
            cj(j) = cj(j) + xc(k1)*zz
         enddo
         do k1 = -jb1d, jb1u
	    zz = cmplx(fw(2*(jb1+k1)),fw(2*(jb1+k1)+1))
            cj(j) = cj(j) + xc(k1)*zz
         enddo
         do k1 = jb1u+1, nspread
	    zz = cmplx(fw(2*(jb1+k1-nf1)),fw(2*(jb1+k1-nf1)+1))
            cj(j) = cj(j) + xc(k1)*zz
         enddo
      enddo
      deallocate(fw)
      return
      end
c
c
c
c
c
c
************************************************************************
      subroutine nufft1d3ff90_ffte(nj,xj,cj, iflag,eps, nk,sk,fk,ier)
C		use moudle to load global vars
      USE NUFFTModule
C	TYPE(C_PTR), i.e., "void *", variable needs "ISO_C_BINDING" to define the TYPE(C_PTR)
      USE , intrinsic ::ISO_C_BINDING
      implicit none

		interface
C		fortran always uses call-by-address, so it passes pointers, not value to c functions.
C		ffts_plan_t *ffts_init_1d(size_t n, int sign);
C 		size_t for 64 bits platform is 8 bytes, fortran uses "call-by-addr",
C		so wrong size passing corrupts the memory
		function ffts_init_1df ( n, sign ) result(ptr) bind ( c )
			use, intrinsic :: iso_c_binding
			implicit none
			TYPE(C_PTR) 	:: ptr
			integer ( c_long ) :: n
			integer ( c_int ) :: sign
		end function ffts_init_1df

C		void ffts_execute (ffts_plan_t *plan , const void *input, void *output )
		subroutine ffts_executef (plan , input, output ) bind ( c )
			use, intrinsic :: iso_c_binding
			implicit none
			TYPE(C_PTR) :: plan
			REAL(C_FLOAT) :: input(*)
			REAL(C_FLOAT) :: output(*)
		end subroutine ffts_executef

C		void ffts_free(ffts_plan_t *plan);
		subroutine ffts_freef ( plan ) bind ( c )
			use, intrinsic :: iso_c_binding
			implicit none
			TYPE(C_PTR) :: plan
		end subroutine ffts_freef

C		unsigned find_best_pow2f(size_t *n);the nearest power of 2 next to n
		function find_best_pow2f ( n ) bind ( c )
			use, intrinsic :: iso_c_binding
			implicit none
			integer ( c_long ) :: n
			integer ( c_int ) :: find_best_pow2f
		end function find_best_pow2f
		end interface
C	special vars passed to c functions
C	https://people.sc.fsu.edu/~jburkardt/f_src/f90_calls_c/f90_calls_c.html
		integer (c_int) :: iflag
		integer (c_long) :: nf1

      integer ier,iw1,iwsave,iwtot,j,jb1,k1,kb1,kmax,nj,nk
      integer next235,nspread
      real*4 ang,cross,cross1,diff1,eps,hx,hs,rat,pi,r2lamb1
      real*4 sm,sb,t1,t2,xm,xb,max_t,nf1_t
      real*4 xc(-147:147), xj(nj), sk(nk)
      parameter (pi=3.14159265358979323846e0)
      complex*8 cj(nj), fk(nk), zz, cs
C 		void *
      TYPE(C_PTR) :: fftp

c ----------------------------------------------------------------------
      integer nw, istart
      real*4, allocatable :: fw(:)
!dir$ attributes align:64
c ----------------------------------------------------------------------
c     if (iflag .ge. 0) then
c
c                  nj
c     fk(sk(k)) = SUM cj(j) exp(+i sk(k) xj(j))  for k = 1,..., nk
c                 j=1
c
c     else
c
c                  nj
c     fk(sk(k)) = SUM cj(j) exp(-i sk(k) xj(j))  for k = 1,..., nk
c                 j=1
c ----------------------------------------------------------------------
c     INPUT:
c
c     nj     number of sources   (integer)
c     xj     location of sources (double array)
c     cj     strengths of sources (double complex array)
c     iflag  determines sign of FFT (see above)
c     eps    precision request  (between 1.0d-33 and 1.0d-1)
c               recomended value is 1d-15 for double precision calculations
c     nk     number of (noninteger) Fourier modes computed
c     sk     k-values (locations) of desired Fourier modes
c
c     OUTPUT:
c
c     fk     Fourier transform values (double complex array)
c     ier    error return code
c
c            ier = 0  => normal execution.
c            ier = 1  => precision eps requested is out of range.
c
c
c     References:
c
c     [DR] Fast Fourier transforms for nonequispaced data,
c          A. Dutt and V. Rokhlin, SIAM J. Sci. Comput. 14,
c          1368-1383, 1993.
c
c     [LG] The type 3 nonuniform FFT and its applications
c          J.-Y. Lee and L. Greengard, J. Comput. Phys. 206, 1-5 (2005).
c
c     The algorithm is essentially a concatenation of the
c     type 1 and 2 transforms.
c
c     1) Gaussian gridding of strengths cj(j) centered at xj
c        to create f_\tau(n \Delta_x)  (in notation of [LG])
c     2) Deconvolve each regular grid mode
c        to create f^{-\sigma}_\tau(n \Delta_x)  (in notation of [LG])
c     3) compute FFT on uniform mesh
c        to create F^{-\sigma}_\tau(m \Delta_s)  (in notation of [LG])
c     4) Gaussian gridding to irregular frequency points
c        to create F_\tau(s_k)  (in notation of [LG])
c     5) Deconvolution of  result
c        to create F(s_k)  (in notation of [LG])
c
c***********************************************************************
      ier = 0
      if ((eps.lt.1d-33).or.(eps.gt.1d-1)) then
         ier = 1
         return
      endif
c
c --- Check the ranges of {xj}/{sk} and the workspace size
c
      t1 = xj(1)
      t2 = xj(1)
      do j = 2, nj
         if (xj(j).gt.t2) then
             t2=xj(j)
         else if (xj(j).lt.t1) then
             t1=xj(j)
         endif
      enddo
      xb = (t1+t2) / 2e0
      xm = max(t2-xb,-t1+xb)  ! max(abs(t2-xb),abs(t1-xb))
c
      t1 = sk(1)
      t2 = sk(1)
      do k1 = 2, nk
         if (sk(k1).gt.t2) then
             t2=sk(k1)
         else if (sk(k1).lt.t1) then
             t1=sk(k1)
         endif
      enddo
      sb = (t1+t2) / 2e0
      sm = max(t2-sb,-t1+sb)
c
c     -------------------------------
c     Precision dependent parameters
c
c     rat is oversampling parameter
c     nspread is number of neighbors to which Gaussian gridding is
c     carried out.
c     -------------------------------
      if (eps.le.1d-11) then
         rat = sqrt(3.0e0)
      else
         rat = sqrt(2.0e0)
      endif
c
      nspread = int(-log(eps)/(pi*(rat-1e0)/(rat-.5e0)) + .5e0)
      t1 = 2e0/pi * xm*sm

C      WRITE(6,*) "rat=",rat, "t1=", t1, "nspread=", nspread
C      WRITE(6,*) "max=", max(rat*t1+2*nspread,2*nspread/(rat-1))

      nf1 = next235(DBLE(rat*max(rat*t1+2*nspread,2*nspread/(rat-1))))
C	ffts supports only radix2
	if(ffte .eq. E_FFTS) then
		nf1=find_best_pow2f(nf1)
	endif
C      WRITE(6,*) "nf1=", nf1

      rat = (sqrt(nf1*t1+nspread**2)-nspread)/t1
c
      r2lamb1 = rat*rat * nspread / (rat*(rat-.5e0))
      hx = pi/(rat*sm)
      hs = 2e0*pi/real(nf1)/hx            ! hx hs = 2.pi/nf1
c
c     -------------------------------
c     Compute workspace size and allocate
c     -------------------------------
c
      kmax = int(nf1*(r2lamb1-nspread)/r2lamb1+.1e0)
      iw1 = 2*nf1
      iwsave = iw1 + nspread+1
C	iwsav should be 16 bytes aligned
      if(mod(iwsave,16).ne.0) then
		iwsave = iwsave + 16 - mod(iwsave,16)
      endif
	if(ffte .eq. E_FFTS) then
C	ffts has its own workspace, but it needs input/output buffer
		iwtot = iwsave + 16+2*nf1
	else
C	ffte needs this extra space
        iwtot = iwsave + 16+4*nf1
	endif

C      allocate ( fw(0:iwtot-1) )
      allocate ( fw(0:iwtot) )
c
c     ---------------------------------------------------------------
c     Precompute spreading constants and initialize fw
c     to hold one term needed for fast Gaussian gridding
c     ---------------------------------------------------------------
c
      t1 = pi/r2lamb1
      do k1 = 1, nspread
         fw(iw1+k1) = exp(-t1*k1**2)
      enddo
c

C Thomas init : CALL ZFFT1F(fw(0), nf1, 0, fw(iwsave))
C Thomas backward/inverse FFT : dcfftb(nf1,fw(0),fw(iwsave)) ==> CALL ZFFT1F(fw(0), nf1, 1, fw(iwsave))
C Thomas forward FFT : dcfftf(nf1,fw(0),fw(iwsave)) ==> CALL ZFFT1F(fw(0), nf1, -1, fw(iwsave))

	if(ffte .eq. E_FFTS) then
C	FFTS
C	iflag, < 0, forward, >=0 back
		fftp = ffts_init_1df(nf1, iflag)
	else if(ffte .eq. E_FFTE) then
C	FFTE
		CALL		ZFFT1F(fw(0), nf1, 0, fw(iwsave))
	else
C	netlib fftpack
		call dcffti(nf1,fw(iwsave))
	endif
	WRITE(6,*) "NUFFT type3: init fftp=",fftp, nf1,iwsave,iflag
C		CALL		DUMPF(fw(0), nf1)

c
c     ---------------------------------------------------------------
c     Initialize fine grid data to zero.
c     ---------------------------------------------------------------
      do k1 = 0, 2*nf1-1
         fw(k1) = cmplx(0e0,0e0)
      enddo
c
c     ---------------------------------------------------------------
c     Step 1/5  - gridding as in type 1 transform.
c     ---------------------------------------------------------------
c
      t1 = pi/r2lamb1
      if (iflag .lt. 0) sb = -sb
      do j = 1, nj
         jb1 = int(real(nf1/2) + (xj(j)-xb)/hx)
         diff1 = real(nf1/2) + (xj(j)-xb)/hx - jb1
         ang = sb*xj(j)
         cs = cmplx(cos(ang),sin(ang)) * cj(j)

         xc(0) = exp(-t1*diff1**2)
         cross = xc(0)
         cross1 = exp(2e0*t1 * diff1)
         do k1 = 1, nspread
            cross = cross * cross1
            xc(k1) = fw(iw1+k1)*cross
         enddo
         cross = xc(0)
         cross1 = 1e0/cross1
         do k1 = 1, nspread-1
            cross = cross * cross1
            xc(-k1) = fw(iw1+k1)*cross
         enddo

         do k1 = -nspread+1, nspread
	    istart = 2*(jb1+k1)
	    zz = xc(k1)*cs
            fw(istart) = fw(istart) + real(zz)
            fw(istart+1) = fw(istart+1) + imag(zz)
         enddo
      enddo
      if (iflag .lt. 0) sb = -sb
c
c ---------------------------------------------------------------
c     Step 2: Deconvolve (amplify) as in Type 2 transform
c     Step 3: Compute FFT with shift
c             (-1)^k F_(k+M/2) = Sum (-1)^j F_(j+M/2) e(2pi ijk/M)
c ---------------------------------------------------------------
c
      t1 = pi * r2lamb1 / real(nf1)**2
      cross1 = (1e0-2e0*mod(nf1/2,2))/r2lamb1
      zz = cmplx(fw(nf1),fw(nf1+1))
      zz = cross1*zz
      fw(nf1) = real(zz)
      fw(nf1+1) = imag(zz)
      do k1 = 1, kmax
         cross1 = -cross1
         cross = cross1*exp(t1*real(k1)**2)
         zz = cmplx(fw(nf1-2*k1),fw(nf1-2*k1+1))
         zz = cross*zz
         fw(nf1-2*k1) = real(zz)
         fw(nf1-2*k1+1) = imag(zz)
         zz = cmplx(fw(nf1+2*k1),fw(nf1+2*k1+1))
         zz = cross*zz
         fw(nf1+2*k1) = real(zz)
         fw(nf1+2*k1+1) = imag(zz)
      enddo
c

C Thomas init : CALL ZFFT1F(fw(0), nf1, 0, fw(iwsave))
C Thomas backward/inverse FFT : dcfftb(nf1,fw(0),fw(iwsave)) ==> CALL ZFFT1F(fw(0), nf1, 1, fw(iwsave))
C Thomas forward FFT : dcfftf(nf1,fw(0),fw(iwsave)) ==> CALL ZFFT1F(fw(0), nf1, -1, fw(iwsave))

	if(ffte .eq. E_FFTS) then
C	output=[real, imaginary, real, imaginary ... interleaving format]
		call ffts_executef(fftp, fw(0), fw(iwsave))
C	copy the output buffer back to the input buffer
C	output : complex, complex.... => real,image,real,image....
C	memcpy(fw(0), fw(iwsave), 2*nf1)
		fw(0:2*nf1-1) = fw(iwsave:iwsave+2*nf1-1)
		call ffts_freef(fftp)
	else
		if (iflag .ge. 0) then
			if(ffte .eq. E_FFTE) then
				CALL		ZFFT1F(fw(0), nf1, 1, fw(iwsave))
			else
				call dcfftb(nf1,fw(0),fw(iwsave))
			endif
		else
			if(ffte .eq. E_FFTE) then
				CALL		ZFFT1F(fw(0), nf1, -1, fw(iwsave))
			else
				call dcfftf(nf1,fw(0),fw(iwsave))
			endif
		endif
	endif
	WRITE(6,*) "NUFFT type3:",nf1,iwsave

      do k1 = 1, kmax+nspread, 2
         fw(nf1+2*k1) = -fw(nf1+2*k1)
         fw(nf1+2*k1+1) = -fw(nf1+2*k1+1)
         fw(nf1-2*k1) = -fw(nf1-2*k1)
         fw(nf1-2*k1+1) = -fw(nf1-2*k1+1)
      enddo
c
c     ---------------------------------------------------------------
c     Step 4 Gaussian gridding to irregular points
c     Step 5 Final deconvolution
c     ---------------------------------------------------------------
      t1 = pi/r2lamb1
      do j = 1, nk
         kb1 = int(real(nf1/2) + (sk(j)-sb)/hs)
         diff1 = real(nf1/2) + (sk(j)-sb)/hs - kb1

         ! exp(-t1*(diff1-k1)**2) = xc(k1)
         xc(0) = exp(-t1*diff1**2)
         cross = xc(0)
         cross1 = exp(2e0*t1 * diff1)
         do k1 = 1, nspread
            cross = cross * cross1
            xc(k1) = fw(iw1+k1)*cross
         enddo
         cross = xc(0)
         cross1 = 1e0/cross1
         do k1 = 1, nspread-1
            cross = cross * cross1
            xc(-k1) = fw(iw1+k1)*cross
         enddo
c
         fk(j) = cmplx(0e0,0e0)
         do k1 = -nspread+1, nspread
	    zz = cmplx(fw(2*(kb1+k1)),fw(2*(kb1+k1)+1))
            fk(j) = fk(j) + xc(k1)*zz
         enddo
      enddo
c
      if (iflag .lt. 0) xb = -xb
      t1 = r2lamb1/(4e0*pi) * hx**2
      do j = 1, nk
         fk(j) = (exp(t1*(sk(j)-sb)**2))*fk(j)
         ang = (sk(j)-sb)*xb
         fk(j) = cmplx(cos(ang),sin(ang)) * fk(j)
      enddo
      deallocate(fw)
      return
      end
C dump array
      SUBROUTINE DUMPF(A,N)
      IMPLICIT REAL*4 (A-H,O-Z)
      COMPLEX*8 A(*)
C
      DO 10 I=1,N
C normalize the output by dividing N
C multiplies N back
		A(I) = DCMPLX(REAL(A(I)), IMAG(A(I))) * N
        WRITE(6,*) I,A(I)
   10 CONTINUE
      RETURN
      END