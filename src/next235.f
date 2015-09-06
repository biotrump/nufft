************************************************************************
C engine of fft
	MODULE NUFFTModule
		implicit none
C constant
		INTEGER, PARAMETER :: E_FFTE = 0, E_FFTS = 1, E_FFTW = 2
		INTEGER :: ffte=E_FFTE
	END MODULE NUFFTModule
C setup fft engine type
C default is ffte
	subroutine nufft_ffte(id)
		USE NUFFTModule
		implicit none
		INTEGER id
		if ((id.le.E_FFTW).and.(id.ge.E_FFTE) ) then
			ffte=id
		else
			ffte=E_FFTE
		endif
	end

      function next235(base)
      implicit none
      integer next235, numdiv
      real*8 base
c ----------------------------------------------------------------------
c     integer function next235 returns a multiple of 2, 3, and 5
c
c     next235 = 2^p 3^q 5^r >= base  where p>=1, q>=0, r>=0
************************************************************************
      next235 = 2 * int(base/2d0+.9999d0)
      if (next235.le.0) next235 = 2

100   numdiv = next235
      do while (numdiv/2*2 .eq. numdiv)
         numdiv = numdiv /2
      enddo
      do while (numdiv/3*3 .eq. numdiv)
         numdiv = numdiv /3
      enddo
      do while (numdiv/5*5 .eq. numdiv)
         numdiv = numdiv /5
      enddo
      if (numdiv .eq. 1) return
      next235 = next235 + 2
      goto 100
      end

