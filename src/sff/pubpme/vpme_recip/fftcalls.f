      subroutine get_fftdims(nfft1,nfft2,nfft3,
     $       nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork,
     $       sizfftab,sizffwrk)
      implicit none
      integer nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3,
     $       nfftable,nffwork,sizfftab,sizffwrk
      integer n,nfftmax

      nfftmax = max(nfft1,nfft2,nfft3)
      nfftdim1 = nfft1
      n = nfft1/2
      if ( nfft1 .eq. 2*n )nfftdim1 = nfft1+1
      nfftdim2 = nfft2
      n = nfft2/2
      if ( nfft2 .eq. 2*n )nfftdim2 = nfft2+1
      nfftdim3 = nfft3
      n = nfft3/2
      if ( nfft3 .eq. 2*n )nfftdim3 = nfft3+1
#ifdef SGIFFT
      nfftable = 2*(nfftdim1+nfftdim2+nfftdim3+50)
      nffwork = 0
      sizfftab = nfftable
      sizffwrk  = nffwork
#endif
#ifdef CRAYFFT
      nfftable = 2*(nfftdim1+nfftdim2+nfftdim3+50)
      nffwork = 4*nfftdim1*nfftdim2*nfftdim3
      sizfftab = nfftable
      sizffwrk  = nffwork
#endif
#ifdef PUBFFT
      nfftable = 4*nfftmax + 15
      nffwork = nfftmax
      sizfftab = 3*nfftable
      sizffwrk  = 2*nfftmax
#endif
      return
      end
c---------------------------------------------------------------
      subroutine fft_setup(array,fftable,ffwork,
     $      nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3,
     $      nfftable,nffwork)
      implicit none

      double precision array(*),fftable(*),ffwork(*)
      integer nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3
      integer nfftable,nffwork

      integer isign,inc1,inc2,inc3
      double precision scale

#ifdef SGIFFT
      call ZFFT3DI(nfft1,nfft2,nfft3,fftable)
#endif
#ifdef CRAYFFT
      write(6,*)'using cray fft code'
      isign = 0
      inc1 = 1
      inc2 = NFFTDIM1
      inc3 = NFFTDIM1*NFFTDIM2
      scale = 1.d0
      call CFFT3D(isign,nfft1,nfft2,nfft3,scale,array,
     $      inc1,inc2,inc3,array,inc1,inc2,inc3,fftable,nfftable,
     $      ffwork,nffwork)
#endif
#ifdef PUBFFT
      write(6,*)'using public domain fft code'
      call pubz3di(nfft1,nfft2,nfft3,fftable,nfftable)
#endif
      return
      end
c-----------------------------------------------------------
      subroutine fft_forward(array,fftable,ffwork,
     $      nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3,
     $      nfftable,nffwork)
      implicit none

      double precision array(*),fftable(*),ffwork(*)
      integer nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3

      integer isign,inc1,inc2,inc3
      double precision scale
      integer nfftable,nffwork

      isign = 1

#ifdef SGIFFT
      call ZFFT3D(isign,nfft1,nfft2,nfft3,array,
     $   nfftdim1,nfftdim2,fftable)
#endif
#ifdef CRAYFFT
      inc1 = 1
      inc2 = nfftdim1
      inc3 = nfftdim1*nfftdim2
      scale = 1.d0
      call CFFT3D(isign,nfft1,nfft2,nfft3,scale,array,
     $      inc1,inc2,inc3,array,inc1,inc2,inc3,fftable,nfftable,
     $      ffwork,nffwork)
#endif
#ifdef PUBFFT
      call pubz3d(isign,nfft1,nfft2,nfft3,array,
     $   nfftdim1,nfftdim2,fftable,nfftable,ffwork,nffwork)
#endif
      return
      end
c-----------------------------------------------------------
      subroutine fft_back(array,fftable,ffwork,
     $      nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3,
     $      nfftable,nffwork)
      implicit none

      double precision array(*),fftable(*),ffwork(*)
      integer nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3
      integer nfftable,nffwork

      integer isign,inc1,inc2,inc3
      double precision scale

      isign = -1

#ifdef SGIFFT
      call ZFFT3D(isign,nfft1,nfft2,nfft3,array,
     $   nfftdim1,nfftdim2,fftable)
#endif
#ifdef CRAYFFT
      inc1 = 1
      inc2 = nfftdim1
      inc3 = nfftdim1*nfftdim2
      scale = 1.d0
      call CFFT3D(isign,nfft1,nfft2,nfft3,scale,array,
     $      inc1,inc2,inc3,array,inc1,inc2,inc3,fftable,nfftable,
     $      ffwork,nffwork)
#endif
#ifdef PUBFFT
      call pubz3d(isign,nfft1,nfft2,nfft3,array,
     $   nfftdim1,nfftdim2,fftable,nfftable,ffwork,nffwork)
#endif
      return
      end
c-----------------------------------------------------------
