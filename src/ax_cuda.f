#ifdef _CUDA

      attributes(global) subroutine ax_cuf2(w,u,ur,us,ut,
     &                gxyz,dxm1,dxtm1)

      include 'SIZE'

      real, intent(out) :: w(lx1,ly1,lz1,lelt)
!      real, intent(in)  :: u(lx1,ly1,lz1,lelt)
      real u(lx1,ly1,lz1,lelt)
      real ur  (lx1,ly1,lz1,lelt)
      real us  (lx1,ly1,lz1,lelt)
      real ut  (lx1,ly1,lz1,lelt)

      real gxyz(lx1,ly1,lz1,2*ldim,lelt)

      real, intent(in) :: dxm1(lx1,lx1)
      real, intent(in) :: dxtm1(lx1,lx1)

      real rtmp,stmp,ttmp,wijke
      real, shared :: shdxm1(lx1,ly1)
      real, shared :: shdxtm1(lx1,ly1)
      integer l,e,i,j,k,kk

      e = blockIdx%x
      k = threadIdx%z
      j = threadIdx%y
      i = threadIdx%x

      if (k.eq.1) then
         shdxm1(i,j) = dxm1(i,j)
         shdxtm1(i,j) = dxtm1(i,j)
      end if
      call syncthreads()

      rtmp = 0.0
      stmp = 0.0
      ttmp = 0.0
      do l=1,lx1
        rtmp = rtmp + shdxm1(i,l)    * u(l,j,k+12,e)
        stmp = stmp + shdxm1(j,l)    * u(i,l,k+12,e)
        ttmp = ttmp + shdxm1(k+12,l) * u(i,j,l,e)
      enddo
      ur(i,j,k+12,e) = gxyz(i,j,k+12,1,e)*rtmp
     $               + gxyz(i,j,k+12,2,e)*stmp
     $               + gxyz(i,j,k+12,3,e)*ttmp
      us(i,j,k+12,e) = gxyz(i,j,k+12,2,e)*rtmp
     $               + gxyz(i,j,k+12,4,e)*stmp
     $               + gxyz(i,j,k+12,5,e)*ttmp
      ut(i,j,k+12,e) = gxyz(i,j,k+12,3,e)*rtmp
     $               + gxyz(i,j,k+12,5,e)*stmp
     $               + gxyz(i,j,k+12,6,e)*ttmp
      call syncthreads()


      rtmp = 0.0
      stmp = 0.0
      ttmp = 0.0
      do l=1,lx1
        rtmp = rtmp + shdxm1(i,l)   * u(l,j,k+8,e)
        stmp = stmp + shdxm1(j,l)   * u(i,l,k+8,e)
        ttmp = ttmp + shdxm1(k+8,l) * u(i,j,l,e)
      enddo
      ur(i,j,k+8,e) = gxyz(i,j,k+8,1,e)*rtmp
     $              + gxyz(i,j,k+8,2,e)*stmp
     $              + gxyz(i,j,k+8,3,e)*ttmp
      us(i,j,k+8,e) = gxyz(i,j,k+8,2,e)*rtmp
     $              + gxyz(i,j,k+8,4,e)*stmp
     $              + gxyz(i,j,k+8,5,e)*ttmp
      ut(i,j,k+8,e) = gxyz(i,j,k+8,3,e)*rtmp
     $              + gxyz(i,j,k+8,5,e)*stmp
     $              + gxyz(i,j,k+8,6,e)*ttmp
      call syncthreads()


      rtmp = 0.0
      stmp = 0.0
      ttmp = 0.0
      do l=1,lx1
        rtmp = rtmp + shdxm1(i,l)   * u(l,j,k+4,e)
        stmp = stmp + shdxm1(j,l)   * u(i,l,k+4,e)
        ttmp = ttmp + shdxm1(k+4,l) * u(i,j,l,e)
      enddo
      ur(i,j,k+4,e) = gxyz(i,j,k+4,1,e)*rtmp
     $              + gxyz(i,j,k+4,2,e)*stmp
     $              + gxyz(i,j,k+4,3,e)*ttmp
      us(i,j,k+4,e) = gxyz(i,j,k+4,2,e)*rtmp
     $              + gxyz(i,j,k+4,4,e)*stmp
     $              + gxyz(i,j,k+4,5,e)*ttmp
      ut(i,j,k+4,e) = gxyz(i,j,k+4,3,e)*rtmp
     $              + gxyz(i,j,k+4,5,e)*stmp
     $              + gxyz(i,j,k+4,6,e)*ttmp
      call syncthreads()


      rtmp = 0.0
      stmp = 0.0
      ttmp = 0.0
      do l=1,lx1
        rtmp = rtmp + shdxm1(i,l) * u(l,j,k,e)
        stmp = stmp + shdxm1(j,l) * u(i,l,k,e)
        ttmp = ttmp + shdxm1(k,l) * u(i,j,l,e)
      enddo

      ur(i,j,k,e) = gxyz(i,j,k,1,e)*rtmp
     $            + gxyz(i,j,k,2,e)*stmp
     $            + gxyz(i,j,k,3,e)*ttmp
      us(i,j,k,e) = gxyz(i,j,k,2,e)*rtmp
     $            + gxyz(i,j,k,4,e)*stmp
     $            + gxyz(i,j,k,5,e)*ttmp
      ut(i,j,k,e) = gxyz(i,j,k,3,e)*rtmp
     $            + gxyz(i,j,k,5,e)*stmp
     $            + gxyz(i,j,k,6,e)*ttmp

      call syncthreads()
      wijke = 0.0
      do l=1,lx1
        wijke = wijke + shdxtm1(i,l)*ur(l,j,k,e) 
     $                + shdxtm1(j,l)*us(i,l,k,e)
     $                + shdxtm1(k,l)*ut(i,j,l,e)
      enddo
      w(i,j,k,e) = wijke


      call syncthreads()
      wijke = 0.0
      do l=1,lx1
        wijke = wijke + shdxtm1(i,l)*ur(l,j,k+4,e) 
     $                + shdxtm1(j,l)*us(i,l,k+4,e)
     $                + shdxtm1(k+4,l)*ut(i,j,l,e)
      enddo
      w(i,j,k+4,e) = wijke

      call syncthreads()
      wijke = 0.0
      do l=1,lx1
        wijke = wijke + shdxtm1(i,l)*ur(l,j,k+8,e) 
     $                + shdxtm1(j,l)*us(i,l,k+8,e)
     $                + shdxtm1(k+8,l)*ut(i,j,l,e)
      enddo
      w(i,j,k+8,e) = wijke

      call syncthreads()
      wijke = 0.0
      do l=1,lx1
        wijke = wijke + shdxtm1(i,l)*ur(l,j,k+12,e) 
     $                + shdxtm1(j,l)*us(i,l,k+12,e)
     $                + shdxtm1(k+12,l)*ut(i,j,l,e)
      enddo
      w(i,j,k+12,e) = wijke

      return
      end

#else

      subroutine ax_cuf2(w,u,ur,us,ut,gxyz,dxm1,dxtm1)
        call err_chk(
     $ 'ERROR: Called ax_cuf2 but did not compile with CUDA')
      return
      end

#endif
