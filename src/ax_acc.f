c-----------------------------------------------------------------------
      subroutine ax_acc(w,u,gxyz,ur,us,ut,wk,n) ! Matrix-vector product: w=A*u

#ifdef TUNED_CUF_KERNEL
      use cudafor
#endif

      include 'SIZE'
      include 'TOTAL'

#ifdef TUNED_CUF_KERNEL
      interface
      attributes(global) subroutine ax_cuf2(w,u,ur,us,ut,
     &                gxyz,dxm1,dxtm1)

      real, intent(out) :: w(nx1,ny1,nz1,nelt)
      real, intent(in)  :: u(nx1,ny1,nz1,nelt)
      real ur  (nx1,ny1,nz1,lelt)
      real us  (nx1,ny1,nz1,lelt)
      real ut  (nx1,ny1,nz1,lelt)

      real gxyz(nx1,ny1,nz1,2*ldim,lelt)

      real, intent(in) :: dxm1(nx1,nx1)
      real, intent(in) :: dxtm1(nx1,nx1)
      end subroutine
      end interface
#endif

      common /mymask/cmask(-1:lx1*ly1*lz1*lelt)

      real w(nx1,ny1,nz1,nelt)
      real u(nx1,ny1,nz1,nelt)
      real gxyz(nx1,ny1,nz1,2*ldim,lelt)

      real ur(nx1,ny1,nz1,lelt)
      real us(nx1,ny1,nz1,lelt)
      real ut(nx1,ny1,nz1,lelt)
      real wk(nx1,ny1,nz1,lelt)

      real wr,ws,wt,tmp
      integer i,j,k,l,e,n

      integer lt
      
      integer cuda_err

      lt = nx1*ny1*nz1*nelt

!$ACC DATA PRESENT(w,u(:,:,:,:),gxyz,ur,us,ut,wk,dxm1,dxtm1)

#ifdef TUNED_CUF_KERNEL

!$acc host_data use_device(w,u(:,:,:,:),ur,us,ut,gxyz,dxm1,dxtm1)
       if (nx1.eq.10) then
         call ax_cuf2<<<nelt,dim3(nx1,ny1,nz1)>>>(w,u,
     $                ur,us,ut,gxyz,dxm1,dxtm1)
       else if (nx1.eq.12) then
         call ax_cuf2<<<nelt,dim3(nx1,ny1,nz1/2)>>>(w,u,
     $                ur,us,ut,gxyz,dxm1,dxtm1)
       else
c         call ax_cuf2<<<nelt,dim3(nx1,ny1,nz1/4)>>>(w,u,
c     $                ur,us,ut,gxyz,dxm1,dxtm1)
         call ax_cuf2<<<nelt,dim3(nx1,ny1,nz1/4)>>>(w,u,
     $         ur,us,ut,gxyz,dxm1,dxtm1) 
       endif
       
       cuda_err = cudaGetLastError()
       if (cuda_err /= cudaSuccess) then
         write(6, 815) cuda_err, cudaGetErrorString(cuda_err)
  815    format('CUDA ERROR', I3, ': ', A)
         call exitt
       endif

       istat = cudaDeviceSynchronize()
       
       cuda_err = cudaGetLastError()
       if (cuda_err /= cudaSuccess) then
         write(6, 815) cuda_err, cudaGetErrorString(cuda_err)
         call exitt
       endif

!$acc end host_data

#else
c ifndef TUNED_CUF_KERNEL
            
!$ACC PARALLEL LOOP COLLAPSE(4) GANG WORKER VECTOR PRIVATE(wr,ws,wt)
!DIR NOBLOCKING
      do e = 1,nelt
         do k=1,nz1
         do j=1,ny1
         do i=1,nx1
            wr = 0
            ws = 0
            wt = 0
!$ACC LOOP SEQ
            do l=1,nx1    ! serial loop, no reduction needed
               wr = wr + dxm1(i,l)*u(l,j,k,e)
               ws = ws + dxm1(j,l)*u(i,l,k,e)
               wt = wt + dxm1(k,l)*u(i,j,l,e)
            enddo
            ur(i,j,k,e) = gxyz(i,j,k,1,e)*wr
     $                  + gxyz(i,j,k,2,e)*ws
     $                  + gxyz(i,j,k,3,e)*wt
            us(i,j,k,e) = gxyz(i,j,k,2,e)*wr
     $                  + gxyz(i,j,k,4,e)*ws
     $                  + gxyz(i,j,k,5,e)*wt
            ut(i,j,k,e) = gxyz(i,j,k,3,e)*wr
     $                  + gxyz(i,j,k,5,e)*ws
     $                  + gxyz(i,j,k,6,e)*wt
         enddo
         enddo
         enddo
      enddo
!$ACC END PARALLEL LOOP

!$ACC PARALLEL LOOP COLLAPSE(4) GANG WORKER VECTOR 
      do e=1,nelt
         do k=1,nz1
         do j=1,ny1
         do i=1,nx1
            w(i,j,k,e) = 0.0
!$ACC LOOP SEQ
            do l=1,nx1    ! serial loop, no reduction needed
               w(i,j,k,e) = w(i,j,k,e) + dxtm1(i,l)*ur(l,j,k,e)
     $                                 + dxtm1(j,l)*us(i,l,k,e)
     $                                 + dxtm1(k,l)*ut(i,j,l,e)
            enddo
         enddo
         enddo
         enddo
      enddo
!$ACC END PARALLEL LOOP

#endif
c endif TUNED_CUF_KERNEL

#ifdef GPUDIRECT
      call dssum(w)         ! Gather-scatter operation  ! w   = QQ  w
#else
      call dssum_acc(w)         ! Gather-scatter operation  ! w   = QQ  w
#endif

      call add2s2_acc(w,u,.1,n)   !2n
      call maskit_acc(w,cmask,nx1,ny1,nz1)  ! Zero out Dirichlet conditions

!$ACC END DATA

      nxyz=nx1*ny1*nz1
      flop_a = flop_a + (19*nxyz+12*nx1*nxyz)*nelt

      return
      end
c-----------------------------------------------------------------------
