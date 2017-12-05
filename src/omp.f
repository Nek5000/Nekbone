#ifdef TIMERS
#define NBTIMER(a) a = dnekclock()
#define STIMER(a) a = dnekclock_sync()
#define ACCUMTIMER(b,a) b = b + (dnekclock()- a)
#else
#define NBTIMER(a)
#define STIMER(a)
#define ACCUMTIMER(a,b)
#endif


      subroutine rzeroi(a,n,start,fin)
        implicit none
  
        real a(n)
        integer n, i, start, fin

        do i = start, fin
          a(i) = 0.0
        end do 

        return
      end subroutine

c----------------------------------------------------------

      subroutine copyi(a,b,n, start, fin)
        implicit none

        real a(n),b(n)
        integer n, i, start, fin

        do i=start,fin
          a(i)=b(i)
        enddo

        return
      end subroutine

c----------------------------------------------------------

      subroutine glsc3i(val,a,b,mult,n,find,lind)
      implicit none

      include 'TIMER'

      real val,a(n),b(n),mult(n)
      real tsum,psum,work(1)
      integer n,find,lind
      integer i, tmt, thread
      integer omp_get_thread_num
    
      save psum
      data psum /0.0/

      thread = 0
#ifdef _OPENMP
      thread = omp_get_thread_num()
#endif
      tmt = thread + 1

      tsum = 0.0
      do i=find, lind
         tsum = tsum + a(i)*b(i)*mult(i)
      end do

c$OMP ATOMIC update
      psum = psum + tsum
c$OMP END ATOMIC

c$OMP BARRIER
      NBTIMER(ttemp4)
c$OMP MASTER
      call gop(psum,work,'+  ',1)
      val = psum
      psum = 0.0
c$OMP END MASTER
c$OMP BARRIER
      ACCUMTIMER(tgop(gopi(tmt),tmt), ttemp4)


      return
      end subroutine

c----------------------------------------------------------

      subroutine solveMi(z,r,n,start,fin)
      implicit none

      real z(n),r(n)
      integer n,start,fin

      call copyi(z,r,n,start,fin) 

      return
      end

c----------------------------------------------------------

      subroutine add2s1i(a,b,c1,n,start,fin)
      implicit none

      real a(n),b(n),c1
      integer n,start,fin
      integer i

      do i= start, fin
        a(i)=c1*a(i)+b(i)
      end do

      return
      end subroutine

c----------------------------------------------------------

      subroutine add2s2i(a,b,c1,n,start,fin)
      implicit none
 
      real a(n),b(n),c1
      integer n,start,fin
      integer i

      do i= start,fin
        a(i)=a(i)+c1*b(i)
      end do

      return
      end subroutine
