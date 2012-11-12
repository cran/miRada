c  Part of R package KernSmooth
c  Copyright (C) 1995  M. P. Wand
c
c  Unlimited use and distribution (see LICENCE).

cccccccccc FORTRAN subroutine linbin.f cccccccccc

c Obtains bin counts for univariate data
c via the linear binning strategy. If "trun=0" then
c weight from end observations is given to corresponding
c end grid points. If "trun=1" then end observations
c are truncated.

c Last changed: 20 MAR 2009

      subroutine gaussbin(X,sigma,n,a,b,M,trun,gcnts)
      double precision X(*),sigma, a,b,gcnts(*),lxi,delta,rem
      integer n,M,i,li,trun

c     Initialize grid counts to zero

      do 10 i=1,M
         gcnts(i) = dble(0)
10    continue

      delta = (b-a)/(M-1)
      do 20 i=1,n
         lxi = ((X(i)-a)/delta) + 1

c        Find integer part of "lxi"

         li = int(lxi) 

         rem = lxi - li
         if (li.ge.1.and.li.lt.M) then
C            gcnts(li) = gcnts(li) + (1-rem)
C            gcnts(li+1) = gcnts(li+1) + rem
            call redistribute(X(i), sigma, a, delta, M, gcnts)
         endif

         if (li.lt.1.and.trun.eq.0) then
C            gcnts(1) = gcnts(1) + 1
            call redistribute(X(i), sigma, a, delta, M, gcnts)
         endif

         if (li.ge.M.and.trun.eq.0) then
C            gcnts(M) = gcnts(M) + 1
            call redistribute(X(i), sigma, a, delta, M, gcnts)
         endif

20    continue

      return
      end

cccccccccc End of linbin.f cccccccccc

CCCCCCCCCCCCCCCCCCCCC  Subroutines to be called 

      SUBROUTINE redistribute(y, sigma, a, delta, M, gcnts)
      double precision ptotal, pb, pa, y, sigma, a, b0, delta, 
     *     gcnts(*), delta2, x0 
      integer i,M
      pa = 0.0
      pb = 0.0
      delta2 = 0.5 * delta
      x0 = a - delta2
      call pnorm(pa,x0,y,sigma)
      b0 = a + (M-1)*delta
      call pnorm(pb,b0+delta2,y,sigma)
      ptotal = pb - pa
      
      do 20 i=1,M
         x0 = x0 + delta
         call pnorm(pb,x0,y,sigma)
         gcnts(i) = gcnts(i) + (pb-pa)/ptotal
         pa = pb
 20   continue
      return
      end

      SUBROUTINE pnorm(Fx,z,mu,sigma)
      double precision Fx, z, x2, mu, sigma
      DATA PI/3.141592653589793238462643d0/
      Fx = (z-mu)/sigma
      x2 = Fx/sqrt(2.)
      Fx = 0.5 + 0.5*erf(x2)
      return
      end
