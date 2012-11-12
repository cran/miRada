c  Wang, B (2012/03/09)
c  Unlimited use and distribution (see LICENCE).

cccccccccc FORTRAN subroutine npqlevel.f cccccccccc

c This program is to find the probability or quantile level for a given
c value on a distribution that has been nonparametrically estimated.

c Last changed: 9 MAR 2012

      subroutine npqlevel(x,n,x0, y0,m)
      double precision x(*), x0(*),y0(*),y1,x1,p
      integer n,m,i,j,k

      y1 = 0
      x1 = x0(1)
      k = 1
      do 100 i=1,n
         do 80 j=k,m
            IF(x(i) .GT. x0(j)) THEN
               y1 = y0(j)
               x1 = x0(j)
            ELSE
               p = (x(i)-x1)/(x0(j)-x1)
               x(i) = p*y0(j)+(1-p)*y1
               k = j
               go to 100
            ENDIF
 80      continue
 100  continue

      return
      end

cccccccccc End of linbin.f cccccccccc
