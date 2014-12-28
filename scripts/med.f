      program med
c     The program to form a quadratic 
c     distribution of the refractive index;
c     nc is the grid  sampling
c     an0 is the axis value of refractive index
c     an1 is n1 in the formula for refractive index: 
c     n=sqrt(an0**2-an1*an0*r^2) where 
c     r is the radial coordinate
c     xmax is the grid size

      read*, nc, an0, an1, xmax
      dx=xmax/(nc-1.)
      n2=nc/2+1
      do i=1,nc
         ii=i-n2
         x=dx*ii
         do j=1,nc
            jj=j-n2
            y=jj*dx
            r2=x*x+y*y
            print*, sqrt(an0*an0-an1*an0*r2)
            end do
            print*
            end do
            stop
            end
