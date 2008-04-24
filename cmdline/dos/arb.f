c     Copyright Gleb Vdovin, 1995
c
c     This Fortran program forms some arbitrary
c     intensity and phase filters in files in and pha, 
c     to use with LightPipes
c     This file may be distributed only with
c     full LighPipes source code
c-----------------------------------------------------
      open(1, file='in', form='formatted')
      open(2, file='pha', form='formatted')

      do i=1,256
c     ii corresponds to X coordinate
         ii=i-129

         do j=1,256
c     jj corresponds to Y coordinate
            jj=j-129

c     ai and af are  some arbitrary function 
c     of the coordinates ii and jj,
c     ai corresponds to intensity and af is the phase:

            ai=abs(sin(i/10.)*cos(j/5.))
            af=cos(i/10.)*sin(j/5.)

c     with this "if" a circular aperture is cutted:
c     writing intensity:
            if (ii**2+jj**2 .le. 64**2) then
            write(1,*) ai
            else
            write(1,*) 0.
            end if

c     writing phase:
            write(2,*) af
         end do

c     empty line as a row separator
c     to make Gnuplot happy:
        write(1,*)
        write(2,*)
       end do

c     all done:
       close(1)
       close(2)
       stop
       end
