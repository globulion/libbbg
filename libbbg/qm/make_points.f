      subroutine make_points(nx,ny,dx,dy,xmin,ymin,points)
C
C     GENERATE THE POINT GRID FOR USE FOR CLEMTP.POTDMA SUBROUTINE
C
C
C                             Bartosz Błasiak, 12 Jan 2015
C
      double precision dx, dy, xmin, ymin
      double precision points(nx*ny*8)
      integer nx, ny, i, j, n, ijx, ijy, ijz
      data zero/0.0D+00/
Cf2py intent(in,out) points
      n = nx*ny
      do i=1, nx
      do j=1, ny
         ijx = ny*8*(i-1) + 8*(j-1) + 1
         ijy = ijx + 1
         ijz = ijy + 1
         points(ijx) = xmin + (i-1)*dx
         points(ijy) = ymin + (j-1)*dy
         points(ijz) = zero
      end do
      end do
      return
      end
