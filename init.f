      subroutine init(f,ux,uy,uz,t,dx,dy,dz)
      !include 'param.h'
      use param
      real*8 f(NX,NY,NZ)
      real*8 ux(NX,NY,NZ),uy(NX,NY,NZ),uz(NX,NY,NZ)
      real*8 t,dx,dy,dz

      real*8  ls(NX,NY,NZ)
      integer i,j,k,is,ii
      real*8  x, y, z, x0, y0, z0, R, tmp 

      t  = 0.0d0
      dx = 1.0/NX !0.01d0
      dy = 1.0/NY !0.01d0
      dz = 1.0/NZ !0.01d0 
      
      R  = 0.15
      x0 = 0.35
      y0 = 0.35
      z0 = 0.35

      do k=1,NZ
      do j=1,NY
      do i=1,NX

         f(i,j,k)  = 0.0d0

         x = (i-1)*dx
         y = (j-1)*dy
         z = (k-1)*dz

c ***
c initializing of level set function 
c ***
         ls(i,j,k) = R - sqrt( (x-x0)**2 + (y-y0)**2 + (z-z0)**2 ) 
      end do
      end do
      end do

      call levelset2vof(nx,ny,nz,ls,f)

      
c ***

      return
      end

