c ****** 1 ******* 2 ******* 3 ******* 4 ******* 5 ******* 6 ******* 7 *
c Split advection of the interface along the x (d=1), y (d=2) and z (d=3)
c direction
c ****** 1 ******* 2 ******* 3 ******* 4 ******* 5 ******* 6 ******* 7 *
      SUBROUTINE swp(us,c,d)
c***
      !include 'param.h'
      use param
      INTEGER i,j,k,invx,invy,invz,d
      DOUBLE PRECISION us(nx,ny,nz),c(nx,ny,nz),mx,my,mz,mm1,mm2
      DOUBLE PRECISION a1,a2,alpha,AL3D,FL3D
      DOUBLE PRECISION vof1(nx,ny,nz),vof2(nx,ny,nz),vof3(nx,ny,nz)
      INTRINSIC DMAX1,DMIN1
      EXTERNAL AL3D,FL3D
c***
      do k=2,nz-1
         do j=2,ny-1
            do i=2,nx-1
              a1 = us(i,j,k) *dt/delta
              if (d.eq.1) then
                a2 = us(i+1,j,k) *dt/delta
              elseif (d.eq.2) then
                a2 = us(i,j+1,k) *dt/delta
              elseif (d.eq.3) then
                a2 = us(i,j,k+1) *dt/delta
              endif
               
c***
c     3 cases: 1: DEFAULT (c=0. and fluxes=0.); 2: c=1.; 3:c>0.
c***
               vof1(i,j,k) = 0.0d0
               vof2(i,j,k) = 0.0d0
               vof3(i,j,k) = 0.0d0

               if (c(i,j,k) .EQ. 1.0d0) then
                  vof1(i,j,k) = DMAX1(-a1,0.d0)
                  vof2(i,j,k) = 1.d0 - DMAX1(a1,0.d0) + DMIN1(a2,0.d0)
                  vof3(i,j,k) = DMAX1(a2,0.d0)

               else if (c(i,j,k) .GT. 0.d0) then
c***
c     (1) normal vector: mx,my,mz; (2) mx,my,mz>0. and mx+my+mz = 1.;
c     (3) get alpha;               (4) back to original plane;
c     (5) lagrangian advection;    (6) get fluxes
c*(1)*
                  mm1 = c(i-1,j-1,k-1)+c(i-1,j-1,k+1)+c(i-1,j+1,k-1)
     %                 +c(i-1,j+1,k+1)+2.0d0*(c(i-1,j-1,k)+c(i-1,j+1,k)
     %                 +c(i-1,j,k-1)+c(i-1,j,k+1))+4.0d0*c(i-1,j,k)
                  mm2 = c(i+1,j-1,k-1)+c(i+1,j-1,k+1)+c(i+1,j+1,k-1)
     %                 +c(i+1,j+1,k+1)+2.0d0*(c(i+1,j-1,k)+c(i+1,j+1,k)
     %                 +c(i+1,j,k-1)+c(i+1,j,k+1))+4.0d0*c(i+1,j,k)
                  mx = mm1 - mm2
                  
                  mm1 = c(i-1,j-1,k-1)+c(i-1,j-1,k+1)+c(i+1,j-1,k-1)
     %                 +c(i+1,j-1,k+1)+2.0d0*(c(i-1,j-1,k)+c(i+1,j-1,k)
     %                 +c(i,j-1,k-1)+c(i,j-1,k+1))+4.0d0*c(i,j-1,k)
                  mm2 = c(i-1,j+1,k-1)+c(i-1,j+1,k+1)+c(i+1,j+1,k-1)
     %                 +c(i+1,j+1,k+1)+2.0d0*(c(i-1,j+1,k)+c(i+1,j+1,k)
     %                 +c(i,j+1,k-1)+c(i,j+1,k+1))+4.0d0*c(i,j+1,k)
                  my = mm1 - mm2
                  
                  mm1 = c(i-1,j-1,k-1)+c(i-1,j+1,k-1)+c(i+1,j-1,k-1)
     %                 +c(i+1,j+1,k-1)+2.0d0*(c(i-1,j,k-1)+c(i+1,j,k-1)
     %                 +c(i,j-1,k-1)+c(i,j+1,k-1))+4.0d0*c(i,j,k-1)
                  mm2 = c(i-1,j-1,k+1)+c(i-1,j+1,k+1)+c(i+1,j-1,k+1)
     %                 +c(i+1,j+1,k+1)+2.0d0*(c(i-1,j,k+1)+c(i+1,j,k+1)
     %                 +c(i,j-1,k+1)+c(i,j+1,k+1))+4.0d0*c(i,j,k+1)
                  mz = mm1 - mm2
c*(2)*  
                  invx = 1
                  invy = 1
                  invz = 1
                  if (mx .LT. 0.0d0) then
                     mx = -mx
                     invx = -1
                  endif
                  if (my .LT. 0.0d0) then
                     my = -my
                     invy = -1
                  endif
                  if (mz .LT. 0.0d0) then
                     mz = -mz
                     invz = -1
                  endif
                  mm2 = mx+my+mz
                  mx = mx/mm2
                  my = my/mm2
                  mz = mz/mm2
c*(3)*  
                  alpha = AL3D(mx,my,mz,c(i,j,k))
c*(4)*  
                  mx = invx*mx
                  my = invy*my
                  mz = invz*mz
                  alpha = alpha + DMIN1(0.d0,mx) + DMIN1(0.d0,my) +
     %                 DMIN1(0.d0,mz)
c*(5)*  

                  mm1 = DMAX1(a1,0.0d0)
                  mm2 = 1.d0 - mm1 + DMIN1(0.d0,a2)
                  
                  if (d.eq.1) then
                    mx = mx/(1.0d0 - a1 + a2)
                    alpha = alpha + mx*a1
                    if (a1 .LT. 0.d0) 
     %                 vof1(i,j,k) = FL3D(mx,my,mz,alpha,a1  ,-a1)
                    if (a2 .GT. 0.d0) 
     %                 vof3(i,j,k) = FL3D(mx,my,mz,alpha,1.d0,a2)
                       vof2(i,j,k) = FL3D(mx,my,mz,alpha,mm1,mm2)
                  elseif (d.eq.2) then
                    my = my/(1.0d0 - a1 + a2)
                    alpha = alpha + my*a1
                    if (a1 .LT. 0.d0) 
     %                 vof1(i,j,k) = FL3D(my,mz,mx,alpha,a1  ,-a1)
                    if (a2 .GT. 0.d0) 
     %                 vof3(i,j,k) = FL3D(my,mz,mx,alpha,1.d0,a2)
                       vof2(i,j,k) = FL3D(my,mz,mx,alpha,mm1,mm2)
                  elseif (d.eq.3) then
                    mz = mz/(1.0d0 - a1 + a2)
                    alpha = alpha + mz*a1
                    if (a1 .LT. 0.d0) 
     %                 vof1(i,j,k) = FL3D(mz,mx,my,alpha,a1  ,-a1)
                    if (a2 .GT. 0.d0) 
     %                 vof3(i,j,k) = FL3D(mz,mx,my,alpha,1.d0,a2)
                       vof2(i,j,k) = FL3D(mz,mx,my,alpha,mm1,mm2)
                  endif
               endif
            enddo
         enddo
      enddo
c***  
c     (1) apply proper boundary conditions to fluxes
c     (2) new values of c and  clip it: 0. <= c <= 1.
c     (3) apply proper boundary conditions to c
c*(1)*  
      call bc_flux(vof1,vof3,us,nx,ny,nz,d)
c*(2)*  
      do k=2,nz-1
         do j=2,ny-1
            do i=2,nx-1
              if (d.eq.1) then
                c(i,j,k) = vof3(i-1,j,k) + vof2(i,j,k) + vof1(i+1,j,k)
              elseif (d.eq.2) then
                c(i,j,k) = vof3(i,j-1,k) + vof2(i,j,k) + vof1(i,j+1,k)
              elseif (d.eq.3) then
                c(i,j,k) = vof3(i,j,k-1) + vof2(i,j,k) + vof1(i,j,k+1)
              endif
              c(i,j,k) = DMAX1(0.0d0,DMIN1(1.0d0,c(i,j,k)))
           enddo
         enddo
      enddo
c*(3)*
      call bc_c(c,nx,ny,nz)
c***
      return
      end





c***********************************************************************
      SUBROUTINE bc_c(scal,nx,ny,nz)
c***
 
      INTEGER nx,ny,nz, i, j, k
      DOUBLE PRECISION scal(nx,ny,nz),y0,z0,zz,yy,r,h
      DOUBLE PRECISION GP_DFETCH
      EXTERNAL GP_DFETCH
      h=1.d0/(nx-2)
c***  
c     periodicity in the x direction
c***  
!      IF(NEIGHB(1).LT.0) THEN
c      write(*,*)'BCC2 PPRANK',PPRANK,NEIGHB(1),NEIGHB(2)
      do k=2,nz-1
         do j=2,ny-1
            scal(1 ,j,k) = 0 !scal(2,j,k)
         enddo 
      enddo
!      ENDIF
c***  
!      IF(NEIGHB(2).LT.0) THEN
      do k=2,nz-1
         do j=2,ny-1
            scal(nx,j,k) = 0 !scal(nx-1,j,k)
         enddo 
      enddo
!      ENDIF
c***  
c     gradient equal to zero in the y and z directions
c***  
      do k=2,nz-1
         do i=1,nx
            scal(i,1,k ) = 0 !scal(i,2   ,k)
            scal(i,ny,k) = 0 !scal(i,ny-1,k)
         enddo
      enddo 
      do j=1,ny
         do i=1,nx
            scal(i,j,1 ) = 0 !scal(i,j,2)
            scal(i,j,nz) = 0 !scal(i,j,nz-1)
         enddo
      enddo 
c***
!      CALL BUPDAT3D(scal,NX,NY,NZ)
      return
      end




c***********************************************************************
      SUBROUTINE bc_flux(vof1,vof3,vel,nx,ny,nz,indx)
c***
 
      INTEGER nx, ny, nz, indx, i, j, k
      DOUBLE PRECISION vof1(nx,ny,nz), vof3(nx,ny,nz), vel(nx,ny,nz)
c***
c     along x direction
c***
      if (indx .eq. 1) then
!        IF(NEIGHB(1).LT.0) THEN
          do k=2,nz-1
            do j=2,ny-1
               vof3(1 ,j,k) = 0.d0
            enddo
          enddo
!        ENDIF
!        IF(NEIGHB(2).LT.0) THEN
          do k=2,nz-1
            do j=2,ny-1
               vof1(nx,j,k) = 0.d0
            enddo
          enddo
!        ENDIF
c***
c     along y direction
c***
      elseif (indx .eq. 2) then
         do k=2,nz-1
            do i=2,nx-1
               vof3(i,1 ,k) = 0 !DMAX1(0.0d0,vel(i,2,k))
               vof1(i,ny,k) = 0.d0
            enddo
         enddo
c***
c     along z direction
c***
      elseif (indx .eq. 3) then
         do j=2,ny-1
            do i=2,nx-1
               vof3(i,j,1) = 0 !DMAX1(0.0d0,vel(i,j,2))
               vof1(i,j,nz) = 0.0d0
            enddo
         enddo
c***  
c     wrong value
c***  
      else 
         stop 'wrong value for indx in one of the swaps'
      endif 
c***
!      CALL BUPDAT3D(vof1,nx,ny,nz)
!      CALL BUPDAT3D(vof3,nx,ny,nz)
      return
      end






















  

