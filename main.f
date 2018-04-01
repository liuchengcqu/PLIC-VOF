
      !include 'param.h'
      use param
      
      character*30 tec, ch
      data tec /'out-t0000000000.vtk'/
      
      real*8,allocatable :: f(:,:,:)
      real*8,allocatable :: ux(:,:,:),uy(:,:,:),uz(:,:,:)
      real*8 t, dx,dy,dz,x,y,z,TT
      integer i,j,k,step,tswap

      real*8,allocatable :: fn(:,:,:)
      real*8 tmp,tmp1,init_f
      integer  nnx
      
      PI = 3.141592657589793d0
      nx = 100
      ny = 100
      nz = 100
      dt = 0.001
      delta = 1.0d0/NX
      
      allocate(ux(nx,ny,nz))
      allocate(uy(nx,ny,nz))
      allocate(uz(nx,ny,nz))
      allocate(f(nx,ny,nz))
      allocate(fn(nx,ny,nz))
      ux = 0.0
      uy = 0.0
      uz = 0.0
      f = 0.0
      fn = 0.0
      
      call init(f,ux,uy,uz,t,dx,dy,dz)
      step = 0 
      
      init_f = 0.0d0
      do k=1,NZ
      do j=1,NY
      do i=1,NX
         fn(i,j,k)=f(i,j,k)
         init_f = init_f+f(i,j,k)
      end do
      end do
      end do
      

      call PARAVIEW(tec,NX,NY,NZ,f) 

      
      open (11,file='mass.txt')
      
      TT = 3
      tswap = 0
      t = 0.0

      do step=1,500000 

c ***
c initializing of velocity field 
c ***

         do k=1,NZ
         do j=1,NY
         do i=1,NX
           x = (i-1)*dx-dx*0.5d0
           y = (j-1)*dy
           z = (k-1)*dz
        ux(i,j,k)=2*(sin(pi*x)**2)*sin(2*pi*y)*sin(2*pi*z)*cos(pi*t/TT) 
         enddo
         enddo
         enddo
         do k=1,NZ
         do j=1,NY
         do i=1,NX
           x = (i-1)*dx
           y = (j-1)*dy-dy*0.5d0
           z = (k-1)*dz 
        uy(i,j,k)= -sin(2*pi*x)*(sin(pi*y)**2)*sin(2*pi*z)*cos(pi*t/TT) 
         enddo
         enddo
         enddo
         do k=1,NZ
         do j=1,NY
         do i=1,NX
           x = (i-1)*dx
           y = (j-1)*dy
           z = (k-1)*dz-dz*0.5d0 
        uz(i,j,k)= -sin(2*pi*x)*sin(2*pi*y)*(sin(pi*z)**2)*cos(pi*t/TT)
         enddo
         enddo
         enddo
          
          tswap = tswap + 1
          if (tswap .gt. 3) tswap = 1
          
          if (MOD(tswap,3) .eq. 0) then
            call swp(uz,f,3)
            call swp(ux,f,1)
            call swp(uy,f,2)
          elseif (MOD(tswap,2) .eq. 0) then
            call swp(uy,f,2)
            call swp(uz,f,3)
            call swp(ux,f,1)
          else 
            call swp(ux,f,1)
            call swp(uy,f,2)
            call swp(uz,f,3)
          endif  
          
          call bc_c(f,nx,ny,nz)


         tmp = 0.0d0
         do k=1,NZ
         do j=1,NY
         do i=1,NX
            tmp = tmp +f(i,j,k)
         end do
         end do
         end do

         write(6,1001) step,t , 'Mass=',(tmp-init_f)*dx*dy*dz
 1001    format(i5,f,2x,a,d30.23)
         
         write(11,*)  step,t/TT, 'Mass=',(tmp-init_f)*dx*dy*dz
         

         if(mod(step,400*1).eq.0 .or. t/TT == 1.0 .or. t/TT == 2.0)then 
              call PARAVIEW(tec,NX,NY,NZ,f)
              !pause
         end if
         
         if(t/TT >= 1.0) then
         call PARAVIEW(tec,NX,NY,NZ,f)
         stop
         endif
         
         t = t+dt
         
      end do
      close(11)
      

      tmp=0.0d0
      tmp1=0.0d0
      do k=1,NZ
      do j=1,NY
      do i=1,NX
         tmp1 = tmp1 + f(i,j,k)
      end do
      end do
      end do
      tmp1 = (tmp1-init_f)/init_f
      write(6,*) '# Mass error=',tmp1
      
      deallocate(ux,uy,uz,f,fn) 

 1002 format(i5,2x,i5,2x,i5,2x,d30.23)

      stop
      end

