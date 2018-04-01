c *** *** *** *** *** *** *** *** *** *** *** *** 
c    convert level set phi to VOF function 
c *** *** *** *** *** *** *** *** *** *** *** ***
      subroutine levelset2vof(nx,ny,nz,ls,cc)

      double precision cc(nx,ny,nz),ls(nx,ny,nz)
      double precision zero, one, normL1
      double precision mm1,mm2,mm3,mx,my,mz,alpha
      double precision fl3d
      integer i,j,k

      zero=0.d0
      one=1.d0
      do k=2,nz-1
         do j=2,ny-1
            do i=2,nx-1
c***
c     (1) normal vector: mx,my,mz; (2) mx,my,mz>0. and mx+my+mz = 1.;
c     (3) shift alpha to origin    (4) get volume from alpha;  
c***

               mm1 = ls(i-1,j-1,k-1)+ls(i-1,j-1,k+1)+ls(i-1,j+1,k-1)
     %              +ls(i-1,j+1,k+1)+2.0d0*(ls(i-1,j-1,k)+ls(i-1,j+1,k)
     %              +ls(i-1,j,k-1)+ls(i-1,j,k+1))+4.0d0*ls(i-1,j,k)
               mm2 = ls(i+1,j-1,k-1)+ls(i+1,j-1,k+1)+ls(i+1,j+1,k-1)
     %              +ls(i+1,j+1,k+1)+2.0d0*(ls(i+1,j-1,k)+ls(i+1,j+1,k)
     %              +ls(i+1,j,k-1)+ls(i+1,j,k+1))+4.0d0*ls(i+1,j,k)
               mx = (mm1 - mm2)/32.d0
               
               mm1 = ls(i-1,j-1,k-1)+ls(i-1,j-1,k+1)+ls(i+1,j-1,k-1)
     %              +ls(i+1,j-1,k+1)+2.0d0*(ls(i-1,j-1,k)+ls(i+1,j-1,k)
     %              +ls(i,j-1,k-1)+ls(i,j-1,k+1))+4.0d0*ls(i,j-1,k)
               mm2 = ls(i-1,j+1,k-1)+ls(i-1,j+1,k+1)+ls(i+1,j+1,k-1)
     %              +ls(i+1,j+1,k+1)+2.0d0*(ls(i-1,j+1,k)+ls(i+1,j+1,k)
     %              +ls(i,j+1,k-1)+ls(i,j+1,k+1))+4.0d0*ls(i,j+1,k)
               my = (mm1 - mm2)/32.d0
               
               mm1 = ls(i-1,j-1,k-1)+ls(i-1,j+1,k-1)+ls(i+1,j-1,k-1)
     %              +ls(i+1,j+1,k-1)+2.0d0*(ls(i-1,j,k-1)+ls(i+1,j,k-1)
     %              +ls(i,j-1,k-1)+ls(i,j+1,k-1))+4.0d0*ls(i,j,k-1)
               mm2 = ls(i-1,j-1,k+1)+ls(i-1,j+1,k+1)+ls(i+1,j-1,k+1)
     %              +ls(i+1,j+1,k+1)+2.0d0*(ls(i-1,j,k+1)+ls(i+1,j,k+1)
     %              +ls(i,j-1,k+1)+ls(i,j+1,k+1))+4.0d0*ls(i,j,k+1)
               mz = (mm1 - mm2)/32.D0
c     *(2)*  
               mx = DABS(MX)
               MY = DABS(MY)
               MZ = DABS(MZ)
               normL1 = mx+my+mz
               mx = mx/normL1
               my = my/normL1
               mz = mz/normL1

               alpha = ls(i,j,k)/normL1
c     *(3)*  
               alpha = alpha + 0.5d0
c     *(4)*  
               if(alpha.ge.1.d0) then 
                  cc(i,j,k) = 1.d0
               else if (alpha.le.0.d0) then
                  cc(i,j,k) = 0.d0
               else 
                  cc(i,j,k) = FL3D(mx,my,mz,alpha,zero,one)
c                  write(*,*) cc(i,j,k)
               end if
            enddo
         enddo
      enddo 
      return
      end












!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE PARAVIEW(tec,NX,NY,NZ,P)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!----------------------------------------------------------------------- 
!----------------------------------------------------------------------- 

      character*30 tec, ch
      real*8  ::  P(NX,NY,NZ)
      
      real*8,allocatable ::   tp(:)
      

      allocate(tp(nx*ny*nz))
      

      n = 0
      do k = 1,nz
      do j = 1,ny
      do i = 1,nx
        n = n + 1
        tp(n) = P(i,j,k)
      enddo
      enddo
      enddo
      
      open (unit=1, file=tec)
      write(1,'(a)') '# vtk DataFile Version 3.0'
      write(1,'(a)') 'Non-uniform Rectilinear - Rectilinear Grid' 
      write(1,'(a)') 'ASCII' 
      write(1,'(a)') 'DATASET RECTILINEAR_GRID' 
      write(1,'(a,3i)') 'DIMENSIONS',nx,ny,nz
      write(1,'(a,i,a)') 'X_COORDINATES',nx,'  float'
      write(1,'(f10.3)') (i*1.0,i=1,nx)
      write(1,'(a,i,a)') 'Y_COORDINATES',ny,'  float'
      write(1,'(f10.3)') (i*1.0,i=1,ny)
      write(1,'(a,i,a)') 'Z_COORDINATES',nz,'  float'
      write(1,'(f10.3)') (i*1.0,i=1,nz)

      write(1,'(a,i)') 'POINT_DATA',nx*ny*nz
      write(1,'(a)') 'SCALARS P float 1'
      write(1,'(a)') 'LOOKUP_TABLE default'
      write(1,'(F)')((tp(i)),i=1,n)
      close(1)
      

      deallocate (tp)
      ENDSUBROUTINE
!-----------------------------------------------------------------------