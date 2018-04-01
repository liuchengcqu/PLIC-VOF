
c ****** 1 ******* 2 ******* 3 ******* 4 ******* 5 ******* 6 ******* 7 *
c PROGRAM TO FIND THE "CUT VOLUME" V0 GIVEN r0, dr0 AND
c m1 x1 + m2 x2 + m3 x3 = alpha
c ****** 1 ******* 2 ******* 3 ******* 4 ******* 5 ******* 6 ******* 7 *
      DOUBLE PRECISION FUNCTION FL3D(m1,m2,m3,alpha,r0,dr0)
c***
      DOUBLE PRECISION m1,m2,m3,alpha,r0,dr0, vm1,vm2,vm3,vm12,a,v
      DOUBLE PRECISION al,al0,n1,n2,n3,b1,b2,b3,b12,bm,tmp,pr,CONST_TINY
      INTRINSIC DMAX1,DMIN1,DABS
      DOUBLE PRECISION, parameter :: ONE = 1.0d0, PB = 1.49d0, 
     &        PC2 = 0.239d0, PC1 = 0.132d0, 
     &        PC0 = (PB * (PB * PC2 + 4d0 * PC1 - 8d0) / 16d0),  
     &        PA = (PB * PB * (PB - 1d0))
      
      CONST_TINY = 1D-25
c***
c     (1) move origin to r0 along r ;  (2) reflect parallelepiped;
c     (3) limit alpha (0<= al0 <=0.5); (4) order coefficients: b1<b2<b3;
c     (5) calculate volume (NOTE: it is assumed:s0=t0=0; ds0=dt0=1.)
c*(1)*
      al = alpha - m1*r0
c*(2)*
      al = al + DMAX1(0.d0,-m1*dr0)+DMAX1(0.d0,-m2)+DMAX1(0.d0,-m3)
      tmp = DABS(m1)*dr0 + DABS(m2) + DABS(m3)
      n1 = DABS(m1)/tmp
      n2 = DABS(m2)/tmp
      n3 = DABS(m3)/tmp
      al = DMAX1(0.d0,DMIN1(1.d0,al/tmp))
c*(3)*
      al0 = DMIN1(al,1.d0-al)
c*(4)*
      b1 = DMIN1(n1*dr0,n2)
      b3 = DMAX1(n1*dr0,n2)
      b2 = n3
      if (b2 .LT. b1) then
         tmp = b1
         b1 = b2
         b2 = tmp
      else if (b2 .GT. b3) then
         tmp = b3
         b3 = b2
         b2 = tmp
      endif
      b12 = b1 + b2
      bm = DMIN1(b12,b3)
      pr = DMAX1(6.d0*b1*b2*b3,1.0d-50)
      
c*5*     

      
      
! Aoki PLIC method
        vm1 = b1 !min(vma, vmb, vmc) ! Sort vm1, vm2, vm3
        vm3 = b3 !max(vma, vmb, vmc) ! such that vm1 ¡Ü vm2 ¡Ü vm3 .
        vm2 = b2 !abs(1d0 - vm3 - vm1)
        vm12 = b12 !vm1 + vm2
        a = al0
        v = 0d0
        if (a > 0.0d0) then
            if (a < vm1) then
                v = a ** 3 / (6d0 * vm1 * vm2 * vm3)
            else if (a < vm2) then
                v = a * (a - vm1) / (2d0 * vm2 * vm3) +  
     &              vm1 ** 2 / (6d0 * vm2 * vm3 + CONST_TINY)
            else if (a < min(vm12, vm3)) then
                v = (a ** 2 * (3d0 * vm12 - a) + 
     &              vm1 ** 2 * (vm1 - 3d0 * a) + 
     &              vm2 ** 2 * (vm2 - 3d0 * a)) / 
     &              (6d0 * vm1 * vm2 * vm3)
            else if (vm3 < vm12) then
                v = (a ** 2 * (3d0 - 2d0 * a) + 
     &              vm1 ** 2 * (vm1 - 3d0 * a) +  
     &              vm2 ** 2 * (vm2 - 3d0 * a) + 
     &              vm3 ** 2 * (vm3 - 3d0 * a)) / 
     &              (6d0 * vm1 * vm2 * vm3)
            else
                v = (a - 0.5d0 * vm12) / vm3
            end if
        endif
        tmp = v
      
        
        
        
        
      
      if (al .LE. 0.5d0) then
         FL3D = tmp*dr0
      else
         FL3D = (1.d0-tmp)*dr0
      endif
c***  
      return
      end

