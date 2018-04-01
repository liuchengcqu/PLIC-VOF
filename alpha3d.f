c ****** 1 ******* 2 ******* 3 ******* 4 ******* 5 ******* 6 ******* 7 *
c PROGRAM TO FIND alpha IN: m1 x1 + m2 x2 + m3 x3 = alpha,
c GIVEN m1+m2+m3=1 (all > 0) AND THE VOLUMETRIC FRACTION cc
c ****** 1 ******* 2 ******* 3 ******* 4 ******* 5 ******* 6 ******* 7 *
      DOUBLE PRECISION FUNCTION AL3D(b1,b2,b3,cc)
c***
C      INCLUDE   'param.h'
      use param
      DOUBLE PRECISION m1,m2,m3,cc,b1,b2,b3,tmp,pr,ch,mm,m12
      DOUBLE PRECISION p,p12,q,teta,cs
      DOUBLE PRECISION UNTIER,V1,V2,V3
      DOUBLE PRECISION alpha, w, vm1, vm2, vm3, vm12 
      DOUBLE PRECISION a0, a1, a2, q0, sp, th ,CONST_TINY,CONST_PI
      
      DOUBLE PRECISION xi, invp, vma, vmb, vmc
      DOUBLE PRECISION, parameter :: ONE = 1.0d0, PB = 1.49d0, 
     &   PC2 = 0.239d0, PC1 = 0.132d0, 
     &   PC0 = (PB * (PB * PC2 + 4d0 * PC1 - 8d0) / 16d0), 
     &   PA = (PB * PB * (PB - 1d0))
      
      PARAMETER (UNTIER=1.d0/3.d0)
      INTRINSIC DMAX1,DMIN1,DSQRT,DACOS,DCOS
      
      CONST_TINY = 1D-25
      CONST_PI = 3.14159265358979323846d0
c***  
c     (1) order coefficients: m1<m2<m3; (2) get ranges: V1<V2<v3;
c     (3) limit ch (0.d0 < ch < 0.5d0); (4) calculate alpha
c*(1)* 
      m1 = DMIN1(b1,b2)
      m3 = DMAX1(b1,b2)
      m2 = b3
      if (m2 .LT. m1) then
         tmp = m1
         m1 = m2
         m2 = tmp
      else if (m2 .GT. m3) then
         tmp = m3
         m3 = m2
         m2 = tmp
      endif
c*(2)*
      m12 = m1 + m2 
      pr  = DMAX1(6.d0*m1*m2*m3,1.d-50)
      V1  = m1*m1*m1/pr
      V2  = V1 + 0.5d0*(m2-m1)/m3
      if (m3 .LT. m12) then
         mm = m3
         V3 = (m3*m3*(3.d0*m12-m3) + m1*m1*(m1-3.d0*m3) +
     %        m2*m2*(m2-3.d0*m3))/pr
      else
         mm = m12
         V3 = 0.5d0*mm/m3
      endif
c*(3)*
      ch = DMIN1(cc,1.d0-cc)
c*(4)*      
      
!original method      
      if (ch .LT. V1) then
c***         AL3D = cbrt(pr*ch)
         AL3D = (pr*ch)**UNTIER
      else if (ch .LT. V2) then
         AL3D = 0.5d0*(m1 + DSQRT(m1*m1 + 8.d0*m2*m3*(ch-V1)))
      else if (ch .LT. V3) then
         p = 2.d0*m1*m2
         q = 1.5d0*m1*m2*(m12 - 2.d0*m3*ch)
         p12 = DSQRT(p)
         teta = DACOS(q/(p*p12))/3.d0
         cs = DCOS(teta)
         AL3D = p12*(DSQRT(3.d0*(1.d0-cs*cs)) - cs) + m12
      else if (m12 .LT. m3) then
         AL3D = m3*ch + 0.5d0*mm
      else 
         p = m1*(m2+m3) + m2*m3 - 0.25d0
         q = 1.5d0*m1*m2*m3*(0.5d0-ch)
         p12 = DSQRT(p)
         teta = DACOS(q/(p*p12))/3.0
         cs = DCOS(teta)
         AL3D = p12*(DSQRT(3.d0*(1.d0-cs*cs)) - cs) + 0.5d0
      endif

      if (cc .GT. 0.5d0)  AL3D = 1.d0 - AL3D
      

      
c***
      return
      end
