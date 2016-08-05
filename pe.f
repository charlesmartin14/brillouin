      implicit real*8(a-b,d-h,o-z)
      implicit complex*16(c)
      dimension ca(4,4),cb(4,4),cgamma(0:3,4,4),
     1     gamma5(4,4),gmunu(0:3,0:3)
      common /dirac/ ci,cgamma,cgamma5,gmunu
      common /con/ pi,hbarc,rc,alpha,GF,GV,GA,dMp,dMn,dme,dmnu,dMd,gnpd,
     1     dkapp,dkapn,dufac,dN0,e


      call const
      call cdirac
 
      
      PRINT '(A)','PNLL code'
      call copyg(cgamma5,ca)
     
      
      STOP
      END


      


      subroutine const
      implicit real*8(a-h,o-z)
      common /con/ pi,hbarc,rc,alpha,GF,GV,GA,dMp,dMn,dme,dmnu,dMd,gnpd
     1    ,dkapp,dkapn,dufac,dN0,e

      pi = 4.d0*datan(1.d0)
      hbarc = 197.3269631d0    ! MeV-fm
      rc = 299792458.d0        ! m/s
      alpha = 7.2973525376d-3
      e = 1.602176487d-13      ! J/MeV
 
      dN0 = 6.02214179d23      ! molˆ-1
      dMp = 938.272013d0       ! MeV
      dMn = 939.565346d0       ! MeV
      dme = 0.510998910d0      ! MeV
      dmnu = 1.d-8              ! MeV
      dMd = 1875.612793d0      ! MeV
      GF = 1.16637d-11         ! 1/MeVˆ2
      GV = 1.d0                ! 1.013d0
      GA = -1.285d0*GV          ! -1.267
      gnpd = 10.6617553d0      ! 11.3 in 9704031 v5.pdf(pg 8-9)
      dkapp= 1.793d0
      dkapn= -1.913d0
      dufac= 0.d0
 
      RETURN
      end


      subroutine cdirac
      implicit real*8(a-b,d-h,o-z)
      implicit complex*16(c)
      dimension cgamma(0:3,4,4),cgamma5(4,4),gmunu(0:3,0:3)
      common /dirac/ ci,cgamma,cgamma5,gmunu
      ci =(0.d0,1.d0)

      do i1 =1,4
         do i2 =1,4
            cgamma5(i1,i2) =(0.d0,0.d0)
            do imu = 0,3
               cgamma(imu,i1,i2) =(0.d0,0.d0)
            enddo
         enddo
      enddo
      cgamma(0,1,1) = 1.d0
      cgamma(0,2,2) = 1.d0
      cgamma(0,3,3) = -1.d0
      cgamma(0,4,4) = -1.d0

      cgamma(1,1,4) = 1.d0
      cgamma(1,2,3) = 1.d0
      cgamma(1,3,2) = -1.d0
      cgamma(1,4,1) = -1.d0
      
      cgamma(2,1,4) = -ci
      cgamma(2,2,3) = ci
      cgamma(2,3,2) = ci
      cgamma(2,4,1) = -ci

      cgamma(3,1,3) = 1.d0
      cgamma(3,2,4) = -1.d0
      cgamma(3,3,1) = -1.d0
      cgamma(3,4,2) = 1.d0
      
      cgamma5(1,3) = 1.d0
      cgamma5(2,4) = 1.d0
      cgamma5(3,1) = 1.d0
      cgamma5(4,2) = 1.d0

      RETURN
      end

      

      subroutine copyg(ca,cb)
      implicit real*8(a-b,d-h,o-z)
      implicit complex*16(c)
      dimension ca(4,4),cb(4,4)

      do i1 =1,4
         do i2 =1,4
            cb(i1,i2) = ca(i1,i2)
         enddo
      enddo

      RETURN
      end



