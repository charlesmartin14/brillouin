      implicit real*8(a-b,d-h,o-z)
      implicit complex*16(c)
      dimension ca(4,4),cb(4,4),cgamma(0:3,4,4),
     1     cgamma5(4,4),gmunu(0:3,0:3),cu1bar(4),cu2(4),cukbar(4),
     2     cue(4),p1(0:3),p2(0:3),pk(0:3),cs1(2),cs2(2),csk(2),
     3     cleft(4,4),cGVGA(4,4),xg(100),wxg(100),pg(100),wpg(100),
     4     pe(0:3),cse(2)
      common /dirac/ ci,cgamma,cgamma5,gmunu
      common /con/ pi,hbarc,rc,alpha,GF,GV,GA,dMp,dMn,dme,dmnu,dMd,gnpd,
     1     dkapp,dkapn,dufac,dN0,e


      call const
      call cdirac
      nx = 8
      call gaus(xg,wxg,nx)
      np = 8
      call gaus(pg,wpg,np)

      do ip =1,np
         pg(ip) =(pg(ip)+1.d0)*pi
         wpg(ip) = wpg(ip)*pi
      enddo
      

      Ekemin =((dMn + dmnu - dme)**2 - dMp**2)/(2.d0*(dMn+dmnu))
      pmin = dsqrt((Ekemin+dme)**2 - dme**2)
      E2check =((dMn + dmnu)**2 -(dMp+dme)**2)/(2.d0*dMp)


      do iep =1,500
         pep =(pmin + dexp(dlog(1.d-11)+ dfloat(iep)*
     1        (dlog(100.d0)-dlog(1.d-11))/500.d0))/ dsqrt(3.d0)
         
         rlconf = 1.d-5*hbarc*pi/pep
         call pekin0(pep,pe,p1,p2,pk)
         
         EKp = p2(0) - dMp
         EKn = p1(0) - dMn
         EKe = pe(0) - dme
         EKnu = pk(0) - dmnu
         write(*,*) rlconf,EKe,EKn,EKnu,EKp
         
      enddo

     
      
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

      
      subroutine spinu(p,rm,cs,cu)
      implicit real*8(a-b,d-h,o-z)
      implicit complex*16(c)
      dimension cgamma(0:3,4,4),cgamma5(4,4),cwu(4,2)
      dimension cu(4),gmunu(0:3,0:3),cs(2),cphi(2),p(0:3)
      common /dirac/ ci,cgamma,cgamma5,gmunu
      fac = dsqrt(p(0) + rm)
      cphi(1) = cs(1)*fac
      cphi(2) = cs(2)*fac
 
      cwu(1,1) = 1.d0
      cwu(1,2) = 0.d0
      cwu(2,1) = 0.d0
      cwu(2,2) = 1.d0
      cwu(3,1) =(p(3))/(p(0) + rm)
      cwu(3,2) =(p(1)- ci*p(2))/(p(0) + rm)
      cwu(4,1) =(p(1)+ ci*p(2))/(p(0) + rm)
      cwu(4,2) =(-p(3))/(p(0) + rm)
      do i1 =1,4
        cu(i1) = 0.d0
         doi2 = 1,2
            cu(i1) = cu(i1) + cwu(i1,i2)*cphi(i2)
         enddo
      enddo
 
      RETURN
      end

      

      subroutine spinubar(p,rm,cs,cubar)
      implicit real*8(a-b,d-h,o-z)
      implicit complex*16(c)

      dimension cgamma(0:3,4,4),cgamma5(4,4),gmunu(0:3,0:3)
      dimension cwu(4,2),cu(4),cubar(4),cs(2),cphi(2),p(0:3)
      common /dirac/ ci,cgamma,cgamma5,gmunu

      fac = dsqrt(p(0) + rm)
      cphi(1) = cs(1)*fac
      cphi(2) = cs(2)*fac
 
      cwu(1,1) = 1.d0
      cwu(1,2) = 0.d0
      cwu(2,1) = 0.d0
      cwu(2,2) = 1.d0
      cwu(3,1) =(p(3))/(p(0) + rm)
      cwu(3,2) =(p(1)- ci*p(2))/(p(0) + rm)
      cwu(4,1) =(p(1)+ ci*p(2))/(p(0) + rm)
      cwu(4,2) =(-p(3))/(p(0) + rm)
      do i1 =1,4
         cu(i1) = 0.d0
         do i2 =1,2
            cu(i1) = cu(i1) + cwu(i1,i2)*cphi(i2)
         enddo
      enddo

      cubar(1) = dconjg(cu(1))
      cubar(2) = dconjg(cu(2))
      cubar(3) = -dconjg(cu(3))
      cubar(4) = -dconjg(cu(4))
      RETURN
      end
      
      

      subroutine ubaru(cubar,ca,cu,cscal)
      implicit real*8(a-b,d-h,o-z)
      implicit complex*16(c)
      dimension cubar(4),ca(4,4),cu(4)
      cscal = 0.d0
      do i1 =1,4
         do i2 =1,4
          cscal = cscal + cubar(i1)*ca(i1,i2)*cu(i2)
        enddo
      enddo

      return
      end
      

      subroutine chsgng(ca)
      implicit real*8(a-b,d-h,o-z)
      implicit complex*16(c)

      dimension ca(4,4)

      do i1 =1,4
        do i2 =1,4
            ca(i1,i2) = -ca(i1,i2)
        enddo
      enddo
      
      RETURN
      end
      

      subroutine vcopyg(mu,caa,cb)
      implicit real*8(a-b,d-h,o-z)
      implicit complex*16(c)
      dimension caa(0:3,4,4),cb(4,4)

      do i1 =1,4
         do i2 =1,4
            cb(i1,i2) = caa(mu,i1,i2)
         enddo
      enddo

      RETURN
      end


      subroutine subg(ca,cb,cc)
      implicit real*8(a-b,d-h,o-z)
      implicit complex*16(c)
      dimension ca(4,4),cb(4,4),cc(4,4)

      do i1 =1,4
         do i2 =1,4
            cc(i1,i2) = ca(i1,i2) - cb(i1,i2)
         enddo
      enddo

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



      subroutine saddg(cscal,ca,cb)
      implicit real*8(a-b,d-h,o-z)
      implicit complex*16(c)
      dimension ca(4,4),cb(4,4)

      do i1 =1,4
        do i2 =1,4
            cb(i1,i2) = ca(i1,i2)
        enddo
      enddo

      do i1 =1,4
         cb(i1,i1) = cb(i1,i1) + cscal
      enddo

      RETURN
      end



      subroutine addg(ca,cb,cc)
      implicit real*8(a-b,d-h,o-z)
      implicit complex*16(c)

      dimension ca(4,4),cb(4,4),cc(4,4)

      do i1 =1,4
        do i2 =1,4
            cc(i1,i2) = ca(i1,i2) + cb(i1,i2)
        enddo
      enddo
      
      RETURN
      end
      


      

      subroutine smultg(cscal,ca,cb)
      implicit real*8(a-b,d-h,o-z)
      implicit complex*16(c)

      dimension ca(4,4),cb(4,4)

      do i1 =1,4
         do i2 =1,4
            cb(i1,i2) = ca(i1,i2)*cscal
         enddo
      enddo

      RETURN
      end

      


      
      subroutine multg(ca,cb,cc)
      implicit real*8(a-b,d-h,o-z)
      implicit complex*16(c)

      dimension ca(4,4),cb(4,4),cc(4,4),cd(4,4)
      do i1 =1,4
         do i2 =1,4
            cd(i1,i2) =(0.d0,0.d0)
            do j =1,4
               cd(i1,i2) = cd(i1,i2) + ca(i1,j)*cb(j,i2)
            enddo
         enddo
      enddo

      do i1 =1,4
         do i2 =1,4
            cc(i1,i2) = cd(i1,i2)
         enddo
      enddo

      RETURN
      end

      
      subroutine pekin0(pep,pe,p1,p2,pk)
      implicit real*8(a-b,d-h,o-z)
      implicit complex*16(c)

      dimension p1(0:3),p2(0:3),pk(0:3),pe(0:3)
      common /con/ pi,hbarc,rc,alpha,GF,GV,GA,dMp,dMn,dme,dmnu,dMd,gnpd
     1  ,dkapp,dkapn,dufac,dN0,e

      p2(0) = dsqrt(dMp**2 + 3*(pep**2))
      pe(0) = dsqrt(dme**2 + 3*(pep**2))

      Etot = p2(0) + pe(0)

      Eknu =((Etot-dmnu)**2 - dMn**2)/(2.d0*Etot)
      Ekn =((Etot-dMn)**2 - dmnu**2)/(2.d0*Etot)
      pnu = dsqrt((Eknu+dmnu)**2 - dmnu**2)

      pk(0) = Eknu + dmnu
      p1(0) = p2(0) + pe(0) - pk(0)

      RETURN
      end

      

      
      SUBROUTINE GAUS(Y,WY,N)
      IMPLICIT REAL*8(A-H,O-Z)                                                  
C     THIS ROUTINE DETERMINES THE N/2 POSITIVE ZEROS(X(N/2)) OF THE             
C     LEGENDRE POLYNOMIAL OF ORDER N.  ALSO THE N/2 DISTINCT WEIGHTS            
C   FOR GAUSS QUADRATURE ARE CALCULATED.  N IS EVEN AND LESS THAN 101.          
      DOUBLE PRECISION V,Z,Z1,P,P1,DP,POL                                       
      DIMENSION X(50),W(50),POL(101)                                            
      DIMENSION Y(*),WY(*)                                                      
      NP1=N+1                                                                   
      M = N/2                                                                   
      MP1=M+1                                                                   
      Z = 3.14159265358979D0/(2.0D0*FLOAT(N))                                   
      DO 20 L=1,M                                                               
      DO 5 K=1,100                                                              
      CALL PL(Z,N,POL)                                                          
      P=POL(NP1)                                                                
      DP=N*(Z*POL(NP1)-POL(N))/(Z*Z-1.D0)                                       
      Z = Z - P/DP                                                              
      IF(DABS(P).LT..1D-11) GO TO 15                                            
    5 CONTINUE                                                                  
      WRITE(6,1000) L                                                           
      RETURN                                                                    
   15 X(L) = Z                                                                  
      V = 2.0D0/((1.0D0 - Z*Z)*DP*DP)                                           
      W(L) = V                                                                  
      IF(L.EQ.M)  GO TO 30                                                      
      DZ=Z                                                                      
      IF (L.GT.1) DZ=(X(L)-X(L-1))*.5D0                                         
      DO 17 K=1,200                                                             
      Z = Z + DZ                                                                
      CALL PL(Z,N,POL)                                                          
      P=POL(NP1)                                                                
      Z1 = Z + DZ                                                               
      CALL PL(Z1,N,POL)                                                         
      P1=POL(NP1)                                                               
      IF(P*P1.LT.0.0D0) GO TO 18                                                
   17 CONTINUE                                                                  
   18 Z = (P1*Z - P*Z1)/(P1 - P)                                                
   20 CONTINUE                                                                  
   30 DO 40 NEG=1,N                                                             
      IF(NEG.LE.M)Y(NEG)=-X(MP1-NEG)                                            
      IF(NEG.GT.M)Y(NEG)=X(NEG-M)                                               
      IF(NEG.LE.M)WY(NEG)=W(MP1-NEG)                                            
      IF(NEG.GT.M)WY(NEG)=W(NEG-M)                                              
   40 CONTINUE                                                                  
      RETURN                                                                    
 1000 FORMAT(5X,40HFAILURE TO GET ZEROS OF LEGENDRES FOR L=,I3)                 
      END                                                      
      SUBROUTINE PL(X,L,POL)                                                    
      IMPLICIT REAL*8(A-H,O-Z)                                                  
C     COMPUTES LEGENDRE POLYNOMINALS OF ORDER ZERO TO ORDER                     
C     (L) FOR A GIVEN ARGUMENT (X=COS(THETA))                                   
      DOUBLE PRECISION S,X,POL                                                  
      DIMENSION POL(*)                                                          
      POL(1)=1.0D0                                                              
      POL(2)=X                                                                  
      IF(L.LE.1)GOTO20                                                          
      DO 10 JJ=2,L                                                              
      S=JJ                                                                      
   10 POL(JJ+1)=((2.D0*S-1.D0)*X*POL(JJ)-(S-1.D0)*POL(JJ-1))/S                  
   20 RETURN                                                                    
      END                                                                 
