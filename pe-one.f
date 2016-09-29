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
      
      PRINT '(A)','PNLL code'
      PRINT '(A)','PG'
      write(*,*) pg(1),pg(2),pg(3),pg(4),pg(5),pg(6),pg(7),pg(8)
      PRINT '(A)','WPG'
      write(*,*) wpg(1),wpg(2),wpg(3),wpg(4),wpg(5),wpg(6),wpg(7),wpg(8)


      Ekemin =((dMn + dmnu - dme)**2 - dMp**2)/(2.d0*(dMn+dmnu))
      pmin = dsqrt((Ekemin+dme)**2 - dme**2)
      E2check =((dMn + dmnu)**2 -(dMp+dme)**2)/(2.d0*dMp)
      PRINT '(A)','EKemin,pmin,E2check'
      write(*,*) Ekemin,pmin,E2check


      call copyg(cgamma5,ca)
      call chsgng(ca)

      cscal = 1.d0
      call saddg(cscal,ca,cleft)
      call copyg(cgamma5,ca)
      cscal = GA
      call smultg(cscal,ca,ca)
      call chsgng(ca)
 
      cscal = GV
      call saddg(cscal,ca,cGVGA)
 
      do is1 =1,2
         cs1(is1) = 0.d0
      enddo
      do is2 =1,2
        cs2(is2) = 0.d0
      enddo
      do ise =1,2
         cse(ise) = 0.d0
      enddo
      do isk =1,2
       csk(isk) = 0.d0
      enddo

      tran2 = 0.d0
      trac = 0.d0
      densep = 7.d0 ! g/cc NiH density ??
      dNiMW = 28.d0 ! g/mol NiH MW ??
      dndens = dN0*densep/dNiMW

      open(9,file='plot.dat',status='UNKNOWN')

      do iep =1,500
      pep =(pmin + dexp(dlog(1.d-11)+ dfloat(iep)*
     1        (dlog(100.d0)-dlog(1.d-11))/500.d0))/ dsqrt(3.d0)
      
      rlconf = 1.d-5*hbarc*pi/pep

      test = 0.d0
      tran2 = 0.d0
      do ix =1,nx
         do ip =1,np
            
      do is1 =1,2
         cs1(is1) = 1.d0
      do is2 =1,2
         cs2(is2) = 1.d0
      do ise =1,2
         cse(ise) = 1.d0
      do isk =1,2
         csk(isk) = 1.d0
         
C     do isgn1 = 1,2
C     do isgn2 = 1,2
C     do isgn3 = 1,2
         isgn1=1
         isgn2=1
         isgn3=1
            
      camp2 = 0.d0
      pep1 = pep*((-1.d0)**isgn1)
      pep2 = pep*((-1.d0)**isgn2)
      pep3 = pep*((-1.d0)**isgn3)
      call pekin(xg(ix),pg(ip),pep1,pep2,pep3,pe,p1,p2,pk)

      call spinubar(p1,dMn,cs1,cu1bar)
      call spinu(p2,dMp,cs2,cu2)
      call spinu(pe,dme,cse,cue)
      call spinubar(pk,dmnu,csk,cukbar)
    
      call vcopyg(0,cgamma,cb)
      call multg(cGVGA,cb,cb)
      call ubaru(cu1bar,cb,cu2,chad0)
      call vcopyg(0,cgamma,cb)
      call multg(cb,clef t,cb)
      call ubaru(cukbar,cb,cue,clep0)
    
    
      call vcopyg(1,cgamma,cb)
      call multg(cGVGA,cb,cb)
      call ubaru(cu1bar,cb,cu2,chad1)
      call vcopyg(1,cgamma,cb)
      call multg(cb,clef t,cb)
      call ubaru(cukbar,cb,cue,clep1)
      call vcopyg(2,cgamma,cb)

      call multg(cGVGA,cb,cb)
      call ubaru(cu1bar,cb,cu2,chad2)
      call vcopyg(2,cgamma,cb)
      call multg(cb,clef t,cb)
      call ubaru(cukbar,cb,cue,clep2)
    
      call vcopyg(3,cgamma,cb)
      call multg(cGVGA,cb,cb)
      call ubaru(cu1bar,cb,cu2,chad3)
      call vcopyg(3,cgamma,cb)
      call multg(cb,clef t,cb)
      call ubaru(cukbar,cb,cue,clep3)
      camp2 = camp2 +
     1     (chad0* clep0 - chad1* clep1 - chad2* clep2 - chad3* clep3)

      write(6,*) "chad0",chad0
      write(6,*) "chad1",chad1
      write(6,*) "chad2",chad2
      write(6,*) "chad3",chad3

      write(6,*) "clep0",clep0
      write(6,*) "clep1",clep1
      write(6,*) "clep2",clep2
      write(6,*) "clep3",clep3

      write(6,*) "camp2",camp2
    
C      enddo
C      enddo
C      enddo

      ftop =rc*((pep/(2.d0*pi))**3)*(GF**2)*(dsqrt(pk(0)**2-dmnu**2))**3
      fbot = 
     1     (512.d0 *(pi **2)* p2(0)* pe(0)* hbarc * 1.d-15*
     1     dabs((pk(1) *(p1(0)* pk(1)-pk(0)* p1(1))) +
     1     (pk(2)*(p1(0)* pk(2)-pk(0)* p1(2))) +
     1     (pk(3)*(p1(0)* pk(3)-pk(0)* p1(3)))))

      fact =rc*((pep/(2.d0*pi))**3)*(GF**2)*(dsqrt(pk(0)**2-dmnu**2))**3
     1     /(512.d0 *(pi **2)* p2(0)* pe(0)* hbarc * 1.d-15*
     1     dabs((pk(1) *(p1(0)* pk(1)-pk(0)* p1(1))) +
     1     (pk(2)*(p1(0)* pk(2)-pk(0)* p1(2))) +
     1     (pk(3)*(p1(0)* pk(3)-pk(0)* p1(3)))))
    
      EKn = p1(0) - dMn
      EKe = pe(0) - dme

      trac = trac +  dreal(camp2* dconjg(camp2)) 
            
      tran2 = tran2 + wxg(ix)*wpg(ip)*
     1    dreal(camp2* dconjg(camp2)) * fact *dndens

      write(6,*) 'ftop', ftop
      write(6,*) 'fbot', fbot
      write(6,*) 'fact', fact
      stop
    
      csk(isk) = 0.d0
      enddo
      cse(ise) = 0.d0
      enddo
      cs2(is2) = 0.d0
      enddo
      cs1(is1) = 0.d0
      enddo

      
      enddo
      enddo
    
      EKp = p2(0) - dMp
      EKn = p1(0) - dMn
      EKe = pe(0) - dme
      EKnu = pk(0) - dmnu
      write(*,*) rlconf,tran2,EKe,EKn

      write(9,*) rlconf,
     1 tran2 *((pe(0) - dme)+(p2(0) - dMp)) * e * 1.d-7,
     1 tran2 *(2.2d0+(p1(0) - dMn)) * e * 1.d-7,EKe,EKn,EKnu,Ekp
      enddo
      close(9)

     
      
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

      
      subroutine pekin(x,phi,pep1,pep2,pep3,pe,p1,p2,pk)
      implicit real*8(a-b,d-h,o-z)
      implicit complex*16(c)

      dimension p1(0:3),p2(0:3),pk(0:3),pe(0:3)
      common /con/ pi,hbarc,rc,alpha,GF,GV,GA,dMp,dMn,dme,dmnu,dMd,gnpd
     1     ,dkapp,dkapn,dufac,dN0,e

      write(6,*) 'pekin',x,phi,pep1,pep2,pep3

      p2(0) = dsqrt(dMp**2 + pep1**2 + pep2**2 + pep3**2)
      p2(1) = -pep1
      p2(2) = -pep2
      p2(3) = -pep3
      pe(0) = dsqrt(dme**2 + pep1**2 + pep2**2 + pep3**2)
      pe(1) = pep1
      pe(2) = pep2
      pe(3) = pep3

      write(6,*) 'dMp',dMp
      write(6,*) 'pre p2',p2
      write(6,*) 'pre pe',pe

      Etot = p2(0) + pe(0)

      Eknu =((Etot-dmnu)**2 - dMn**2)/(2.d0*Etot)
      Ekn =((Etot-dMn)**2 - dmnu**2)/(2.d0*Etot)
      pnu = dsqrt((Eknu+dmnu)**2 - dmnu**2)

      write(6,*) 'Ees',Etot, Eknu, Ekn, pnu

      pk(0) = Eknu + dmnu
      pk(1) = pnu*dsqrt(1.d0 - x**2)*dcos(phi)
      pk(2) = pnu*dsqrt(1.d0 - x**2)*dsin(phi)
      pk(3) = pnu*x
      p1(0) = p2(0) + pe(0) - pk(0)
      p1(1) = p2(1) + pe(1) - pk(1)
      p1(2) = p2(2) + pe(2) - pk(2)
      p1(3) = p2(3) + pe(3) - pk(3)

      write(6,*) 'p1',p1
      write(6,*) 'p2',p2
      write(6,*) 'pe',pe
      write(6,*) 'pk',pk
      
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
