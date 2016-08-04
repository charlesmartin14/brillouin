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
