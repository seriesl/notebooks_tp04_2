      INTEGER FUNCTION IDMAX(NN, SX, IX)                               
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                              
      DIMENSION SX(NN)                                                 
      IDMAX = 1                                                        
      DO 10 J = 2, NN                                                  
         IF(SX(1+(IDMAX-1)*IX) .LT. SX(1+(J-1)*IX)) IDMAX=J            
10    CONTINUE                                                         
      RETURN                                                           
      END                                                              
C----------------------------------------------------------------------
      INTEGER FUNCTION IDMIN(NN, SX, IX)                               
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                              
      DIMENSION SX(NN)                                                 
      IDMIN = 1                                                        
      DO 10 J = 2, NN                                                  
         IF(SX(1+(IDMIN-1)*IX) .GT. SX(1+(J-1)*IX)) IDMIN=J            
10    CONTINUE                                                         
      RETURN                                                           
      END                                                              
C----------------------------------------------------------------------
      INTEGER FUNCTION IDAMAX(NN, SX, IX)                              
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                              
      DIMENSION SX(NN)                                                 
      IDAMAX = 1                                                       
      DO 10 J = 2, NN                                                  
         IF(DABS(SX(1+(IDAMAX-1)*IX)) .LT. DABS(SX(1+(J-1)*IX)))       
     &  IDAMAX=J                                                      
10    CONTINUE                                                         
      RETURN                                                           
      END                                                              
C----------------------------------------------------------------------
      INTEGER FUNCTION IDAMIN(NN, SX, IX)                              
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                              
      DIMENSION SX(NN)                                                 
      IDAMIN = 1                                                       
      DO 10 J = 2, NN                                                  
         IF(DABS(SX(1+(IDAMIN-1)*IX)) .GT. DABS(SX(1+(J-1)*IX)))       
     &  IDAMIN=J                                                      
10    CONTINUE                                                         
      RETURN                                                           
      END                                                              
C----------------------------------------------------------------------
      SUBROUTINE  DCOPY(NN, SX, INCX, SY, INCY)                        
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                              
      DIMENSION SX(NN), SY(NN)                                         
      DO 10 I = 1, NN                                                  
         SY(1+(I-1)*INCY) = SX(1+(I-1)*INCX)                           
10    CONTINUE                                                         
      RETURN                                                           
      END                                                              
C----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION DNRM2(NN, SX, INCX)                    
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                              
      DIMENSION SX(NN)                                                 
      DNRM2 = 0.0D0                                                    
      DO 10 I=1, NN                                                    
         DNRM2 = DNRM2 + SX(1+(I-1)*INCX)**2                           
c         print *,'SX(1+(I-1)*INCX)', SX(1+(I-1)*INCX),I
10    CONTINUE                                                         
      DNRM2 = DSQRT(DNRM2)                                             
      RETURN                                                           
      END                                                              
C----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION DDOT(NN, SX, INCX, SY, INCY)                    
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                              
      DIMENSION SX(NN), SY(NN)                                                
      DDOT = 0.0D0                                                    
      DO 10 I=1, NN                                                    
         DDOT = DDOT + SX(1+(I-1)*INCX)*SY(1+(I-1)*INCY)
10    CONTINUE                                                         
      RETURN                                                           
      END                                                              
C----------------------------------------------------------------------
      SUBROUTINE DAXPY(NN, SA, SX, INCX, SY, INCY)                     
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                              
      DIMENSION SX(NN), SY(NN)                                         
      DO 10 I = 1, NN                                                  
         SY(1+(I-1)*INCY) = SY(1+(I-1)*INCY) + SA*SX(1+(I-1)*INCX)     
10    CONTINUE                                                         
      RETURN                                                           
      END                                                              
C----------------------------------------------------------------------
      SUBROUTINE DSCAL(NN, SA, SX, INCX)                               
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                              
      DIMENSION SX(NN)                                                 
      DO 10 I = 1, NN                                                  
         SX(1+(I-1)*INCX) = SA*SX(1+(I-1)*INCX)                        
10    CONTINUE                                                         
      RETURN                                                           
      END                                                              
