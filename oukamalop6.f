C            ***********************************************************
C               STUDENT NAME - ***
C               EMAIL - ***
C               Electrical and Computer Engineering
C            ***********************************************************
C               PURPOSE - USE GUASSIAN ELIMINATION TO SOLVE 2*NUM+2
C               SYSTEM TO DETERMINE THE DISTANCES,EXTENSIONS AND
C               TENSIONS FOR A FRAMEWORK
C            ***********************************************************
C               DATA FILES
C               cp60950.txt -- input data file
C               oukamalop6.out -- output data file to which
C                 results are written
C            ***********************************************************
      PROGRAM MAIN
C     ******************************************************************
C        DECLARE NECESSARY VARIALBES:
C     ******************************************************************
      INTEGER NUM,I,M,N,FLAG,INPUTSET
      REAL F1,F2,LENGTH(10),AREA(10),ANGLE(10),COEF(22,22),RHS(22)
      REAL AUGM(22,23),SCLF(22),SOL(22)
      OPEN(UNIT=5,FILE='cp60950.txt',STATUS='OLD')
      OPEN(UNIT=6,FILE='oukamalop6.out')
C                  ****************************************************
C                    INPUT EACH DATA SET
C                  ****************************************************
      WRITE(6,*)'ECE 3331 COMPUTER PROBLEM 6'
      WRITE(6,*)'---------------------------'
      WRITE(6,*)
      INPUTSET=0
30    READ(5,*,END=99)NUM,F1,F2
      INPUTSET=INPUTSET+1
      WRITE(6,*)'Framework # ',INPUTSET
      WRITE (*,*)
      WRITE(6,*)'Applied force in horizontal direction(N) =',F1
      WRITE(6,*)'Applied force in vertical direction(N)   =',F2

C                  ****************************************************
C                  WRITE HEADINGS FOR OUTPUT
C                  ****************************************************
      WRITE(6,*)
      WRITE(6,*)' Pin  Angle(degrees)  Length(m)      Area(m^2)'
      WRITE(6,*)'---------------------------------------------------'

C                  ****************************************************
C                     ECHO CHECK INPUT PARAMETERS
C                  ****************************************************
      I=0
60    I=I+1
      READ(5,*)LENGTH(I),AREA(I),ANGLE(I)
      WRITE(6,*)' ',I,'   ',ANGLE(I),'        ',LENGTH(I),'  ',AREA(I)
      IF(I.LT.NUM)go to 60
      WRITE(6,*)
      CALL BUILD(NUM,LENGTH,AREA,ANGLE,F1,F2,COEF,RHS)
      WRITE(6,*)'Coefficient matrix & RHS :'
      WRITE(6,*)'Order of Unknowns is tensions, then extensions'
      WRITE(6,*)' then distances and last is RHS'
      DO 33 M=1,(num+1)*2
      WRITE(6,*)'EQ',M,':',(COEF(M,N),N=1,NUM*2+2),RHS(M)
      WRITE(6,*)
33    CONTINUE
      CALL DUMGAS(NUM*2+2,COEF,RHS,SOL,FLAG,0,22,AUGM,SCLF)
      DO 22 M=1,NUM
      WRITE(6,*)'Tension of member(N)',M,'= ',SOL(M)
22    CONTINUE
      WRITE(6,*)
      DO 44 M=1,NUM
      WRITE(6,*)'Extension of member(m)',M,'= ',SOL(NUM+M)
44    CONTINUE
      WRITE(6,*)
      WRITE(6,*)'Distance moved horizontally(m) =',SOL(NUM*2+1)
      WRITE(6,*)'Distance moved vertically(m)   =',SOL(NUM*2+2)
      WRITE(6,*)
      GO TO 30
99    END

      REAL FUNCTION flex(len,area, E)
      REAL len,area,E
      flex=len/(area*E)
      RETURN
      END

      SUBROUTINE build(NUM,L,A,THETA,F1,F2,COEF,RHS)
      REAL  F1,F2
      REAL  COEF(22,22),RHS(22),L(10),A(10),THETA(10)
      INTEGER NUM,J

      DO 10 J=1,2*(NUM+1)
      RHS(J)=0
      DO 13 K=1,2*(NUM+1)
      COEF(J,K)=0
13    CONTINUE
10    CONTINUE
      RHS(1)=F1
      RHS(2)=F2
      I=0
67    I=I+1

      COEF(1,i)=-cos(THETA(i)*3.1415927/180)
      COEF(2,i)=sin(THETA(i)*3.1415927/180)

      COEF(i+2,num+i)=1
      COEF(i+2,num*2+1)=cos(THETA(i)*3.1415927/180)
      COEF(i+2,num*2+2)=-sin(THETA(i)*3.1415927/180)

      COEF(2+num+i,i)=-flex(L(I),A(I),20.4E10)
      COEF(2+num+i,num+i)=1

      IF(I.LT.NUM)go to 67

      RETURN
      END
C
      SUBROUTINE DUMGAS(N,A,B,X,IER,IWRT,LDA,UL,S)
C     ***********************************************************
C  PURPOSE:  TO SOLVE A SYSTEM OF N LINEAR EQUATIONS IN N UNKNOWNS
C      USING GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING AND ROW SCALING
C  IDENTIFICATION OF PARAMETERS:
C      N - INPUT SCALAR, NUMBER OF EQUATIONS = NUMBER OF UNKNOWNS
C      A - INPUT DOUBLY SUBSCRIPTED REAL ARRAY, COEFFICIENT MATRIX
C      B - INPUT SINGLY SUBSCRIPTED REAL ARRAY, RIGHT HAND SIDE
C      X - OUTPUT SINGLY SUBSCRIPTED REAL ARRAY, SOLUTION VECTOR
C   IER - OUTPUT INTEGER USED AS FLAG:
C          IER=0 MEANS NO UNIQUE SOLUTION
C          IER=1 MEANS UNIQUE SOLUTION FOUND
C    IWRT - INPUT INTEGER USED AS A FLAG FOR WRITING
C         IWRT=0 MEANS ALL WRITING IS SUPPRESSED
C         IWRT = 1 MEANS INTERMEDIATE REDUCED MATRICES ARE WRITTEN
C   LDA - ROW DIMENSION OF A EXACTLY AS DECLARED IN THE CALLING PROGRAM
C      UL - AUGMENTED MATRIX IN WHICH ALL REDUCTION WILL TAKE PLACE
C      S - SCALE VECTOR
C     **************************************************************
        REAL A(LDA,LDA),B(LDA),X(LDA),UL(LDA,LDA+1),S(LDA)

C     *****************************************************************
       IF(IWRT.EQ.1)THEN
C     *****************************************************************
C    THIS SUBROUTINE CONTAINS WRITE STATEMENTS SO THAT INTERMEDIATE STEPS
C      IN THE REDUCTION CAN BE WRITTEN IF THE WRITE FLAG IS SET TO 1
C    THIS IS PARTICULARLY HELPFUL IN FINDING ERRORS
C     ******************************************************************
        WRITE(6,200)
200     FORMAT(' ','ORIGINAL AUGMENTED MATRIX')
        DO 1 I = 1,N
        WRITE(6,201)I,(A(I,J),J=1,N),B(I)
201     FORMAT(' Row ', I2,(1x,10(E11.4,1x)))
1       CONTINUE
        ENDIF
C     **********************************************************************
C      ZERO WILL BE USED AS A CRITERION FOR SAYING A REAL NUMBER IS 0; ITS
C      SIZE MAY HAVE TO BE ADJUSTED FOR SPECIFIC DATA
C      *********************************************************************
        zero = 1.E-5
C        *******************************************************
C        initialize flag for solution & set up augmented matrix
C        ********************************************************
        IER = 1
        DO 5 I = 1,N
        DO 4 J = 1,N
        UL(I,J) = A(I,J)
4       CONTINUE
C         ********************************************************
C         determine scale factors for each row
C         ******************************************************
        s(i)=abs(ul(i,1))
        do 3 j=2,n
        big=abs(ul(i,j))
        if(big.gt.s(i))s(i)=big
3       continue
C         **************************************************
C         test for all zero row in which case no unique
C         solution can exist since system is square
C         ***********************************************
        if(s(i).lt.zero)then
           ier=0
             return
         endif
         UL(I,N+1) = B(I)
5        CONTINUE
C             ********************************************
C             begin elimination loop
C             ************************************************
        DO 17 K = 1,N-1
C         ***********************************************
C         find pivot row for column K
C         ***********************************************
        BIG = 0.0
        DO 11 I = K,N
        SIZE = ABS(UL(I,K))/s(i)
        IF(SIZE.LE.BIG)GO TO 11
        BIG = SIZE
        IDXPIV = I
11      CONTINUE
        IF(BIG.LT.ZERO)then
C            **************************************
C            no nonzero choice-cannot get a one
C             so no unique solution
C             **************************************
              IER = 0
        RETURN
        endif
C          ***************************************
C          see if row interchange is needed
C          ************************************
        IF(IDXPIV.ne.K)then
C          ***********************************
C          yes, interchange row k & row idxpiv
C          **************************************
        DO 7 J = 1,N+1
        TEP = UL(K,J)
        UL(K,J) = UL(IDXPIV,J)
7       UL(IDXPIV,J) = TEP
C          *************************************
C           also must interchange scale factors
C          *************************************
        tep=s(k)
        s(k)=s(idxpiv)
        s(idxpiv)=tep
        IF(IWRT.EQ.1)then
        WRITE(6,317)
317     FORMAT(' ','NEW MATRIX')
        DO 50 I = 1,N
        WRITE(6,201)I,(UL(I,J),J=1,N+1)
50      CONTINUE
        endif
        endif
C        *****************************************
C       now eliminate below
C       ****************************************
        DO 16 I = K+1,N
        P = -UL(I,K)/UL(K,K)
        DO 16 J = K,N+1
        UL(I,J) = UL(I,J) +P*UL(K,J)
16      CONTINUE
        IF(IWRT.EQ.1)then
        WRITE(6,317)
        DO 703 IR =1,N
        WRITE(6,201)IR,(UL(IR,J),J=1,N+1)
703     CONTINUE
        endif
17      continue
C        ********************************************
C        elimination is over,check for last pivot being zero
C        *******************************************
        IF(ABS(UL(N,N)).LT.ZERO)then
        ier = 0
        return
         endif
C        *****************************************
C         begin back substitution
C        ******************************************
        X(N) = UL(N,N+1)/UL(N,N)
        DO 8 IBACK = 2,N
        M = N+1-IBACK
        BACSOL = UL(M,N+1)
        DO 9 JBACK = M+1,N
9       BACSOL = BACSOL-UL(M,JBACK)*X(JBACK)
8       X(M) = BACSOL/UL(M,M)
        IER = 1
        RETURN
        END


