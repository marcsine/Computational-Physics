      IMPLICIT NONE
!Variables de l'EXERCICI 1:

      INTEGER NDAT,K,NBOX1,ISEED
      PARAMETER(NDAT=20000,NBOX1=30)
      DOUBLE PRECISION XHIS(NBOX1),VHIS(NBOX1),ERRHIS(NBOX1),
     + FUN,A,B,XNUMS(NDAT),XDATA(NDAT),M,BOXSIZE,PI
      EXTERNAL FUN

!Variables de l'Exercici 2:

      INTEGER NDAT2,NBOX2
      DOUBLE PRECISION DESVEST,MITJANA,VAR,SUMA
      PARAMETER(NDAT2=14000,NBOX2=120)
      DOUBLE PRECISION XEXPO(NDAT2),XHIS2(NBOX2),VHIS2(NBOX2),
     + ERRHIS2(NBOX2),BOXSIZE2

      double PRECISION xx(1000),yy(1000)

C______________________________________________________________________
C______________________________________________________________________



CCCCCCCCCCC          SUPER IMPORTANT!              CCCCCCCCCCCCCC
CCCCCCCCCCCCCC                             CCCCCCCCCCCCCCC
CCCCCCCCCCCCCC                           CCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCC                          CCCCCCCCCCCCCCCCCCCCCC

C SRAND(ISEED) NOMÉS SHA DE POSAR AL PROGRAMA PRINCIPAL I MAI S'HA
C DE TORNAR A MENCIONAR ENLLOC (NI AL PROGRAMA NI A LES SUBROUTINES)
C LA FUNCIÓ SRAND() ÉS UNA FUNCIO INTRINSECA DE FORTRAN, COM SQRT(), 
C SIN() I COS(), QUE GENERA UN NOMBRE ALEATORI EN (0,1) DISTRIBUIT 
c UNIFORMEMENT. SI AFEGIM L'ARGUMENT ISEED A SRAND, SRAND(ISEED),
C LLAVORS ELS NOMBRES ALEATORIS ES GENEREN D'UNA MANERA DETERMINADA A 
C PARTIR DE ISEED. AL DEFINIR X1=RAND(),X2=RAND(),...,Xn=RAND()
C AQUESTS NOMBRES X1,X2,...,Xn SON ALEATORIS ENTRE ELLS PERÒ ELS SEUS 
C VALORS SERAN ELS MATEIXOS SEMPRE, INDEPENDENTMENT DE QUANTES VEGADES O
C EN QUIN ORDINADOR COMPILEM I EXECUTEM EL PROGRAMA. AIXÒ AJUDA A VEURE
C PATRONS I TROBAR POSSIBLES ERRORS.

C PER TANT, SI CRIDEM A SRAND(ISEED) AL PROGRAMA PRINCIPAL I DEFINIM 
C X1=RAND(),X2=RAND(),X3=RAND(),... I DESPRÉS TORNEM A CRIDAR A 
C SRAND(ISEED) A UNA SUBROUTINE I DEFINIM Y1=RAND(),Y2=RAND(),Y3=RAND(),
C LLAVORS X1 = Y1, X2 = Y2, X3 = Y3,... PER TANT X1,Y1,X2,Y2,X3,Y3... 
C NI DE CONYA ESTAN DISTRIBUITS ALEATORIAMENT ENTRE ELLS!
C NOMÉS HEM DE CRIDAR SRAND(ISEED) UN SOL COP AL PRINCIPI DEL PROGRAMA.
      
      pi=DACOS(-1D0)


      ISEED=17777777
      CALL SRAND(ISEED)




      OPEN(1,file="res.dat")
      OPEN(2,file="his.dat")

!Ex 1--------------------------------------------------

       !cota superior aproximada de fun que ha sigut calculada
       !aproximadament
      M=0.8d0
      print*,M

      CALL acceptrebuig(NDAT,XNUMS,0D0,pi,M,FUN)
      CALL HISTOGRAMA(NDAT,XNUMS,0D0,pi,NBOX1,XHIS,VHIS,ERRHIS,BOXSIZE)

      DO K=1,NBOX1
      	WRITE (2,*) XHIS(K),VHIS(K),ERRHIS(K)
      ENDDO
      
      WRITE(2,*)''
      WRITE(2,*)''
!Ex 2 -----------------------------------------------------

      CALL SUBEXPO(NDAT2,pi,XEXPO)

c APARTAT B) Calculem la variancia,mitjana etc. i les escrivim

      SUMA=0.D0
      DO K=1,NDAT2
      	SUMA=SUMA+XEXPO(K)
      ENDDO
      MITJANA=SUMA/REAL(NDAT2)
      SUMA=0.D0
      DO K=1,NDAT2
      	SUMA=SUMA+(XEXPO(K)-MITJANA)**2.D0
      ENDDO
      VAR=SUMA/REAL(NDAT2)
      DESVEST=DSQRT(VAR)

 30   FORMAT(3(E20.14,2X))
      WRITE(1,*)"#Mijana,variància i desv. est. dels punts distribuits"
     + ,"#segons l'exponencial"
      WRITE(1,*) MITJANA,VAR,DESVEST

C APARTAT C) Fem l'histograma

      CALL HISTOGRAMA(NDAT2,XEXPO,0d0,3d0,NBOX2,XHIS2,VHIS2,ERRHIS2,
     + BOXSIZE2)
      DO K=1,NBOX2
      	WRITE (2,*) XHIS2(K),VHIS2(K),ERRHIS2(K)
      ENDDO
      CLOSE(1)
      CLOSE(2)
      END

c_______________________________________________________________________
c_______________________________________________________________________


       SUBROUTINE HISTOGRAMA(NDAT,XDATA,XA,XB,NBOX,XHIS,VHIS,ERRHIS
     + ,BOXSIZE)
       IMPLICIT NONE

C INPUT/OUTPUT VARIABLES

       INTEGER NDAT,NBOX
       DOUBLE PRECISION XDATA(NDAT),XA,XB
       DOUBLE PRECISION XHIS(NBOX),VHIS(NBOX),ERRHIS(NBOX)
       INTEGER IERR
C 
       INTEGER I,IBOX,ICOUNT
       DOUBLE PRECISION BOXSIZE

       IF (XA.GE.XB) THEN 
          IERR=1
          RETURN
       ENDIF
C BOX SIZE
         BOXSIZE=(XB-XA)/NBOX

C COUNTS NUMBER OF POINTS WITHIN THE INTERVAL XA,XB
       ICOUNT=0

C SETS ALL TO ZERO
       DO I=1,NBOX
          VHIS(I)=0
          ERRHIS(I)=0
       ENDDO

C WE RUN THROUGH THE DATASET
       DO I=1,NDAT
C CHECKS IF DATA LIES WITHIN XA,XB
         IF (XDATA(I).GE.XA.AND.XDATA(I).LE.XB) THEN 
            IBOX=INT((XDATA(I)-XA)/BOXSIZE)+1
C PUTS XB INTO THE LAST BOX, IF NEEDED
            IF (IBOX.EQ.NBOX+1) IBOX=NBOX 

            VHIS(IBOX)=VHIS(IBOX)+1
            ICOUNT=ICOUNT+1
         ENDIF
        ENDDO

        IF (ICOUNT.EQ.0) THEN 
           IERR=2
           RETURN
        ENDIF

        IERR=0
        PRINT*,"ACCEPTED:",ICOUNT,
     1  " OUT OF:",NDAT

        DO I=1,NBOX
c CENTRAL VALUE OF THE BAR
           XHIS(I)=XA+BOXSIZE/2.D0+(I-1)*BOXSIZE
C  ERROBAR, STANDARD DEVIATION OF CORRESPONDING BINOMIAL
           ERRHIS(I)=SQRT(VHIS(I)/ICOUNT*(1.D0-VHIS(I)/ICOUNT))
     1      /BOXSIZE / SQRT(DBLE(ICOUNT))
C NORMALIZED VALUE OF THE BAR
           VHIS(I)=VHIS(I)/ICOUNT/BOXSIZE
        ENDDO
        END
       



      SUBROUTINE SUBEXPO(NDAT,XLAMBDA,XEXPO)
      IMPLICIT NONE
      INTEGER NDAT,K
      DOUBLE PRECISION XEXPO(1:NDAT),XLAMBDA,Y

	  	DO K=1,NDAT
	  		Y=RAND()
	  		XEXPO(K)=(-1D0/XLAMBDA)*DLOG(1d0-Y)
	  ENDDO
	  END




	  DOUBLE PRECISION FUNCTION FUN(X)
	  IMPLICIT NONE
	  DOUBLE PRECISION X,PI
	  PI=4*DATAN(1D0)
	  FUN=(12/(PI*(2*PI**2-3)))*X**2*(DSIN(X))**2
	  END 


C Genera nombres aleatoris distribuïts segons la funció externa 'FUN'
c en l'interval A B pel mètode d'acceptacio i rebuig.

      SUBROUTINE acceptrebuig(NDAT,XNUMS,A,B,M,FUN)
	  IMPLICIT NONE
	  INTEGER NDAT,K,NBOX
	  PARAMETER(NBOX=30)
	  DOUBLE PRECISION XNUMS(NDAT),FUN,A,B,M
	  DOUBLE PRECISION X1,X2,X,P,SUMA,VAR,DESVEST,MITJANA
	  DOUBLE PRECISION XHIS(NBOX),VHIS(NBOX)
	  DOUBLE PRECISION ERRHIS(NBOX)
	  DO K=1,NDAT
 5	    X1=RAND()
	  	X2=RAND()
	  	X=(B-A)*X1+A
	  	P=M*X2
      

      IF (FUN(X).Ge.P) THEN
      XNUMS(K)=X

      ELSE

      GOTO 5
      ENDIF

	  ENDDO


c A partir d'aqui el mètode ha acabat. Com ens demanen al guió fem que
c la subroutine ens calculi mitjana etc. i ho escrigui al fitxer de 
c resultats. 

	  SUMA=0.D0
	  DO K=1,NDAT
	  	SUMA=SUMA+XNUMS(K)
	  ENDDO
	  MITJANA=SUMA/REAL(NDAT)
	  SUMA=0.D0
	  DO K=1,NDAT
	  	SUMA=SUMA+(XNUMS(K)-MITJANA)**2.D0
	  ENDDO
	  VAR=SUMA/REAL(NDAT)
	  DESVEST=DSQRT(VAR)

 30   FORMAT(3(E20.14,2X))
      WRITE(1,*) '#Mitjana, var. i desv. est. dels punts distribuïts'
     + ,"#segons la funció d'entrada"
	  WRITE(1,*) MITJANA,VAR,DESVEST
	  WRITE(1,*)''
	  WRITE(1,*)''
      END