
	  PROGRAM PrePractica6
	  IMPLICIT NONE
	  double precision F1,F2,F3,F4,F5,G,fun
	  External F1,F2,F3,F4,F5,G,fun
	  CALL SRAND(17777777)
	  OPEN(15,FILE="res.dat")
      call MontecarloP6
      call MultidMcP6

      CLOSE(15)

	  END

C 	APARTAT 1) Mètode de Montecarlo cru
	
	  DOUBLE PRECISION FUNCTION F1(x)
	  	IMPLICIT NONE
	  	double precision x,pi
	  	pi=dacos(-1d0)
		F1=DSQRT(pi**2+x**2)	
	 	RETURN
	  END

	  DOUBLE PRECISION FUNCTION F2(x)
	 	IMPLICIT NONE
	 	double precision x
		F2= (x+3*x**2*dsin(x)-x**3)*(dcos(x))**2*dsin(x)
	 	RETURN
	  END 

	  DOUBLE PRECISION FUNCTION F3(x)
	  	IMPLICIT NONE
	  	double precision x,pi
	  	pi=dacos(-1.d0)
	  	F3=4d0/5d0*x**2*(1d0-dexp(-1d0*pi))
	  	RETURN
	  END

	  DOUBLE PRECISION FUNCTION F4(x)
	  	IMPLICIT NONE
	  	double precision x,pi
	  	pi =dacos(-1d0)
	  	F4=4d0/5d0*dexp(dabs(x))*(dsin(x))**(-2d0)*(1d0-dexp(-1d0*pi))*
     + dexp(-x**2/2d0)*(dcos(x))**2d0*(1d0+x**2)
	  	RETURN
	  END

	  DOUBLE PRECISION FUNCTION F5(x)
	  	IMPLICIT NONE
	  	double precision x,pi
	  	pi=dacos(-1d0)
	  	F5=dsqrt(2d0*pi)*dexp(-(x**2)/2d0)*(dsin(x))**4.d0*x**2
	  	RETURN
	  END

	  

	  SUBROUTINE MontecarloP6

	  IMPLICIT NONE
	  double precision a1,a2,b1,b2,pi,a,b,m,x,Nr,e,h
	  double precision sumaF1,suma2F1,sumaF2,suma2F2
	  double precision sumaF3,suma2F3,sumaF4,suma2F4,sumaF5,suma2F5
	  double precision I1,I2,errorI1,errorI2,I1exacte,I2exacte
	  double precision I3,I4,I5,errorI3,errorI4,errorI5
	  double precision F1,F2,F3,F4,F5,fun
	  EXTERNAL F1,F2,F3,F4,F5,fun
	  INTEGER*4 N,ndat,iseed,k
	  PARAMETER(ndat=1050000)
	  double precision xgauss(ndat),xnums(ndat)
      COMMON/dades/xgauss
	  pi=dacos(-1.d0)
      e=dexp(1d0)
      m=0d0
      do k=1,100000
      	h=(2d0*pi)*k/100000
      	x=fun(-1d0*pi+h)
      	if (x.gt.m) m=x
      	ENDDO

C   1a) Càlcul de I1 i I2. Escriurem:  N  I1  σI1  I2   σI2 
	  a1=-e
	  b1=e
	  a2=-1d0*pi
	  b2=pi
	  
	  I1exacte=e*dsqrt(e**2+pi**2)+pi**2*dlog(e/pi+dsqrt((e/pi)**2+1))
	  I2exacte=1061d0*pi/288d0-5d0*pi**3/12

	  sumaF1=0.d0
	  suma2F1=0.d0
	  sumaF2=0.d0
	  suma2F1=0.d0

	  WRITE (15,*)
	  WRITE (15,*)
	  WRITE (15,*) "#APARTAT 1a"
	  WRITE (15,*) "#Numero iteracio, Valor integral, Error integral"
      WRITE (15,*) 'Valor real I1',I1exacte,'Valor real I2',I2exacte
	  DO N=1,ndat     
		x=rand()	
		sumaF1=sumaF1+(b1-a1)*F1((b1-a1)*x+a1)
		suma2F1=suma2F1+((b1-a1)*F1((b1-a1)*x+a1))**2.d0
		sumaF2=sumaF2+(b2-a2)*F2((b2-a2)*x+a2)
		suma2F2=suma2F2+((b2-a2)*F2((b2-a2)*x+a2))**2.d0
		IF(mod(N,2500).EQ.0) THEN
			Nr=real(N)
			I1=sumaF1/Nr
	        I2=sumaF2/Nr
		 	errorI1=(1.d0/DSQRT(Nr))*DSQRT(suma2F1/Nr-I1**2.d0)
		 	errorI2=(1.d0/DSQRT(Nr))*DSQRT(suma2F2/Nr-I2**2.d0)
		 	WRITE(15,30) N,I1-I1exacte,errorI1,I2-I2exacte,errorI2
		ENDIF
	  ENDDO
 30   FORMAT (I10,2x,6(E20.14,2x))

C   1b) i 1c)

c	    iseed=17777777
c	    CALL SRAND(iseed)
	    a=-1.d0*pi
	    b=pi
	    CALL SubGauss(ndat,xgauss)
      print*,'La cota superior de que utilitza subAiR és',m
	    CALL Subair(ndat,xnums,fun,a,b,m)


C   1d) Càlcul de I3,I4 i I5

        WRITE(15,*)
        WRITE(15,*)
        WRITE(15,*) "#APARTAT 1d"
        WRITE(15,*) "#Numero iteracio, Valor integral, Error integral"
        sumaF3=0.d0
        suma2F3=0.d0
        sumaF4=0.d0
        suma2F4=0.d0
        sumaF5=0.d0
        suma2F5=0.d0

        DO N=1,ndat
          sumaF3=sumaF3+F3(xnums(N))
          suma2F3=suma2F3+(F3(xnums(N)))**2.d0
          sumaF4=sumaF4+F4(xnums(N))
          suma2F4=suma2F4+(F4(xnums(N)))**2.d0
          sumaF5=sumaF5+F5(xgauss(N))
          suma2F5=suma2F5+(F5(xgauss(N)))**2.d0
          IF (mod(N,5000).EQ.0) THEN
          	Nr=dble(N)
            I3=sumaF3/Nr
            I4=sumaF4/Nr
            I5=sumaF5/Nr
          	errorI3=(1.d0/DSQRT(Nr))*(DSQRT((suma2F3/Nr)-(I3**2.d0)))
            errorI4=(1.d0/DSQRT(Nr))*(DSQRT((suma2F4/Nr)-(I4**2.d0)))
            errorI5=(1.d0/DSQRT(Nr))*(DSQRT((suma2F5/Nr)-(I5**2.d0)))
            WRITE(15,*) N,I3,errorI3,I4,errorI4,I5,errorI5,'quo',
     * errorI3/I3,errorI4/I4,errorI5/I5
          ENDIF
        ENDDO

	  RETURN
	  END SUBROUTINE


C   APARTAT 2)
C funció de l'integral de l'apartat 2, que utilitzarem a la subrutina multidmcP6

      DOUBLE PRECISION FUNCTION G(x1,x2,x3,x4,x5)
      IMPLICIT NONE 
      double precision x1,x2,x3,x4,x5,pi
      pi=dacos(-1d0)
      G=dexp(x1*dcos(x2+x3))*(x3**2*x4**2*x5**2+(dcos(x3+x4))**2*
     + x5*dsin(x5))*dexp(-(x1**2+x2**2+x3**2+x4**2+x5**2)/2)*(dsqrt(2d0
     + *pi))**5
      RETURN
      END




C   Subrutina del Mètode de Montecarlo multidimensional

	  SUBROUTINE MultidMcP6

        IMPLICIT NONE
C ndat2 es el nombre de particules que volem agafar de ndat totals
        INTEGER*4 ndat2,N,numero,ndat
        double precision sumaF6,suma2F6,G,I6,errorI6,Nr
        double precision x1,x2,x3,x4,x5
        parameter(ndat2=210000)
        parameter(ndat=1000000)
        double precision xgauss(ndat)
        External G
        COMMON/dades/xgauss
        WRITE(15,*)
        WRITE(15,*)
        WRITE(15,*) "#APARTAT 2"
        WRITE(15,*) "#Numero iteracio, Valor integral, Error integral"

        sumaF6=0.d0
        suma2F6=0.d0
        numero=1
        DO N=1,ndat2
          x1=xgauss(numero)
          numero=numero+1
          x2=xgauss(numero)
          numero=numero+1
          x3=xgauss(numero)
          numero=numero+1
          x4=xgauss(numero)
          numero=numero+1
          x5=xgauss(numero)
          numero=numero+1
          sumaF6=sumaF6+G(x1,x2,x3,x4,x5)
          suma2F6=suma2F6+(G(x1,x2,x3,x4,x5))**2.d0
          IF (mod(N,1500).EQ.0) THEN
              Nr=dble(N)
              I6=sumaF6/Nr
              errorI6=(1.d0/DSQRT(Nr))*(DSQRT((suma2F6/Nr)-(I6**2.d0)))   
           	  WRITE(15,31) N,I6,errorI6
          ENDIF
        ENDDO
 31     FORMAT(I10,2X,2(e24.14,2X))       

        RETURN
        END SUBROUTINE





C   Funció que fem servir a SubAiR: A i R de Acceptacio i Rebuig
	
	  DOUBLE PRECISION FUNCTION fun(x)
	   IMPLICIT NONE
	   double precision x, pi
		  pi=dacos(-1.d0)
	     IF ((x.GE.(-1.d0*pi)).AND.(x.LE.(pi))) THEN
	  fun=(5d0/4d0)*dexp(-dabs(x))*((dsin(x))**2)/(1-dexp(-1d0*pi))
	     ENDIF
	     IF((x.GT.pi).OR.(x.LE.(-1.d0*pi))) THEN	
		  fun=0.d0
	     ENDIF	
	  RETURN
	  END 



C   Subrutina que genera 'NDAT' números gaussians de valor mitjà 0 i variança 1

	  SUBROUTINE SUBGAUSS(NDAT,XGAUSS)

	  	IMPLICIT NONE
	  	INTEGER*4 NDAT,ISEED,RESIDU,NCAIXES,K
	  	PARAMETER(NCAIXES=100)
	  	double precision XGAUSS(NDAT),SUMA,XX,PI,ALPHA
	  	PI=DACOS(-1.D0)

	  	SUMA=0.D0
	  	DO K=1,NDAT
	  		XX=RAND()
	  		ALPHA=RAND()
	  		XGAUSS(K)=((-2.D0*DLOG(1.D0-XX))**0.5D0)*DSIN(2*PI*ALPHA)
	  		SUMA=SUMA+XGAUSS(K)
	  	ENDDO

		RETURN
		END SUBROUTINE



C   Subrutina que genera nombres aleatoris distribuïts segons la funció externa 'FUN' en un interval donat.

	  SUBROUTINE SUBAIR(NDAT,XNUMS,FUN,A,B,M)

	  IMPLICIT NONE
	  INTEGER*4 NDAT,K,NCAIXES,ISEED
	  PARAMETER(NCAIXES=50)
	  double precision XNUMS(NDAT),FUN,A,B,M
	  double precision X1,X2,X,P,SUMA

c	  CALL SRAND(ISEED)
	  DO K=1,NDAT
 5	    X1=RAND()
	  	X2=RAND()
	  	X=(B-A)*X1+A
	  	P=M*X2
	  	IF(FUN(X).GE.P) THEN
	 	 	XNUMS(K)=X
	  	ELSE
	  		GOTO 5
	  	ENDIF
	  ENDDO
      END SUBROUTINE

