      IMPLICIT NONE
      DOUBLE PRECISION V,VARREL,EPS,A,B,H,PI
      INTEGER K,NITER,NUMARRELS,NPUNTS,MPUNTS,J

      PARAMETER(MPUNTS=420)
      parameter(npunts=34)

      DOUBLE PRECISION X(1:2*NPUNTS),FVEC(2*NPUNTS) 
      !No és necessari posar 2*NPUNTS. Podriem posar NPUNTS, però no menys.
      !NPUNTS és el tamany mínim que necessitem per fer l'exercici
       
      DOUBLE PRECISION DFVEC(2*NPUNTS),DFVEC2(2*MPUNTS)
      DOUBLE PRECISION X2(2*MPUNTS),FVEC2(2*MPUNTS)
      DOUBLE PRECISION DFUN2(2*MPUNTS),DFUN(2*NPUNTS)

      double precision F,DF !SI AMB UBUNTU NO US COMPILA I 
      !WINDOWS SI, POT SER PERQUE LES FUNCIONS NO ESTAN DEFINIDES COM A 
      !DOUBLE PRECISION AL PROGRAMA PRINCIPAL. A MI, PER EXEMPLE, SENSE 
      !AQUEST double precision F,DF,fun EM COMPILA.
      EXTERNAL F,DF,FUN

      PI=datan(1.d0)*4.d0



C   2a) Càlcul de la funció P(v) i la seva derivada en v=[0,2*pi]

      OPEN (15,FILE="p3-1920-res.dat")
      WRITE (15,*)'#	valors de E, F(E), dF(E)'
      DO K=0,400
      	V=dble(K)*2d0*pi/400d0

      	WRITE (15,21) V,F(V),DF(V)
21    FORMAT (E20.12,2X,E20.12,2X,E20.12)
      ENDDO
      WRITE (15,*)''
      WRITE (15,*)''
  

C   2b) Càlcul de les 2 arrels de P(v) amb la subrutina de bisecció i escriptura en un fitxer.
      WRITE (15,*)"#	arrels de F(E)"
      EPS=10**(-12)

      NUMARRELS=2
      DO K=1,NUMARRELS
      	IF (K.EQ.1) THEN
      	!Els valors de A i B a continuacio els coneixem perque se suposa que
      	!hem realitzat una gràfica i a ull podem determinar-los. Ojo que convé
      	!posar "set yrange[  :  ]" al script de gnuplot (s.txt) amb valors adequats 
      	!per poder veure-la be, com ara [-10:10]
      		A=1D0
            B=2D0
        ENDIF
        IF (K.EQ.2) THEN
        	A=2D0
        	B=3D0
        ENDIF


      	CALL BISECTION(A,B,EPS,FUN,NITER,VARREL)


      	WRITE (15,*) NITER,VARREL
 37   	FORMAT (I4,2X,E20.12)
      ENDDO
      WRITE (15,*)''
      WRITE (15,*)''

C   2c) Càlcul de les arrels de F amb la subrutina de Newton-Raphson 
C      començant des de 10 punts diferents.
      WRITE (15,*) "#	Arrels de F(E): #iteracions, valor arrel "

      EPS=1.D-12

      OPEN(101,FILE="conv_NR.dat")

      CALL NEWTONRAPHSON (0.1D0,EPS,FUN,NITER,VARREL) 

      WRITE (15,60) NITER,VARREL
60    FORMAT (I8,2X,E20.12)
      CALL NEWTONRAPHSON (0.2D0,EPS,FUN,NITER,VARREL)
      WRITE (15,60) NITER,VARREL
      CALL NEWTONRAPHSON (0.67D0,EPS,FUN,NITER,VARREL)
      WRITE (15,60) NITER,VARREL
      CALL NEWTONRAPHSON (0.7D0,EPS,FUN,NITER,VARREL)
      WRITE (15,60) NITER,VARREL
      CALL NEWTONRAPHSON (1.5D0,EPS,FUN,NITER,VARREL)
      WRITE (15,60) NITER,VARREL
      CALL NEWTONRAPHSON (2.5D0,EPS,FUN,NITER,VARREL)
      WRITE (15,60) NITER,VARREL
      CALL NEWTONRAPHSON (2.6D0,EPS,FUN,NITER,VARREL)
      WRITE (15,60) NITER,VARREL
      CALL NEWTONRAPHSON (4.0D0,EPS,FUN,NITER,VARREL)
      WRITE (15,60) NITER,VARREL
      CALL NEWTONRAPHSON (5.4D0,EPS,FUN,NITER,VARREL)
      WRITE (15,60) NITER,VARREL

      CLOSE(101)


C   4) Càlcul numèric de la derivada de F amb la subrutina 'DERIFUN'.
c        Cas 34 punts:____________________________________________

      OPEN (18,FILE="P3-1920-res3_n34.dat")
      A=0.D0
      B=2.D0*PI
      H=(B-A)/(dble(NPUNTS-1))
      DO J=0,(NPUNTS-1)
      	X(J+1)=A+dble(J)*H
      	FVEC(J+1)=F(X(J+1))
      	DFVEC(J+1)=DF(X(J+1)) !DERIVADA EXACTE
      ENDDO

      CALL DERIFUN (NPUNTS,X,FVEC,DFUN)

      DO J=1,NPUNTS
      	WRITE (18,90) X(J),FVEC(J),DFUN(J),DFVEC(J)
 90     FORMAT (4(E20.12,2X))
      ENDDO
      CLOSE(18)

c        Cas 420 punts:____________________________________________
      OPEN (19,FILE="P3-1920-res3_n420.dat")
      A=0.D0
      B=2.D0*PI
      H=(B-A)/(REAL(MPUNTS-1))
      DO J=0,(MPUNTS-1)
      	X2(J+1)=A+J*H
      	FVEC2(J+1)=F(X2(J+1))
      	DFVEC2(J+1)=DF(X2(J+1))
      ENDDO

      	CALL DERIFUN (MPUNTS,X2,FVEC2,DFUN2)

      DO J=1,(MPUNTS)
      	WRITE (19,91) X2(J),FVEC2(J),DFUN2(J),DFVEC2(J)
91      FORMAT (4(E20.12,2X))
      ENDDO
      CLOSE(19)
      CLOSE(15)

      END !podem posar END PROGRAM si volem


C   Funció externa que retorna el valor del polinomi considerat P(v), per poder-lo fer servir com a argument en les subrutines d'interpolació.

      DOUBLE PRECISION FUNCTION F(V)
      	IMPLICIT NONE
      	DOUBLE PRECISION P,V
      	P=35D0/16D0+0.5D0*V-(61D0/20D0)*V**2+V**3
        F=P*DSINh(V)

      END !PODEM POSAR END FUNCTION SI VOLEM

C   Funció externa que retorna el valor de la derivada, per poder-la fer servir com a argument en la subrutina de Newton-Raphson.

      DOUBLE PRECISION FUNCTION DF(V)
      	IMPLICIT NONE
      	DOUBLE PRECISION DP,V,P
      	P=35D0/16D0+0.5D0*V-(61D0/20D0)*V**2+V**3
      	DP=0.5D0-2*(61D0/20D0)*V+3*V**2
      	DF=DP*DSINh(V)+P*DCOSh(V)

      END FUNCTION

c   Subroutina que retorna en FU, DFU els valors de les dues funcions externes
c   F i DF evaluades en l'input X
      SUBROUTINE FUN(X,FU,DFU) 
	      IMPLICIT NONE
	      DOUBLE PRECISION X,F,DF,FU,DFU
	      external F,DF
	      FU=F(X)
	      DFU=DF(X)
      END SUBROUTINE


C   Subrutina que calcula zeros de la funció 'FUN' amb derivada 'DFUN' donat un punt inicial 'X0' amb una precisió 'EPS'.
C   Retorna el número d'iteracions necessàries 'NITER' i el valor del zero 'XARREL'

      SUBROUTINE NEWTONRAPHSON (X0,EPS,FUN,NITER,XARREL)
	      IMPLICIT NONE
	      EXTERNAL FUN
	      DOUBLE PRECISION X0,EPS,XARREL,X1,DELTA,FX0,DFX0,PUNT
	      INTEGER NITER,MAXITER,IERR,I
	      MAXITER=10**5
 

C--------------- PART només PER ESCRIURE LA CONVERGENCIA COM DEMANEN A C)
	      IF (((X0.GT.0.1999D0).AND.(X0.LT.0.2001D0)).OR.
     + ((X0.GT.0.6999D0).AND.(X0.LT.0.70001D0)).OR.
     + ((X0.GT.1.4999D0).AND.(X0.LT.1.50001D0))) THEN !!!!!!!!!!!!

		      WRITE(101,*) '# X0 =',X0
		      PUNT=X0
		      niter=1
		      delta=eps+1d0
		      DO while ((DELTA.gt.EPS).and.(niter.lt.MAXITER))
		      	CALL FUN(PUNT,FX0,DFX0)

		      	X1=PUNT-(FX0/DFX0)

		      	WRITE(101,*) niter,X1

		      	DELTA=ABS(FX0/DFX0)
		      	PUNT=X1
		      niter=niter+1
		      ENDDO
		      xarrel=x1
		      WRITE(101,*) ''
		      WRITE(101,*) ''

C--------------------------------------------------------------
	      ELSE

		      PUNT=X0
		      niter=1
		      delta=eps+1d0      
		      DO while ((DELTA.gt.EPS).and.(niter.lt.MAXITER))
		      	CALL FUN(PUNT,FX0,DFX0)
		      	X1=PUNT-(FX0/DFX0) !TEORIA
		      	DELTA=ABS(FX0/DFX0) !TEORIA
		      	PUNT=X1
		      	niter=niter+1
		      ENDDO
		      	xarrel=x1
      	ENDIF
      END



C   Subrutina que calcula zeros de la funció 'FUN' en l'interval (A,B) amb una precisió 'EPS'. 
C   Retorna el número d'iteracions necessàries per obtenir la precisió desitjada 'NITER', 
C   el valor del zero 'XARREL'. És necesari que només hi hagi un zero en [A,B]
c   per garantir que funcioni sempre i que B >= A.

      SUBROUTINE BISECTION (A,B,EPS,FUN,NITER,XARREL)
	      IMPLICIT NONE
	      EXTERNAL FUN
	      DOUBLE PRECISION A,B,EPS,XARREL
	      DOUBLE PRECISION C,FA,FB,FC
	      INTEGER NITER,I,MAXITER,acabat
	      DOUBLE PRECISION NADA !Variable que no s'utilitza per res pero que
	      !necesitem per posar tots els arguments al cridar a FUN
	      MAXITER=10**9
	      acabat=0
	      xarrel=0
	      niter=1
	 !__________________________________________________
	      DO while ((acabat.eq.0).and.(niter.lt.MAXITER))

	      	C=(A+B)/2.D0
	      	CALL FUN(A,FA,NADA)
	      	CALL FUN(B,FB,NADA)
	      	CALL FUN(C,FC,NADA)
	      	

	      	IF ((FA*FB).GE.(0.D0)) THEN !Aquesta condició és per evitar errors
	      	! moooolt improvables quan FA=0.0000000000000 o FB=0.000000000000000

	      		IF (FA.EQ.0d0) THEN
	      			XARREL = A
	      			ACABAT=1


	      		ELSE IF (FB.EQ.0d0) THEN
	      			XARREL = B
	      			ACABAT =1

	      		ELSE

	      		acabat = -1
	      		endif

	      	else IF ((FA*FC).LT.(0.D0)) THEN !Bisecció de l'interval
	      		B=C
	      	ELSE
	      		A=C
	      ENDIF

	      	IF ((B-A).LE.EPS) THEN
	      		XARREL = C
	      		acabat=1

	      	endif
	      	niter=niter+1
	      ENDDO !__________________________________________________

	      IF (acabat.eq.(-1)) then
	      print*, 'la biseccio no ha trobat en A = ', A 
	      print*,'B = ', B

	      ELSE if (acabat.eq.0) then
	      	XARREL=C 
	    print*,'El valor no ha convergit amb la epsilon que es demana'
	      endif

      END


C   Subrutina que calcula la derivada primera d'una funció dins de l'interval [A,B].
c   Els inputs son dos vectors X,FU que corresponen a valors d'x i a les seves
c   imatges f(x). Ndat és el tamany que necessitem del vector. 
c   Retorna el vector DFU que és la derivada de FU.

      SUBROUTINE DERIFUN (NDAT,X,FU,DFU)
	      IMPLICIT NONE
	      INTEGER NDAT,K
	      DOUBLE PRECISION FU(NDAT),X(NDAT),DFU(NDAT),H

	      DO K=1,NDAT
	      	if (k.ne.ndat) THEN
	      	h=x(k+1)-x(k)
	      	endif

	      	IF (K.EQ.1) THEN
	      		DFU(K)=(FU(K+1)-FU(K))/H
	      	ELSE IF (K.EQ.NDAT) THEN
	      		DFU(K)=(FU(K)-FU(K-1))/H
	      	ELSE
	      		DFU(K)=(FU(K+1)-FU(K-1))/(2.D0*H)
	      	ENDIF

	      ENDDO
      END