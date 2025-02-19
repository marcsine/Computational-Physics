      IMPLICIT NONE
      double precision theta0,dtheta0,pi,longitud
     + ,massa,g,OmegaN,TN,Ecine,Epoten
      INTEGER k,pas
      EXTERNAL Ecine,Epoten
      common/dades_g_l_m/g,longitud,massa
      g = 8.87d0
      longitud = 0.45d0
      massa = 0.510d0


      pi = DACOS(-1.d0)
      OmegaN = DSQRT(g/longitud)
      TN = (2.d0*pi)/OmegaN

      OPEN(15,FILE="res.dat")

C   Apartat a) OScil·lacions grans. Cridem les dues subrutines.

      theta0=pi-0.02d0
      dtheta0=0.d0

      WRITE(15,*)
      WRITE(15,*)
      WRITE(15,*) "#Apartat a) Oscil·lacions grans"
      WRITE(15,*) "#Mètode d'Euler: Temps, Theta, dTheta"
      CALL Euler(theta0,dtheta0,0.d0,5*TN,1800,'n')

      WRITE(15,*)
      WRITE(15,*)
      WRITE(15,*) "#Mètode d'Euler millorat: Temps, Theta, dTheta"
      CALL EulerM(theta0,dtheta0,0.d0,5*TN,1800,'n')

C   Apartat b) Petites oscil·lacions. Cridem les dues subrutines.

      theta0 = 0.04d0
      dtheta0 = 0.d0

      WRITE(15,*)
      WRITE(15,*)
      WRITE(15,*) "#Apartat b) Petites oscil·lacions"
      WRITE(15,*) "#Mètode d'Euler: Temps, Theta, dTheta"
      CALL Euler(theta0,dtheta0,0.d0,5*TN,1300,'s')

      WRITE(15,*)
      WRITE(15,*)
      WRITE(15,*) "#Mètode d'Euler millorat: Temps, Theta, dTheta"
      CALL EulerM(theta0,dtheta0,0.d0,5*TN,1300,'s')

C   Apartat c) Energia. 

      theta0 = 1.d0
      dtheta0=0.d0

      WRITE(15,*)
      WRITE(15,*)
      WRITE(15,*) "#Apartat c) Energia. Theta0 = 1 rad"
      WRITE(15,*) "#Euler: Temps, Theta, dTheta, Ecine, Epoten, ETot"
      CALL Euler(theta0,dtheta0,0.d0,5*TN,2000,'n')

      WRITE(15,*)
      WRITE(15,*)
      WRITE(15,*) "#Euler millorat:Temps, Theta, dTheta,Ecine,"
     + ,"Epoten,ETot"
      CALL EulerM(theta0,dtheta0,0.d0,5*TN,2000,'n')

      theta0=pi-0.02d0
      dtheta0=0.d0

      WRITE(15,*)
      WRITE(15,*)
      WRITE(15,*) "#Theta0 = PI-0.02 rad"
      WRITE(15,*) "#Euler: Temps, Theta, dTheta, Ecine, Epoten, ETot"
      CALL Euler(theta0,dtheta0,0.d0,5*TN,2000,'n')

      WRITE(15,*)
      WRITE(15,*)
      WRITE(15,*) "#Euler millorat:Temps, Theta, dTheta,Ecine,Epoten,"
     + ,"ETot"
      CALL EulerM(theta0,dtheta0,0.d0,5*TN,2000,'n')


C   Apartat d) Transició.

      theta0=0.d0
      dtheta0=2.d0*DSQRT(g/longitud)+0.03d0

      WRITE(15,*)
      WRITE(15,*)
      WRITE(15,*) "#Apartat d) Transició."
      WRITE(15,*) "#Euler millorat: Temps, Theta, dTheta, Ecine,"
     + ,"Epoten,ETot"
      WRITE(15,*) "#dTheta0=2sqrt(g/l)+0.03 rad/s"
      CALL EulerM(theta0,dtheta0,0.d0,11*TN,2100,'n')

      theta0=0.d0
      dtheta0=2.d0*DSQRT(g/longitud)-0.03d0

      WRITE(15,*)
      WRITE(15,*)
      WRITE(15,*) "#Euler millorat:Temps, Theta, dTheta,Ecine,"
     + ,"Epoten,ETot"
      WRITE(15,*) "#dTheta0=2sqrt(g/l)-0.03 rad/s"
      CALL EulerM(theta0,dtheta0,0.d0,11*TN,2100,'n')


C   Apartat e) Convergència del mètode

      theta0=2.8d0
      dtheta0=0.d0

      WRITE(15,*)
      WRITE(15,*)
      WRITE(15,*) "#Apartat e) Convergència del mètode"
      WRITE(15,*) "#Euler millorat:Temps, Theta, dTheta,"
     + ,"Ecine,Epoten,ETot"
      WRITE(15,*) "#600 passos"
      CALL EulerM(theta0,dtheta0,0.d0,12*TN,600,'n')

      WRITE(15,*)
      WRITE(15,*)
      WRITE(15,*) "#Euler millorat:Temps, Theta, dTheta,"
     + , "Ecine,Epoten,ETot"
      WRITE(15,*) "#1300 passos"
      CALL EulerM(theta0,dtheta0,0.d0,12*TN,1300,'n')

      WRITE(15,*)
      WRITE(15,*)
      WRITE(15,*) "#Euler millorat:Temps, Theta, dTheta,Ecine,"
     + ,"Epoten,ETot"
      WRITE(15,*) "#2600 passos"
      CALL EulerM(theta0,dtheta0,0.d0,12*TN,2600,'n')

      WRITE(15,*)
      WRITE(15,*)
      WRITE(15,*)"#Euler millorat:Temps, Theta, dTheta,Ecine,"
     + ,"Epoten,ETot"
      WRITE(15,*) "#15000 passos"
      CALL EulerM(theta0,dtheta0,0.d0,12*TN,15000,'n')

      CLOSE(15)

      END



C   Subrutina que calcula l'angle 'theta1' i la velocitat angular 'dtheta1' 
c      d'un pèndol a mesura que avança el temps 'temps', pel mètode d'Euler, 
c      i escriu el resultat en un fitxer.

      SUBROUTINE Euler(theta0,dtheta0,ti,tf,npassos,aprox)
      IMPLICIT NONE
      double precision theta0,dtheta0,theta1,dtheta1,ti,tf
     + ,dt,temps,g,longitud
     + ,Ecine,Epoten,eK,eP,eT,In1,In2,massa
      INTEGER npassos,k
      character aprox
      EXTERNAL Ecine,Epoten
      common/dades_g_l_m/g,longitud,massa

       !Aquestes dues variables a continuacio nomes serveixen just al final
       ! de la subroutine. Com els valors a l'input detheta0 i 
       ! dtheta0 es canviaran durant l'algoritme, els guardem a In1 i In2.
       ! Al final de la subroutine els tornem a recuperar per fer que
       ! theta0 i dtheta0 valguien el que valien inicialment, i aixi
       !evitar errors no desitjats. Si no ho fessim hauriem de redefinir
       ! al programa principal els valors theta0 i dtheta0 cada cop que
       !volguem tornar a trucar a la subroutina.

      In1=theta0      
      In2=dtheta0

      dt=(tf-ti)/real(npassos)
      DO k=1,npassos
      	temps=k*dt+ti
      	theta1=theta0+dt*dtheta0
      	if (aprox.eq.'s') then

      		dtheta1=dtheta0-dt*(g/longitud)*(theta0)
      	else
      		dtheta1=dtheta0-dt*(g/longitud)*DSIN(theta0)
      	endif

      		eK=Ecine(dtheta1)
      		eP=Epoten(theta1)
      		eT=eK+eP
c      	IF (mod(k,100).EQ.0) THEN
      		WRITE(15,30) temps,theta1,dtheta1,eK,eP,eT
c      	ENDIF
      	theta0=theta1
      	dtheta0=dtheta1
      ENDDO
 30		FORMAT(6(E20.14,2x))
      theta0 = In1
      dtheta0 = In2
      END SUBROUTINE



      SUBROUTINE EulerM(theta0,dtheta0,ti,tf,npassos,a)
      IMPLICIT NONE
      double precision theta0,dtheta0,theta1,dtheta1,ti,
     + tf,theta2,dtheta2,dt,temps,g,longitud,
     + Ecine,Epoten,eK,eP,eT,In1,In2,massa

      INTEGER npassos,k
      character a
      EXTERNAL Ecine,Epoten
      common/dades_g_l_m/g,longitud,massa

       !Aquestes dues variables a continuacio nomes serveixen just al final
       ! de la subroutine. Com els valors a l'input detheta0 i 
       ! dtheta0 es canviaran durant l'algoritme, els guardem a In1 i In2.
       ! Al final de la subroutine els tornem a recuperar per fer que
       ! theta0 i dtheta0 valguien el que valien inicialment, i aixi
       !evitar errors no desitjats. Si no ho fessim hauriem de redefinir
       ! al programa principal els valors theta0 i dtheta0 cada cop que
       !volguem tornar a trucar a la subroutina.

      In1=theta0      
      In2=dtheta0

      dt=(tf-ti)/real(npassos)
            theta1=theta0+dt*dtheta0
      	if (a.eq.'s') then
      		dtheta1=dtheta0-dt*(g/longitud)*(theta0)
      	else
      		dtheta1=dtheta0-dt*(g/longitud)*DSIN(theta0)
      	endif

      DO k=2,npassos
      	temps=k*dt+ti
      	theta2=theta0+2.d0*dt*dtheta1
      	if (a.eq.'s') then
      		dtheta2=dtheta0-2.d0*dt*(g/longitud)*(theta1)
      	else
      		dtheta2=dtheta0-2.d0*dt*(g/longitud)*DSIN(theta1)
      	endif

      		eK = Ecine (dtheta2)
      		eP = Epoten (theta2)
      		eT=eK+eP
c      	IF (mod(k,100).EQ.0) THEN
      		WRITE(15,31) temps,theta2,dtheta2,eK,eP,eT
c      		ENDIF
      	theta0=theta1
      	dtheta0=dtheta1
      	theta1=theta2
      	dtheta1=dtheta2
      ENDDO
 31		FORMAT(6(E20.14,2x))
      theta0 = In1
      dtheta0 = In2
      END SUBROUTINE



      DOUBLE PRECISION FUNCTION Ecine (dtheta)
      IMPLICIT NONE
      double precision dtheta,massa,longitud,g
      common/dades_g_l_m/g,longitud,massa
      Ecine=(1.d0/2.d0)*massa*(dtheta**2.d0)*(longitud**2.d0)

      END FUNCTION


      DOUBLE PRECISION FUNCTION Epoten (theta)
      IMPLICIT NONE
      double precision theta,massa,longitud,g
      common/dades_g_l_m/g,longitud,massa
      Epoten=-1.d0*massa*g*longitud*DCOS(theta)
      
      END FUNCTION

