      IMPLICIT NONE
      integer icontrol,apartat
      double precision Tint

C   APARTAT 3: Convergència dels mètodes
     
      apartat=3 !Aquesta variable val 1 pel primer apartat i 2 pel segon
           ! La subroutine temperatura realitzarà el primer apartat o segon
           ! segons el valor que tingui "apartat" quan la truquem

      OPEN(15,FILE="res.dat")

C   Mètode de Gauss-Seidel "Viejo" per diferents 'Tint'
      
      icontrol=1
      WRITE(15,*)
      WRITE(15,*)
      WRITE(15,*) "#Mètode de Gauss-Seidel (viejo)"
      WRITE(15,*) "#Numero iteracio      Temperatura al punt"
      Tint=6d0
      CALL Temperatura(apartat,icontrol,Tint)
      WRITE(15,*)
      WRITE(15,*)
      Tint=19d0
      CALL Temperatura(apartat,icontrol,Tint)
      WRITE(15,*)
      WRITE(15,*)
      Tint=320d0
      CALL Temperatura(apartat,icontrol,Tint)
  

C   Mètode de Jacobi (nuevo) per diferents 'Tint'
  
      icontrol=2
      WRITE(15,*)
      WRITE(15,*)
      WRITE (15,*) "#Mètode de Jacobi (nuevo)"
      WRITE(15,*) "#Numero iteracio      Temperatura al punt"
      Tint=6d0
      CALL Temperatura(apartat,icontrol,Tint)
      WRITE(15,*)
      WRITE(15,*)
      Tint=19d0
      CALL Temperatura(apartat,icontrol,Tint)
      WRITE(15,*)
      WRITE(15,*)
      Tint=320d0
      CALL Temperatura(apartat,icontrol,Tint)
      WRITE(15,*)
      WRITE(15,*)

C   Mètode de sobrerelaxació per diferents 'Tint'
  
      icontrol=3
      WRITE(15,*)
      WRITE(15,*)
      WRITE (15,*) "#Mètode de sobrerelaxació"
      WRITE(15,*) "#Numero iteracio      Temperatura al punt"
      Tint=6d0
      CALL Temperatura(apartat,icontrol,Tint)
      WRITE(15,*)
      WRITE(15,*)
      Tint=19d0
      CALL Temperatura(apartat,icontrol,Tint)
      WRITE(15,*)
      WRITE(15,*)
      Tint=320d0
      CALL Temperatura(apartat,icontrol,Tint)
      WRITE(15,*)
      WRITE(15,*)

      CLOSE(15)


C   APARTAT 4: Càlcul del mapa de temperatures final

      apartat=4
      icontrol=1
      Tint=0d0
      OPEN(16,FILE="res2.dat") !Guardarem les dades del grafic 3d a res2.dat
      WRITE(16,*) "#Mapa final de temperatures"
      WRITE(16,*) "#Posició x       posició y      temperatura"
      CALL Temperatura(apartat,icontrol,Tint)
      CLOSE(16)
 
 	  END

c_______________________________________________________________________
c_______________________________________________________________________

c Aquesta funció calcula el valor de la densitat de temperatura segons
c l'enunciat

      DOUBLE PRECISION FUNCTION FuncioRho(x,y)  
	      IMPLICIT NONE
	      double precision x,y,rho1,rho2,r,rho10,rho20
	      double precision cx,cy,dx,dy,lx,ly !(cx,cy) es el punt central del rec
	      PARAMETER(rho10=1.33d0)
	      PARAMETER(rho20=1.3d0)
	      PARAMETER(cx=22d0)
	      PARAMETER(cy=14d0)
	      PARAMETER(lx=4d0)
	      PARAMETER(ly=2d0)
	      r=DSQRT(((x-7.d0)**2.d0)+((y-8.d0)**2.d0))
	      rho1=rho10*DEXP(-1.d0*((r-4d0)**2.d0)/(0.4d0**2.d0))
	      dx=dABS(cx-x)
	      dy=dABS(cy-y)
	      IF ((dx.LE.(lx/2.d0)).AND.(dy.LE.(ly/2.d0))) THEN
	      	rho2=rho20
	      ELSE
	      	rho2=0.d0
	      ENDIF
	      FuncioRho=rho1+rho2
      END FUNCTION

c_______________________________________________________________________
c_______________________________________________________________________

      SUBROUTINE Temperatura(apartat,icontrol,Tint)

	      IMPLICIT NONE
	      double precision h,lx,ly
	      double precision suma,DELTA,FuncioRho,rho,omega,Tint
	      double precision puntx,punty
	      integer nx,ny,nkmax,niter,icontrol,ix,iy,apartat
	      PARAMETER(h=0.5d0) ! distancia entre nodes de la malla
	      PARAMETER(lx=32.5d0) ! tamany del domini de x
	      PARAMETER(ly=16.5d0) ! tamany del domini de y
	      PARAMETER(nx=INT(lx/h)) ! nombre de nodes al llarg de x
	      PARAMETER(ny=INT(ly/h)) ! nombre de nodes al llarg de y
	      PARAMETER(nkmax=10000) ! nombre maxim d'iteracions
	      PARAMETER(omega=1.73d0) ! paràmetre omega de la sobrerelaxació successiva
	      PARAMETER(puntx=18d0) ! coordenada x del punt del que mirem la convergència
	      PARAMETER(punty=12.5d0) ! coordenada y del punt del que mirem la convergència
	      PARAMETER(ix=INT(puntx/h)) ! index de la coordenada x del punt
	      PARAMETER(iy=INT(punty/h)) ! index de la coordenada y del punt 
	      double precision tnew(0:nx,0:ny),told(0:nx,0:ny),error,tol
	      integer i,j,k



	      tol=0d0 ! En el nostre cas no hem de definir cap tolerancia.
	       ! La fem 0 per que no jugui cap paper

C     Inicialitzem la temperatura a 'Tint' a tots els punts que no són extrems
	      DO i=1,nx-1
	      	DO j=1,ny-1
	      		told(i,j)=Tint
	      	ENDDO
	      ENDDO
C   Imposem condicions de contorn de Dirichlet (Dirichlet equival a condicions
	       ! a la frontera del domini de la Funció Fi. Si
	       ! les condicions a la frontera fossin a la primera derivada de Fi serien
	       ! condicions de contorn de Neumann)

	      DO i=0,nx
	      	told(i,0)=3.36d0
	      	told(i,ny)=23.1d0
	      ENDDO
	      DO j=0,ny
	      	told(0,j)=4d0
	      	told(nx,j)=25d0
	      ENDDO
	      tnew=told
	      niter=1

	      error=tol+1d0

C   Programació dels 3 mètodes; Gauss-Seidel (VIEJO) per icontrol=1, Jacobi 
C (NUEVO) per icontrol=2 i sobrerelaxació per icontrol=3
c   El criteri per finalitzar el mètode es el que segueix:
c si l'error es mes petit o igual que la tolerancia o el nombre de iteracions
c arriba al límit nkmax el programa acaba. 

C Evidentment, si fixem tol=0d0, el programa acabarà només si es fan totes
c les iteracions fins nkmax, ja que mai es donarà error < tol .

	      do while((dabs(error).ge.tol).and.(niter.lt.nkmax)) 

	      	DO i=1,nx-1
	      	  DO j=1,ny-1
	      		rho=FuncioRho(h*i,h*j)

	      		IF (icontrol.EQ.1) THEN

	      			suma=told(i+1,j)+told(i-1,j)+told(i,j+1)+told(i,j-1)
	      			tnew(i,j)=(suma+rho*h**2.d0)/4.d0
	      			 ! nomes serveix si voleu saber
	      			 ! una estimacio del error del valor final


	      		ELSE IF (icontrol.EQ.2) THEN

c COMENTARI SENSE IMPORTÀNCIA:
c Nomes em refereixo a les parts DRETES de les igualtats seguents:

	  ! Els (i+1,j) (i,j+1) poden ser de la new o old. Valen el mateix.
	  ! Els (i-1,j) (i,j-1) han de ser new perque es el metode Jacobi "NUEVO"

	      			suma=told(i+1,j)+tnew(i-1,j)+told(i,j+1)+tnew(i,j-1)
	      			tnew(i,j)=(suma+rho*h**2.d0)/4.d0
	      			

	      		ELSE IF (icontrol.EQ.3) THEN

c COMENTARI SENSE IMPORTÀNCIA:
c Nomes em refereixo a les parts DRETES de les igualtats seguents:

	! Els (i+1,j) (i,j+1) poden ser de la new o old. Valen el mateix.
	! Els (i-1,j) (i,j-1) han de ser new perque aixi es el metode de sobrerelaxació
	! El (i,j) evidentment nomes pot ser old

	      			suma=told(i+1,j)+tnew(i-1,j)+told(i,j+1)+tnew(i,j-1)
	      			DELTA=(suma+rho*h**2.d0-4.d0*told(i,j))*omega/4.d0
	      			tnew(i,j)=told(i,j)+DELTA
	      			
	      		ENDIF
	      	  ENDDO

	      	ENDDO



C Per l'apartat 3, escrivim el nombre d'iteracions i la temperatura
c corresponent a la temperatura al punt (18. , 12.5)

	      	IF ((apartat.EQ.3).AND.(MOD(niter,100).EQ.0)) THEN
	      		WRITE(15,30) niter,tnew(ix,iy)
30    		FORMAT(I6,2X,2(E20.14,2X))
	      	ENDIF

	      error=-100000d0

	      DO i=1,nx-1
	      	  DO j=1,ny-1
	      	  	if (dabs(tnew(i,j)-told(i,j)).gt.error) then
	      	  		error=dabs(tnew(i,j)-told(i,j))
	      	  	endif
	      	 enddo
	      enddo

C   Li donem els valors nous a la matriu vella abans de tornar a recórrer tota la malla de punts
	      told=tnew
	      niter=niter+1

	      enddo !END DO WHILE

C  Un cop fetes totes les iteracions, per l'apartat 4,
c  escrivim cada punt de la malla i la temperatura corresponent per poder
c  graficar

	      IF (apartat.EQ.4) THEN
	      	DO i=0,nx
	 			DO j=0,ny
	      			WRITE(16,31) i*h,j*h,tnew(i,j)
31                FORMAT(2(E12.4,2X),E20.14)
	      		ENDDO
	      		WRITE (16,*)
	      	ENDDO

	      ENDIF
      END SUBROUTINE


