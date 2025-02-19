! La sintaxi en la que està escrit es Fortran Modern enlloc de 
! Fortran Fixed Form


program Pr8
	implicit none
	integer npassos,i,nmax,k,j
	double precision E,V,phi0,dphi0,e1,e2,e3,phi3(5000),dx,x
	double precision integrand(5000),Norm !ojo.
	              !el tamany de phi3 i integrand és suficient que sigui el mateix
	              !que nmax. Aixi que si es vol canviar nmax, cal nomes 
	  !donarli el valor desitjat a nmax i canviar el tamany de phi3
	  ! i integrand a aquest mateix.
	 
 !(Aixo ho fem perque fortran no feixa ficar un parametre en un common, i la idea es 
	   !tenirho tot en un common per passar el valor nmax a les subroutines.)

	common/dades/E,V,nmax
	nmax=5000 !tamany de les llistes. Esta pensat per que sgui mes gran
	        !que el que es demana. Millor que sobri que que falti
	V= -1.2d0

	open(1,file='res.dat')
	open(101,file='secant.dat')


	do k=1,2
	phi0=0.d0
	dphi0=0.25d0

	! k controla si fem per 50 passos o 640
		if (k.eq.1) npassos = 50
		if (k.eq.2) npassos = 640
	dx=1d0/real(npassos)
!Provem el mètode del tir per a diferents energies:
!________________________________________________________________________
!________________________________________________________________________
!________________________________________________________________________
! primer autovalor

	x=0.d0
	e1=7.5d0  
	e2=8.6d0
	call tirenergies(npassos,phi0,dphi0,e1,e2,e3,phi3)
	print*,'autoval1',e3
	! ara vole, integrar. COm el que sintegra es el valor abs 
	! de la func dona al quadrat, creem aquesta llista a integrand
	do j=1,nmax
		integrand(j)=dabs(phi3(j)*phi3(j))
	enddo
	!calculem la ctant de normalitzacio Norm
	if (k.eq.1) then
	call trapezoids(0d0,1d0,50,phi3*phi3,Norm)
	endif

	if (k.eq.2) then
	call trapezoids(0d0,1d0,640,phi3*phi3,Norm)
	endif

	do i=1,npassos
		x=x+dx
		write(1,*)x,phi3(i)/dsqrt(Norm)
	enddo
	write(1,*)
	write(1,*)


!________________________________________________________________________
!________________________________________________________________________
!________________________________________________________________________

!EL MATEIX PEL SEGON AUTOVALOR (lunic diferent respecte el primer autoval
                                   ! es e1 i e2)
	x=0.d0
	e1=18.5d0  
	e2=19d0

	call tirenergies(npassos,phi0,dphi0,e1,e2,e3,phi3)
	print*,'autoval2',e3
	! ara vole, integrar. COm el que sintegra es el valor abs 
	! de la func dona al quadrat, creem aquesta llista a integrand
	do j=1,nmax
		integrand(j)=dabs(phi3(j)*phi3(j))
	enddo
	!calculem la ctant de normalitzacio Norm
	if (k.eq.1) then
	call trapezoids(0d0,1d0,50,phi3*phi3,Norm)
	endif

	if (k.eq.2) then
	call trapezoids(0d0,1d0,640,phi3*phi3,Norm)
	endif

	do i=1,npassos
		x=x+dx
		write(1,*)x,phi3(i)/dsqrt(Norm)
	enddo
	write(1,*)
	write(1,*)

!________________________________________________________________________
!________________________________________________________________________
!________________________________________________________________________

!EL MATEIX PEL TERCER AUTOVALOR (lunic diferent respecte el primer autoval
                                   ! es e1 i e2)
	x=0.d0
	e1=43.5d0  
	e2=44d0

	call tirenergies(npassos,phi0,dphi0,e1,e2,e3,phi3)
	print*,'autoval3',e3
	! ara vole, integrar. COm el que sintegra es el valor abs 
	! de la func dona al quadrat, creem aquesta llista a integrand
	do j=1,nmax
		integrand(j)=dabs(phi3(j)*phi3(j))
	enddo
	!calculem la ctant de normalitzacio Norm
	if (k.eq.1) then
	call trapezoids(0d0,1d0,50,phi3*phi3,Norm)
	endif

	if (k.eq.2) then
	call trapezoids(0d0,1d0,640,phi3*phi3,Norm)
	endif

	do i=1,npassos
		x=x+dx
		write(1,*)x,phi3(i)/dsqrt(Norm)
	enddo
	write(1,*)
	write(1,*)
!________________________________________________________________________
!________________________________________________________________________
!________________________________________________________________________

!EL MATEIX PEL quart AUTOVALOR (lunic diferent respecte el primer autoval
                                   ! es e1 i e2)
	x=0.d0
	e1=77.5d0  
	e2=78d0

	call tirenergies(npassos,phi0,dphi0,e1,e2,e3,phi3)
	print*,'autoval4',e3
	! ara vole, integrar. COm el que sintegra es el valor abs 
	! de la func dona al quadrat, creem aquesta llista a integrand
	do j=1,nmax
		integrand(j)=dabs(phi3(j)*phi3(j))
	enddo
	!calculem la ctant de normalitzacio Norm
	if (k.eq.1) then
	call trapezoids(0d0,1d0,50,phi3*phi3,Norm)
	endif

	if (k.eq.2) then
	call trapezoids(0d0,1d0,640,phi3*phi3,Norm)
	endif

	do i=1,npassos
		x=x+dx
		write(1,*)x,phi3(i)/dsqrt(Norm)
	enddo

	write(1,*)
	write(1,*)

	enddo
	close(1)
	close(101)
	end program Pr8


!------------------------------------------------------------------
!			FUNCIONS I SUBRUTINES
!------------------------------------------------------------------

! Fem una subrutina que calculi un pas del mètode de Ralston 3er ordre.

subroutine RLSTN3(x,dx,nequs,yyin,yyout)
	implicit none
	integer nequs
	double precision yyin(nequs),yyout(nequs),x,x2,x3,dx
	double precision k1(nequs),k2(nequs),k3(nequs)
	call derivades(nequs,x,yyin,k1)

 !Usem x2 per si necessitessim una variable x (li dic x2) per evaluar 
 ! derivades. En aquesta prepractica, la subroutine derivades no necessita
 ! saber el valor de x, ja que la segona derivada, Fi'' es igual a 
 ! 2*Fi*(V-E). Pero si fos per exemple Fi''=2 * Fi * (V-E) * sin(x), SÍ
 ! necessitariem un valor de x. En aquest cas, per calcular k2, hem devaluar
 ! x a una distancia dx/2 respecte lanterior, segons el metode de Ralston. Per 
 ! aixo definim x2=x+dx/2. 

 	x2=x+dx/2.d0 
	call derivades(nequs,x2,yyin + dx*k1/2.d0,k2)

! La derivada es fa ara a (3/4)*dx del punt inicial x. Usem x3 = x+(3.d0/4.d0)*dx
	x3=x+(3.d0/4.d0)*dx
	call derivades(nequs,x3,yyin + dx*k2*3.d0/4.d0,k3)

!El resultat es
	yyout=yyin+dx/9d0*(2.d0*k1+3.d0*k2+4.d0*k3)
	end subroutine

! Fem una subrutina que ens calculará la primer i segona derivada de la nostre funció	
subroutine derivades(nequ,x,yin,dyout)
	implicit none
	integer nequ,nmax
	double precision x,yin(nequ),dyout(nequ),E,V
	common/dades/E,V,nmax
	dyout(1)=yin(2)
	dyout(2)= 2d0*yin(1)*(V-E)
	end subroutine

subroutine tirenergies(npassos,phi0,dphi0,e1,e2,e3,phi3)
	implicit none
	double precision yyin(2),yyout(2),x,dx,phi0,dphi0
	double precision phi1(nmax),phi2(nmax),phi3(nmax)
	double precision e1,e2,e3,E,V
	integer npassos,i,nmax
	common/dades/E,V,nmax
	dx=1.d0/real(npassos)
	phi3(npassos)=1.d0
	do while (dabs(phi3(npassos)).gt.(10.d0)**(-5d0))
	
!!! Calculem phi per a la primera energia
		E=e1
		yyin(1)=phi0
		yyin(2)=dphi0
		x=0.d0
		do i=1,npassos
			call RLSTN3(x,dx,2,yyin,yyout)
			x=x+dx
			yyin(1)=yyout(1)
			yyin(2)=yyout(2)
			phi1(i)=yyout(1)
		enddo

!!! El mateix per a la segona energia
		E=e2
		yyin(1)=phi0
		yyin(2)=dphi0
		x=0.d0
		do i=1,npassos
			call RLSTN3(x,dx,2,yyin,yyout)
			x=x+dx
			yyin(1)=yyout(1)
			yyin(2)=yyout(2)
			phi2(i)=yyout(1)
		enddo
		
!!! Calculem la tercera energia amb el mètode de la secant.
		e3=(e1*phi2(npassos)-e2*phi1(npassos))/(phi2(npassos)-phi1(npassos))
!!! I tornem a calcular phi per a la nova energia
		E=e3
		yyin(1)=phi0
		yyin(2)=dphi0
		x=0.d0
		do i=1,npassos
			call RLSTN3(x,dx,2,yyin,yyout)
			x=x+dx
			yyin(1)=yyout(1)
			yyin(2)=yyout(2)
			phi3(i)=yyout(1)
			enddo

		e1=e2
		e2=e3
		write(101,*) e3
		write(101,*) ''
		write(101,*) ''

		enddo

	end subroutine


!!! Utilitzarem una subrutina d trapezoides per poder normalitzar les nostres autofuncions.
subroutine trapezoids(x1,x2,k,funci,integral)
	implicit none
	double precision x1,x2,funci(nmax),h,F1,F2,integral,x,E,V
	integer k,i,nmax
	common/dades/E,V,nmax
	h=(x2-x1)/dble(k)
	F1=funci(1)
	F2=funci(k)
	x=x1
	integral=(h/2.d0)*(F1+F2)
	do i=1,(k-1)
		x=x+h
		integral=integral+funci(i)*h
		enddo
	end subroutine

