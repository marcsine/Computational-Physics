  	  implicit none
      double precision x0,dx,x1,h
      double precision pi,a0,densitat0,L
      double precision aux1,aux2,longitud,integral1,integral2
      double precision longi,densi,densi2
      real x2
      external longi, densi,densi2
  	  integer k,npunts1,i
      pi=4.D0*datan(1.d0)


  		 npunts1=2000
  		 x0=0.d0
  		 dx=0.02d0
 

  		 open(1,file='metode1.dat')
  		 open(2,file='metode2.dat')
  		 x1=x0
  		 x2=x0
  		 do i=0,200000000
  			if (mod(i,100000).eq.0) write(1,*) x1,x2
11		 format(e14.8) !(1P,e14.8)
  			x1=x1+dx
  			x2=x2+dx
  		 enddo
  		 print*, 'acabat 1' 
  		 close(1)
  		 x1=x0
  		 h=1000.d0
  		 do i=0,npunts1
        x1=2*h*i
        x2=2*h*i
  			write(2,*) x1,x2
  		 enddo
  		 print*,'acabat 2'
  		 close(2)
c      -----------------------APARTAT 2-------------------

       open(3,file='P4-1920-res1.dat')
       a0 = 0.35d0
       densitat0=0.42d0
       L=126.32d0/2d0
       k=18 
       call trapezoids(-pi,pi,k,longi,integral1)
       call simpson(-pi,pi,k,longi,integral2)
 12    FORMAT(a,E14.8,2X,E14.8,2X)
       write(3,12) 'longitud= ',a0*integral1,a0*integral2 !CAL POSAR FORMAT 12
       write(3,*)''
       write(3,*)''
       call trapezoids(-L,L,k,densi,integral1)
       call simpson(-L,L,k,densi,integral2)
       write(3,12) 'densitat lineal = ',integral1,integral2 !CAL POSAR FORMAT 12
       close(3)




C APARTAT 3 - PRIMERA PART

C Càlculs de l'apartat a) per trapezis i simpson en funció de h. 
c Fem un bucle variant k=4,5,...,22 tal que el nombre d'intervals es
c 2**4,2**5,...,2**22

 	  OPEN(13,FILE="p4-1920-res2.dat")
 	  DO k=4,22

 	  	CALL trapezoids(-pi,pi,k,longi,integral1)
        CALL simpson(-pi,pi,k,longi,integral2)
 	  	WRITE(13,*) dble(2**k),a0*integral1,a0*integral2 !CAL POSAR FORMAT 30
c 30     FORMAT(3(E14.8,2X))
      ENDDO
 		CLOSE(13)

C APARTAT 3 - SEGONA PART

 	  OPEN(14,FILE="p4-1920-res3.dat")
 	  DO k=4,22

 	  	CALL trapezoids(-L,L,k,densi,integral1)
        CALL simpson(-L,L,k,densi,integral2)
 	  	WRITE(14,*) dble(2**k),integral1,integral2 !CAL POSAR FORMAT 30
      ENDDO
 		CLOSE(14)

C APARTAT 4

      OPEN(15, FILE="p4-1920-res4.dat")
      do k=4,20
       call trapezoids(-pi/2d0,pi/2d0,k,densi2,integral1)
       call simpson(-pi/2d0,pi/2d0,k,densi2,integral2)
 	  	WRITE(15,*) dble(2**k),integral1,integral2 !CAL POSAR FORMAT 30
      ENDDO
 		CLOSE(15)
  	   end



c 		subrutina que em calcula el valor de integral entre x2 i x1 de la funcio funci
c 		fent servir la regla trapezoidal composta amb 2**k intervals
c      la funcio funci no es una llista, el que posare sera funci, es a dir el valor de la funcio en un punt
       subroutine trapezoids(x1,x2,k,funci,integral)
	       implicit none
	       double precision x1,x2,integral,h,F1,F2,x,funci
	       integer k,i
	       external funci

	       h = (x2-x1)/(2.d0**k)!interval h = (b-a)/(n)

	       F1= funci(x1)
	       
	       X=X1
	       	integral= (h/2.d0)*F1

	       DO i=1,int(2d0**k-1)
	       		X=X1+h*i
	       		integral=integral+funci(X)*h

	       ENDDO   

	       F2= funci(x2)

	       integral = integral + (h/2.d0)*F2
       end

       subroutine simpson(x1,x2,k,funci,integral)
	       implicit none
	       double precision x1,x2,funci,integral,h,x12,F1,F2,F12,x
	       integer k,i
	       external funci
	       h= (x2-x1)/(2.d0**k)!interval h = (b-a)/(n)
	       integral=0.d0
       	   do i=1,int(2**k)
       		 X=X1+i*h
       		 integral=h/6d0*(funci(x-h)+4d0*funci(x-h/2d0)+funci(x))+
     + 				integral
           ENDDO
       end
       

c_______________________________________________________________________
c_______________________________________________________________________

c     la funcio em retorna el valor de la funcio longitud en el valor x demanat
       double precision function longi(x)
	       implicit none
	       double precision x,aux1,aux2,pi
	       pi=4.D0*datan(1.d0)
	       aux1 = (((dcos(x-2.d0))**2)*dexp(-(x)**2+dsin(x)))**2
	       aux2 = dsqrt(pi-x)
	       longi = aux1*aux2
       end function
c     lo mateix amb la densitat




       double precision function densi(x)
	       implicit none
	       double precision x,densitat0,aux1,aux2,L
	       L=126.32d0/2.d0
	       densitat0=0.42d0
	       aux1=densitat0*(dsqrt(1d0-(x/L)**2d0))*(1d0-(x/L))
	       aux2=((x/L)**2d0)+(x/L)+ 1d0
	       densi=aux1*aux2
       end



       double precision function densi2(t)
	       implicit none
	       double precision x,densitat0,aux1,aux2,L,t
	       L=126.32d0/2.d0
	       
	       densitat0=0.42d0
	       aux1=densitat0*(dsqrt(1d0-(dsin(t))**2d0))*(1d0-(dsin(t)))
	       aux2=((dsin(t))**2d0)+(dsin(t))+ 1d0
	       densi2=aux1*aux2*L*dcos(t) !EL L*COS SURT DE dx=dx/dt * dt. dx/dt=Lcos(t)
       end











